//
// C++ Implementation: conjugategradient
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//



#include "conjugategradient.h"
#include "inversediagonal.h"
#include "tridiagonal.h"
#include "ssor.h"
#include "gaussseidellstep.h"
#include "incompletecholeskidecomposition.h"
#include "eigenvalues.h"
#include <limits>
#include <sys/time.h>
#include <omp.h>

using namespace Amie ;

ConjugateGradient::ConjugateGradient( Assembly* a ) :LinearSolver(a), r(x.size()),z(x.size()),p(x.size()) ,q(x.size()), xmin(x.size()), cleanup(false), P(nullptr), nit(0) { }

bool ConjugateGradient::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{

    double realeps = std::max(1e-24, eps) ;
    size_t Maxit ;
    if(maxit != -1)
        Maxit = maxit ;
    else
        Maxit = x.size() ;
    if(x0.size() == assembly->getForces().size())
    {
        x.resize(x0.size());
        x = x0 ;
    }
    else
    {

        x = 0 ;
        for(size_t i = 0 ; i < std::min(assembly->getForces().size(), x0.size()) ; i++)
            x[i] = x0[i] ;
    }
// 	x = 0 ;


    if(rowstart)
    {
        for(size_t i = 0 ; i < rowstart ; i++)
            x[i] = assembly->getForces()[i] ;
    }
    InverseDiagonal P0(assembly->getMatrix()) ;

    
    if(precond == nullptr && !cleanup)
    {
        cleanup = true ;
// 		P = new InCompleteCholesky(A) ;
        P = new InverseDiagonal(assembly->getMatrix()) ;
//		P = new InverseDiagonalSquared(assembly->getMatrix()) ;
//		P = new Inverse2x2Diagonal(assembly->getMatrix()) ;
        //   0.1      0.2   0.3   0.4   0.5     0.6   0.7   0.8     0.9  1.0  1.1   1.2   1.3   1.4   1.5   1.6  1.9
        //   505     16    15    16    10.6    15    14    10.6    15   14   10    11    10.3  10.2  10.6  10.7
// 		P = new Ssor(A, 1.5) ;
//  		P = new InverseLumpedDiagonal(assembly->getMatrix()) ;
// 		P = new TriDiagonal(A) ;
// 		P = new NullPreconditionner() ;
// 		P = new GaussSeidellStep(A) ;
    }
    else if (precond != nullptr)
    {
        delete P ;
        cleanup = false ;
        P = precond ;
    }

    assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
    int vsize = r.size() ;
    double err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
    xmin = x ;
    errmin = err0 ;
    r*=-1 ;

    if (err0 < realeps)
    {
        if(verbose)
            std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        return true ;
    }
    //*************************************


    z = r ;
    P->precondition(r,z) ;

    p = z ;
    q = assembly->getMatrix()*p ;

    double last_rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;
    double pq = parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
    double alpha = last_rho/pq ;

    #pragma omp parallel for schedule(static) if (vsize > 10000)
    for(int i = rowstart ; i < vsize ; i++)
    {
        r[i] -= q[i]*alpha ;
        x[i] += p[i]*alpha ;
    }
    nit++ ;
    //****************************************

    assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
    #pragma omp parallel for schedule(static) if (vsize > 10000)
    for(int i = rowstart ; i < vsize ; i++)
        r[i] *= -1 ;

    err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
    if(err0 < errmin)
    {
        errmin = err0 ;
        xmin = x ;
    }

    if (err0 < realeps)
    {
        if(verbose)
            std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

        return true ;
    }
#ifdef HAVE_OMP
    double t0 = omp_get_wtime() ;
#else

#endif
// 	double neps = /*std::min(*/realeps*realeps/*, err0*realeps)*/ ; //std::max(err0*realeps, realeps*realeps) ;
    double rho = 0 ;
    double beta = 0 ;
    double lastReset = rho ;
    int resetIncreaseCount = 0 ;
    while((last_rho*last_rho > std::max(realeps*realeps*err0, realeps*realeps) && nit < Maxit ) || nit < 16)
    {
//             if(nit < 256)
        P->precondition(r, z) ;
//             else
//                 P0.precondition(r, z) ;

        rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;


        beta = rho/last_rho ;


        #pragma omp parallel for schedule(static) if (vsize > 10000)
        for(int i = rowstart ; i < vsize ; i++)
            p[i] = p[i]*beta+z[i] ;

        assign(q, assembly->getMatrix()*p, rowstart, colstart) ;
        pq =  parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
        alpha = rho/pq;

        if(std::abs(pq) < POINT_TOLERANCE*POINT_TOLERANCE)
        {
            last_rho = 0 ;
            break ;
        }
        #pragma omp parallel for schedule(static) if (vsize > 10000)
        for(int i = rowstart ; i < vsize ; i++)
        {
            r[i] -= q[i]*alpha ;
            x[i] += p[i]*alpha ;
        }

        if(sqrt(rho) < errmin)
        {
            errmin = sqrt(rho) ;
            xmin = x ;
        }

        if(nit%256 == 0)
        {
            assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
            #pragma omp parallel for schedule(static) if (vsize > 10000)
            for(int i = rowstart ; i < vsize ; i++)
                r[i] *= -1 ;
        }
        if( verbose && nit%256 == 0 )
        {
            if(rho > lastReset)
                resetIncreaseCount++ ;
            else
                resetIncreaseCount = 0 ;

            lastReset = rho ;

            std::cerr << resetIncreaseCount << "\t" << sqrt(rho) << std::endl  ;

            if(resetIncreaseCount > maxIncreaseReset)
            {
                x = (xmin+x)*0.5 ;
                assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
                #pragma omp parallel for schedule(static) if (vsize > 10000)
                for(int i = rowstart ; i < vsize ; i++)
                    r[i] *= -1 ;
                resetIncreaseCount = 0 ;
            }
        }


        last_rho = rho ;
        nit++ ;
    }
#ifdef HAVE_OMP
    double delta = std::max((omp_get_wtime() - t0)*1e6, 1e-14) ;
#else
    double delta = 1 ;
#endif

    std::cerr << "mflops: "<< nit*((2.+2./256.)*assembly->getMatrix().array.size()+(4+1./256.)*p.size())/delta << std::endl ;

    assign(r,assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
    double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;

    if(verbose)
    {
        if(nit <= Maxit && last_rho*last_rho< std::max(realeps*realeps*err0, realeps*realeps))
            std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        else
        {
            x = xmin ;
            std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        }
    }

    return nit <= Maxit && last_rho*last_rho < std::max(realeps*realeps*err0, realeps*realeps);
}

