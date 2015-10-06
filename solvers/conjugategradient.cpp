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
    nit = 0 ;
    if(std::abs(assembly->getForces()).max() < eps)
    {
        x = 0 ;
        std::cerr << "\n CG "<< x.size() << " homogeneous." << std::endl ;
        return true ;
    }
    
    double realeps = std::max(1e-12, eps) ;
    size_t Maxit ;
    if(maxit != -1)
        Maxit = maxit ;
    else
        Maxit = assembly->getForces().size()*4 ;
    if(x0.size() == assembly->getForces().size())
    {
        x.resize(assembly->getForces().size());
        xmin.resize(assembly->getForces().size(), 0.);
        x = x0 ;
    }
    else
    {
        x.resize(assembly->getForces().size(), 0.);
        xmin.resize(assembly->getForces().size(), 0.);
        for(size_t i = 0 ; i < std::min(assembly->getForces().size(), x0.size()) ; i++)
            x[i] = x0[i] ;
    }


    if(rowstart)
    {
        for(size_t i = 0 ; i < rowstart ; i++)
            x[i] = assembly->getForces()[i] ;
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

    if(!precond)
    {
        cleanup = true ;
//      P = new InCompleteCholesky(A) ;
//         if(err0 < 1e-4)   
//             P = new Ssor(assembly->getMatrix(), 1.5) ;
//         else
            P = new InverseDiagonal(assembly->getMatrix()) ;
//      P = new InverseDiagonalSquared(assembly->getMatrix()) ;
//      P = new Inverse2x2Diagonal(assembly->getMatrix()) ;
        //   0.1      0.2   0.3   0.4   0.5     0.6   0.7   0.8     0.9  1.0  1.1   1.2   1.3   1.4   1.5   1.6  1.9
        //   505     16    15    16    10.6    15    14    10.6    15   14   10    11    10.3  10.2  10.6  10.7
        
//          P = new InverseLumpedDiagonal(assembly->getMatrix()) ;
//      P = new TriDiagonal(A) ;
//      P = new NullPreconditionner() ;
//      P = new GaussSeidellStep(A) ;
    }
    else
    {
        delete P ;
        cleanup = false ;
        P = precond ;
    }
    
    
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
    //****************************************

    assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
    #pragma omp parallel for schedule(static) if (vsize > 10000)
    for(int i = rowstart ; i < vsize ; i++)
        r[i] *= -1 ;

    err0 = std::max(1., sqrt( std::abs(parallel_inner_product(&r[rowstart], &z[rowstart], vsize-rowstart)))) ;
    if(verbose)
        std::cerr << 0 << "\t" << sqrt(err0) << std::endl  ;
    if(err0 < errmin)
    {
        errmin = err0 ;
        xmin = x ;
    }

#ifdef HAVE_OMP
    double t0 = omp_get_wtime() ;
#else

#endif
    double rho = 0 ;
    double beta = 0 ;
    double lastReset = rho ;
    int resetIncreaseCount = 0 ;
    bool pqconvergence = false ;
    Vector rcompensate(0., r.size()) ;
    Vector xcompensate(0., x.size()) ;
    while((sqrt(std::abs(last_rho)) > realeps*err0 && nit < Maxit ) || nit < 100)
    {
//         if(nit < 16)
//             err0 = std::max(sqrt(last_rho), err0) ;
        if(nit%512 == 0)
        {
             assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
             #pragma omp parallel for schedule(static) if (vsize > 10000)
             for(int i = rowstart ; i < vsize ; i++)
                 r[i] *= -1 ;
            rcompensate= 0 ;  
            std::cerr << resetIncreaseCount << "\t" << sqrt(last_rho) << std::endl  ;
        }

        P->precondition(r, z) ;

        rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;

        beta = rho/last_rho ;


        #pragma omp parallel for schedule(static) if (vsize > 10000)
        for(int i = rowstart ; i < vsize ; i++)
            p[i] = p[i]*beta+z[i] ;

        assign(q, assembly->getMatrix()*p, rowstart, colstart) ;
        pq =  parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
        if(std::abs(pq) < realeps*realeps*.25)
        {
            last_rho = rho ;
            pqconvergence = true ;
            break ;
        }
        alpha = rho/pq;


        #pragma omp parallel for schedule(static) if (vsize > 10000)
        for(int i = rowstart ; i < vsize ; i++)
        {
            double yr = -q[i]*alpha - rcompensate[i];
            double yx =  p[i]*alpha - xcompensate[i];
            double rtot = r[i]+yr ;
            double xtot = x[i]+yx ;
            rcompensate[i] = (rtot-r[i])-yr ;
            xcompensate[i] = (xtot-x[i])-yx ;
            r[i] = rtot ;
            x[i] = xtot ;
//             r[i] -= q[i]*alpha ;
//             x[i] += p[i]*alpha ;
        }

        if( sqrt(rho) < errmin )
        {
            errmin = sqrt(rho) ;
            xmin = x ;
        }


        if( verbose && nit%128 == 0 )
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
		rcompensate = 0 ;
		xcompensate= 0 ;
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
        if((nit <= Maxit && sqrt(last_rho)< realeps*err0) || pqconvergence)
            std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        else
        {
//             x = xmin ;
            std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        }
    }

    return (nit <= Maxit && sqrt(last_rho) < realeps*err0) || pqconvergence;
}

