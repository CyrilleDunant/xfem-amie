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
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace Amie {

ConjugateGradient::ConjugateGradient( Assembly* a ) :LinearSolver(a), r(x.size()),z(x.size()),p(x.size()) ,q(x.size()), xmin(x.size()), cleanup(false), P(nullptr), nit(0) 
{ 
//     if(assembly->getMatrix().stride > 3)
//         nssor = 0 ;
}
/*
double ConjugateGradient::ssorSmooth(Ssor * ssor, Vector * xcompensate, size_t maxIt, double epsilon)
{
    if(nssor > 0)
    {
        *xcompensate = 0 ;
        assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
        double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], xcompensate->size()-rowstart)) ;
        double perr = err ;
        size_t iter = 0 ;
        while(iter++ < maxIt && err > epsilon)
        {    
            ssor->precondition(r,r);
            for(size_t i = rowstart ; i < xcompensate->size() ; i++)
            {
                double yx =  -r[i] - (*xcompensate)[i];
                double xtot = x[i]+yx ;
                (*xcompensate)[i] = (xtot-x[i])-yx ;
                x[i] = xtot ;
            }
            assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
            err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], xcompensate->size()-rowstart)) ;
            std::cout << err << std::endl ;
            perr = std::min(err, perr) ;
            if(err > perr*1.01)
                break ;
        }
    }
    assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
    r*=-1 ;
    
    return sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], xcompensate->size()-rowstart)) ;
}*/

bool ConjugateGradient::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    nssor = 64 ;
    nit = 0 ;
    
    if(std::abs(assembly->getForces()).max() < eps)
    {
        std::cerr << "\n CG "<< x.size() << " homogeneous." << std::endl ;
        return true ;
    }
    
    double realeps = std::max(1e-12, eps) ;
    size_t Maxit ;
    if(maxit != -1)
        Maxit = maxit ;
    else
        Maxit = assembly->getForces().size()*10. ;
    
    if(x0.size() == assembly->getForces().size())
    {
        x = x0 ;
    }
    else
    {
        for(size_t i = 0 ; i < std::min(assembly->getForces().size(), x0.size()) ; i++)
            x[i] = x0[i] ;
    }

    if(rowstart)
    {
        for(size_t i = 0 ; i < rowstart ; i++)
            x[i] = assembly->getForces()[i] ;
    }
    int vsize = r.size() ;
     
    //smooth the initial guess
    Ssor ssor (assembly->getMatrix(), rowstart, colstart) ;
    Vector rcompensate(0., r.size()) ;
    Vector xcompensate(0., x.size()) ; 
    if(nssor > 0)
    {
        double err = 2 ;
        double perr = 0 ;
        size_t iter = 0 ;
        while(iter++ < nssor && err > realeps)
        {

            assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
            perr = err ;
            err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
            if(err > perr)
                break ;
            ssor.precondition(r,r);
            for(int i = rowstart ; i < vsize ; i++)
            {
                double yx =  -r[i] - xcompensate[i];
                double xtot = x[i]+yx ;
                xcompensate[i] = (xtot-x[i])-yx ;
                x[i] = xtot ;
            }
        }
    }
    assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
    r*=-1 ;
   
    double err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
    xmin = x ;
    errmin = err0 ;
    
    //try an alternative start
    if(err0 > 1e12)
    {
        x = 0 ;
        if(nssor > 0)
        {
            double err = 2 ;
            double perr = 0 ;
            size_t iter = 0 ;
            while(iter++ < nssor && err > realeps)
            {

                assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
                perr = err ;
                err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
                if(perr > err)
                    break ;
                ssor.precondition(r,r);
                for(int i = rowstart ; i < vsize ; i++)
                {
                    double yx =  -r[i] - xcompensate[i];
                    double xtot = x[i]+yx ;
                    xcompensate[i] = (xtot-x[i])-yx ;
                    x[i] = xtot ;
                }
            }
        }
        assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
        r*=-1 ;
    
        err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
        xmin = x ;
        errmin = err0 ;
    }

    if (err0 < realeps)
    {
        if(verbose)
            std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", last rho = " << 0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        return true ;
    }
    if(verbose)
        std::cerr << "p" << "\t" << err0 << std::endl  ;

    if(!precond)
    {
        cleanup = true ;
         P = new InverseDiagonal(assembly->getMatrix()) ;
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
    if(std::abs(pq) < realeps*realeps*last_rho)
    {
        if(verbose)
            std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        return true ;
    }
    double alpha = last_rho/pq ;

    #pragma omp parallel for schedule(static) if (vsize > 10000)
    for(int i = rowstart ; i < vsize ; i++)
    {
        r[i] -= q[i]*alpha ;
        x[i] += p[i]*alpha ;
    }

    err0 = 1 ;

#ifdef HAVE_OMP
    double t0 = omp_get_wtime() ;
#else

#endif
    double rho = 0 ;
    double beta = 0 ;
    double lastReset = rho ;
    int resetIncreaseCount = 0 ;

    while((sqrt(std::abs(last_rho)) > realeps*err0*5. && nit < Maxit ))
    {
        if(nit && nit%256 == 0)
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
        if(std::abs(pq) < realeps*realeps*last_rho)
        {
            last_rho = rho ;
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
        if(sqrt(rho) > 1./realeps && nit > Maxit/10)
        {
            x = xmin ;
            break ;
        }

        last_rho = rho ;
        nit++ ;
    }
#ifdef HAVE_OMP
    double delta = std::max((omp_get_wtime() - t0)*1e6, 1e-14) ;
#else
    double delta = 1 ;
#endif

    std::cerr << "mflops: "<< 1e-6*nit*((2.+2./256.)*assembly->getMatrix().array.size()+(4+1./256.)*p.size())/delta << std::endl ;
    
    assign(r,assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
    double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;

    xcompensate= 0 ;
    if(nssor > 0)
    {
        double minerr = err*2. ;
        while(err > realeps*err0)
        {
            assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
            ssor.precondition(r,r);
            for(int i = rowstart ; i < vsize ; i++)
            {
                double yx =  -r[i] - xcompensate[i];
                double xtot = x[i]+yx ;
                xcompensate[i] = (xtot-x[i])-yx ;
                x[i] = xtot ;
            }
            
            minerr = std::min(err, minerr) ;
            err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
            
            if(err > minerr)
                break ;
        }
    }
    
    if(verbose)
    {
        if(err < sqrt(realeps*err0))
            std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        else
        {
            std::cerr << "\n CG " << p.size() << " did not converge after " << nit << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        }
    }

    return err <  sqrt(realeps*err0);
}

}
