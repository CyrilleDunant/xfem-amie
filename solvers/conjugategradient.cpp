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
//     nssor = 32 ;
    size_t attempts = 0 ;
    
    if(std::abs(assembly->getForces()).max() < eps*eps)
    {
        std::cerr << "\n CG "<< x.size() << " homogeneous. " << std::abs(assembly->getForces()).max() << std::endl ;
        return true ;
    }
    
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
    
    double realeps = std::max(1e-12, eps) ;
    size_t Maxit = (maxit != -1) ? maxit : assembly->getForces().size()/2 ;
    
    if(x0.size() == assembly->getForces().size())
    {
        x = x0 ;
    }
    else
    {
        #pragma omp parallel for 
        for(size_t i = 0 ; i < std::min(assembly->getForces().size(), x0.size()) ; i++)
            x[i] = x0[i] ;
    }

    if(rowstart)
    {
        #pragma omp parallel for 
        for(size_t i = 0 ; i < rowstart ; i++)
            x[i] = assembly->getForces()[i] ;
    }
    
    int vsize = r.size() ;
//     Ssor ssor (assembly->getMatrix(), rowstart, colstart) ;
    Vector rcompensate(0., r.size()) ;
    Vector xcompensate(0., x.size()) ; 
    nit = 0 ;
    double perr = 1 ;
    int multi = 1 ;
    xmin = x ;
    for( ; nit < Maxit ; )
    {
        size_t localnit = 0 ;
        rcompensate = 0 ;
        xcompensate = 0 ;
        
        //smooth the initial guess
        if(nssor)
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
                P->precondition(r,r);
                #pragma omp parallel for 
                for(int i = rowstart ; i < vsize ; i++)
                {
                    double yx =  -r[i] - xcompensate[i];
                    double xtot = x[i]+yx ;
                    xcompensate[i] = (xtot-x[i])-yx ;
                    x[i] = xtot ;
                }
            }
            xcompensate = 0 ;
        }
        assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
        
        #pragma omp parallel for 
        for(int i = rowstart ; i < vsize ; i++)
            r[i] *= -1 ;
    
        double err0 = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
        if(std::isnan(err0))
        {
            assembly->print() ;
            std::cout << "NaN error" << std::endl ;
            exit(0) ;
        }
        
        if(nit == 0)
            errmin = err0 ;
        
        if (err0 < realeps)
        {
            std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", last rho = " << 0 << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

            return true ;
        }
        
        std::cerr << "p\t" << err0 << std::endl  ;
        if(nit && std::abs(err0-perr) < 1e-4*std::max(err0, perr))
            multi *=2 ;
        else if(nit)
            multi = 1 ;
        perr = err0 ;
        z = r ;
        P->precondition(r,z) ;

        p = z ;
        q = assembly->getMatrix()*p ;

        double last_rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;
        double pq = parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
        if(std::abs(pq) < 1e-12*last_rho)
        {
             std::cerr << "\n CG "<< p.size() << " converged after " << nit << " iterations. Error : " << err0 << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;

            return true ;
        }
        double alpha = last_rho/pq ;

        #pragma omp parallel for 
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

    #ifdef HAVE_OPENMP
        double t0 = omp_get_wtime() ;
    #endif
        double rho = 0 ;
        double beta = 0 ;

        while(sqrt(std::abs(last_rho)) > realeps && localnit < assembly->getForces().size())
        {
            localnit++ ;

            P->precondition(r, z) ;
            rho = parallel_inner_product_restricted(&r[rowstart], &z[rowstart], vsize-rowstart) ;

            beta = rho/last_rho ;

            #pragma omp parallel for 
            for(int i = rowstart ; i < vsize ; i++)
	    {
                p[i] = p[i]*beta+z[i] ;
	    }

            assign(q, assembly->getMatrix()*p, rowstart, colstart) ;
            pq =  parallel_inner_product_restricted(&q[rowstart], &p[rowstart], vsize-rowstart);
            if(std::abs(pq) < 1e-24*rho)
            {
                last_rho = rho ;
                break ;
            }
            alpha = rho/pq;

            #pragma omp parallel for
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

            last_rho = rho ;
            nit++ ;
        }
    #ifdef HAVE_OPENMP
        double delta = std::max((omp_get_wtime() - t0)*1e6, 1e-32) ;
    #else
        double delta = 1 ;
    #endif

        std::cerr << "mflops: "<< nit*(2.*assembly->getMatrix().array.size()+4.*p.size())/delta << std::endl ;

        assign(r,assembly->getMatrix()*x-assembly->getForces(), rowstart, rowstart) ;
        double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
        if(err < errmin )
        {
            errmin = sqrt(std::abs(rho)) ;
            xmin = x ;
        }

        if(nssor)
        {
            size_t iters = 0 ;
            while( err > realeps && iters++ < assembly->getForces().size() )
            {
                assign(r, assembly->getMatrix()*x-assembly->getForces(), rowstart, colstart) ;
                P->precondition(r,r);
                #pragma omp parallel for 
                for(int i = rowstart ; i < vsize ; i++)
                {
                    double yx =  -r[i] - xcompensate[i];
                    double xtot = x[i]+yx ;
                    xcompensate[i] = (xtot-x[i])-yx ;
                    x[i] = xtot ;
                }
                
                double perr = err;
                err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
                if(iters%256 == 0 )
                    std::cerr << "m" << "\t" << err << std::endl  ;
                if( err > perr )
                    break ;
                else
                    xmin = x ;
            }
            if( iters > 2 )
                x = xmin ;

        }

        
        if( std::min(err,sqrt(std::abs(last_rho))) < sqrt(realeps) )
        {
            std::cerr << "\n CG " << p.size() << " converged after " << nit << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
            return true ;
        }


    }
    
    double err = sqrt( parallel_inner_product(&r[rowstart], &r[rowstart], vsize-rowstart)) ;
    double last_rho = err ;
    std::cerr << "\n CG " << p.size() << " did not converge after " << nit*attempts << " iterations. Error : " << err << ", last rho = " << last_rho << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
    return false ;
}

}
