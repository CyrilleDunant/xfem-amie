//
// C++ Implementation: gausseidell
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gausseidell.h"

namespace Amie {

GaussSeidel::GaussSeidel(Assembly * a) :LinearSolver(a) { }

bool GaussSeidel::solve(const Vector &x0, Preconditionner * precond, const double eps, const int maxit, bool verbose)
{
    double err=10 ;
    x.resize(assembly->getForces().size()) ;
    if(x0.size() == assembly->getForces().size())
    {
        x = x0 ;
    }
    else
        x = 0 ;
// 	return true ;
    size_t nit=0 ;
    size_t Maxit ;
    if(maxit > 0)
        Maxit = maxit ;
    else
        Maxit = 2000*assembly->getForces().size() ;

// 	double omega = 1 ;

    Vector inverseDiagonal = assembly->getMatrix().inverseDiagonal() ;
    int stride =assembly->getMatrix().stride ;
    Vector xprev = x ;
    while((err > eps) && nit<Maxit)
    {
        for(size_t i = 0 ; i < x.size() ; i+=stride)
        {
            Vector delta = assembly->getMatrix()[i]*x ;
            for(size_t j = 0 ; j < delta.size() ; j++)
            {
                delta[j] -= x[i+j]/inverseDiagonal[i] ;
                x[i+j] = (assembly->getForces()[i+j] - delta[j]) * inverseDiagonal[i+j] ;
            }

        }

        err = std::abs(x-xprev).max() ;
        xprev = x ;
        if(nit%1000 == 0 && verbose)
        {
            std::cout << "error :"<< err <<", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
        }

        nit++ ;

    }

    if(nit< Maxit)
        return true ;
    return false ;
    if(verbose)
        std::cout << "converged after " << nit << " iterations. Error : " << err << ", max : "  << x.max() << ", min : "  << x.min() <<std::endl ;
}



}
