//
// C++ Interface: conjugategradient
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "solver.h"

namespace Amie
{
class LinearSolver ;
/** \brief preconditionned Conjugate Gradient for symmetric systems*/
struct ConjugateGradient : public LinearSolver
{
    Vector r ;
    Vector z ;
    Vector p ;
    Vector q ;
    Vector xmin ;
    double errmin = 1e9 ;
    int maxIncreaseReset = 16 ;
    bool cleanup ;
    Preconditionner * P ;
    size_t nit ;
    virtual ~ConjugateGradient() {
        if(cleanup) delete P ;
    } ;
    ConjugateGradient(Assembly * a) ;
    virtual bool solve(const Vector &x0, Preconditionner * precond = nullptr, const double eps = 1e-10, const int maxit = -1, bool verbose = false)  ;
} ;

}

#endif
