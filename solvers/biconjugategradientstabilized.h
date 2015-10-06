//
// C++ Interface: biconjugategradientstabilized
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef BICONJUGATE_GRADIENT_H
#define BICONJUGATE_GRADIENT_H

#include "solver.h"

namespace Amie
{
struct BiConjugateGradientStabilized : public LinearSolver
{
    virtual ~BiConjugateGradientStabilized() { } ;
    BiConjugateGradientStabilized(Assembly * a) ;
    virtual bool solve(const Vector &x0, Preconditionner * precond = nullptr, const double eps = default_solver_precision, const int maxit = -1, bool verbose = false)  ;
} ;

}

#endif
