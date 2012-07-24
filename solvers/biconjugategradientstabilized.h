//
// C++ Interface: biconjugategradientstabilized
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef BICONJUGATE_GRADIENT_H
#define BICONJUGATE_GRADIENT_H

#include "solver.h"

namespace Mu 
{
struct BiConjugateGradientStabilized : public LinearSolver
{
	virtual ~BiConjugateGradientStabilized() { } ;
	BiConjugateGradientStabilized(const Mu::CoordinateIndexedSparseMatrix& A_, Vector& b_) ;
	virtual bool solve(const Vector &x0, Preconditionner * precond = nullptr, const double eps = 1e-10, const int maxit = -1, bool verbose = false)  ;
} ;

} ;

#endif
