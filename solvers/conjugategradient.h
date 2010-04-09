 //
// C++ Interface: conjugategradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "solver.h"

namespace Mu 
{
	class LinearSolver ;
/** \brief preconditionned Conjugate Gradient for symmetric systems*/
	struct ConjugateGradient : public LinearSolver
	{
		Vector r ;
		Vector z ;
		Vector p ;
		Vector q ;
		bool cleanup ;
		Preconditionner * P ;
		size_t nit ;
		virtual ~ConjugateGradient() { if(cleanup) delete P ;} ;
		ConjugateGradient(const Mu::CoordinateIndexedSparseMatrix& A_, Vector& b_) ;
		virtual bool solve(const Vector &x0, Preconditionner * precond = NULL, const double eps = 1e-7, const int maxit = -1, bool verbose = false)  ;
	} ;

} ;

#endif
