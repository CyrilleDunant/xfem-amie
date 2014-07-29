//
// C++ Interface: gausseidell
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUGAUSSEIDELL_H
#define MUGAUSSEIDELL_H

#include "solver.h"

namespace Amie 
{

/** \brief Gauß-Seidell solver for symmetric systems
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	This is an implementation of the classic Gauß-Seidell iteration
*/
struct GaussSeidel : public LinearSolver
{
	virtual ~GaussSeidel() { } ;
	GaussSeidel(const Amie::CoordinateIndexedSparseMatrix& A_, Vector& b_) ;
	virtual bool solve(const Vector &x0, Preconditionner * precond= nullptr, const double eps = 1e-20, const int maxit = -1, bool verbose = false)  ;
} ;

} ;

#endif
