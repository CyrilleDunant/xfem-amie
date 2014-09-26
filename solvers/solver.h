//
// C++ Interface: solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUSOLVER_H
#define MUSOLVER_H

#include "../sparse/sparse_matrix.h"
#include "assembly.h"
#include "preconditionners.h"

namespace Amie {

class Assembly ;

	/** \brief Generic interface for solvers
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
		Generic interface for solvers
	*/
	class Solver
	{
	public:
		Solver();
		unsigned long colstart ;
		unsigned long rowstart ;
	
		virtual ~Solver();
		
		virtual bool solve(const Vector &x0 = Vector(0), Preconditionner * precond = nullptr, const double eps = 1e-9, const int maxit = -1, bool verbose = false)  = 0 ;
	
	};

	/**  \brief Generic interface for linear solvers
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	 * Generic interface for linear solvers. linear solvers 
	 * do not require knowledge about assemblies, only the matrices.
	*/
	struct LinearSolver : public Solver
	{
		virtual ~LinearSolver() { } ;
		const Vector & b ;
		const CoordinateIndexedSparseMatrix & A ;
		Vector x ;
		LinearSolver(const CoordinateIndexedSparseMatrix &, const Vector &) ;
		virtual bool solve(const Vector &x0, Preconditionner * precond = nullptr, const double eps = 1e-10, const int maxit = -1, bool verbose = false)  = 0 ;
	};
	
	/**   \brief Generic interface for non-linear solvers
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	 * Generic interface for non-linear solvers. NL solvers 
	 * do require knowledge about assemblies, not only the matrices.
	*/
	struct NonLinearSolver : public Solver
	{
		Assembly * assembly ;
		Vector x ;
		virtual ~NonLinearSolver() { } ;
		NonLinearSolver(Assembly * a) ;
		virtual bool solve(const Vector &x0 = Vector(0), Preconditionner * precond = nullptr, const double eps = 1e-10, const int maxit = -1, bool verbose = false)  = 0 ;
	};

} ;

#endif
