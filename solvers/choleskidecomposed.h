//
// C++ Interface: choleskidecomposed
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

//
// C++ Interface: solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef CHOLESKI_DECOMPOSED_H
#define CHOLESKI_DECOMPOSED_H

#include "solver.h"

namespace Mu 
{

/** \brief direct solver for lower-triangular system*/
	struct LowerTriangular : public LinearSolver 
	{
		Vector d ;
		virtual ~LowerTriangular() { } ;
		LowerTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) ;
		virtual bool solve(const Vector &x0, const Preconditionner * precond= NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	} ;
	
/** \brief direct solver for upper-triangular system*/
	struct UpperTriangular : public LinearSolver 
	{
		Vector d ;
		virtual ~UpperTriangular() { } ;
		UpperTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) ;
		virtual bool solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	} ;
	
/** \brief Direct Solver for Symmetric Systems. The Matrix is assumed to have been Cholesky-decomposed*/
	struct CholeskiDecomposed : public LinearSolver
	{
		const Vector &d ;
		Vector y ;
		virtual ~CholeskiDecomposed() { } ;
		CholeskiDecomposed(const CoordinateIndexedSparseMatrix &A_, const Vector &b_, const Vector &d_) ;
		virtual bool solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	};

} ;

#endif
