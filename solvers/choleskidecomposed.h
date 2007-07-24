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

	struct LowerTriangular : public LinearSolver 
	{
		Vector d ;
		virtual ~LowerTriangular() { } ;
		LowerTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) ;
		virtual Vector & solve(const Vector &x0, const Preconditionner * precond= NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	} ;
	
	struct UpperTriangular : public LinearSolver 
	{
		Vector d ;
		virtual ~UpperTriangular() { } ;
		UpperTriangular(const CoordinateIndexedSparseMatrix &A_, const Vector &b_) ;
		virtual Vector & solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	} ;
	
	struct CholeskiDecomposed : public LinearSolver
	{
		const Vector &d ;
		Vector y ;
		virtual ~CholeskiDecomposed() { } ;
		CholeskiDecomposed(const CoordinateIndexedSparseMatrix &A_, const Vector &b_, const Vector &d_) ;
		virtual Vector & solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-12 , const int maxit = -1, bool verbose = true)  ;
	};

} ;

#endif
