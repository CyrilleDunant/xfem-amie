//
// C++ Interface: inclompletecholeskidecomposition
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCOMPLETE_CHOLESKI_H
#define __INCOMPLETE_CHOLESKI_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

namespace Mu
{

	struct InCompleteCholesky  : public Preconditionner
	{
		bool stable ;
		Vector d ;
		CoordinateIndexedSparseMatrix A ;
		virtual ~InCompleteCholesky() { } ;
		InCompleteCholesky(const CoordinateIndexedSparseMatrix &A) ;
		virtual void precondition(const Vector &v,Vector &)  ;
	} ;


} ;

#endif 
