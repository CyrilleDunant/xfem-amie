//
// C++ Interface: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __TRI_DIAG_H
#define __TRI_DIAG_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

namespace Mu
{

	/** \brief Tridiagonal preconditionner*/
	struct TriDiagonal  : public Preconditionner
	{
		Vector diagonal ;
		Vector upper ;
		virtual ~TriDiagonal() {  ;}
		TriDiagonal(const CoordinateIndexedSparseMatrix &A) ;
		virtual void precondition(const Vector &v,Vector &) const ;
	} ;
	

} ;

#endif
