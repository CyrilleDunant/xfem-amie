//
// C++ Interface: inversediagonal
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INV_DIAG_H
#define __INV_DIAG_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"

namespace Mu
{

	/** \brief Preconditionner: inverse diagonal*/
	struct InverseDiagonal  : public Preconditionner
	{
		Vector diagonal ;
		virtual ~InverseDiagonal() { }
		InverseDiagonal(const CoordinateIndexedSparseMatrix &A) ;
		virtual void precondition(const Vector &v,Vector &)  ;
	} ;
	
	/** \brief Preconditionner: inverse lumped diagonal*/
	struct InverseLumpedDiagonal  : public Preconditionner
	{
		Vector diagonal ;
		virtual ~InverseLumpedDiagonal() { }
		InverseLumpedDiagonal(const CoordinateIndexedSparseMatrix &A) ;
		virtual void precondition(const Vector &v,Vector &)  ;
	} ;
	
	/** \brief Preconditionner: inverse diagonal squared*/
	struct InverseDiagonalSquared  : public Preconditionner
	{
		Vector * diagonal ;
		virtual ~InverseDiagonalSquared() { delete diagonal ; }
		InverseDiagonalSquared(const CoordinateIndexedSparseMatrix &A) ;
		virtual void precondition(const Vector &v,Vector &)  ;
	} ;

} ;

#endif
