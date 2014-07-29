//
// C++ Interface: gaussseidellstep
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUGAUSSSEIDELLSTEP_H
#define MUGAUSSSEIDELLSTEP_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"
#include "gausseidell.h"


namespace Amie {

/** \brief Preconditionner, perform a GS step
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
struct GaussSeidellStep : public Preconditionner
{
	Vector b ;
	GaussSeidel gs ;
	virtual ~GaussSeidellStep() { } ;
	GaussSeidellStep(const CoordinateIndexedSparseMatrix &A) ;
	virtual void precondition(const Vector& v, Vector& t)  ;

};

}

#endif
