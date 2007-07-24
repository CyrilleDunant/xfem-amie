//
// C++ Interface: gaussseidellstep
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUGAUSSSEIDELLSTEP_H
#define MUGAUSSSEIDELLSTEP_H

#include "preconditionners.h"
#include "../sparse/sparse_matrix.h"


namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
struct GaussSeidellStep : public Preconditionner
{

	const CoordinateIndexedSparseMatrix & A ;
	virtual ~GaussSeidellStep() { } ;
	GaussSeidellStep(const CoordinateIndexedSparseMatrix &A) ;
	virtual void precondition(const Vector &v,Vector &) const ;

};

}

#endif
