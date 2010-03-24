//
// C++ Implementation: gaussseidellstep
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gaussseidellstep.h"
#include "gausseidell.h"

namespace Mu {

GaussSeidellStep::GaussSeidellStep(const CoordinateIndexedSparseMatrix &A_) : A(A_) { };

void GaussSeidellStep::precondition(Vector& v, Vector& t) const
{
	GaussSeidel gs(A, v) ;
	gs.solve(v, NULL, 0, 4, false) ;
	t = gs.x ;
}


}
