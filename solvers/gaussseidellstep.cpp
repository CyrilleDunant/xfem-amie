//
// C++ Implementation: gaussseidellstep
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gaussseidellstep.h"

namespace Mu {

GaussSeidellStep::GaussSeidellStep(const CoordinateIndexedSparseMatrix &A_) : b(0., A_.row_size.size()*(A_.stride+A_.stride%2)), gs(A_, b) { };

void GaussSeidellStep::precondition(const Vector& v, Vector& t) 
{

	gs.b = v ;
	gs.solve(v, NULL, 0, 1, false) ;
	t = gs.x ;
}


}
