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

void GaussSeidellStep::precondition(const Vector &v,Vector &t) const
{
	t = GaussSeidel(A, v).solve(v, NULL, 1e-16, 1, false) ;
}

}
