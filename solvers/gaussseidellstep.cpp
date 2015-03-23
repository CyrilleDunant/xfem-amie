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

namespace Amie {

GaussSeidellStep::GaussSeidellStep(Assembly * a) : b(0., a->getMatrix().row_size.size()*(a->getMatrix().stride+a->getMatrix().stride%2)), gs(a) { }

void GaussSeidellStep::precondition(const Vector& v, Vector& t)
{
    GaussSeidel gs_(gs.assembly);
    gs_.solve(v, nullptr, 0, 1, false) ;
    t = gs_.x ;
}


}
