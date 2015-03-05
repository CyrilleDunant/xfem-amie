//
// C++ Implementation: solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "solver.h"

namespace Amie {

Solver::Solver(Assembly * a) : colstart(0), rowstart(0), x(0., a->getMaxDofID()*a->getMatrix().stride), assembly(a)
{
}


Solver::~Solver()
{
}

LinearSolver::LinearSolver(Assembly * a) : Solver(a) {} ;


NonLinearSolver::NonLinearSolver(Assembly * a) : Solver(a) { } ;


}
