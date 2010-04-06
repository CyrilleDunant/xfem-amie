//
// C++ Implementation: solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "solver.h"

namespace Mu {

Solver::Solver()
{
}


Solver::~Solver()
{
}

LinearSolver::LinearSolver(const CoordinateIndexedSparseMatrix &A_,Vector &b_) : b(b_), A(A_), x(0., b.size()) {} ;


NonLinearSolver::NonLinearSolver(Assembly * a) : assembly(a) { } ;


}
