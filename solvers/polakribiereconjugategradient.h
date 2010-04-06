//
// C++ Interface: polakribiereconjugategradient
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef POLAK_RIBIERE_GRADIENT_H
#define POLAK_RIBIERE_GRADIENT_H

#include "solver.h"
#include "conjugategradient.h"

namespace Mu
{

	struct Assembly ;

/** \brief Non-linear solver for symmetric systems*/
	struct ConjugateGradientWithSecant : public NonLinearSolver
	{
		virtual ~ConjugateGradientWithSecant() { } ;
		ConjugateGradientWithSecant ( Assembly * a ) ;
		virtual bool solve ( const Vector &x0= Vector ( 0 ), const Preconditionner * precond = NULL, const double eps = 5e-8, const int maxit = -1, bool verbose = false )  ;
	} ;


} ;

#endif

