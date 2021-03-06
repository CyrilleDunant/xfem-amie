//
// C++ Interface: polakribiereconjugategradient
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef POLAK_RIBIERE_GRADIENT_H
#define POLAK_RIBIERE_GRADIENT_H

#include "solver.h"
#include "conjugategradient.h"

namespace Amie
{

class Assembly ;

/** \brief Non-linear solver for symmetric systems*/
struct ConjugateGradientWithSecant : public NonLinearSolver
{
    size_t nssor = 20 ;

    virtual ~ConjugateGradientWithSecant() { } ;
    ConjugateGradientWithSecant ( Assembly * a, size_t n = 20 ) ;
    virtual bool solve ( const Vector &x0= Vector ( 0 ), Preconditionner * precond = nullptr, const double eps = 5e-8, const int maxit = -1, bool verbose = false )  ;
} ;


}

#endif

