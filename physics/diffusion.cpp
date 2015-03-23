//
// C++ Implementation: stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "diffusion.h"

using namespace Amie ;

Diffusion::Diffusion(const Matrix & rig) : LinearForm(rig, false, false, 1)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 9)
        v.push_back(ZETA);

    v.push_back(TIME_VARIABLE);

}

Diffusion::~Diffusion() { }

void Diffusion::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{

    ret[0][0] = (vm->ieval(VectorGradientDot(p_i) * param * VectorGradient(p_j, true),  gp, Jinv, v)
                 + vm->ieval(VectorGradient(p_i) * param * VectorGradientDot(p_j, true),  gp, Jinv, v)) ;
}

void Diffusion::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{

    ret[0][0] =  (vm->ieval(Differential(p_i, TIME_VARIABLE) * Differential(p_j, TIME_VARIABLE),  gp, Jinv, v)
                  + vm->ieval(DoubleDifferential(p_j, TIME_VARIABLE, TIME_VARIABLE) * p_i,  gp, Jinv, v)) ;
}

bool Diffusion::fractured() const
{
    return false ;
}

Form * Diffusion::getCopy() const
{
    Diffusion * copy = new Diffusion(*this) ;

    return copy ;
}

