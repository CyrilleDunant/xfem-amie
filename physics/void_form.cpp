//
// C++ Implementation: void_form
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "void_form.h"

namespace Amie {


VoidForm::VoidForm() : LinearForm(Matrix(1,1),false, false, 0 )
{
    num_dof = 0 ;
    type = VOID_BEHAVIOUR ;
    time_d = false ;
}

void VoidForm::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
}

void VoidForm::step(double timestep, ElementState & currentState, double maxscore)
{

}

void VoidForm::updateElementState(double timestep, ElementState & s) const
{

}

bool VoidForm::fractured() const
{
    return false ;
}

VoidForm::~VoidForm() { }

Form * VoidForm::getCopy() const
{
    return new VoidForm(*this) ;
}

}
