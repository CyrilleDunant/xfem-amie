//
// C++ Implementation: void_form
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "void_form.h"

using namespace Mu ;


VoidForm::VoidForm() : LinearForm(Matrix(1,1),false, false, 0 )
{
	this->num_dof = 0 ;
	this->type = VOID_BEHAVIOUR ;
	this->time_d = false ;
}

void VoidForm::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
}

void VoidForm::step(double timestep, const ElementState & currentState) 
{
	
}

void VoidForm::updateElementState(double timestep, ElementState & s) const
{
	
}

bool VoidForm::fractured() const 
{
	return false ;
}

VoidForm::~VoidForm() { } ;

Form * VoidForm::getCopy() const 
{
	return new VoidForm(*this) ;
}

