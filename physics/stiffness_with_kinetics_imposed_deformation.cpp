//
// C++ Implementation: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_kinetics_imposed_deformation.h"

using namespace Mu ;

StiffnessWithKineticsImposedDeformation::StiffnessWithKineticsImposedDeformation(const Matrix & rig, Vector imposedDef, Function k) : StiffnessWithImposedDeformation(rig, imposedDef), kinetics(k)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	this->time_d = true ;
} ;

StiffnessWithKineticsImposedDeformation::~StiffnessWithKineticsImposedDeformation() { } ;

Matrix StiffnessWithKineticsImposedDeformation::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void StiffnessWithKineticsImposedDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithKineticsImposedDeformation::fractured() const
{
	return false ;
}

Form * StiffnessWithKineticsImposedDeformation::getCopy() const 
{
	return new StiffnessWithKineticsImposedDeformation(*this) ;
}

bool StiffnessWithKineticsImposedDeformation::hasInducedForces() const 
{
	return true ; 
} 

Vector StiffnessWithKineticsImposedDeformation::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

void StiffnessWithKineticsImposedDeformation::step(double timestep, ElementState & currentState)
{
	double fprevious = VirtualMachine().eval(kinetics, currentState.getTime()) ;
	double fnext = VirtualMachine().eval(kinetics,currentState.getTime()+timestep) ;
	fnext /= fprevious ;

	imposed *= fnext ;
}


void StiffnessWithKineticsImposedDeformation::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
}

