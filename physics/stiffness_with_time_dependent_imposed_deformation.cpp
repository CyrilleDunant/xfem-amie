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

#include "stiffness_with_time_dependent_imposed_deformation.h"

using namespace Mu ;

StiffnessWithTimeDependentImposedDeformation::StiffnessWithTimeDependentImposedDeformation(const Matrix & rig, Vector imposedDef, Function k) : StiffnessWithImposedDeformation(rig, imposedDef), kinetics(k)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	this->time_d = false ;
} ;

StiffnessWithTimeDependentImposedDeformation::~StiffnessWithTimeDependentImposedDeformation() { } ;

Matrix StiffnessWithTimeDependentImposedDeformation::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void StiffnessWithTimeDependentImposedDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithTimeDependentImposedDeformation::fractured() const
{
	return false ;
}

Form * StiffnessWithTimeDependentImposedDeformation::getCopy() const 
{
	return new StiffnessWithTimeDependentImposedDeformation(*this) ;
}

bool StiffnessWithTimeDependentImposedDeformation::hasInducedForces() const 
{
	return true ; 
} 

Vector StiffnessWithTimeDependentImposedDeformation::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

void StiffnessWithTimeDependentImposedDeformation::step(double timestep, ElementState & currentState)
{
	double fprevious = VirtualMachine().eval(kinetics, currentState.getTime()) ;
	double fnext = VirtualMachine().eval(kinetics,currentState.getTime()+timestep) ;
	fnext /= fprevious ;

	imposed[0] *= fnext ;
	imposed[1] *= fnext ;
	imposed[2] *= fnext ;
}


void StiffnessWithTimeDependentImposedDeformation::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
}

