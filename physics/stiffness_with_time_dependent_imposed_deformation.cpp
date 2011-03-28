//
// C++ Implementation: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
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


std::vector<BoundaryCondition * > StiffnessWithTimeDependentImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
	
	std::vector<BoundaryCondition * > ret ;
	ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, s.getParent(), id, f[0]);
	ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, s.getParent(), id, f[1]);
	if(f.size() == 3)
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, s.getParent(), id, f[2]);

	return ret ;
}


