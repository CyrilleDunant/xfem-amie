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

#include "stiffness_with_kinetics_imposed_deformation.h"
#include "../features/boundarycondition.h"

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

Vector StiffnessWithKineticsImposedDeformation::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

Vector StiffnessWithKineticsImposedDeformation::getImposedStrain(const Point & p) const
{
	return imposed ;
}

void StiffnessWithKineticsImposedDeformation::step(double timestep, ElementState & currentState)
{
	double fprevious = VirtualMachine().eval(kinetics, currentState.getTime()) ;
	double fnext = VirtualMachine().eval(kinetics,currentState.getTime()+timestep) ;
	fnext /= fprevious ;

	imposed *= fnext ;
}


std::vector<BoundaryCondition * > StiffnessWithKineticsImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
	
	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}


