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

#include "stiffness_with_variable_imposed_deformation_and_fracture.h"
#include <limits>

using namespace Mu ;

StiffnessWithVariableImposedDeformationAndFracture::StiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * c) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef), dfunc(rig.numRows()-1),criterion(c)
{
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previousDamage = 0 ;
	damage = 0 ;

	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

StiffnessWithVariableImposedDeformationAndFracture::~StiffnessWithVariableImposedDeformationAndFracture() { } ;

Matrix StiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{

	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void StiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithVariableImposedDeformationAndFracture::hasInducedForces() const 
{
	return true ; 
} 

void StiffnessWithVariableImposedDeformationAndFracture::step(double timestep, ElementState & currentState)
{
	if(!frac && criterion->met(currentState) )
	{
		dfunc.step(currentState) ;
		previousDamage = damage ;
		
		Vector state = dfunc.damageState() ;
		damage = 0 ;
		for(size_t i = 0 ; i < state.size() ; i++)
			damage += state[i] ;
		currentState.getParent()->behaviourUpdated = true ;
		if(damage > .9)
		{
			frac = true ;
// 			damage = .9999 ;
// 			param[0][1] = 0 ;param[0][1] = 0 ;
// 			param[2][2] *= 0.0001 ;
// 			this->type = VOID_BEHAVIOUR ;
		}
	}

	if(!currentState.getParent()->behaviourUpdated && timestep > std::numeric_limits<double>::epsilon())
	{
		double randomVar = (double)rand()/(double)RAND_MAX ;
		imposed[0] += timestep*randomVar ;
		imposed[1] += timestep*randomVar ;
		currentState.getParent()->behaviourUpdated = true ;
	}

	if(frac)
	{
		imposed = 0 ;
	}

	change = currentState.getParent()->behaviourUpdated  ;
}

bool StiffnessWithVariableImposedDeformationAndFracture::changed() const
{
	return change ;
} 

bool StiffnessWithVariableImposedDeformationAndFracture::fractured() const
{
	return frac;
}

Vector StiffnessWithVariableImposedDeformationAndFracture::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

Form * StiffnessWithVariableImposedDeformationAndFracture::getCopy() const 
{
	StiffnessWithVariableImposedDeformationAndFracture * copy = new StiffnessWithVariableImposedDeformationAndFracture(param, imposed, criterion->getCopy()) ;
	copy->damage = damage ;
	return copy ;
}

void StiffnessWithVariableImposedDeformationAndFracture::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	f = VirtualMachine().ieval(Gradient(p_i,true) * (param * imposed), gp, Jinv,v) ;
}
