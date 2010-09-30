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
#include "../features/boundarycondition.h"

using namespace Mu ;

StiffnessWithVariableImposedDeformationAndFracture::StiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * c) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef),criterion(c)
{
	dfunc = new IsotropicLinearDamage(rig.numRows()-1, .01) ;
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

StiffnessWithVariableImposedDeformationAndFracture::~StiffnessWithVariableImposedDeformationAndFracture() 
{ 
	delete dfunc ;
} ;

void StiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

void StiffnessWithVariableImposedDeformationAndFracture::step(double timestep, ElementState & currentState)
{
	if(!frac && criterion->met(currentState) )
	{
		dfunc->step(currentState) ;
		previousDamage = damage ;
		
		Vector state = dfunc->damageState() ;
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

Vector StiffnessWithVariableImposedDeformationAndFracture::getPreviousDamage()
{
	Vector previous(1) ;
	previous[0] = previousDamage ;
	return previous ;
}


void StiffnessWithVariableImposedDeformationAndFracture::artificialDamageStep(double d)
{
	previousDamage = damage ;
	dfunc->artificialDamageStep(d) ;
	Vector state = dfunc->damageState() ;
	damage = 0 ;
	for(size_t i = 0 ; i < state.size() ; i++)
		damage += state[i] ;
	if(damage > .9)
		frac = true ;

	if(frac)
		imposed = 0 ;

	change = true ;
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

std::vector<BoundaryCondition * > StiffnessWithVariableImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
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
