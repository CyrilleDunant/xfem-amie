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

#include "weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.h"

using namespace Mu ;

WeibullStiffnessWithVariableImposedDeformationAndFracture::WeibullStiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * c) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef), dfunc(rig.numRows()-1, .01),criterion(c)
{
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previousDamage = 0 ;
	damage = 0 ;
	variability = .2 ;	
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

WeibullStiffnessWithVariableImposedDeformationAndFracture::~WeibullStiffnessWithVariableImposedDeformationAndFracture() { } ;

Matrix WeibullStiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

void WeibullStiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool WeibullStiffnessWithVariableImposedDeformationAndFracture::hasInducedForces() const 
{
	return true ; 
} 

Vector WeibullStiffnessWithVariableImposedDeformationAndFracture::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

void WeibullStiffnessWithVariableImposedDeformationAndFracture::step(double timestep, ElementState & currentState)
{
	change = false ;

	if(!frac && criterion->met(currentState) )
	{
		dfunc.step(currentState) ;
		previousDamage = damage ;
		
		Vector state = dfunc.damageState() ;
		damage = 0 ;
		for(size_t i = 0 ; i < state.size() ; i++)
			damage += state[i] ;
		change = true ;//std::abs(damage-previousDamage) > 1e-12 ;
		if(damage > .9)
		{
			frac = true ;
			damage = .9999 ;
// 			param[0][1] = 0 ;param[0][1] = 0 ;
// 			param[2][2] *= 0.0001 ;
// 			this->type = VOID_BEHAVIOUR ;
		}
	}

	if(!change && !frac)
	{
		double randomVar = (double)rand()/(double)RAND_MAX ;
		imposed[0] += timestep*randomVar ;
		imposed[1] += timestep*randomVar ;
	}
	
	if(frac)
	{
		imposed = 0 ;
	}
}

void WeibullStiffnessWithVariableImposedDeformationAndFracture::artificialDamageStep(double d)
{
	dfunc.artificialDamageStep(d) ;
	previousDamage = damage ;
		
	Vector state = dfunc.damageState() ;
	damage = 0 ;
	for(size_t i = 0 ; i < state.size() ; i++)
		damage += state[i] ;
	change = true ;//std::abs(damage-previousDamage) > 1e-12 ;
	if(damage > .9)
		frac = true ;

	if(frac)
		imposed = 0 ;

}

Vector WeibullStiffnessWithVariableImposedDeformationAndFracture::getPreviousDamage()
{
	Vector previous(1) ;
	previous[0] = previousDamage ;
	return previous ;
}



bool WeibullStiffnessWithVariableImposedDeformationAndFracture::changed() const
{
	return change ;
} 

bool WeibullStiffnessWithVariableImposedDeformationAndFracture::fractured() const
{
	return frac;
}

Form * WeibullStiffnessWithVariableImposedDeformationAndFracture::getCopy() const 
{
	double randomVar = (double)rand()/(double)RAND_MAX ;
	randomVar = 1.*pow(-log(randomVar),1./2.) ;
	Matrix newTensor(param*(1.-variability)+param*randomVar*variability) ;
	
	StiffnessWithVariableImposedDeformationAndFracture * copy = new StiffnessWithVariableImposedDeformationAndFracture(newTensor, imposed, criterion->getCopy()) ;
	copy->damage = damage ;

	return copy ;
}

void WeibullStiffnessWithVariableImposedDeformationAndFracture::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	f = VirtualMachine().ieval(Gradient(p_i) * ( (param * damage) * imposed), gp, Jinv,v) ;
}

