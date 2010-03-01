//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mohrcoulomb.h"

using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit) : LinearForm(rig, false, true, rig.numRows()/3+1),dfunc(rig.numRows()-1)/*dfunc(rig.numRows()-1)*/, eps(0.02)
{
	criterion = crit ;
	crit->setNeighbourhoodRadius(eps) ;
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previousDamage = 0 ;
	intermediateDamage = 0 ;
	count = 0 ;
	previousPreviousDamage = 0 ;
	damage = 0 ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
} ;

FractureCriterion * StiffnessAndFracture::getFractureCriterion() const
{
	return criterion ;
}

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{

	VirtualMachine vm ;
	return vm.ieval(Gradient(p_i) * dfunc.apply(param) * Gradient(p_j, true), e,v) ;
}

void StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc.apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}


void StiffnessAndFracture::stepBack()
{
	damage = previousDamage ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessAndFracture::step(double timestep, ElementState & currentState) 
{
	change = false ;
	currentState.getParent()->behaviourUpdated = false ;
	if(!frac && criterion->met(currentState) )
	{
		
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal *=.9 ;
		dfunc.step(currentState) ;
		if(timestep > 0)
		{
			previousPreviousDamage = previousDamage ;
			previousDamage = damage ;
		}
		Vector state = dfunc.damageState() ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal /=1.-damage ;
		damage = 0 ;
		for(size_t i = 0 ; i < state.size() ; i++)
			damage += state[i] ;

// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal *=1.-damage ;
// 		double deltadamage = damage-previousDamage ;
// 		
// 		if(damage < .1)
// 			dynamic_cast<MohrCoulomb *>(criterion)->upVal += (previousDamage/damage)*.1*23. ;
// 		else if (damage > .3)
// 			dynamic_cast<MohrCoulomb *>(criterion)->upVal *=0.9 ;

		change = true ;
		currentState.getParent()->behaviourUpdated = true ;
		frac = dfunc.fractured() ;
	}

}

bool StiffnessAndFracture::changed() const
{
	return change ;
} 

bool StiffnessAndFracture::fractured() const
{
	return frac;
}

Form * StiffnessAndFracture::getCopy() const 
{
	StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy()) ;
	copy->damage = damage ;
	return copy ;
}

void StiffnessAndFracture::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

Matrix StiffnessAndFracture::getTensor(const Point & p) const
{
	return dfunc.apply(param) ;
}

