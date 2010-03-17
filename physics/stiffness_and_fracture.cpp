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
	previouschange = false ;
	previousDamage.resize(dfunc.damageState().size()) ; previousDamage =0 ;
	intermediateDamage.resize(dfunc.damageState().size()) ;intermediateDamage = 0 ;
	count = 0 ;
	previousPreviousDamage.resize(dfunc.damageState().size()) ;previousPreviousDamage = 0 ;
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
	change = previouschange ;
	damage.resize(previousDamage.size()) ;
	damage = previousDamage ;
	dfunc.damageState() = damage ;
	frac = dfunc.fractured() ;
	previousDamage.resize(previousPreviousDamage.size()) ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessAndFracture::step(double timestep, ElementState & currentState) 
{
	previouschange = change ;
	change = false ;
	currentState.getParent()->behaviourUpdated = false ;
	if(!frac && criterion->met(currentState) )
	{
		
		dfunc.step(currentState) ;
// 		if(timestep > 0)
// 		{
			previousPreviousDamage.resize(previousDamage.size()) ;
			previousPreviousDamage = previousDamage ;
			previousDamage.resize(damage.size()) ;
			previousDamage = damage ;
// 		}
		Vector d = dfunc.damageState() ;
		damage.resize(d.size()) ;
		damage = d ;

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
	return dfunc.fractured() ;
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

