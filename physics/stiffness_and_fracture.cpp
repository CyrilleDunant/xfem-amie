//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"

using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit) : LinearForm(rig, false, true, rig.numRows()/3+1)
{
	dfunc = new /*NonLocal*/IsotropicLinearDamage() ;
	criterion = crit ;
	init = param[0][0] ;
	previouschange = false ;
	previousDamage.resize(dfunc->getState().size()) ; 
	previousDamage =0 ;
	
	intermediateDamage.resize(dfunc->getState().size()) ;
	intermediateDamage = 0 ;
	
	previousPreviousDamage.resize(dfunc->getState().size()) ;
	previousPreviousDamage = 0 ;
	
	count = 0 ;
	damage = 0 ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

void StiffnessAndFracture::setNeighbourhoodRadius(double d)
{
	criterion->setNeighbourhoodRadius(d);
}

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
	delete dfunc ;
} ;

FractureCriterion * StiffnessAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessAndFracture::getDamageModel() const
{
	return dfunc ;
}

void StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}


void StiffnessAndFracture::stepBack()
{

	damage.resize(previousDamage.size()) ;
	damage = previousDamage ;
	dfunc->getState(true) = damage ;

	previousDamage.resize(previousPreviousDamage.size()) ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessAndFracture::step(double timestep, ElementState & currentState) 
{
	previouschange = dfunc->changed() ;
	dfunc->step(currentState) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
	
	if(dfunc->changed())
	{
		
		previousPreviousDamage.resize(previousDamage.size()) ;
		previousPreviousDamage = previousDamage ;
		previousDamage.resize(damage.size()) ;
		previousDamage = damage ;
		
		Vector d = dfunc->getState() ;
		damage.resize(d.size()) ;
		damage = d ;
	}
}

void StiffnessAndFracture::artificialDamageStep(double d)
{
	previouschange = dfunc->changed() ;

	dfunc->artificialDamageStep(d) ;
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d_ = dfunc->getState() ;
	damage.resize(d_.size()) ;
	damage = d ;
}

void StiffnessAndFracture::artificialPreviousDamage(Vector previous, Vector previousprevious)
{
	previousDamage.resize(damage.size()) ;
	if(previous.size() < previousDamage.size())
	{
		for(size_t i = 0 ; i < previous.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
		for(size_t j = previous.size() ; j < previousDamage.size() ; j++)
			previousDamage[j] = std::min(damage[j],previous[previous.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousDamage.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
	}
	previousPreviousDamage.resize(damage.size()) ;
	if(previousprevious.size() < previousPreviousDamage.size())
	{
		for(size_t i = 0 ; i < previousprevious.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
		for(size_t j = previous.size() ; j < previousPreviousDamage.size() ; j++)
			previousPreviousDamage[j] = std::min(previousDamage[j],previousprevious[previousprevious.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousPreviousDamage.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
	}
}



bool StiffnessAndFracture::changed() const
{
	return dfunc->changed() ;
} 

bool StiffnessAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

Form * StiffnessAndFracture::getCopy() const 
{
	StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy()) ;
	copy->damage = damage ;
	copy->dfunc->getState(true).resize(dfunc->getState().size());
	copy->dfunc->getState(true) = dfunc->getState() ;
	copy->dfunc->getPreviousState().resize(dfunc->getPreviousState().size());
	copy->dfunc->getPreviousState() = dfunc->getPreviousState() ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->criterion->setNeighbourhoodRadius(criterion->getNeighbourhoodRadius()) ;
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	return copy ;
}

Matrix StiffnessAndFracture::getTensor(const Point & p) const
{
	return dfunc->apply(param) ;
}

Material StiffnessAndFracture::toMaterial()
{
	Material mat(getTensor(Point(0,0))) ;
	mat.setProperties(criterion->toMaterial()) ;
	return mat ;
}



