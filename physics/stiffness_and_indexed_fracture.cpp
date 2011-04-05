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

#include "stiffness_and_indexed_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"

using namespace Mu ;


StiffnessAndIndexedFracture::StiffnessAndIndexedFracture(const Matrix & rig, FractureCriterion * crit, double eps) : LinearForm(rig, false, true, rig.numRows()/3+1), /*dfunc(rig.numRows()-1)*/damagedAtStep(true), eps(eps)
{
	dfunc = new IndexedLinearDamage(rig.numRows()-1,1., crit) ; 
	criterion = crit ;
	crit->setNeighbourhoodRadius(eps) ;
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previouschange = false ;
	previousDamage.resize(dfunc->getState().size()) ; previousDamage =0 ;
	intermediateDamage.resize(dfunc->getState().size()) ;intermediateDamage = 0 ;
	count = 0 ;
	previousPreviousDamage.resize(dfunc->getState().size()) ;previousPreviousDamage = 0 ;
	damage = 0 ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

void StiffnessAndIndexedFracture::setNeighbourhoodRadius(double d)
{
	criterion->setNeighbourhoodRadius(d);
	eps = d ;
}

StiffnessAndIndexedFracture::~StiffnessAndIndexedFracture() 
{ 
	delete criterion ;
	delete dfunc ;
} ;

FractureCriterion * StiffnessAndIndexedFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessAndIndexedFracture::getDamageModel() const
{
	return dfunc ;
}

void StiffnessAndIndexedFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}


void StiffnessAndIndexedFracture::stepBack()
{
// 	if(change)
// 	{
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal /= .95 ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->downVal /= .95 ;
// 	}
	change = previouschange ;
	damage.resize(previousDamage.size()) ;
	damage = previousDamage ;
	dfunc->getState() = damage ;
	frac = dfunc->fractured() ;

	previousDamage.resize(previousPreviousDamage.size()) ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessAndIndexedFracture::step(double timestep, ElementState & currentState) 
{
	previouschange = change ;
	change = false ;
	currentState.getParent()->behaviourUpdated = false ;
	if(currentState.getDeltaTime() > POINT_TOLERANCE_2D)
		damagedAtStep = false ;
		
	if( currentState.getDeltaTime() <= POINT_TOLERANCE_2D  /*criterion->getScoreAtState() > 0*/|| damagedAtStep)
	{
		damagedAtStep = true ;
		double pdamage = dfunc->getState()[0] ;
		dfunc->step(currentState) ;
		change = dfunc->changed() ;
		currentState.getParent()->behaviourUpdated = change ;
		frac = dfunc->fractured() ;
	}
	if(currentState.getDeltaTime() > POINT_TOLERANCE_2D)
		change = true ;
	
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d = dfunc->getState() ;
	damage.resize(d.size()) ;
	damage = d ;
}

void StiffnessAndIndexedFracture::artificialDamageStep(double d)
{
	previouschange = change ;
	change = false ;

	dfunc->artificialDamageStep(d) ;
	change = true ;
	frac = dfunc->fractured() ;
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d_ = dfunc->getState() ;
	damage.resize(d_.size()) ;
	damage = d ;
}

void StiffnessAndIndexedFracture::artificialPreviousDamage(Vector previous, Vector previousprevious)
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



bool StiffnessAndIndexedFracture::changed() const
{
	return change ;
} 

bool StiffnessAndIndexedFracture::fractured() const
{
	return dfunc->fractured() ;
}

Form * StiffnessAndIndexedFracture::getCopy() const 
{
	StiffnessAndIndexedFracture * copy = new StiffnessAndIndexedFracture(param, criterion->getCopy(), criterion->getMaterialCharacteristicRadius()) ;
	copy->damage = damage ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->criterion->setNeighbourhoodRadius(criterion->getNeighbourhoodRadius()) ;
	copy->dfunc->setMaterialCharacteristicRadius(dfunc->getCharacteristicRadius());
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	return copy ;
}

Matrix StiffnessAndIndexedFracture::getTensor(const Point & p) const
{
	return dfunc->apply(param) ;
}


Material StiffnessAndIndexedFracture::toMaterial()
{
	Material mat(getTensor(Point(0,0))) ;
	mat.setProperties(criterion->toMaterial()) ;
	return mat ;
}

