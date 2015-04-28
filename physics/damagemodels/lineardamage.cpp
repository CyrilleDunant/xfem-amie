//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

LinearDamage::LinearDamage()
{
	getState(true).resize(2, 0.) ;
	isNull = false ;
	state = 0 ;

	inCompression = false ;
	inTension = false ;
}

std::pair< Vector, Vector > LinearDamage::computeDamageIncrement(ElementState &s)
{
	inCompression = false ;
	inTension = false ;
	Vector ret = state ; 
	double compressionDamage = 0 ;
	double tensionDamage = 0 ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0))
	{
		inCompression = true ;
		compressionDamage = 1 ; 
	}
	else
	{
		inCompression = false ;
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0))
	{
		inTension = true ;
		
		tensionDamage = 1 ; 
	}
	else
	{
		inTension = false ;
	}
	
	ret[0] = compressionDamage ;
	ret[1] = tensionDamage ;
	
	return std::make_pair(state, ret) ;
}

void LinearDamage::computeDelta(ElementState & s)
{
	Vector ret = state ; 
	double compressionDamage = 0 ;
	double tensionDamage = 0 ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0))
	{
		compressionDamage = 1 ; 
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0))
	{
		tensionDamage = 1 ; 
	}
	
	ret[0] = compressionDamage ;
	ret[1] = tensionDamage ;
	
	delta = (ret-state).max() ;
	
}

Matrix LinearDamage::apply(const Matrix & m, const Point & p , const IntegrableEntity * e , int g) const
{
// 	std::cout << damageDensityIncrement<< "   "<< tensionDamage << "  " << compressionDamage << std::endl ;
	
	if(fractured())
		return m*0.;

	if(inTension && !inCompression)
	{
		return m*(1.-getState()[1]) ;
	}
	else if(inCompression && !inTension)
	{
		return m*(1.-getState()[0]) ;
	}
	
	return m*(1.-getState().max()) ;
}



bool LinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
		if(inCompression)
			return getState()[0] >= thresholdDamageDensity ;
		if(inTension)
			return getState()[1] >= secondaryThresholdDamageDensity ;
		
	return  getState()[1] >= secondaryThresholdDamageDensity || getState()[0] >= thresholdDamageDensity ;
}

LinearDamage::~LinearDamage()
{
}


}
