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

namespace Mu {

LinearDamage::LinearDamage()
{
	getState(true).resize(2, 0.) ;
	getPreviousState().resize(2, 0.);
	isNull = false ;
	state = 0 ;

	inCompression = false ;
	inTension = false ;
}

std::pair< Vector, Vector > LinearDamage::computeDamageIncrement(ElementState &s)
{
	inCompression = false ;
	inTension = false ;
	Vector ret(0., 2) ; 
	double compressionDamage = 0 ;
	double tensionDamage = 0 ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		compressionDamage = 1 ; 
		tensionDamage = 1 ; 
	}
	else
	{
		inCompression = false ;
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
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

void LinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < getState().size() -1 ; i++)
		getState(true)[i] = std::min(getState()[i]+d,0.9999999) ;
}

Matrix LinearDamage::apply(const Matrix & m) const
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
	
	return m ;
}

Matrix LinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return ret*0.;
	//this is a silly way of distinguishing between 2D and 3D
	for(size_t i = 0 ; i < (m.numRows()+1)/2 ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(getPreviousState()[0]) ;
		}
	}

	for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(getPreviousState()[i]) ;
		}
	}

	return ret ;
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
