//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

LinearDamage::LinearDamage(int numDof, double characteristicRadius) : DamageModel(characteristicRadius)
{
	state.resize(2, 0.) ;
	previousstate.resize(2, 0.);
	isNull = false ;
	state = 0 ;
	tensionDamage = 0 ;
	compressionDamage = 0 ;
	inCompression = false ;
	inTension = false ;
}

const Vector & LinearDamage::damageState() const
{
	return state ;
}


Vector & LinearDamage::damageState() 
{
	return state ;
}

Vector LinearDamage::computeDamageIncrement(ElementState & s)
{
	inCompression = false ;
	inTension = false ;
	Vector ret(2) ; ret = 0 ;
	compressionDamage = 0 ;
	tensionDamage = 0 ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		
		compressionDamage = 0.5 ; 
// 		compressionDamage = std::min(thresholdDamageDensity+POINT_TOLERANCE, compressionDamage) ;
// 		compressionDamage = std::min(.9999999, compressionDamage) ;
// 		compressionDamage = std::max(0., compressionDamage) ;
		
		tensionDamage = 0.5 ; 
// 		tensionDamage = std::min(secondaryThresholdDamageDensity+POINT_TOLERANCE, tensionDamage) ;
// 		tensionDamage = std::min(.9999999, tensionDamage) ;
// 		tensionDamage = std::max(0., tensionDamage) ;
		
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
	{
		inTension = true ;

		tensionDamage = 0.5 ; 
// 		tensionDamage = std::min(secondaryThresholdDamageDensity+POINT_TOLERANCE, tensionDamage) ;
// 		tensionDamage = std::min(.9999999, tensionDamage) ;
// 		tensionDamage = std::max(0., tensionDamage) ;
	}
	ret[0] = compressionDamage ;
	ret[1] = tensionDamage ;
// 	std::cout << state.sum() << std::flush ;
}

void LinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < state.size() -1 ; i++)
		state[i] = std::min(state[i]+d,0.9999999) ;
}

Matrix LinearDamage::apply(const Matrix & m) const
{
// 	std::cout << damageDensityIncrement<< "   "<< tensionDamage << "  " << compressionDamage << std::endl ;
	
	if(fractured())
		return m*0.;

	if(inTension && !inCompression)
	{
		return m*(1.-state[1]) ;
	}
	else if(inCompression && !inTension)
	{
		return m*(1.-state[0]) ;
	}
	
	return m*(1.-std::max(state[1], state[0])) ;
}

Matrix LinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return ret*0.0001	;
	//this is a silly way of distinguishing between 2D and 3D
	for(size_t i = 0 ; i < (m.numRows()+1)/2 ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(previousstate[0]) ;
		}
	}

	for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(previousstate[i]) ;
		}
	}

	return ret ;
}

bool LinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
	return state[1] >= secondaryThresholdDamageDensity || state[0] >= thresholdDamageDensity ;
}

LinearDamage::~LinearDamage()
{
}


}
