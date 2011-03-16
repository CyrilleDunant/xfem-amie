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
#include "fractiondamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

FractionLinearDamage::FractionLinearDamage(int numDof, double characteristicRadius, Matrix remnant, double phi) : DamageModel(characteristicRadius), remnant(remnant), phi(phi)
{
	state.resize(2, 0.) ;
	previousstate.resize(2, 0.) ;
	isNull = false ;
	state = 0 ;

	inCompression = false ;
	inTension = false ;
}

Vector FractionLinearDamage::computeDamageIncrement(ElementState & s)
{
	Vector ret(state.size()) ; ret = 0;

	inCompression = false ;
	inTension = false ;

	double compressionDamage = 0 ;
	double tensionDamage = 0 ;

	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		
		compressionDamage = 1.-state[0] ; 
// 		compressionDamage = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, compressionDamage) ;
// 		compressionDamage = std::min(.99999, compressionDamage) ;
// 		compressionDamage = std::max(0., compressionDamage) ;
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
	{
		inTension = true ;

		tensionDamage = 1.-state[1] ; 
// 		tensionDamage = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamage) ;
// 		tensionDamage = std::min(.99999, tensionDamage) ;
// 		tensionDamage = std::max(0., tensionDamage) ;
	}
	ret[0] = compressionDamage ;
	ret[1] = tensionDamage ;
	return ret ;
// 	std::cout << state.sum() << std::flush ;
}

void FractionLinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < state.size() -1 ; i++)
		state[i] = std::min(state[i]+d,0.9999) ;
}

Matrix FractionLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m*(1.-phi)+remnant*phi) ;
	
	if(fractured())
		return remnant*phi;

	if(inTension)
	{
		return m*(1.-state[0])*(1.-phi)+remnant*phi ;
	}
	if(inCompression)
	{
		return m*(1.-state[1])*(1.-phi)+remnant*phi ;
	}
	return ret*(1.-std::max(state[0], state[1]))*(1.-phi)+remnant*phi ;
}

Matrix FractionLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return remnant*phi;
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

bool FractionLinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
	return (state[0] >= secondaryThresholdDamageDensity/fraction) || (state[1] >= thresholdDamageDensity/fraction)  ;
}

FractionLinearDamage::~FractionLinearDamage()
{
}


}
