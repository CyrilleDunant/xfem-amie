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
	isNull = false ;
	state = 0 ;
	tensionDamage = 0 ;
	compressionDamage = 0 ;
	inCompression = false ;
	inTension = false ;
}

const Vector & FractionLinearDamage::damageState() const
{
	return state ;
}


Vector & FractionLinearDamage::damageState() 
{
	return state ;
}
void FractionLinearDamage::step(ElementState & s)
{
	previousstate.resize(state.size());
	previousstate = state ;
	inCompression = false ;
	inTension = false ;
	if(fraction < 0)
	{
		double volume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			volume = sqrt(s.getParent()->area()) ;
		else
			volume = pow(s.getParent()->volume(), 2./3.) ;
		
		double charVolume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			charVolume = sqrt(M_PI*characteristicRadius*characteristicRadius) ;
		else
			charVolume = pow(4./3.*M_PI*characteristicRadius*characteristicRadius*characteristicRadius, 2./3.) ;
		fraction = volume/charVolume ;
		if(fraction > 1)
			std::cout << "elements too large for damage characteristic radius!" << std::endl ;
		fraction = std::min(fraction, 1.) ;
	}
	

	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		
		compressionDamage += damageDensityIncrement*fraction ; 
		compressionDamage = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, compressionDamage) ;
		compressionDamage = std::min(.99999, compressionDamage) ;
		compressionDamage = std::max(0., compressionDamage) ;
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
	{
		inTension = true ;

		tensionDamage += damageDensityIncrement*fraction ; 
		tensionDamage = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamage) ;
		tensionDamage = std::min(.99999, tensionDamage) ;
		tensionDamage = std::max(0., tensionDamage) ;
	}
	state[0] = compressionDamage ;
	state[1] = tensionDamage ;
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
		return m*(1.-tensionDamage)*(1.-phi)+remnant*phi ;
	}
	if(inCompression)
	{
		return m*(1.-compressionDamage)*(1.-phi)+remnant*phi ;
	}
	return ret*(1.-std::max(tensionDamage, compressionDamage))*(1.-phi)+remnant*phi ;
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
	
	return (tensionDamage >= secondaryThresholdDamageDensity) || (compressionDamage >= thresholdDamageDensity)  ;
}

FractionLinearDamage::~FractionLinearDamage()
{
}


}
