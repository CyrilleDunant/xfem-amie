//
// C++ Interface: damagemodel
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "damagemodel.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu
{

/** \brief Damage model interface */


	void DamageModel::step(ElementState & s )
	{
		if(fraction < 0)
		{
			damageIncrement.resize(state.size(), 0.);
			previousDamageIncrement.resize(state.size(), 0.);
			
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
		
		int testRank = s.getParent()->getBehaviour()->getFractureCriterion()->getRank(200,s) ;
		if(s.getDeltaTime() < POINT_TOLERANCE && lastRank == 1) // we are within the two iteration (bissection and damage)
		{
			change = true ;
			damageIncrement = computeDamageIncrement(s) ;
			getState() -= previousDamageIncrement ;
			if(testRank != 1)
			{
				
				upFactor = currentFactor ;
				currentFactor = downFactor ;
				
				getState() += damageIncrement*currentFactor ; 
				previousDamageIncrement = damageIncrement*currentFactor ;
				
				if(std::abs(currentFactor-upFactor) < damageDensityTolerance || (wasBroken && fractured())) // we have converged
				{
					if(fractured())
						wasBroken = true ;
					change = false ;
					lastRank = 2 ;
					upFactor = 1 ;
					downFactor = 0 ;
					currentFactor = 0.5 ;
				}
			}
			else
			{
				if(currentFactor < 1)
					downFactor = currentFactor ;
					
				currentFactor = (upFactor + downFactor)*.5 ;
				getState() += damageIncrement*currentFactor ; 
				previousDamageIncrement = damageIncrement*currentFactor ;
				
				if(std::abs(currentFactor-upFactor) < damageDensityTolerance|| wasBroken && fractured()) // we have converged
				{
					if(fractured())
						wasBroken = true ;
					change = false ;
					lastRank = 2 ;
					upFactor = 1 ;
					downFactor = 0 ;
					currentFactor = 1 ;
				}
			}
		}
		else if (s.getDeltaTime() < POINT_TOLERANCE && lastRank != 1) // we can enter the damage iteration
		{
			if( testRank == 1)
			{
				change = true ;
				damageIncrement = 0 ;
				getState() += damageIncrement*.5 ;
				previousDamageIncrement = 0 ;
				lastRank = 1 ;
				upFactor = 1 ;
				downFactor = 0 ;
				currentFactor = 0.5 ;
			}
			else
			{
				change = false ;
				lastRank = 2 ;
			}
		}
		else if(s.getDeltaTime() > POINT_TOLERANCE)
		{
			getPreviousState() = getState() ;
			if(fractured())
				wasBroken = true ;
			
			if(testRank == 1 && !fractured())
			{
				change = true ;
				damageIncrement = computeDamageIncrement(s) ;
				getState() += damageIncrement*.5 ;
				previousDamageIncrement = damageIncrement*.5 ;
				lastRank = 1 ;
				upFactor = 1 ;
				downFactor = 0 ;
				currentFactor = 0.5 ;
			}
			else
			{
				change = false ;
				previousDamageIncrement = 0 ;
				lastRank = 2 ;
				upFactor = 1 ;
				downFactor = 0 ;
				currentFactor = 0.5 ;
			}
		}
		
		for(size_t i = 0 ; i < getState().size() ; i++)
			getState()[i] = std::min(getState()[i], 1.) ;
	}
	
	DamageModel::DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		wasBroken = false ;
		change = false ;
		isNull = true ; 
		thresholdDamageDensity = .999999999999 ;
		secondaryThresholdDamageDensity = .99999999999999 ;
		damageDensityTolerance = .00001 ;
		fraction = -1 ;
		
		upFactor = 1;
		downFactor = 0;
		currentFactor = 0.5;
		
		lastRank = 2 ;
	} ;
	
	double DamageModel::getThresholdDamageDensity() const
	{
		return thresholdDamageDensity ;
	}
	
	double DamageModel::getSecondaryThresholdDamageDensity() const
	{
		return secondaryThresholdDamageDensity ;
	}
	
	double DamageModel::getCharacteristicRadius() const
	{
		return characteristicRadius ;
	}
	
	void DamageModel::setThresholdDamageDensity(double d)
	{
		thresholdDamageDensity = d ;
	}
	
	void DamageModel::setSecondaryThresholdDamageDensity(double d)
	{
		secondaryThresholdDamageDensity = d ;
	}
	
	void DamageModel::setMaterialCharacteristicRadius(double d)
	{
		characteristicRadius = d ;
	}
	
	void DamageModel::setDamageDensityTolerance(double d)
	{
		damageDensityTolerance = d ;
	}
	
	bool DamageModel::changed() const
	{
		return change ;
	}
} ;	
