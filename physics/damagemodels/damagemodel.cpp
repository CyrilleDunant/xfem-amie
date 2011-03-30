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
#include "../../mesher/delaunay.h"

namespace Mu
{

/** \brief Damage model interface */


	void DamageModel::step(ElementState & s )
	{
		
		double fractionIncrease = 0.5;
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
		
		change = false ;
		
		if(wasBroken)
			return ;
		
		std::pair<bool, bool> inSetAndSetChanged = s.getParent()->getBehaviour()->getFractureCriterion()->inSetAndSetChanged(200,s) ;
		if(!inSetAndSetChanged.first)
			return ;
		
		if(s.getParent()->getBehaviour()->getFractureCriterion()->met(s) 
			&& !wasBroken 
			&& inSetAndSetChanged.first
			&& !inSetAndSetChanged.second
			&& converged
			&& s.getDeltaTime() < POINT_TOLERANCE_2D
		) //we are at a checkpoint during an iteration
		{
// 			std::cout << "checkpoint! " << dynamic_cast<DelaunayTriangle *>(s.getParent())->index << "  "<< getState().max() <<std::endl; 
			if(fractured())
			{
				wasBroken ;
				return ;
			}
			
			converged = false ;
			damageIncrement = computeDamageIncrement(s) ;
			damageIncrement /= damageIncrement.max() ;
// 			damageIncrement = 1 ;
			
			getState() += damageIncrement*fractionIncrease ;
			previousDamageIncrement = damageIncrement*fractionIncrease ;
			upFactor = 1 ;
			downFactor = 0 ;
			currentFactor = fractionIncrease ;
			change = true ;
			for(size_t i = 0 ; i < getState().size() ; i++)
				getState()[i] = std::max(0., std::min(getState()[i], 1.)) ;
			return ;
		}
		
		if(s.getDeltaTime() > POINT_TOLERANCE_2D ) // initiate iteration
		{
			getPreviousState() = getState() ;
			if(fractured())
				wasBroken = true ;
			
			if(inSetAndSetChanged.first && !wasBroken)
			{
				converged = false ;
				damageIncrement = computeDamageIncrement(s) ;
				damageIncrement /= damageIncrement.max() ;
// 				damageIncrement = 1 ;
				
				getState() += damageIncrement*fractionIncrease ;
				previousDamageIncrement = damageIncrement*fractionIncrease ;
				upFactor = 1 ;
				downFactor = 0 ;
				currentFactor = fractionIncrease ;
			}
			else
			{
				converged = true ;
				previousDamageIncrement = 0 ;
				upFactor = 1 ;
				downFactor = 0 ;
				currentFactor = fractionIncrease ;
			}
		}
		else if(s.getDeltaTime() < POINT_TOLERANCE_2D 
			&& inSetAndSetChanged.first && !converged)
		{
			damageIncrement = computeDamageIncrement(s) ;
			damageIncrement /= damageIncrement.max() ;
// 			damageIncrement = 1 ;
			bool frac  = fractured() ;
			getState() -= previousDamageIncrement ;
				
			if(inSetAndSetChanged.second || frac) //the damage was increased too much
			{
				upFactor = currentFactor ;
				currentFactor = (upFactor + downFactor)*.5 ;
				
				getState() += damageIncrement*currentFactor ; 
				previousDamageIncrement = damageIncrement*currentFactor ;
				
				if(std::abs(downFactor-upFactor) < damageDensityTolerance ) // we have converged
				{
					converged = true ;
					upFactor = 1 ;
					downFactor = 0 ;
					currentFactor = fractionIncrease ;
					previousDamageIncrement = 0 ;
				}
			}
			else 
			{
				downFactor = currentFactor ;
				currentFactor = (upFactor + downFactor)*.5 ;
				getState() += damageIncrement*currentFactor ; 
				previousDamageIncrement = damageIncrement*currentFactor ;
				
				if(std::abs(downFactor-upFactor) < damageDensityTolerance ) // we have converged
				{
					converged = true ;
					upFactor = 1 ;
					downFactor = 0 ;
					currentFactor = fractionIncrease ;
					previousDamageIncrement = 0 ;
				}
			}

		}

		change = !converged ;
		
		for(size_t i = 0 ; i < getState().size() ; i++)
			getState()[i] = std::max(0., std::min(getState()[i], 1.)) ;
	}
	
	DamageModel::DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		wasBroken = false ;
		change = false ;
		isNull = true ; 
		thresholdDamageDensity = .999999999999 ;
		secondaryThresholdDamageDensity = .99999999999999 ;
		damageDensityTolerance = 1e-7 ;
		fraction = -1 ;
		
		upFactor = 1;
		downFactor = 0;
		currentFactor = 0.5;
		
		converged = true ;
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
