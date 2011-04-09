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
		
		double phi = (1. + sqrt(5.)) *.5 ;
		double resphi = 2. - phi ;   //goldensearch
// 		resphi = .5 ;              //bisection
// 		resphi = .1 ;                //down bias
		
		if(fraction < 0)
		{
			upState.resize(state.size(), 0.);
			downState.resize(state.size(), 0.);
			
			if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
				fraction = s.getParent()->area() ;
			else
				fraction = s.getParent()->volume();
			
		}
		
		change = false ;
		
		if(wasBroken)
		{
			converged = true ;
			return ;
		}
		
		int setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange(s) ;
		
		if(!s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
		{
			s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint(false);
			converged = true ;
			if(fractured())
				wasBroken = true ;
			return ;
		}
		
		
		bool checkpoint = s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ;

		Vector startPosition = getState() ;
		if(checkpoint) // initiate iteration
		{
			s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint(false);

			getPreviousState() = getState() ;
			if(!fractured())
			{
				converged = false ;
				change = true ;
				Vector damageIncrement = computeDamageIncrement(s) ;
				if(damageIncrement.max() > POINT_TOLERANCE_2D)
					damageIncrement /= damageIncrement.max() ;
				double factor = 0.5 ;
				double down = 0 ;
				double up = 1 ;
				//convergence requires consistent errors
				while(std::abs(up-down) > damageDensityTolerance*.25) 
				{
					getState() = getPreviousState()+damageIncrement*factor ;
					if(getState().max() >=1 || fractured())
					{
						up = factor ;
						factor = (up+down)/2 ;
						Vector startingState = getState() ;
						getState() = getPreviousState()+damageIncrement*up ;
						while(!fractured())
						{
							up += 2.*damageDensityTolerance ;
							getState() = getPreviousState()+damageIncrement*up ;
						}
						getState() = startingState ;
					}
					else
					{
						down = factor ;
						factor = (up+down)/2 ;
					}
				}
				downState = getPreviousState() ;
				upState = getPreviousState()+damageIncrement*up ;
				getState() = downState+(upState-downState)*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				
			}
			else
			{
				wasBroken = true ;
				converged = true ;
			}
		}
		else if(!converged)
		{
			change = true ;

				
			if(setChange != 0 || fractured())
			{
				upState = getState() ;
				getState() = downState + resphi * (upState - downState) ;
			}
			else
			{
				downState = getState() ;
				getState() = downState + resphi * (upState - downState) ;
			}
		}

		if(!converged && std::abs(startPosition-getState()).max()*fraction < damageDensityTolerance )
		{

					
			if(setChange == 0)
				std::cerr << "0" << std::flush ;
			else
				std::cerr << "1" <<std::flush ;
			
// 			getState() = upState ;
			converged = true ;
		}
		
	}
	
	DamageModel::DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		wasBroken = false ;
		change = false ;
		isNull = true ; 
		thresholdDamageDensity = 1 ;
		secondaryThresholdDamageDensity = 1 ;
		
		fraction = -1 ;
		converged = true ;
		iterationcount = 0 ;
		
		// The exploration increment is crucial for finding 
		// the correct distribution of damage: the effect
		// of damage increment on the distribution of
		// fracture criterion scores is non-monotonic.
		explorationIncrement = 4./pow(2., 8) ;
		damageDensityTolerance = 1./pow(2., 8) ; // about 1e-4
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
