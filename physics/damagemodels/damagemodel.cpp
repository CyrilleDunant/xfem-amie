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
		double resphi = 2. - phi ; //goldensearch
		resphi = .5 ;              //bisection
		
		if(fraction < 0)
		{
			upState.resize(state.size(), 0.);
			downState.resize(state.size(), 0.);
			
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
				std::cerr << "elements too large for damage characteristic radius!" << std::endl ;
			fraction = std::min(fraction, 1.) ;
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
		
		Vector originalState = getState() ;
		bool checkpoint = s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ;

		if( checkpoint ) // initiate iteration
		{
			s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint(false);

			getPreviousState() = getState() ;
			if(!fractured())
			{
				converged = false ;
				exploring = true ;
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
					getState() = originalState+damageIncrement*factor ;
					if(getState().max() >=1 || fractured())
					{
						up = factor ;
						factor = (up+down)/2 ;
						Vector startingState = getState() ;
						getState() = originalState+damageIncrement*up ;
						while(!fractured())
						{
							up += 2.*damageDensityTolerance ;
							getState() = originalState+damageIncrement*up ;
						}
						getState() = startingState ;
					}
					else
					{
						down = factor ;
						factor = (up+down)/2 ;
					}
				}
				
				downState = originalState ;
				upState = originalState+damageIncrement*(up) ;
				getState() = originalState + damageIncrement*explorationIncrement ;
				
				//no point in performing an iteration. Also, the element should be broken
				if(std::abs(upState-downState).max() < damageDensityTolerance) 
				{
					std::cerr << ":" << std::flush ;
					getState() = upState ;
					converged = true ;
					exploring = false ;
					return ;
				}
				
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
			if(exploring)
			{
				if(setChange == 0)
				{
					Vector damageIncrement = computeDamageIncrement(s) ;
					if(damageIncrement.max() > POINT_TOLERANCE_2D)
						damageIncrement /= damageIncrement.max() ;
					downState = getState() ;
					getState() = downState + damageIncrement*explorationIncrement ;
					if(getState().max() > 1)
					{
						getState() /= getState().max() ;
						converged = true;
						exploring = false ;
					}
				}
				else //the range for damage has been found
				{
					std::cerr << "!" << std::flush ;
					Vector delta = explorationIncrement*(upState-downState) ;
					upState = downState + delta  ;
					getState() = downState + resphi * (upState - downState) ;
					exploring = false ;

				}
				return ;
			}
			else
			{
				
				if(setChange != 0)
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
			
			if(std::abs(upState-downState).max() < damageDensityTolerance)
			{
				
				if(setChange == 0)
					std::cerr << "0" << std::flush ;
				else
					std::cerr << "1" <<std::flush ;
				
				getState() = upState ;
				exploring = false ;
				converged = true ;
			}
		}
		
	}
	
	DamageModel::DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		wasBroken = false ;
		change = false ;
		isNull = true ; 
		exploring = true ;
		thresholdDamageDensity = 1 ;
		secondaryThresholdDamageDensity = 1 ;
		
		fraction = -1 ;
		converged = true ;
		explorationIncrement = 0.01 ;
		damageDensityTolerance = 1e-6; //1./pow(2., 8) ; // about 1e-4
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
