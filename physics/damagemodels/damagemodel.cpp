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
				std::cout << "elements too large for damage characteristic radius!" << std::endl ;
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
				lastDirectionUp = true ;
				damageAndSetInPhase = true ;
				converged = false ;
				exploring = true ;
				change = true ;
				Vector damageIncrement = computeDamageIncrement(s) ;
				if(damageIncrement.max() > POINT_TOLERANCE_2D)
					damageIncrement /= damageIncrement.max() ;
				double factor = 0.5 ;
				double down = 0 ;
				double up = 1 ;
				while(std::abs(up-down) > 1e-12)
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
							up +=2e-12 ;
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
				getState() = originalState + (upState-originalState)*explorationIncrement ;
				iterationcount = 1 ;
				
				//no point in performing an iteration. Also, the element should be broken
				if(std::abs(upState-downState).max() < damageDensityTolerance) 
				{
					std::cout << ":" << std::flush ;
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
			if(exploring)
			{
				// There can be two outcomes: either a change occurs or not.
				// When no change occurs, either the element is fractured --
				// and its fracturation does not affect the rest of the mesh
				// or it is not damaged enough.
				//
				// When the change occurs, we can determine the phase of damage
				// increment and set size.
				if(setChange == 0)
				{
					getState() += explorationIncrement*iterationcount*(upState-downState) ; //harmonic search
					if(getState().max() >= upState.max() || fractured())
					{
						for(size_t i = 0 ; i < getState().size() ; i++)
						{
							if(getState()[i] > 1 || getState()[i] >= upState.max()) //this ensures we are properly fractured
								getState()[i] = 1 ;
						}
						exploring = true ;
						converged = true ;
						change = true ;
						std::cout << "¡" << std::flush ;
						return ;
					}
				}
				else //the range for damage has been found
				{
					std::cout << "!" << std::flush ;
					Vector delta = explorationIncrement*(iterationcount)*(upState-downState) ;
					upState = getState() ;
					
					downState = upState - delta ;
					getState() = downState + resphi * (upState - downState) ;
					exploring = false ;
					change = true ;
					if(setChange > 0)
					{
						damageAndSetInPhase = true ;
					}
					else
					{
						damageAndSetInPhase = false ;
					}
					return ;
				}
			}
			
			iterationcount++ ;
			change = true ;
			if(!exploring)
			{
				// finding the relation of the direction and damage increment sign
				// is equivalent to determining whether we are in a parallel or
				// series setup.
				
				bool overdamaged = fractured() || std::abs(getState().max()-1.) < POINT_TOLERANCE_2D;
				bool underdamaged = getState().min() <= POINT_TOLERANCE_2D ;
							
				if(lastDirectionUp && setChange > 0 || !lastDirectionUp && setChange < 0)
					damageAndSetInPhase = true ;
				if(lastDirectionUp && setChange < 0 || !lastDirectionUp && setChange > 0)
					damageAndSetInPhase = false ;
				
				if(!underdamaged && (damageAndSetInPhase && setChange > 0 || !damageAndSetInPhase && setChange < 0 || overdamaged)) //the damage was increased too much
				{
					lastDirectionUp = false ;
					
					upState = getState() ;
					getState() = downState + resphi * (upState - downState) ;
				}
				else if(damageAndSetInPhase && setChange < 0 || !damageAndSetInPhase && setChange > 0)
				{
					lastDirectionUp = true ;
					
					downState = getState() ;
					getState() = downState + resphi * (upState - downState) ;
				}
				else if(setChange == 0 && (lastDirectionUp != damageAndSetInPhase))
				{
					lastDirectionUp = false ;
					
					upState = getState() ;
					getState() = downState + resphi * (upState - downState) ;
				}
				else // setChange == 0 && (lastDirectionUp == damageAndSetInPhase)
				{
					lastDirectionUp = true ;
					
					downState = getState() ;
					getState() = downState + resphi * (upState - downState) ;
				}
			}
			if(std::abs(upState-downState).max() < damageDensityTolerance)
			{
				getState() = upState ;
				std::cout << "*" << std::flush ;
				exploring = true ;
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
		damageDensityTolerance = 1./pow(2., 20) ; // about 1e-6
		fraction = -1 ;
		damageAndSetInPhase = false ;
		lastDirectionUp = true ;
		iterationcount = 0 ;
		converged = true ;
		explorationIncrement = 0.1 ;
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
