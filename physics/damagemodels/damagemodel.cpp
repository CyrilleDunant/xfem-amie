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
		
		double setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange(s) ;
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

				states.clear() ;
				states.push_back(PointState(s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange, 0.)) ;
				downState = getPreviousState();
				upState = getPreviousState()+damageIncrement*up ;
				getState() =upState ;
				trialRatio = 1 ;
				if((upState-downState).min() < 0)
				{
					getState() = upState ;
					converged = true ;
					wasBroken = true ;
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
			
			for(auto i = states.begin() ; i != states.end() ; i++)
			{
				if(trialRatio > i->fraction)
				{
					
					states.insert(i, PointState(s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange, trialRatio)) ;
					break;
				}
			}

			//find the most likely midPoint
			RangeState bestRange(states[1], states[0]) ;
			PointState bestState = bestRange.extrapolate() ;
			bool maxDamage = bestState.isMet ;
			for (size_t i = 1 ; i < states.size()-1 ;i++ )
			{
				RangeState trialRange(states[i+1], states[i]) ;
				PointState trialState = trialRange.extrapolate() ;
				if(trialState.isMet == false)
				{	
					maxDamage = true ;
					bestState = trialState ;
					bestRange = trialRange ;
					break ;
				}
				if(trialState.delta < bestState.delta)
				{
					bestState = trialState ;
					bestRange = trialRange ;
				}
			}
			
			if(maxDamage)
			{
				trialRatio = bestState.fraction ;
				getState() = upState*trialRatio + downState*(1.-trialRatio) ;
			}
			else
			{
				double zero = bestRange.zeroLocation() ;
				if(zero > 0)
					trialRatio = zero ;
				else
					trialRatio = bestState.fraction ;
				getState() = upState*trialRatio + downState*(1.-trialRatio) ;
			}
			
			
			if(std::abs(bestRange.up.fraction*upState+(1.-bestRange.up.fraction)*downState - bestRange.down.fraction*upState - (1.-bestRange.down.fraction)*downState).max()/std::abs(upState-downState).max() < damageDensityTolerance)
			{
				converged = true ;
			}
			

//			if(setChange != 0 || fractured())
//			{
// 				if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() >= 0 && false)
// 				{
// 					Vector delta = upState - downState ;
// 					upState = getState() ;
// 					getState() += s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() * delta ;
// 				}
// 				else
// 				{
//					upState = getState() ;
//					getState() = downState + resphi * (upState - downState) ;
// 				}
//			}
//			else
//			{
// 				if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() <= 0 && false)
// 				{
// 					Vector delta = upState - downState ;
// 					downState = getState() ;
// 					getState() += s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() * delta ;
// 				}
// 				else
// 				{
//					downState = getState() ;
//					getState() = downState + resphi * (upState - downState) ;
// 				}
//			}
		}

//		if(!converged && std::abs(startPosition-getState()).max()*fraction < damageDensityTolerance )
//		{
//
//					
//			if(setChange == 0)
//				std::cerr << "0" << std::flush ;
//			else
//				std::cerr << "1" <<std::flush ;
//			
// 			getState() = downState ;
//			converged = true ;
//		}
		
	}
	
	DamageModel::DamageModel() 
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
		explorationIncrement = 4e-4;4./pow(2., 16) ;
		damageDensityTolerance = 1e-4;1./pow(2., 16) ; // about 1e-4
	} ;
	
	double DamageModel::getThresholdDamageDensity() const
	{
		return thresholdDamageDensity ;
	}
	
	double DamageModel::getSecondaryThresholdDamageDensity() const
	{
		return secondaryThresholdDamageDensity ;
	}
	
	void DamageModel::setThresholdDamageDensity(double d)
	{
		thresholdDamageDensity = d ;
	}
	
	void DamageModel::setSecondaryThresholdDamageDensity(double d)
	{
		secondaryThresholdDamageDensity = d ;
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
