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
		std::cerr << "\n h" << std::flush ;
		double setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange(s) ;
		std::cerr << "\n i" << std::flush ;
		double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;
		std::cerr << "\n j" << std::flush ;
		bool isInDamagingSet = s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() ;
		std::cerr << "\n k" << std::flush ;
		if(!isInDamagingSet)
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
				states.push_back(PointState(s.getParent()->getBehaviour()->getFractureCriterion()->met(), -setChange, 0., score)) ;
				downState = getPreviousState();
				upState = getPreviousState()+damageIncrement*up ;
				getState() = upState ;
				for(size_t j = 0 ; j < upState.size() ; j++)
					upState[j] = std::min(1.,std::max(downState[j], upState[j])) ;
				for(size_t j = 0 ; j < downState.size() ; j++)
					downState[j] = std::min(upState[j],std::max(0., downState[j])) ;
				trialRatio = 1 ;
				if((upState-downState).min() < 0)
				{
					getState() = upState ;
					converged = true ;
					wasBroken = true ;
				}
				if( (upState-downState).max() < damageDensityTolerance)
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

			states.push_back(PointState(s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange, trialRatio, score)) ;
			std::stable_sort(states.begin(), states.end()) ;	
			double stol = 0.25*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() ;

			//find the most likely midPoint
			RangeState bestRangeNonMet(states[1], states[0]) ;
			RangeState bestRangeZero(states[1], states[0]) ;
			PointState bestStateNonMet = bestRangeNonMet.extrapolate(0.5) ;
			PointState bestStateZero = bestRangeZero.extrapolate(0.5) ;
			bool equilibrium = (states[0].score*states[1].score < 0) ;
			bool foundzero = (states[0].delta*states[1].delta < 0) ;
			double minDelta = states[0].delta ;
			double minScore = states[0].score ;
			size_t minIndexNonMet = 0 ;
			size_t minIndexZero = 0 ;
			if(equilibrium && std::abs(states[0].fraction -states[1].fraction)*(upState-downState).max()  < damageDensityTolerance)
			{
				getState() = downState + (upState-downState)*std::max(states[1].fraction, explorationIncrement/(upState-downState).max()) ;
				converged = true ;
				for(size_t j = 0 ; j < getState().size() ; j++)
					getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
				return ;
			}
			if(foundzero  && std::abs(states[0].fraction -states[1].fraction)*(upState-downState).max()  < damageDensityTolerance )
		       {
				getState() = downState + (upState-downState)*std::max(states[1].fraction, explorationIncrement/(upState-downState).max()) ;
				converged = true ;
				for(size_t j = 0 ; j < getState().size() ; j++)
					getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
				return ;
			}
			if(foundzero && std::abs(bestStateZero.score) < stol || equilibrium && std::abs(bestStateZero.delta) < stol)
			{
				getState() = downState+(upState-downState)*std::max(bestStateZero.fraction, explorationIncrement/(upState-downState).max()) ;
				for(size_t j = 0 ; j < getState().size() ; j++)
					getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
				converged = true ;
				return ;
			}

		       if(std::abs(states[0].fraction -states[1].fraction)*(upState-downState).max()  < damageDensityTolerance)
		       {
				getState() = downState + (upState-downState)*std::max(states[1].fraction, explorationIncrement/(upState-downState).max()) ;
				for(size_t j = 0 ; j < getState().size() ; j++)
					getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
				converged = true ;
				return ;
		       }

			if(!foundzero && !equilibrium)
			{
				for (size_t i = 1 ; i < states.size()-1 ;i++ )
				{
					if(states[i].delta < minDelta)
					{
						minDelta = states[i].delta ;
						minIndexZero = i ;
					}
					if(states[i].score < minScore)
					{
						minScore = states[i].score ;
						minIndexNonMet = i ;
					}
					RangeState trialRange(states[i+1], states[i]) ;
					PointState trialState = trialRange.extrapolate(0.5) ;
					if(states[i].score*states[i+1].score < 0 )
					{
						if(std::abs(states[i].fraction -states[i+1].fraction)*(upState-downState).max()  < damageDensityTolerance)
						{
							getState() = downState+(upState-downState)*std::max(states[i+1].fraction, explorationIncrement/(upState-downState).max())  ;
							for(size_t j = 0 ; j < getState().size() ; j++)
								getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
							converged = true ;
							return ;
						}
						if(std::abs(trialState.score) < stol)
						{
							getState() = downState+(upState-downState)*std::max(trialState.fraction, explorationIncrement/(upState-downState).max()) ;
							for(size_t j = 0 ; j < getState().size() ; j++)
								getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
							converged = true ;
							return ;
						}	
						equilibrium = true ;
						bestStateNonMet = trialState ;
						bestRangeNonMet = trialRange ;
						break ;
					}
					if(states[i].delta*states[i+1].delta < 0 )
					{
                                                if(std::abs(states[i].fraction -states[i+1].fraction)*(upState-downState).max()  < damageDensityTolerance) 
						{
							getState() = downState+(upState-downState)*std::max(states[i+1].fraction, explorationIncrement/(upState-downState).max()) ;
							for(size_t j = 0 ; j < getState().size() ; j++)
								getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
							converged = true ;
							return ;
						}
						if(std::abs(trialState.delta) < stol)
						{
							getState() = downState+(upState-downState)*std::max(trialState.fraction, explorationIncrement/(upState-downState).max()) ;
							for(size_t j = 0 ; j < getState().size() ; j++)
								getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
							converged = true ;
							return ;
						}
						if(states[i+1].fraction*std::abs(upState-downState).max() > damageDensityTolerance )
						{
							foundzero = true ;
							bestStateZero = trialState ;
							bestRangeZero = trialRange ;
							break ;
						}
					}
					
					if(std::abs(states[i].fraction -states[i+1].fraction)*(upState-downState).max()  < damageDensityTolerance)
					{
						getState() = downState+(upState-downState)*std::max(states[i+1].fraction, explorationIncrement/(upState-downState).max())  ;
						for(size_t j = 0 ; j < getState().size() ; j++)
							getState()[j] = std::min(1.,std::max(downState[j], getState()[j])) ;
						converged = true ;
						return ;
					}
					
				}
			}
			RangeState bestRange =  bestRangeNonMet;
			double zeroloc = bestRangeZero.zeroLocation() ;
			double equilibriumloc = bestRangeNonMet.equilibriumLocation() ;
			if(equilibrium)
			{
				trialRatio = bestStateNonMet.fraction ;
// 				if(equilibriumloc > 0)
// 					trialRatio = equilibriumloc ;
				getState() = downState + (upState-downState)*trialRatio ;

			}
			else if(foundzero)
			{
				bestRange = bestRangeZero ;
				trialRatio = bestStateZero.fraction ;
// 				if(zeroloc > 0)
// 					trialRatio = zeroloc ;
				getState() = downState + (upState-downState)*trialRatio ;
			}
			else if(states.size() < 4)
			{
				size_t idx = minIndexNonMet ;
				if(std::abs(minDelta) < std::abs(minScore))
					idx = minIndexZero ;
				bestRange = RangeState(states[idx], states[idx+1]) ;
				trialRatio = bestRange.extrapolate(0.5).fraction ;
				getState() = downState + (upState-downState)*trialRatio ;
			}
			else // there probably is no solution
			{
				getState() = upState ;
				converged = true ;
			}
		}
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
		explorationIncrement = 0.01;
		damageDensityTolerance = 0.5e-3;
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
