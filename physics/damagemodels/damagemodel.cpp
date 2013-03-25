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


void DamageModel::step( ElementState &s , double maxscore)
{
	elementState = &s ;
	
	double iterationNumber = 32 ;
	double phi = ( 1. + sqrt( 5. ) ) * .5 ;
	double resphi = 2. - phi ;   //goldensearch
// 		resphi = .5 ;              //bisection
// 		resphi = .1 ;                //down bias

	if( fraction < 0 )
	{
		upState.resize( state.size(), 0. );
		downState.resize( state.size(), 0. );

		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			fraction = s.getParent()->area() ;
		else
			fraction = s.getParent()->volume();

	}
	
	change = false ;
	if(!s.getParent()->getBehaviour()->getFractureCriterion())
	{
		converged = true ;
		return ;
	}
	double maxScoreInNeighbourhood = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood() ;
	std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , maxScoreInNeighbourhood) ;
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;//maxscore ;
	bool isInDamagingSet = s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() ;
	if( !isInDamagingSet )
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		
		// this is necessary because we want to trigger side-effects
		//for example, plasticstrain gets a pointer to s
		computeDamageIncrement( s ) ;
		converged = true ;

		return ;
	}
	
	std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ; 
	
	if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ) // initiate iteration
	{
		error = score ;
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		states.clear() ;
		if(!fractured())
		{
			effectiveDeltaFraction = s.getParent()->getBehaviour()->getFractureCriterion()->getMinDeltaInNeighbourhood()/getDelta() ;
			converged = false ;
			
			change = true ;
			
			downState = damageIncrement.first;
			upState = damageIncrement.second;

// 			for(size_t i = 0 ; i < upState.size() ; i++)
// 				upState[i] = std::min(upState[i], 1.-2.*thresholdDamageDensity) ;
// 			std::cout << "0'   "<< score << std::endl ;
// 			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.01, -1 ) ) ;
			trialRatio = 0. ;
			getState( true ) = downState + ( upState - downState ) * trialRatio*effectiveDeltaFraction  ;
			
		}
		else
		{
			converged = true ;
			alternate = true ;
		}

	}
	else if( !converged )
	{
		double scoreTolerance = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() ;
		double globalAngleShift = s.getParent()->getBehaviour()->getFractureCriterion()->maxAngleShiftInNeighbourhood ;
		int globalMode = s.getParent()->getBehaviour()->getFractureCriterion()->maxModeInNeighbourhood ;
		change = true ;
		
		if(alternating && !alternate)
		{
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.001, globalMode ) ) ;
			alternate = !alternate ;
		}
		else if(alternating && alternate)
		{
			alternate = !alternate ;
			return ;
		}
		else
		{
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.001, globalMode ) ) ;
		}
		
		double n = 2. ;
		if(states.size() <= n)
		{
			if(states.size() == 1)
			{
					getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-3*effectiveDeltaFraction ;
					trialRatio = 0 ;
					return ;
			}
			trialRatio = (double)(states.size()-1)/n ;
			getState( true ) = downState + ( upState - downState ) * trialRatio * effectiveDeltaFraction ;
			for(size_t i = 0 ; i < upState.size() ; i++)
					getState( true )[i] = std::min(getState( true )[i], thresholdDamageDensity-damageDensityTolerance) ;
			return ;
		}
// 		else
// 		{
// 			for(size_t i = 0 ; i < states.size() ; i++ )
// 			{
// 				states[i].print();
// 			}
// 		}

		
		std::stable_sort( states.begin(), states.end() ) ;

		double minFraction = states[0].fraction ;
		double maxFraction = states[1].fraction ;
		double prevDelta = states[0].delta ;
		double prevScore = states[0].score ;
		double prevProximity = states[0].proximity ;
		double prevShift = states[0].angleShift ;
		double prevMode = states[0].mode ;
		double currentDelta = states[0].delta ;
		double currentScore = states[0].score ;
		double currentProximity = states[0].proximity ;
		double currentShift = states[0].angleShift ;
		double currentMode = states[0].mode ;
		bool deltaRoot = false ;
		bool scoreRoot = false ;
		bool proximityRoot = false ;
		bool shiftRoot = false ;
		bool modeRoot = false ;
		for( int i = 1 ; i < states.size() ; i++ )
		{
			currentDelta = states[i].delta ;
			currentScore = states[i].score ;
			currentProximity = states[i].proximity ;
			currentShift = states[i].angleShift ;
			currentMode = states[i].mode ;
			
			maxFraction = states[i].fraction ;

			if( currentDelta * prevDelta < 0 || 
				currentScore * prevScore < 0 || 
				currentProximity * prevProximity < 0 || 
				currentShift * prevShift < 0 || 
				currentMode * prevMode < 0)
			{
				deltaRoot = currentDelta * prevDelta < 0 ;
				if(deltaRoot)
					error = std::abs(currentDelta-prevDelta) ;
				scoreRoot = currentScore * prevScore < 0 ;
				if(scoreRoot)
					error = std::min(error, std::abs(currentScore-prevScore)) ;
				proximityRoot = currentProximity * prevProximity < 0 ;
				if(proximityRoot)
					error = std::min(error,std::abs(currentProximity-prevProximity)) ;
				shiftRoot = currentShift * prevShift < 0 ;
				if(shiftRoot)
					error = std::min(error, std::abs(currentShift-prevShift)) ;
				modeRoot =  currentMode * prevMode < 0 ;

				break ;
				
			}
			else
			{
				prevDelta = states[i].delta ;
				prevScore = states[i].score ;
				prevProximity = states[i].proximity ;
				prevShift = states[i].angleShift ;
				prevMode = states[i].mode ;
				minFraction = states[i].fraction ;
			}
		}
		
		trialRatio = minFraction*0.5 + maxFraction*0.5  ;
		
		getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction ;
// 		for(size_t i = 0 ; i < states.size() ; i++ )
// 		{
// 			states[i].print();
// 		}
		
		if( /*std::abs( minFraction - maxFraction ) < damageDensityTolerance*effectiveDeltaFraction  */states.size() > iterationNumber+n
			&& (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)
		)
		{
			Vector delta = ( upState - downState ) ;
// 			std::cout << deltaRoot << scoreRoot << proximityRoot << shiftRoot << modeRoot << "  "<< setChange.first << "  "<< score<< std::endl ;
			if(ctype == DISSIPATIVE_CENTER)
			{
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-2*effectiveDeltaFraction ; //+ 2.*damageDensityTolerance;
			}
			else if(ctype == CONSERVATIVE_CENTER)
			{
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-3*effectiveDeltaFraction;
			}
			else if(ctype == DISSIPATIVE_MIN)
			{
				trialRatio = minFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-2*effectiveDeltaFraction; //+ 2.*damageDensityTolerance;
			}
			else if(ctype == DISSIPATIVE_MAX)
			{
				trialRatio = maxFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-2*effectiveDeltaFraction; //+ 2.*damageDensityTolerance;
			}
			else if(ctype == CONSERVATIVE_MAX)
			{
				trialRatio = maxFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction ;
			}
			else if(ctype == CONSERVATIVE_MIN)
			{
				trialRatio = minFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction ;
			}
			converged = true ;
			alternate = true ;
			
			for(size_t i = 0 ; i <  state.size() ; i++)
				state[i] = std::min(state[i], 1.) ;

// 			if(trialRatio < 1e-4)
// 			{
// 				std::cout << "case green"<<error << std::endl ;
// 				for(size_t i = 0 ; i < states.size() ; i++ )
// 				{
// 					states[i].print();
// 				}
// 				if(states.size() < 10)
// 				{
// 					std::cout << "\n" << effectiveDeltaFraction << ", "<<error << std::endl ;
// 					for(size_t i = 0 ; i < states.size() ; i++ )
// 					{
// 						states[i].print();
// 					}
// 					exit(0) ;
// 				}
// 			}
			trialRatio = 0 ;
		}
		else if(states.size() > iterationNumber+n)
		{
// 			std::cout << "case red"<<error << std::endl ;
// 			for(size_t i = 0 ; i < states.size() ; i++ )
// 			{
// 				states[i].print();
// 			}
			
			double minscoreRatio = states[0].fraction ;
			double mindeltaRatio = states[0].fraction ;
			double minscore =  states[0].score ;
			double mindelta =  states[0].delta ;
			for(size_t i = 1 ; i <  states.size() ; i++)
			{
				if(states[i].score < minscore)
				{
					minscoreRatio = states[i].fraction ;
					minscore = states[i].score ;
				}
				if(states[i].delta < mindelta)
				{
					mindeltaRatio = states[i].fraction ;
					mindelta = states[i].delta ;
				}
			}
// 			std::cout << minscoreRatio << std::endl ;
			getState( true ) = downState + ( upState - downState )*.5 +damageDensityTolerance;
// 			getState( true ) = downState + ( upState - downState )*effectiveDeltaFraction*.5 ;
			for(size_t i = 0 ; i <  state.size() ; i++)
				state[i] = std::min(state[i], 1.) ;
// 			
			converged = true ;
			alternate = true ;
			trialRatio = 0 ;
			

// 				std::cout << "\n" << effectiveDeltaFraction << ", "<<error << std::endl ;
// 				for(size_t i = 0 ; i < states.size() ; i++ )
// 				{
// 					states[i].print();
// 				}
// 				if(std::max(minscoreRatio,mindeltaRatio) < 1e-3)
// 					exit(0) ;
		}
	}
}


void DamageModel::postProcess()
{
}

DamageModel::DamageModel(): state(0)
{
	elementState = nullptr ;
	change = false ;
	isNull = true ;
	haslimit = false ;
	error = 1 ;

	ctype = DISSIPATIVE_CENTER ;
	fraction = -1 ;
	converged = true ;
	delta = 1 ;
	effectiveDeltaFraction = 1 ;
	alternating = false ;
	alternate = false ;
	// The exploration increment is crucial for finding
	// the correct distribution of damage: the effect
	// of damage increment on the distribution of
	// fracture criterion scores is non-monotonic.
	damageDensityTolerance =  1e-5 ; //1e-8 ;//1. / pow( 2., 14 );
	thresholdDamageDensity = 1.-2e-5 ;
	secondaryThresholdDamageDensity = 1.-1e-5 ;
} ;

double DamageModel::getThresholdDamageDensity() const
{
	return thresholdDamageDensity ;
}

double DamageModel::getSecondaryThresholdDamageDensity() const
{
	return secondaryThresholdDamageDensity ;
}

Vector &DamageModel::getState( bool )
{
	return state ;
}

void DamageModel::setThresholdDamageDensity( double d )
{
	thresholdDamageDensity = d ;
}

void DamageModel::setSecondaryThresholdDamageDensity( double d )
{
	secondaryThresholdDamageDensity = d ;
}

void DamageModel::setDamageDensityTolerance( double d )
{
	damageDensityTolerance = d ;
}

bool DamageModel::changed() const
{
	return change ;
}
} ;
