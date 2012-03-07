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


void DamageModel::stepBack()
{
	change = previouschange ;
	previousstate = previouspreviousstate ;
	state = previousstate ;
	previousauxiliarystate = previouspreviousauxiliarystate ;
	auxiliarystate = previousauxiliarystate ;
}

void DamageModel::step( ElementState &s )
{
	previouschange = change ;
	previouspreviousstate = previousstate ;
	previousstate = state;
	previouspreviousauxiliarystate = previousauxiliarystate ;
	previousauxiliarystate = auxiliarystate ;
	elementState = &s ;
	
	
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
	std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s ) ;
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood() ;
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

// 			std::cout << "0'   "<< score << std::endl ;
// 			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.01, -1 ) ) ;
			trialRatio = 1 ;
			getState( true ) = downState + ( upState - downState ) * trialRatio*effectiveDeltaFraction  ;
			
		}
		else
		{
			converged = true ;
		}

	}
	else if( !converged )
	{
		double scoreTolerance = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() ;
		double globalAngleShift = s.getParent()->getBehaviour()->getFractureCriterion()->maxAngleShiftInNeighbourhood ;
		int globalMode = s.getParent()->getBehaviour()->getFractureCriterion()->maxModeInNeighbourhood ;
		change = true ;
		
		states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.01, globalMode ) ) ;

		if(states.size() == 1)
		{
			trialRatio = 0 ;
			getState( true ) = downState + ( upState - downState ) * 2.*damageDensityTolerance*effectiveDeltaFraction ;
			for(size_t i = 0 ; i <  state.size() ; i++)
				state[i] = std::min(state[i], 1.) ;
			return ;
		}
		if(states.size() == 2)
		{
			if(states[1].score < 0 /*||
				states[1].proximity < 1e-5 && states[1].score < 1e-5
			*/)
			{
				change = false ;
				converged = true ;
				return ;
			}
			
			trialRatio = .5 ;
			getState( true ) = downState + ( upState - downState ) * .5*effectiveDeltaFraction ;
			return ;
		}
		if(states.size() == 3)
		{
			trialRatio = .75 ;
			getState( true ) = downState + ( upState - downState ) * .75*effectiveDeltaFraction ;
			return ;
		}
		if(states.size() == 4)
		{
			trialRatio = .25 ;
			getState( true ) = downState + ( upState - downState ) * .25*effectiveDeltaFraction ;
			return ;
		}
		
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
		
		if( std::abs( minFraction - maxFraction ) < damageDensityTolerance*effectiveDeltaFraction  
			&& (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)
		)
		{

// 			std::cout << dynamic_cast<DelaunayTriangle *>(s.getParent())->index<< "  "<< getState().max() << "   "<< score << std::endl ;
			if(ctype == DISSIPATIVE_CENTER)
			{
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-4*effectiveDeltaFraction ; //effectiveDeltaFraction*2.*damageDensityTolerance;
			}
			else if(ctype == CONSERVATIVE_CENTER)
			{
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction ;
			}
			else if(ctype == DISSIPATIVE_MIN)
			{
				trialRatio = minFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-4*effectiveDeltaFraction ; //+ effectiveDeltaFraction*2.*damageDensityTolerance;
			}
			else if(ctype == DISSIPATIVE_MAX)
			{
				trialRatio = maxFraction ;
				getState( true ) = downState + ( upState - downState ) *trialRatio*effectiveDeltaFraction + 1e-4*effectiveDeltaFraction; //+ effectiveDeltaFraction*2.*damageDensityTolerance;
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
			
			if(states.size() < 10)
			{
				std::cout << "\n" << effectiveDeltaFraction << ", "<<error << std::endl ;
				for(size_t i = 0 ; i < states.size() ; i++ )
				{
					states[i].print();
				}
				exit(0) ;
			}
		}
		else if(std::abs( minFraction - maxFraction ) < damageDensityTolerance*effectiveDeltaFraction)
		{
			getState( true ) = downState + ( upState - downState ) *1.*effectiveDeltaFraction ;
			
			for(size_t i = 0 ; i <  state.size() ; i++)
				state[i] = std::min(state[i], 1.) ;
			
// 			for(size_t i = 0 ; i <  states.size() ; i++)
// 				states[i].print() ;
// 			
// 			exit(0) ;
			converged = true ;
		}
	}
}


void DamageModel::postProcess()
{
}

DamageModel::DamageModel(): state(0), previousstate(0), previouspreviousstate(0), auxiliarystate(0), previousauxiliarystate(0),previouspreviousauxiliarystate(0)
{
	elementState = NULL ;
	previouschange = false ;
	change = false ;
	isNull = true ;
	haslimit = false ;
	error = 1 ;

	ctype = CONSERVATIVE_CENTER ;
	fraction = -1 ;
	converged = true ;
	delta = 1 ;
	effectiveDeltaFraction = 1 ;
	// The exploration increment is crucial for finding
	// the correct distribution of damage: the effect
	// of damage increment on the distribution of
	// fracture criterion scores is non-monotonic.
	damageDensityTolerance =  1e-4 ; //1e-8 ;//1. / pow( 2., 14 );
	thresholdDamageDensity = 1.-std::max(1e-4, 2.*damageDensityTolerance) ;
	secondaryThresholdDamageDensity = 1.-std::max(1e-4, 2.*damageDensityTolerance) ;
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


Vector DamageModel::smoothedState( const ElementState &s , bool setUpdate) const
{
	if( !s.getParent()->getBehaviour()->getFractureCriterion() )
		return getState() ;


	double fracturedFraction = 0 ;
	auto fiterator =  s.getParent()->getBehaviour()->getFractureCriterion()->factors.begin() ;
	std::vector <unsigned int> cache = s.getParent()->getBehaviour()->getFractureCriterion()->cache ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		Vector stra = getState()*(*fiterator) ;
		fiterator++ ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *s.getParent()->getBehaviour()->getFractureCriterion()->mesh2d )[cache[i]] ) ;
			if(ci->getBehaviour()->fractured())
			{
				fracturedFraction += *fiterator ;
				fiterator++ ;
				continue ;
			}

			stra += ci->getBehaviour()->getDamageModel()->getState()*(*fiterator) ;
			fiterator++ ;
		}
		stra /= s.getParent()->getBehaviour()->getFractureCriterion()->factors.back()-fracturedFraction ;
		return stra ;
	}
	else if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		Vector stra = getState() *(*fiterator) ;
		fiterator++ ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *s.getParent()->getBehaviour()->getFractureCriterion()->mesh3d )[cache[i]] ) ;
			
			if( ci->getBehaviour()->fractured() || *fiterator < POINT_TOLERANCE_2D)
			{
				fracturedFraction += *fiterator ;
				fiterator++ ;
				continue ;
			}
			stra += ci->getState().getStrainAtCenter()*(*fiterator) ;
			fiterator++ ;
		}
		stra /= s.getParent()->getBehaviour()->getFractureCriterion()->factors.back()-fracturedFraction ;
		return stra ;
	}
	
	return getState() ;
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
