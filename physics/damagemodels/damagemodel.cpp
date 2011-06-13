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


void DamageModel::step( ElementState &s )
{
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

	if( wasBroken )
	{
		converged = true ;
		return ;
	}

	Vector damageIncrement = computeDamageIncrement( s ) ;
	std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s ) ;
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;
	bool isInDamagingSet = s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() ;

	if( !isInDamagingSet )
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		converged = true ;

		if( fractured() )
			wasBroken = true ;

		return ;
	}

	bool checkpoint = s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ;

	Vector startPosition = getState() ;

	if( checkpoint ) // initiate iteration
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		states.clear() ;

		getPreviousState() = getState() ;

		if( !fractured() )
		{
			converged = false ;
			change = true ;

			if( damageIncrement.max() > POINT_TOLERANCE_2D )
				damageIncrement /= damageIncrement.max() ;

			double factor = 0.5 ;
			double down = 0 ;
			double up = 1 ;

			//convergence requires consistent errors
			while( std::abs( up - down ) > damageDensityTolerance * .25 )
			{
				getState() = getPreviousState() + damageIncrement * factor ;

				if( getState().max() >= 1 || fractured() )
				{
					up = factor ;
					factor = ( up + down ) / 2 ;
					Vector startingState = getState() ;
					getState() = getPreviousState() + damageIncrement * up ;

					while( !fractured() )
					{
						up += 2.*damageDensityTolerance ;
						getState() = getPreviousState() + damageIncrement * up ;
					}

					getState() = startingState ;
				}
				else
				{
					down = factor ;
					factor = ( up + down ) / 2 ;
				}
			}

			
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), -setChange.first, 0., score, setChange.second ) ) ;
			downState = getPreviousState();
			upState = getPreviousState() + damageIncrement * up ;


			for( size_t j = 0 ; j < upState.size() ; j++ )
				upState[j] = std::min( 1., std::max( downState[j], upState[j] ) ) ;

			for( size_t j = 0 ; j < downState.size() ; j++ )
				downState[j] = std::min( upState[j], std::max( 0., downState[j] ) ) ;

			trialRatio = (1.-0.5*damageDensityTolerance)*up ;
			getState() = getPreviousState() + damageIncrement * up *(1.-0.5*damageDensityTolerance);
// 			std::cout << dynamic_cast<DelaunayTriangle *>(s.getParent())->index << "  " <<  0 << "  "<< score << "  "<< -setChange.first << std::endl ;
			if( ( upState - downState ).min() < 0 )
			{
				getState() = upState ;
				converged = true ;
				wasBroken = true ;
			}

			if( ( upState - downState ).max() < 2.*damageDensityTolerance )
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
	else if( !converged )
	{
		double scoreTolerance = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() ;
		change = true ;

		states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second ) ) ;
		std::stable_sort( states.begin(), states.end() ) ;

// 		std::cout << dynamic_cast<DelaunayTriangle *>(s.getParent())->index << "  " <<  trialRatio << "  "<< score << "  "<< setChange.first << std::endl ;
		double minFraction = states[0].fraction ;
		double maxFraction = states[1].fraction ;
		double prevDelta = states[0].delta ;
		double prevScore = states[0].score ;
		double currentDelta = states[0].delta ;
		double currentScore = states[0].score ;
// 		std::cout << 0 << "  " <<  states[0].fraction << "  "<< states[0].score << "  "<< states[0].delta << std::endl ;
		for(int i = 1 ; i < states.size() ; i++)
		{
// 			std::cout << i << "  " <<  states[i].fraction << "  "<< states[i].score << "  "<< states[i].delta << std::endl ;
			currentDelta = states[i].delta ;
			currentScore = states[i].score ;
			minFraction = states[i-1].fraction ;
			maxFraction = states[i].fraction ;
			if(currentDelta*prevDelta < 0 || currentScore*prevScore < 0 )
			{
				break ;
			}
			else
			{
				prevDelta = states[i].delta ;
				prevScore = states[i].score ;
			}
		}
		trialRatio = (minFraction+maxFraction)*.5 ;
		getState() = downState + ( upState - downState ) * trialRatio ;
		if(std::abs( minFraction - maxFraction ) * ( upState - downState ).max()  < damageDensityTolerance)
		{
			getState() = downState + ( upState - downState ) * std::min(trialRatio+damageDensityTolerance, 1.) ;
			if(states.size() < 4)
				getState() = upState ;
// 			std::cout << std::endl ;
// 			for(int i = 0 ; i < states.size() ; i++)
// 				std::cout << states[i].fraction << "   " << states[i].score << "   " << states[i].delta << std::endl ;
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

	// The exploration increment is crucial for finding
	// the correct distribution of damage: the effect
	// of damage increment on the distribution of
	// fracture criterion scores is non-monotonic.
	damageDensityTolerance = 1e-4;
} ;

double DamageModel::getThresholdDamageDensity() const
{
	return thresholdDamageDensity ;
}

double DamageModel::getSecondaryThresholdDamageDensity() const
{
	return secondaryThresholdDamageDensity ;
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
