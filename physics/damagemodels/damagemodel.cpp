//
// C++ Interface: damagemodel
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "damagemodel.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../../mesher/delaunay.h"

namespace Amie
{

/** \brief Damage model interface */

void DamageModel::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    size_t maxit = iterationNumber ;
    if(alternating)
        maxit /=2 ;

//     double phi = ( 1. + sqrt( 5. ) ) * .5 ;
//     double resphi = 2. - phi ;   //goldensearch
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

        initalState.resize( state.size(), 0. );
    }

    change = false ;
    if(!s.getParent()->getBehaviour()->getFractureCriterion())
    {
        alternate = true ;
        converged = true ;
        return ;
    }
    double max = needGlobalMaximumScore?maxscore:s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s) ;

    std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , max) ;
    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;//maxscore ;
    if( !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
    {
        s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
        alternateCheckpoint = false ;
        // this is necessary because we want to trigger side-effects
        //for example, plasticstrain gets a pointer to s
        computeDamageIncrement( s ) ;
        converged = true ;
        alternate = true ;
        return ;
    }

    std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ;
    
    if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() || (alternateCheckpoint && alternating) ) // initiate iteration
    {

        initalState = state ;
        error = score ;
        s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
        alternateCheckpoint = false ;
        states.clear() ;
        if(!fractured())
        {

            converged = false ;
            change = true ;

            downState = damageIncrement.first ;
            upState = damageIncrement.second ;
            upState= downState+(upState-downState) * s.getParent()->getBehaviour()->getFractureCriterion()->getMinDeltaInNeighbourhood();

            states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.05, -1 ) ) ;
            trialRatio = 1. ;
            getState( true ) = upState ;

        }
        else
        {
            converged = true ;
            alternate = true ;
        }
//        std::cout << "\n-++-"<< damageIncrement.first[0] << "  ;  "<< damageIncrement.second[0] << "-++-" << std::endl ;  
//       std::cout << "\n----"<< score << "  ;  "<< trialRatio << "----" << std::endl ;

    }
    else if( !converged )
    {

//       if(states.size() == 0)
//       {
// 	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.05, -1 ) ) ;
// 	 trialRatio = 1. ;
// 	 getState( true ) = upState ;
// 	 change = true ;
// 	 return ; 
//       }
      
//        [1] 0.0001000000 0.0002137962 0.0004570882 0.0009772372 0.0020892961
//  [6] 0.0044668359 0.0095499259 0.0204173794 0.0436515832 0.0933254301

      double globalAngleShift = std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->maxAngleShiftInNeighbourhood) ;
      int globalMode = s.getParent()->getBehaviour()->getFractureCriterion()->maxModeInNeighbourhood ;
      if(states.size() == 1)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,1., score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = .0001 ;
	 getState( true ) = (upState+downState)*.0001 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 2)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,.0001, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0002137962 ;
	 getState( true ) = downState+(upState-downState)*0.0002137962 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 3)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0002137962, score, setChange.second, globalAngleShift-M_PI*.05, globalMode) ) ;
	 trialRatio = 0.0004570882 ;
	 getState( true ) = downState+(upState-downState)*0.0004570882 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 4)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0004570882, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0009772372 ;
	 getState( true ) = downState+(upState-downState)*0.0009772372 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 5)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0009772372, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0020892961 ;
	 getState( true ) = downState+(upState-downState)*0.0020892961 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 6)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0020892961, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0044668359 ;
	 getState( true ) = downState+(upState-downState)*0.0044668359 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 7)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0044668359, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0095499259 ;
	 getState( true ) = downState+(upState-downState)*0.0095499259 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 8)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0095499259, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0204173794 ;
	 getState( true ) = downState+(upState-downState)*0.0204173794 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 9)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0204173794, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0436515832 ;
	 getState( true ) = downState+(upState-downState)*0.0436515832 ;
	 change = true ;
	 return ; 
      }
      if(states.size() == 10)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0436515832, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = 0.0933254301 ;
	 getState( true ) = downState+(upState-downState)*0.0933254301 ;
	 change = true ;
	 return ; 
      }
      if(allowBackSearch && states.size() == 11)
      {
	 states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0.0933254301, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
	 trialRatio = -.01 ;
	 getState( true ) = downState-(upState-downState)*.01 ;
	 change = true ;
	 return ; 
      }

//       std::cout << "\n----"<< score << "  ;  "<< trialRatio << "----" << std::endl ;
//       if(states.size() == 12)
// 	exit(0) ;

        change = true ;

        states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.05, globalMode ) ) ;
        std::sort( states.begin(), states.end() ) ;

        double minFraction = states[0].fraction ;
        double nextFraction = states[1].fraction ;
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
        error = 10 ;
	double minscore = currentScore ;
	
        for( size_t i = 1 ; i < states.size() ; i++ )
        {
            currentDelta = states[i].delta ;
            currentScore = states[i].score ;
            currentProximity = states[i].proximity ;
            currentShift = states[i].angleShift ;
            currentMode = states[i].mode ;
	    minscore = std::min(minscore, currentScore);
            if(     ((currentDelta > 0    && prevDelta  < 0)     ||
                    (currentDelta < 0     && prevDelta  > 0 ))   ||
                    (std::abs(currentDelta) < POINT_TOLERANCE && std::abs(prevDelta) < POINT_TOLERANCE) ||
                    ((currentScore > 0     && prevScore  < 0  )   ||
                     (currentScore < 0     && prevScore  > 0))    ||
                    (std::abs(currentScore) < POINT_TOLERANCE  && std::abs(prevScore)  < POINT_TOLERANCE) ||
                    ((currentProximity > 0 && prevProximity < 0 )  ||
                     (currentProximity < 0 && prevProximity > 0 )) ||
                    ((currentShift > 0     && prevShift < 0 )     ||
                     (currentShift < 0     && prevShift > 0 ) )   ||
                    currentMode * prevMode < 0
              )
            {
                deltaRoot = (currentDelta > 0     && prevDelta < 0)   ||
                            (currentDelta < 0     && prevDelta  > 0)  ||
                            (std::abs(currentDelta) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5 && 
                            std::abs(prevDelta) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5) ;
                if(deltaRoot)
                    error = std::min(error, std::abs(currentDelta-prevDelta)) ;

                scoreRoot = (currentScore > 0     && prevScore  < 0 ) ||
                            (currentScore < 0     && prevScore  > 0 ) ||
                            (std::abs(currentScore) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5  && 
                            std::abs(prevScore)  < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5) ;
                if(scoreRoot)
                    error = std::min(error, std::abs(currentScore-prevScore)) ;

                proximityRoot = (currentProximity > 0 && prevProximity < 0 )||
                                (currentProximity < 0 && prevProximity > 0 );
                if(proximityRoot)
                    error = std::min(error,std::abs(currentProximity-prevProximity)) ;

                shiftRoot = (currentShift > 0     && prevShift < 0 )    ||
                            (currentShift < 0     && prevShift > 0 );
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
                if(i == states.size()-1 )
                {
                    minFraction = states[i-2].fraction ;
                    nextFraction = states[i-1].fraction ;
                }
            }
        }
        
        if(!(deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)) // we should then minimise the score
	{

	  error = 10 ;
	  for( size_t i = 0 ; i < states.size() ; i++ )
	  {
	    currentScore = states[i].score ;

	    if(std::abs(currentScore-minscore)< 1e-8)
	    {
	      minFraction = states[i].fraction ;
	      if(i == states.size()-1)
	      {
		minFraction = states[i-1].fraction ;
		nextFraction = states[i].fraction ;
		scoreRoot = true ;
		break ;
	      }
	      else if(i == 0)
	      {
		minFraction = states[0].fraction ;
		nextFraction = states[1].fraction ;
		scoreRoot = true ;
		break ;
	      }
	      else
	      {
		if(states[i+1].fraction-states[i].fraction > states[i].fraction-states[i-1].fraction)
		{
		  minFraction = states[i].fraction ;
		  nextFraction = states[i+1].fraction ;
		  scoreRoot = true ;
		}
		else
		{
		  minFraction = states[i-1].fraction ;
		  nextFraction = states[i].fraction ;
		  scoreRoot = true ;
		}
	      }
	    }    
	  }  
	}
        

        trialRatio = (minFraction+nextFraction)*.5/*std::min(minFraction, states[states.size()-2].fraction ) + .27/pow(2, states.size()-1)*/  ;

        getState( true ) = downState + ( upState - downState ) *trialRatio ;

        if( states.size() > maxit-1 && (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot))
        {
// 	  std::cout << "ah " << trialRatio << std::endl ;

            if(ctype == DISSIPATIVE)
            {
                for(size_t i = 0 ; i <  state.size() ; i++)
                {
                    if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE)
                    {
                        state[i] += damageDensityTolerance ;
                        state[i] = std::min(state[i], 1.) ;
                    }
                }

            }
            else if(ctype == CONSERVATIVE)
            {
                for(size_t i = 0 ; i <  state.size() ; i++)
                {
                    if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE)
                    {
                        state[i] += 0.5*damageDensityTolerance ;
                        state[i] = std::min(state[i], 1.) ;
                    }
                }
            }

            converged = true ;
            alternate = true ;

            for(size_t i = 0 ; i <  state.size() ; i++)
                state[i] = std::min(state[i], 1.) ;

            trialRatio = 0 ;
            initalState = state ;
        }
        else if(states.size() > maxit-1)
        {
	    std::cout << "ouch" << std::endl ;
            for(size_t i = 0 ; i <  state.size() ; i++)
            {
                if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE)
                {
                    state[i] = downState[i] + 1e-4;
                    state[i] = std::min(state[i], 1.) ;
                }
            }

            converged = true ;
            alternate = true ;
            trialRatio = 0 ;
            initalState = state ;
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
    iterationNumber = 22 ;

    ctype = DISSIPATIVE ;
    fraction = -1 ;
    converged = true ;
    alternateCheckpoint = false ;
    delta = 1 ;
    effectiveDeltaFraction = 1 ;
    alternate = true ;
    alternating = false ;
    needGlobalMaximumScore = true ;
    // The exploration increment is crucial for finding
    // the correct distribution of damage: the effect
    // of damage increment on the distribution of
    // fracture criterion scores is non-monotonic.
    damageDensityTolerance =  std::max(0.25/pow(2.,iterationNumber), 1e-4) ; //1e-8 ;//1. / pow( 2., 14 );
    thresholdDamageDensity = 1. ;
    secondaryThresholdDamageDensity = 1. ;
    allowBackSearch = false ;
} 

void DamageModel::copyEssentialParameters( const DamageModel * dam ) 
{
    thresholdDamageDensity = dam->getThresholdDamageDensity() ;
    secondaryThresholdDamageDensity = dam->getSecondaryThresholdDamageDensity() ;
    damageDensityTolerance = dam->getDamageDensityTolerance() ;
    residualStiffnessFraction = dam->getResidualStiffnessFraction() ;
}

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
} 
