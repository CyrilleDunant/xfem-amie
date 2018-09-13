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
#include "../collisiondetectors/collisiondetector.h"
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
    double score = maxscore ;
    if(s.getParent()->getBehaviour()->getCollisionDetection())
        score = std::max(score, s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState()) ;
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

            states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.025, -1 ) ) ;
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

      double globalAngleShift = std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->maxAngleShiftInNeighbourhood) ;
      int globalMode = s.getParent()->getBehaviour()->getFractureCriterion()->maxModeInNeighbourhood ;

        change = true ;

    states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.025, globalMode ) ) ;
	
	if(states.size()-2  < ratios.size()/*trialRatio < 1.001*/)
	{
//       std::cout << trialRatio << "  " << states.back().score<< "  " << states.back().delta<< "  " << states.back().proximity<< "  " << states.back().angleShift <<std::endl ;
	  trialRatio = ratios[states.size()-2] ;
	  getState( true ) = downState + ( upState - downState ) *trialRatio ;
      
	  return ;
	}
// 	exit(0) ;
	double initialRatio = trialRatio ;
	
	
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
        
        for(  size_t i = 1 ; i < states.size() ; i++ )
        {
            currentDelta = states[i].delta ;
            currentScore = states[i].score ;
            currentProximity = states[i].proximity ;
            currentShift = states[i].angleShift ;
            currentMode = states[i].mode ;
            minscore = std::min(minscore, currentScore);
            if(     ((currentDelta > 0    && prevDelta  < 0)     ||
                    (currentDelta < 0     && prevDelta  > 0 ))   ||
                    ((currentScore > 0     && prevScore  < 0  )   ||
                     (currentScore < 0     && prevScore  > 0))    ||
                    ((currentProximity > 0 && prevProximity < 0 )  ||
                     (currentProximity < 0 && prevProximity > 0 )) ||
                    ((currentShift > 0     && prevShift < 0 )     ||
                     (currentShift < 0     && prevShift > 0 ) )   ||
                    currentMode * prevMode < 0
              )
            {
                deltaRoot = (currentDelta > 0     && prevDelta <= 0)   ||
                            (currentDelta <= 0     && prevDelta  > 0)  ||
                            (std::abs(currentDelta) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5 && 
                            std::abs(prevDelta) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5) ;
                if(deltaRoot)
                    error = std::min(error, std::abs(currentDelta-prevDelta)) ;

                scoreRoot = (currentScore > 0     && prevScore  <= 0 ) ||
                            (currentScore <= 0     && prevScore  > 0 ) ||
                            (std::abs(currentScore) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5  && 
                            std::abs(prevScore)  < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*.5) ;
                if(scoreRoot)
                    error = std::min(error, std::abs(currentScore-prevScore)) ;

                proximityRoot = (currentProximity > 0 && prevProximity <= 0 )||
                                (currentProximity <= 0 && prevProximity > 0 );
                if(proximityRoot)
                    error = std::min(error,std::abs(currentProximity-prevProximity)) ;

                shiftRoot = (currentShift > 0     && prevShift <= 0 )    ||
                            (currentShift <= 0     && prevShift > 0 );
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
        
        trialRatio = (minFraction+nextFraction)*.5  ;
        if(deltaRoot)
        {
            trialRatio = (minFraction*std::abs(currentDelta)/(std::abs(prevDelta)+std::abs(currentDelta)) +nextFraction*std::abs(prevDelta)/(std::abs(prevDelta)+std::abs(currentDelta))) ;
        }
        else if(scoreRoot)
        {
            trialRatio = (minFraction*std::abs(currentScore)/(std::abs(prevScore)+std::abs(currentScore)) +nextFraction*std::abs(prevScore)/(std::abs(prevScore)+std::abs(currentScore))) ;
        }
        else if(proximityRoot)
        {
            trialRatio = (minFraction*std::abs(currentProximity)/(std::abs(prevProximity)+std::abs(currentProximity)) +nextFraction*std::abs(prevProximity)/(std::abs(prevProximity)+std::abs(currentProximity))) ;
        }
        
        
//         if(!(deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)) // we should then minimise the score or proximity.
//         {
//             double del = 0 ;
//             double fac = std::min(std::max(.75/((std::abs(score) + std::abs(setChange.first) + std::abs(setChange.second) )), 0.25), 1.) ;
//             if(std::max(setChange.first, setChange.second) > 4.*damageDensityTolerance)
//             {
//                 if(setChange.first <= setChange.second && setChange.first <= score)
//                 del = fac*setChange.first ;
//                 else if(setChange.second <= setChange.first && setChange.second <= score)
//                 del = fac*setChange.second ;
//                 else 
//                 del = fac*score ;
//             }
//             else
//                 del = fac*score ;
//                 
//             
//             trialRatio = 1e-3 ; //std::min(std::max(initialRatio + del, 0.), 1.) ;//initialRatio+damageDensityTolerance*.175 ;
//         // 	  std::cout << "! "<< trialRatio << std::endl ;
//             deltaRoot = true ;
//         }
        trialRatio = std::max(std::min(trialRatio, 1.), 0.) ;
    //     std::cout << trialRatio << "  " << score <<std::endl ;
        getState( true ) = downState + ( upState - downState ) *trialRatio /*+ damageDensityTolerance*.25*/;


        if( states.size() > maxit-1 && (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot))
        {
// 	  std::cout << "ah " << trialRatio << std::endl ;

            double initalScore = states[0].score ;

            if(score > initalScore)
            {
                state = downState ;
            }
            
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
                        state[i] += 0.05*damageDensityTolerance ;
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
            std::cout << "\033[1;31m\'\033[0m" << std::flush ;
            for(size_t i = 0 ; i <  state.size() ; i++)
            {
                if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE)
                {
                    state[i] = downState[i] + ( upState[i] - downState[i] ) *.1;
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
    iterationNumber = ratios.size()+6 ;

    ctype = DISSIPATIVE ;
    fraction = -1 ;
    converged = true ;
    alternateCheckpoint = false ;
    delta = 1 ;
    effectiveDeltaFraction = 1 ;
    alternate = true ;
    alternating = false ;
    needGlobalMaximumScore = false ;
    // The exploration increment is crucial for finding
    // the correct distribution of damage: the effect
    // of damage increment on the distribution of
    // fracture criterion scores is non-monotonic.
    damageDensityTolerance =  1e-4 ; //std::max(0.25/pow(2.,iterationNumber), 0.5e-3) ; //1e-8 ;//1. / pow( 2., 14 );
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
    getState(true).resize(dam->getState().size()) ;
    getState(true) = dam->getState() ;
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
