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

    double phi = ( 1. + sqrt( 5. ) ) * .5 ;
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
    double max = -1 ;
    if(needGlobalMaximumScore)
        max = maxscore ;
    else
        max = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s) ;

    std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , max) ;
    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;//maxscore ;
    if( !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
    {
        s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );

        // this is necessary because we want to trigger side-effects
        //for example, plasticstrain gets a pointer to s
        computeDamageIncrement( s ) ;
        converged = true ;
        alternate = true ;
        return ;
    }

    std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ;


    if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && !alternate) // initiate iteration
    {
        initalState = state ;
        error = score ;
        s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
        states.clear() ;
        if(!fractured())
        {
// 			effectiveDeltaFraction = s.getParent()->getBehaviour()->getFractureCriterion()->getMinDeltaInNeighbourhood()/getDelta() ;
//             iterationNumber = round(log2(s.getParent()->getBehaviour()->getFractureCriterion()->getMinDeltaInNeighbourhood()/(getDelta()*damageDensityTolerance))*.25) ;
            converged = false ;

            change = true ;

            downState = damageIncrement.first ;
            upState = damageIncrement.second ;
            upState= downState+(upState-downState) * s.getParent()->getBehaviour()->getFractureCriterion()->getMinDeltaInNeighbourhood();

            states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second, -M_PI*.01, -1 ) ) ;
            trialRatio = 1. ;
            getState( true ) = upState ;

        }
        else
        {
            converged = true ;
            alternate = true ;
        }

    }
    else if( !converged && !alternate)
    {

        double globalAngleShift = s.getParent()->getBehaviour()->getFractureCriterion()->maxAngleShiftInNeighbourhood ;
        int globalMode = s.getParent()->getBehaviour()->getFractureCriterion()->maxModeInNeighbourhood ;
        change = true ;

        states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.001, globalMode ) ) ;
        std::sort( states.begin(), states.end() ) ;

        double minFraction = states[0].fraction ;
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
        for( size_t i = 1 ; i < states.size() ; i++ )
        {
            currentDelta = states[i].delta ;
            currentScore = states[i].score ;
            currentProximity = states[i].proximity ;
            currentShift = states[i].angleShift ;
            currentMode = states[i].mode ;
            if(     ((currentDelta > 0    && prevDelta  < 0)     ||
                    (currentDelta < 0     && prevDelta  > 0 ))   ||
                    (std::abs(currentDelta) < POINT_TOLERANCE_2D && std::abs(prevDelta) < POINT_TOLERANCE_2D) ||
                    ((currentScore > 0     && prevScore  < 0  )   ||
                     (currentScore < 0     && prevScore  > 0))    ||
                    (std::abs(currentScore) < POINT_TOLERANCE_2D  && std::abs(prevScore)  < POINT_TOLERANCE_2D) ||
                    ((currentProximity > 0 && prevProximity < 0 )  ||
                     (currentProximity < 0 && prevProximity > 0 )) ||
                    ((currentShift > 0     && prevShift < 0 )     ||
                     (currentShift < 0     && prevShift > 0 ) )   ||
                    currentMode * prevMode < 0
              )
            {
                deltaRoot = (currentDelta > 0     && prevDelta ) < 0  ||
                            (currentDelta < 0     && prevDelta ) > 0  ||
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
            }
        }

        trialRatio = std::min(minFraction, states[states.size()-2].fraction ) + .27/pow(2, states.size()-1)  ;

        getState( true ) = downState + ( upState - downState ) *trialRatio ;


        if( states.size() > iterationNumber && (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)
          )
        {
// 			std::cout << effectiveDeltaFraction << "  " << trialRatio<< "  "<< getState().max() << std::endl ;

            if(ctype == DISSIPATIVE)
            {
                getState( true ) = downState + ( upState - downState ) * trialRatio;
                for(size_t i = 0 ; i <  state.size() ; i++)
                {
                    if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE_2D)
                    {
                        state[i] += damageDensityTolerance ;
                        state[i] = std::min(state[i], 1.) ;
                    }
                }

            }
            else if(ctype == CONSERVATIVE)
            {
                getState( true ) = downState + ( upState - downState ) * trialRatio  ;
                for(size_t i = 0 ; i <  state.size() ; i++)
                {
                    if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE_2D)
                    {
                        state[i] += 0.25*damageDensityTolerance ;
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
        else if(states.size() > iterationNumber)
        {

            for(size_t i = 0 ; i <  state.size() ; i++)
            {
                if(std::abs( upState[i] - downState [i]) > POINT_TOLERANCE_2D)
                {
                    state[i] = downState[i] + 0.05;
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
    iterationNumber = 8 ;

    ctype = DISSIPATIVE ;
    fraction = -1 ;
    converged = true ;
    delta = 1 ;
    effectiveDeltaFraction = 1 ;
    alternate = false ;
    needGlobalMaximumScore = false ;
    // The exploration increment is crucial for finding
    // the correct distribution of damage: the effect
    // of damage increment on the distribution of
    // fracture criterion scores is non-monotonic.
    damageDensityTolerance =  std::max(0.25/pow(2.,iterationNumber), 1e-6) ; //1e-8 ;//1. / pow( 2., 14 );
    thresholdDamageDensity = 1. ;
    secondaryThresholdDamageDensity = 1. ;
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
