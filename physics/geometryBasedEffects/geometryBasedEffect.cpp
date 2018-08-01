//
// C++ Interface: contact model
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "geometryBasedEffect.h"
#include "../collisiondetectors/collisiondetector.h"
#include "../../mesher/delaunay.h"

namespace Amie
{

/** \brief Damage model interface */

void GeometryBasedEffect::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    size_t maxit = iterationNumber ;
    if(alternating)
        maxit /=2 ;

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
    if(!s.getParent()->getBehaviour()->getCollisionDetection() )
    {
        alternate = true ;
        converged = true ;
        return ;
    }
    

    double max = maxscore ;
//     std::cout << "a "<< converged << std::endl ;
    std::pair<double, double> setChange = s.getParent()->getBehaviour()->getCollisionDetection()->setChange( s , max) ;
    double score = maxscore ;
    if(s.getParent()->getBehaviour()->getFractureCriterion())
        score = std::max(score, s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) ;
    if( !s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() )
    {
        s.getParent()->getBehaviour()->getCollisionDetection()->setCheckpoint( false );
        alternateCheckpoint = false ;
        // this is necessary because we want to trigger side-effects
        //for example, plasticstrain gets a pointer to s
        computeDamageIncrement( s ) ;
        converged = true ;
        alternate = true ;
        return ;
    }

    std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ;
//     std::cout << "b "<< converged << std::endl ;
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint() || (alternateCheckpoint && alternating) ) // initiate iteration
    {
        initalState = state ;
        error = score ;
        s.getParent()->getBehaviour()->getCollisionDetection()->setCheckpoint( false );

        alternateCheckpoint = false ;
        states.clear() ;
        if(!fractured())
        {
            converged = false ;
            change = true ;

            downState = damageIncrement.first ;
            upState = damageIncrement.second ;
            upState = downState+(upState-downState) * s.getParent()->getBehaviour()->getCollisionDetection()->getMinDeltaInNeighbourhood();

            states.push_back( PointState( s.getParent()->getBehaviour()->getCollisionDetection()->met(), setChange.first,0., score, setChange.second, -M_PI*.025, -1 ) ) ;
            trialRatio = 1. ;
            getState( true ) = upState ;
//             std::cout << "c "<< converged << std::endl ;

        }
        else
        {
            converged = true ;
            alternate = true ;
        }

    }
    else if( !converged )
    {
        double globalAngleShift = std::abs(s.getParent()->getBehaviour()->getCollisionDetection()->getMaxAngleShiftInNeighbourhood()) ;
        int globalMode = s.getParent()->getBehaviour()->getCollisionDetection()->getMaxAngleShiftInNeighbourhood() ;

        change = true ;

        states.push_back( PointState( s.getParent()->getBehaviour()->getCollisionDetection()->met(), setChange.first, trialRatio, score, setChange.second, globalAngleShift-M_PI*.025, globalMode ) ) ;

        if(states.size()-2  < ratios.size())
        {
            trialRatio = ratios[states.size()-2] ;
            getState( true ) = downState + ( upState - downState ) *trialRatio ;
            return ;
        }
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
                deltaRoot = (currentDelta > 0     && prevDelta < 0)   ||
                            (currentDelta < 0     && prevDelta  > 0)  ||
                            (std::abs(currentDelta) < s.getParent()->getBehaviour()->getCollisionDetection()->getScoreTolerance()*.5 &&
                             std::abs(prevDelta) < s.getParent()->getBehaviour()->getCollisionDetection()->getScoreTolerance()*.5) ;
                if(deltaRoot)
                    error = std::min(error, std::abs(currentDelta-prevDelta)) ;

                scoreRoot = (currentScore > 0     && prevScore  < 0 ) ||
                            (currentScore < 0     && prevScore  > 0 ) ||
                            (std::abs(currentScore) < s.getParent()->getBehaviour()->getCollisionDetection()->getScoreTolerance()*.5  &&
                             std::abs(prevScore)  < s.getParent()->getBehaviour()->getCollisionDetection()->getScoreTolerance()*.5) ;
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
        

        if(!(deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot)) // we should then minimise the score or proximity.
        {
            double del = 0 ;
            double fac = std::min(std::max(.75/((std::abs(score) + std::abs(setChange.first) + std::abs(setChange.second) )), 0.25), 1.) ;
            if(std::max(setChange.first, setChange.second) > 4.*damageDensityTolerance)
            {
                if(setChange.first <= setChange.second && setChange.first <= score)
                    del = fac*setChange.first ;
                else if(setChange.second <= setChange.first && setChange.second <= score)
                    del = fac*setChange.second ;
                else
                    del = fac*score ;
            }
            else
                del = fac*score ;


            trialRatio = std::min(std::max(initialRatio + del, 0.), 1.) ;//initialRatio+damageDensityTolerance*.175 ;
            
            deltaRoot = true ;
        }
        
        
        trialRatio = std::max(std::min(trialRatio, 1.), 0.) ;
        getState( true ) = downState + ( upState - downState ) *trialRatio /*+ damageDensityTolerance*.25*/;


        if( states.size() > maxit-1 && (deltaRoot || scoreRoot || proximityRoot || shiftRoot || modeRoot))
        {

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


GeometryBasedEffect::GeometryBasedEffect(): DamageModel() {
//     change = false ;
//     iterationNumber = 64 ;
}


}
