//
// C++ Implementation: isotropiclineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "spacetimebadisotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

// WARNING: these damage models were implemented for testing purposes only!

namespace Amie {

DamageModel * SpaceTimeFixedPointIsotropicLinearDamage::getCopy() const
{
    SpaceTimeFixedPointIsotropicLinearDamage * dam = new SpaceTimeFixedPointIsotropicLinearDamage(fibreFraction, thresholdDamageDensity) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

void SpaceTimeFixedPointIsotropicLinearDamage::step( ElementState &s , double maxscore)
{

    elementState = &s ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
    }

    change = false ;
    converged = true ;
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0 || fractured())
    {
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtTimeStepEnd() ;
    lastState = state[0] ;

    if(score > -POINT_TOLERANCE)
        state[0] = std::min(1., std::max( lastDamage, 1. - (1.-lastState)*(1.-score) ) ) ;
    else
        state[0] = std::max(0., std::max( lastDamage, 1. - (1.-lastState)/(1.-score) ) ) ;

//    std::cout << maxscore << "\t" << lastDamage << "\t" << score << "\t"  << state[0] << std::endl ;

    if(std::abs( state[0] - lastState) > 1e-6)
    {
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }
    else if( score > 0 && score < POINT_TOLERANCE )
    {
        state[0] += fibreFraction ;
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }

//    std::cout << std::abs( state[0] - lastState) << " (" << score << ")"  ;
/*    else
    {
        state[0] = std::max( lastState, state[0] ) ;
    }*/


}

DamageModel * SpaceTimeSequentialIsotropicLinearDamage::getCopy() const
{
    SpaceTimeSequentialIsotropicLinearDamage * dam = new SpaceTimeSequentialIsotropicLinearDamage(fibreFraction, timeTolerance, thresholdDamageDensity) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

void SpaceTimeSequentialIsotropicLinearDamage::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    converged = true ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();

//         if(state.size() +1 != s.getParent()->timePlanes())
//         {
//             state.resize(s.getParent()->timePlanes()-1) ;
//             s.getParent()->getBehaviour()->setTimeDependent( s.getParent()->timePlanes() > 2) ;
//         }

    }

    change = false ;
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
    {
        
        return ;
    }

    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;

	if(!fractured() && score >= 1)/* && score == maxScoreInNeighbourhood*/
	{
		double maxScoreInNeighbourhood = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s) ;
		if(score == maxScoreInNeighbourhood)
		{
			state[state.size()-1] += fibreFraction ;
			change = true ;
			converged = true ;
			s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
		}
	}
    else if(!fractured() && score > 0 && (maxscore - score) < timeTolerance)
    {
        double target = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtTimeStepEnd() ;
        double current = state[0] ;

        state[0] = std::min(1., std::max( current+fibreFraction, 1. - (1.-current)*(1.-target) ) ) ;


        for(size_t i = 0 ; i < state.size() ; i++)
        {
            if(state[i] > 1)
                state[i] = 1. ;
        }
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }
    return ;
}



}
