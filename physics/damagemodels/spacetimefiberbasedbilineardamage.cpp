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
#include "spacetimefiberbasedbilineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

Matrix SpaceTimeFiberBasedBilateralLinearDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;

    std::pair<double, double> props = Tensor::getIsotropicMaterialParameters( m, mode, pt ) ;
    if(state[0] < thresholdDamageDensity)
        props.first *= 1.-state[0] ;
    else
        props.first *= residualStiffnessFraction ;

    if(state[1] < secondaryThresholdDamageDensity)
        props.second *= 1.-state[1] ;
    else
        props.second *= residualStiffnessFraction ;

    Matrix test = Tensor::cauchyGreen( props, m.numCols() == 6 ? SPACE_THREE_DIMENSIONAL : SPACE_TWO_DIMENSIONAL, pt, mode ) ;
    if(test[0][0] < 0) { test = m*residualStiffnessFraction ; }
    return test ;
}

Matrix SpaceTimeFiberBasedBilateralLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;

    std::pair<double, double> props = Tensor::getIsotropicMaterialParameters( m, mode, pt ) ;
    if(state[0] < thresholdDamageDensity)
        props.first *= 1.-state[0] ;
    else
        props.first *= residualStiffnessFraction ;

    if(state[1] < secondaryThresholdDamageDensity)
        props.second *= 1.-state[1] ;
    else
        props.second *= residualStiffnessFraction ;

    Matrix test = Tensor::cauchyGreen( props, m.numCols() == 6 ? SPACE_THREE_DIMENSIONAL : SPACE_TWO_DIMENSIONAL, pt, mode ) ;
    if(test[0][0] < 0) { test = m*residualStiffnessFraction ; }
//    test.print() ;
    return test ;
}


DamageModel * SpaceTimeFiberBasedBilateralLinearDamage::getCopy() const
{
    SpaceTimeFiberBasedBilateralLinearDamage * dam = new SpaceTimeFiberBasedBilateralLinearDamage(fibreFraction, timeTolerance, thresholdDamageDensity, secondaryThresholdDamageDensity, mode, pt) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

bool SpaceTimeFiberBasedBilateralLinearDamage::fractured(int direction) const
{
    if(fraction < 0)
        return false ;

    if(direction == 0)
       return getState()[0] >= thresholdDamageDensity ;
    if(direction == 1)
       return getState()[1] >= secondaryThresholdDamageDensity ;

    return (getState().max() >= thresholdDamageDensity && getState()[1] > secondaryThresholdDamageDensity) ;
}


void SpaceTimeFiberBasedBilateralLinearDamage::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    converged = true ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
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
			if( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0, 1.-timeTolerance) )
				state[0] += fibreFraction ;
			if( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1, 1.-timeTolerance) )
				state[1] += fibreFraction ;
			change = true ;
			converged = true ;
			s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
		}
	}
    else if(!fractured() && score > 0 && (maxscore - score) < timeTolerance)
    {
	if( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0, maxscore-timeTolerance) )
		state[0] += fibreFraction ;
	if( s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1, maxscore-timeTolerance) )
		state[1] += fibreFraction ;

        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }

    if(state[0] > thresholdDamageDensity)
        state[0] = 1. ;
    if(state[1] > secondaryThresholdDamageDensity)
        state[1] = 1. ;

//    std::cout << "|||" << state[0] << ";" << state[1] << std::endl ;

        for(size_t i = 0 ; i < state.size() ; i++)
        {
            if(state[i] > 1)
                state[i] = 1. ;
        }

    return ;
}



}
