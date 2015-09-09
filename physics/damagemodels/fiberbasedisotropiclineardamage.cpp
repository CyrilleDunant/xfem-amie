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
#include "fiberbasedisotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

FiberBasedIsotropicLinearDamage::FiberBasedIsotropicLinearDamage(double f, double c)  : fibreFraction(f)
{
    alt = false ;
    thresholdDamageDensity = c ;
    getState(true).resize(1, 0.);
    residualStiffnessFraction = 1e-2 ;
    isNull = false ;
}

std::pair< Vector, Vector > FiberBasedIsotropicLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
    return std::make_pair(state, Vector(1., 1)) ;
}

void FiberBasedIsotropicLinearDamage::computeDelta(ElementState & s)
{
    delta = 1.-getState()[0] ;
}

Matrix FiberBasedIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

    if(fractured())
        return m*residualStiffnessFraction ;

    return m*(1.-state[0]) ;
}



bool FiberBasedIsotropicLinearDamage::fractured() const
{
    if(fraction < 0)
        return false ;
    return getState().max() >= thresholdDamageDensity ;
}


FiberBasedIsotropicLinearDamage::~FiberBasedIsotropicLinearDamage()
{
}

void FiberBasedIsotropicLinearDamage::step( ElementState &s , double maxscore)
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
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || !s.getParent()->getBehaviour()->getFractureCriterion()->met())
    {
        converged = true ;
        return ;
    }
    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;//maxscore ;
    double maxScoreInNeighbourhood = maxscore ; //s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s) ;

    if(!fractured() && s.getParent()->getBehaviour()->getFractureCriterion()->met() && std::abs(score - maxScoreInNeighbourhood) < s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()*maxScoreInNeighbourhood)
    {
        state += fibreFraction ;
        for(size_t i = 0 ; i < state.size() ; i++)
        {
            if(state[i] > 1)
                state[i] = 1 ;
        }
        change = true ;
    }
    alt = !alt ;
    converged = true ;
    return ;
}

DamageModel * FiberBasedIsotropicLinearDamage::getCopy() const
{
    FiberBasedIsotropicLinearDamage * ret = new FiberBasedIsotropicLinearDamage(fibreFraction) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}

