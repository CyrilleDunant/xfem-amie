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
#include "spacetimeisotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../../polynomial/vm_function_extra.h"

namespace Amie {

SpaceTimeIsotropicLinearDamage::SpaceTimeIsotropicLinearDamage(double density) 
{
    thresholdDamageDensity = density ;
    getState(true).resize(1, 0.);
    isNull = false ;
}

std::pair< Vector, Vector > SpaceTimeIsotropicLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
    return std::make_pair(state, Vector(1., 1)) ;
}

void SpaceTimeIsotropicLinearDamage::computeDelta(ElementState & s)
{
    delta = 1.-getState()[0] ;
}

Matrix SpaceTimeIsotropicLinearDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE&& accelerate < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*1e-6 ;
    
    double factor = (p.getT()+1.)*.5 ;
    double d = std::min(state[0]+factor*accelerate, 1.) ;
    return m*(1.-d) ;


}

Matrix SpaceTimeIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
// return m;

    if(state.max() < POINT_TOLERANCE && accelerate < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*1e-6 ;
    
    double factor = (p.getT()+1.)*.5 ;
    double d = std::min(state[0]+factor*accelerate, 1.) ;
    return m*(1.-d) ;

}


bool SpaceTimeIsotropicLinearDamage::fractured() const
{
    if(fraction < 0)
        return false ;
    return getState().max() >= thresholdDamageDensity ;
}


SpaceTimeIsotropicLinearDamage::~SpaceTimeIsotropicLinearDamage()
{
}

void SpaceTimeIsotropicLinearDamage::step( ElementState &s , double maxscore)
{
    elementState = &s ;
    converged = true ;
    double timetol = 1e-3 ;
    s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = false ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
        dt = s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT() ;
    }

    change = false ;
    if(accelerate > 0)
    {
        if(maxscore > 0 && maxscore < 1)
            state[0] = std::min(state[0]+(1.-maxscore)*accelerate, 1.) ;
        else if (maxscore < 0)
            state[0] = std::min(state[0]+accelerate, 1.) ;
        change = true ;
    }
    
    
    dt = std::max(s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT(), timetol) ;
    double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
    
    if(!s.getParent()->getBehaviour()->getFractureCriterion() 
        || score <= maxscore-timetol
        ||!s.getParent()->getBehaviour()->getFractureCriterion()->met()
        || maxscore < 0
    )
    {

        accelerate = 0 ;
        return ;
    }
    else if(!fractured() && score > maxscore-timetol)
    {
        
        double originalState = state[0] ;
        accelerate = 0 ;
        double downDamage = originalState ;
        double upDamage = 1 ;

        while(upDamage-downDamage > 1e-6)
        {
            state[0] = (downDamage+upDamage)*.5 ;
            double scoreAtEnd = s.getParent()->getBehaviour()->getFractureCriterion()->gradeAtTime( s, 1 ) ;
            if(scoreAtEnd > 0)
                downDamage = state[0] ;
            else
                upDamage = state[0] ;
        }
        double maxDamage = (downDamage+upDamage)*.5 ; 
        
        upDamage = maxDamage ;
        downDamage = originalState ;
        double damageInitiationTime = 1. - score*2. ;
        while(upDamage-downDamage > 1e-6)
        {
            state[0] = (downDamage+upDamage)*.5 ;
            double scoreAtEnd = s.getParent()->getBehaviour()->getFractureCriterion()->gradeAtTime( s, damageInitiationTime+timetol ) ;
            if(scoreAtEnd > 0)
                downDamage = state[0] ;
            else
                upDamage = state[0] ;
        }
        state[0] = upDamage ;

        double downAccelerate = 0 ;
        double upAccelerate =  (state[0]-originalState)/(timetol/score) ;
        while(upAccelerate-downAccelerate > 1e-6)
        {
            accelerate = (downAccelerate+upAccelerate)*.5 ;
            double scoreAtEnd = s.getParent()->getBehaviour()->getFractureCriterion()->gradeAtTime( s, 1) ;
            if(scoreAtEnd > 0)
                downAccelerate = accelerate ;
            else
                upAccelerate = accelerate ;
        }
        accelerate = downAccelerate*.5 ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }

}


DamageModel * SpaceTimeIsotropicLinearDamage::getCopy() const
{
    SpaceTimeIsotropicLinearDamage * dam = new SpaceTimeIsotropicLinearDamage(thresholdDamageDensity) ;
    dam->copyEssentialParameters( this ) ;
    return dam ;
}

}

