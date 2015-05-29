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

SpaceTimeIsotropicLinearDamage::SpaceTimeIsotropicLinearDamage(double f, double density)  : fibreFraction(f)
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
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*0 ;

    double factor = (p.getT()+1.)*.5 ;
    double d = std::min(state[0]+factor*std::max(dt, 1e-4)*accelerate*fibreFraction, 1.) ;
    return m*(1.-d) ;


}

Matrix SpaceTimeIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*0 ;
    
    double factor = (p.getT()+1.)*.5 ;
    double d = std::min(state[0]+factor*std::max(dt, 1e-4)*accelerate*fibreFraction, 1.) ;
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
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
        dt = s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT() ;
    }

    change = false ;    
    state[0] = std::min(state[0]+dt*accelerate*fibreFraction, 1.) ;
    dt = s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT() ;
    
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || !s.getParent()->getBehaviour()->getFractureCriterion()->met() || std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() - maxscore) >= 1e-2)
    {
        accelerate = 0 ;
        return ;
    }
    else if(!fractured() && std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() - maxscore) < 1e-2)
    {
        dt *= s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
        accelerate += .5 ;
        double maxAccelerate = 1 ;
        if(fibreFraction*std::max(dt, 1e-4) > POINT_TOLERANCE)
            maxAccelerate = (1.-state[0])/(fibreFraction*std::max(dt, 1e-4)) ;
        accelerate = std::min(accelerate, std::min(maxAccelerate, 10000.)) ;
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }

}


DamageModel * SpaceTimeIsotropicLinearDamage::getCopy() const
{
    SpaceTimeIsotropicLinearDamage * dam = new SpaceTimeIsotropicLinearDamage(fibreFraction, thresholdDamageDensity) ;
    return dam ;
}

}

