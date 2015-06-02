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

SpaceTimeIsotropicLinearDamage::SpaceTimeIsotropicLinearDamage(double f, double density)  : overdamage(f)
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
        return m*residualStiffnessFraction ;

    double factor = (p.getT()+1.)*.5*overdamage ;
    double d = std::min(state[0]+factor*std::max(dt, 1e-4)*accelerate, 1.) ;
    return m*(1.-d) ;


}

Matrix SpaceTimeIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
    if(state.max() < POINT_TOLERANCE)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;
    
    double factor = (p.getT()+1.)*.5*overdamage ;
    double d = std::min(state[0]+factor*std::max(dt, 1e-4)*accelerate, 1.) ;
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
    double timetolerance = 1e-3 ;
    if( fraction < 0 )
    {
        if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
            fraction = s.getParent()->area() ;
        else
            fraction = s.getParent()->volume();
        dt = s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT() ;
    }

    change = false ;    
    
    
    //special case when there is a local snap/back instability
    if(accelerate)
    {
//         std::cout  <<  state[0] << "->" << std::flush ;
        if(!fractured() && s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > 1.-timetolerance)
        {
            state[0] = std::min(state[0]+timetolerance, 1.) ;
        }
        else
        {
            double timeDelta = std::max((1.-maxscore)*dt, timetolerance) ;
            state[0] = std::min(state[0]+timeDelta*overdamage*accelerate, 1.)  ; 
        }
        
//         std::cout  <<  state[0] << std::endl ;
    }
    dt = s.getParent()->getBoundingPoint(s.getParent()->getBoundingPoints().size()-1).getT() - s.getParent()->getBoundingPoint(0).getT() ;
    
    if(!s.getParent()->getBehaviour()->getFractureCriterion() || 
        !s.getParent()->getBehaviour()->getFractureCriterion()->met() || 
        std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() - maxscore) >= timetolerance 
    )
    {
        accelerate = 0 ;
    }
    else if(!fractured() 
            && (std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() - maxscore) <timetolerance 
                ||  s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > 1.-timetolerance ))
    {
//         std::cout << dynamic_cast<DelaunayTriangle *>(s.getParent())->index << std::endl ;
        double maxAccelerate = 1 ;
        if(std::max(dt, timetolerance) > POINT_TOLERANCE)
            maxAccelerate = 1./std::max(dt, timetolerance) ;
        if(accelerate < 1)
        {
            double denominator = (1.-maxscore) ;
            if(denominator < POINT_TOLERANCE)
                denominator = POINT_TOLERANCE ;
            accelerate = 1./denominator ;
        }
        else
            accelerate *= 2. ;
        accelerate = std::min(maxAccelerate, accelerate) ;
        change = true ;
        s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
    }

}


DamageModel * SpaceTimeIsotropicLinearDamage::getCopy() const
{
    SpaceTimeIsotropicLinearDamage * dam = new SpaceTimeIsotropicLinearDamage(overdamage, thresholdDamageDensity) ;
    dam->setResidualStiffnessFraction( residualStiffnessFraction ) ;
    return dam ;
}

}

