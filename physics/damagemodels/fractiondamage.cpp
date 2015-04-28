//
// C++ Implementation: lineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fractiondamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

FractionLinearDamage::FractionLinearDamage( Matrix remnant, double phi) : remnant(remnant), phi(phi)
{
    state.resize(2, 0.) ;
    isNull = false ;
    state = 0 ;

    inCompression = false ;
    inTension = false ;
}

std::pair< Vector, Vector > FractionLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
    Vector ret= state;

    inCompression = false ;
    inTension = false ;

    double compressionDamage = 0 ;
    double tensionDamage = 0 ;

    if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0))
    {
        inCompression = true ;
        compressionDamage = 1 ;
        tensionDamage = 1 ;
// 		compressionDamage = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, compressionDamage) ;
// 		compressionDamage = std::min(.99999, compressionDamage) ;
// 		compressionDamage = std::max(0., compressionDamage) ;
    }

    if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0))
    {
        inTension = true ;

        tensionDamage = 1 ;
// 		tensionDamage = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamage) ;
// 		tensionDamage = std::min(.99999, tensionDamage) ;
// 		tensionDamage = std::max(0., tensionDamage) ;
    }
    ret[0] = compressionDamage ;
    ret[1] = tensionDamage ;
    return std::make_pair(state, ret) ;
// 	std::cout << state.sum() << std::flush ;
}

void FractionLinearDamage::computeDelta(ElementState & s)
{
    Vector ret= state;

    inCompression = false ;
    inTension = false ;

    double compressionDamage = 0 ;
    double tensionDamage = 0 ;

    if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0))
    {
        inCompression = true ;
        compressionDamage = 1 ;
        tensionDamage = 1 ;
        // 		compressionDamage = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, compressionDamage) ;
        // 		compressionDamage = std::min(.99999, compressionDamage) ;
        // 		compressionDamage = std::max(0., compressionDamage) ;
    }

    if(s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0))
    {
        inTension = true ;

        tensionDamage = 1 ;
        // 		tensionDamage = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamage) ;
        // 		tensionDamage = std::min(.99999, tensionDamage) ;
        // 		tensionDamage = std::max(0., tensionDamage) ;
    }
    ret[0] = compressionDamage ;
    ret[1] = tensionDamage ;
    delta = (ret-state).max() ;
}

Matrix FractionLinearDamage::apply(const Matrix & m, const Point & p , const IntegrableEntity * e, int g) const
{
    if(fractured())
        return remnant*phi;

    if(inTension)
    {
        return m*(1.-state[0])*(1.-phi)+remnant*phi ;
    }
    if(inCompression)
    {
        return m*(1.-state[1])*(1.-phi)+remnant*phi ;
    }

    return m*(1.-std::max(state[0], state[1]))*(1.-phi)+remnant*phi ;
}


bool FractionLinearDamage::fractured() const
{
    if (fraction < 0)
        return false ;

    return (state[0] >= secondaryThresholdDamageDensity) || (state[1] >= thresholdDamageDensity)  ;
}

FractionLinearDamage::~FractionLinearDamage()
{
}


}
