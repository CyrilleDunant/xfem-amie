//
// C++ Interface: lineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
// WARNING: these damage models were implemented for testing purposes only!

#ifndef MU_SPACETIME_BAD_ISOTROPIC_DAMAGE_H
#define MU_SPACETIME_BAD_ISOTROPIC_DAMAGE_H

#include "damagemodel.h"
#include "spacetimefiberbasedisotropiclineardamage.h"

namespace Amie {

/*PARSE SpaceTimeFixedPointIsotropic DamageModel
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[maximum_damage] 0.6 // damage above which an element is considered broken
*/
class SpaceTimeFixedPointIsotropicLinearDamage : public SpaceTimeFiberBasedIsotropicLinearDamage
{
    double lastDamage = 0 ;
    double lastState  = 0 ;

public:
    SpaceTimeFixedPointIsotropicLinearDamage(double f = 0.1, double cutoff = 0.6) : SpaceTimeFiberBasedIsotropicLinearDamage(f, 1., cutoff) {  }

    virtual ~SpaceTimeFixedPointIsotropicLinearDamage() { }

    virtual void step(ElementState & s, double maxscore)  ;
    virtual DamageModel * getCopy() const ;
    virtual void prepare() { lastDamage = state[0] ; }

} ;

/*PARSE SpaceTimeSequentialIsotropic DamageModel
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 1e-9 // minimum time between two successive damage events
    @value[maximum_damage] 0.6 // damage above which an element is considered broken
*/
class SpaceTimeSequentialIsotropicLinearDamage : public SpaceTimeFiberBasedIsotropicLinearDamage
{

public:
    SpaceTimeSequentialIsotropicLinearDamage(double f = 0.1, double tol = 1e-9, double cutoff = 0.6) : SpaceTimeFiberBasedIsotropicLinearDamage(f, tol, cutoff) {  }

    virtual ~SpaceTimeSequentialIsotropicLinearDamage() { }

    virtual void step(ElementState & s, double maxscore)  ;
    virtual DamageModel * getCopy() const ;

} ;


}

#endif
