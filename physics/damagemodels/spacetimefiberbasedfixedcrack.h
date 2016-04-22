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
#ifndef MU_SPACETIME_FIBER_FIXEDCRACK_H
#define MU_SPACETIME_FIBER_FIXEDCRACK_H

#include "damagemodel.h"
#include "spacetimefiberbasedisotropiclineardamage.h"

namespace Amie {

/*PARSE SpaceTimeFiberBasedFixedCrack DamageModel --no-suffix
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 1e-5 // time lapsed between two damage events
    @value[maximum_damage] 0.999 // damage above which an element is considered broken
    @string<FieldType>[orientation_field] REAL_STRESS_FIELD // field used to get the orientation of the crack
    @string<bool>[external_orientation] false // find the orientation from the ElementState
*/
class SpaceTimeFiberBasedFixedCrack : public SpaceTimeFiberBasedIsotropicLinearDamage
{
protected:
    double orientation ;
    FieldType orientationField ;
    bool fixedOrientation ;
    std::map<std::string, double> values ;

public:
    SpaceTimeFiberBasedFixedCrack(double f = 0.1, double tol = 0.001, double cutoff = 0.6, FieldType type = REAL_STRESS_FIELD, bool externalOrientation = false) : SpaceTimeFiberBasedIsotropicLinearDamage(f,tol,cutoff), orientation(0.), orientationField(type), fixedOrientation(externalOrientation) { getState(true).resize(2,0.) ; }

    virtual ~SpaceTimeFiberBasedFixedCrack() { } 

    virtual void step(ElementState & s, double maxscore)  ;

    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual DamageModel * getCopy() const ;

    virtual bool fractured(int direction = -1) const  ;
};

}

#endif
