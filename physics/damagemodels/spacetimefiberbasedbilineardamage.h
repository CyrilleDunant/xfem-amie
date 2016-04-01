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
#ifndef MU_SPACETIME_BILINEAR_DAMAGE_H
#define MU_SPACETIME_BILINEAR_DAMAGE_H

#include "damagemodel.h"
#include "spacetimefiberbasedisotropiclineardamage.h"

namespace Amie {

/*PARSE SpaceTimeFiberBasedBilateral DamageModel
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 0.001 // minimum time between two successive damage events
    @value[maximum_damage] 0.6 // damage above which an element is considered broken
    @value[secondary_maximum_damage] -1 // damage above which an element is considered broken 
    @string<IsotropicMaterialParameters> BULK_SHEAR // independent material parameters affected
    @string<planeType> PLANE_STRESS // 2D planar assumption
*/
class SpaceTimeFiberBasedBilateralLinearDamage final: public SpaceTimeFiberBasedIsotropicLinearDamage
{
    planeType pt ;
    IsotropicMaterialParameters mode ;

public:
    SpaceTimeFiberBasedBilateralLinearDamage(double f = 0.1, double tol = 0.001, double cutoff = 0.6, double shearCut = -1, IsotropicMaterialParameters m = BULK_SHEAR, planeType pt_ = PLANE_STRESS) : SpaceTimeFiberBasedIsotropicLinearDamage(f, tol, cutoff), pt(pt_), mode(m)
    { 
        if(shearCut > 0)
            secondaryThresholdDamageDensity = shearCut ;
        else
            secondaryThresholdDamageDensity = cutoff ;
        getState(true).resize(2, 0.);
    }

    virtual ~SpaceTimeFiberBasedBilateralLinearDamage() { }

    virtual void step(ElementState & s, double maxscore)  ;
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual bool fractured(int direction = -1) const  ;

    virtual DamageModel * getCopy() const ;
};


}

#endif
