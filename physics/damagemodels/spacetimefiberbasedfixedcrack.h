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

class SpaceTimeFiberBasedFixedCrack : public SpaceTimeFiberBasedIsotropicLinearDamage
{
protected:
    double orientation ;
    bool fixedOrientation ;
    std::map<std::string, double> values ;

public:
    SpaceTimeFiberBasedFixedCrack(double f = 0.1, double tol = 0.001, double cutoff = 0.6, bool externalOrientation = false) : SpaceTimeFiberBasedIsotropicLinearDamage(f,tol,cutoff), orientation(0.), fixedOrientation(externalOrientation) { getState(true).resize(2,0.) ; }

    virtual ~SpaceTimeFiberBasedFixedCrack() { } 

    virtual void step(ElementState & s, double maxscore)  ;

    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual DamageModel * getCopy() const ;
};

}

#endif
