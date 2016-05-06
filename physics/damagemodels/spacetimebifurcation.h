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
#ifndef MU_SPACETIME_BIFURCATION_H
#define MU_SPACETIME_BIFURCATION_H

#include "damagemodel.h"
#include "spacetimefiberbasedisotropiclineardamage.h"

namespace Amie {

class SpaceTimeBifurcation : public SpaceTimeFiberBasedIsotropicLinearDamage
{
protected:
    Vector residualStrain ;
    Vector residualStress ;
    double stiffnessFactor ;

public:
    SpaceTimeBifurcation(double f = 0.1, double tol = 0.001) : SpaceTimeFiberBasedIsotropicLinearDamage(1., tol, 1.), stiffnessFactor(f) {  }

    virtual ~SpaceTimeBifurcation() { }

    virtual void step(ElementState & s, double maxscore)  ;
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual bool fractured(int direction = -1) const  { return ( direction == 0 ? state.max() > POINT_TOLERANCE : false ) ; }

    virtual DamageModel * getCopy() const ;

    virtual bool hasInducedBoundaryConditions() const { return true ; }
    virtual bool hasInducedForces() const { return true ; }

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
    virtual Vector getImposedStress(const Point & p) const ;
    virtual Vector getImposedStrain(const Point & p) const ;

};

class SpaceTimeBifurcationAndDamage : public SpaceTimeBifurcation
{
protected:
    DamageModel * secondary ;

public:
    SpaceTimeBifurcationAndDamage( double f, DamageModel * dam, double tol = 0.001 ) : SpaceTimeBifurcation(f,tol), secondary(dam) { state.resize(0., dam->getState().size()+1) ; }

    virtual ~SpaceTimeBifurcationAndDamage() { if(secondary) { delete secondary ; } }
    virtual DamageModel * getCopy() const ;

    virtual bool fractured(int direction = -1) const ;
    virtual Vector getImposedStress(const Point & p) const ;
    virtual Vector getImposedStrain(const Point & p) const ;

    virtual void step(ElementState & s, double maxscore)  ;
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const { return this->apply(m,p,e,g) ; }



};



}

#endif
