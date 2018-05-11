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
#ifndef LINEARCONTACT_H
#define LINEARCONTACT_H

#include "contactmodel.h"

namespace Amie {


class LinearContactForce final: public ContactModel
{
public:
    Geometry * geo ;   
    double stiffness ;
    double tangentStiffness ;
    Vector forces ;
    Vector deltaForce ;
    Vector tangentForces ;
    Vector tangentDeltaForce ;
    bool reset = false ;
    
    ElementState * es = nullptr;
public:

    LinearContactForce(Geometry  *geo, double stiffness = 500e9, double tangentStiffness = 50e9) ;

    virtual ~LinearContactForce();

//     virtual void step( ElementState &s , double maxscore) ;

    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;


    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual bool fractured(int direction = -1) const {return false ;};

    virtual DamageModel * getCopy() const ;

    virtual void computeDelta(ElementState & s) ;
    
    virtual bool hasInducedBoundaryConditions() const {
        return true ;
    } ;
             
    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
    
    virtual void postProcess() ;
    virtual bool hasInducedForces() const {
        return true ;
    }


};

}

#endif
