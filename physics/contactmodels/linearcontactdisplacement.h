//
// C++ Interface: linear contact displacement
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LINEARCONTACTDISP_H
#define LINEARCONTACTDISP_H

#include "contactmodel.h"

namespace Amie {


class LinearContactDisplacement final: public ContactModel
{
public:
    Geometry * geo ;   
    double tangentStiffness ;
    Vector deltaPosition ;
    Vector initialPosition ;
    Vector tangentForces ;
    Vector tangentDeltaForces ;
    
    bool vert = false ;
    
    ElementState * es = nullptr;
public:

    LinearContactDisplacement(Geometry  *geo, double tangentStiffness = 500e9) ;

    virtual ~LinearContactDisplacement();

//     virtual void step( ElementState &s , double maxscore) ;

    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;


    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual bool fractured(int) const {return false ;};

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
