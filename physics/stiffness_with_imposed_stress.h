//
// C++ Interface: stiffness_with_imposed_stress
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_WITH_STRESS
#define __STIFFNESS_WITH_STRESS

#include "physics_base.h"

namespace Amie
{
struct StiffnessWithImposedStress : public LinearForm
{
    std::vector<Variable> v ;
    Vector imposed ;

    StiffnessWithImposedStress(const Matrix & rig) ;
    StiffnessWithImposedStress(const Matrix & rig, Vector imposedStress) ;
    StiffnessWithImposedStress(double E, double nu, double beta, SpaceDimensionality dim) ;
    virtual ~StiffnessWithImposedStress() ;

    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual bool fractured() const ;
    virtual Form * getCopy() const ;
    virtual bool hasInducedForces() const {
        return true ;
    }

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

    virtual void step(double timestep, ElementState & currentState, double maxScore) ;

    virtual void resetImposedStress() ;
    virtual void setImposedStress(Vector imp) ;
    virtual void addImposedStress(Vector imp) ;

} ;

}

#endif
