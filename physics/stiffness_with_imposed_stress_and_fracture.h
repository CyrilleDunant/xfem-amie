//
// C++ Interface: stiffness_with_imposed_stress_and_fracture
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_WITH_STRESS_FRACTURE
#define __STIFFNESS_WITH_STRESS_FRACTURE

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/lineardamage.h"
#include "damagemodels/anisotropicdamage.h"
#include "damagemodels/isotropiclineardamage.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"

namespace Amie
{
struct StiffnessWithImposedStressAndFracture : public LinearForm
{
    std::vector<Variable> v ;
    Vector imposed ;
    DamageModel * dfunc ;
    FractureCriterion * criterion ;

    StiffnessWithImposedStressAndFracture(const Matrix & rig, Vector imposedStress, FractureCriterion * c, DamageModel * d = nullptr) ;
    virtual ~StiffnessWithImposedStressAndFracture() ;

    virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

    virtual bool fractured() const ;
    virtual bool changed() const ;
    virtual bool hasInducedForces() const {
        return true ;
    }

    virtual Form * getCopy() const ;

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e, int g = -1) const ;

    std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

    virtual void step(double timestep, ElementState & currentState, double maxScore) ;

    virtual FractureCriterion * getFractureCriterion() const ;
    virtual DamageModel * getDamageModel() const ;

    virtual void resetImposedStress() ;
    virtual void setImposedStress(Vector beta) ;
    virtual void addImposedStress(Vector beta) ;

    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

} ;


}

#endif
