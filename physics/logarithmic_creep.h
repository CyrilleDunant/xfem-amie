// C++ Interface: logarithmic_creep
//
// Description: Logarithmic visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LOGARITHMIC_CREEP_H_
#define __LOGARITHMIC_CREEP_H_

#include "physics_base.h"
#include "viscoelasticity.h"
#include "material_laws/logcreep_accumulator.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

namespace Amie
{

/*PARSE LogarithmicCreep Form 
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @value[creep_modulus] -1 // viscosity of the dashpot
    @value[creep_poisson] -1 // Poisson ratio of the dashpot
    @value[creep_characteristic_time] -1 // creep rate in the log scale
    @value[recoverable_modulus] -1 // modulus of the Kelvin_Voigt unit
    @value[recoverable_poisson] -1 // Poisson ratio of the Kelvin-Voigt unit
    @object<LogCreepAccumulator>[accumulator] RealTimeLogCreepAccumulator() // how the logrithmic creep is obtained
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
 */
struct LogarithmicCreep : public Viscoelasticity
{
    // material parameters
    Matrix C ; //stiffness
    Matrix E ; //stiffness of Maxwell spring
    Matrix R ; //stiffness of Kelvin-Voigt spring
    double tau ; //viscosity of material
    double reducedTimeStep ;
    LogCreepAccumulator * accumulator ;
    bool isPurelyElastic ;
    bool updated ;
    bool timeDependentIntegration ;
    bool fixCreepVariable ;
    Matrix prevParam ;
    Matrix prevEta ;

    // constructor for pure elasticity
    LogarithmicCreep( const Matrix & rig, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator()) ;
    // standard constructor
    LogarithmicCreep( const Matrix & rig, const Matrix & v, const Matrix & r, double t, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator()) ;

    LogarithmicCreep(double young, double poisson, double E_creep = -1, double nu_creep = -1, double tau = -1, double E_recovery = -1, double nu_recovery = -1, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS) ;

    virtual ~LogarithmicCreep() {
        if(accumulator) {
            delete accumulator ;
        }
    } ;

    virtual Matrix getElasticTensor(const Point & p) const { return C ; }

    virtual Form * getCopy() const ;

    virtual void print() const ;

    virtual void preProcess( double timeStep, ElementState & currentState ) ;
    virtual void step(double timestep, ElementState &s, double maxScore) ;

    virtual bool hasInducedForces() const {
        return fixCreepVariable ;
    }

    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;

struct LogarithmicCreepWithImposedDeformation : public LogarithmicCreep
{
    // size = 0 indicates no imposed deformation
    Vector imposed ;
    Vector prevImposed ;

    LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc = new LogCreepAccumulator() ) ;
    LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator() ) ;

    LogarithmicCreepWithImposedDeformation(double young, double poisson, double alpha, double E_creep = -1, double nu_creep = -1, double tau = -1, double E_recovery = -1, double nu_recovery = -1, LogCreepAccumulator * acc = new RealTimeLogCreepAccumulator(), SpaceDimensionality = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS) ;

    virtual ~LogarithmicCreepWithImposedDeformation() { } ;

    virtual Form * getCopy() const ;

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual bool hasInducedForces() const {
        return fixCreepVariable || imposed.size()>0 ;
    }

    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual void preProcess( double timeStep, ElementState & currentState ) ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
} ;

struct LogarithmicCreepWithImposedDeformationAndFracture : public LogarithmicCreepWithImposedDeformation
{

    DamageModel * dfunc ;
    FractureCriterion * criterion ;
    bool noFracture ;
    bool etaSet ;
    bool paramSet ;
    Matrix currentEta ;
    Matrix currentParam ;

    LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc = new LogCreepAccumulator()) ;
    LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator()) ;
    LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc = new LogCreepAccumulator()) ;
    LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator()) ;

    virtual ~LogarithmicCreepWithImposedDeformationAndFracture() ; //{ if(dfunc){delete dfunc ;} if(criterion){delete criterion ;} }

    virtual Form * getCopy() const ;

    virtual bool changed() const {
        return (noFracture ? false : dfunc->changed()) ;
    }
    virtual bool fractured() const {
        return (noFracture ? false : dfunc->fractured()) ;
    }

    virtual void step(double timestep, ElementState &s, double maxScore) ;

    virtual FractureCriterion * getFractureCriterion() const {
        return criterion ;
    }
    virtual DamageModel * getDamageModel() const {
        return dfunc ;
    }
    virtual void setFractureCriterion(FractureCriterion *frac) ;

    virtual void preProcess( double timeStep, ElementState & currentState ) ;

    virtual bool hasInducedForces() const {
        return fixCreepVariable || imposed.size()>0 || (dfunc && dfunc->hasInducedForces() );
    }

    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;


}


#endif
