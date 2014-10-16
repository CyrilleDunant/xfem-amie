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
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

namespace Amie
{

typedef enum
{
    LOGCREEP_CONSTANT,
    LOGCREEP_FORWARD,
    LOGCREEP_PREDICTED,
    LOGCREEP_AGEING,
    LOGCREEP_EXPSTRAIN,
} LogCreepStressAccumulator ;


struct LogarithmicCreep : public Viscoelasticity
{
	// material parameters
	Matrix C ; //stiffness
	Matrix E ; //viscosity
	Matrix Eprev ;
	Matrix Edot ;
	double tau ; //characteristic time
    LogCreepStressAccumulator accumulator ;
	bool fixCreepVariable ;

	double accumulatedStress ;
	double currentStress ;
    double previousTimeStep ;
    double reducedTimeStep ;
	bool sign ;

	Matrix accumulatedStressTensor ;
	Matrix currentStressTensor ;

/*    std::vector<Matrix> E_t ;
    std::vector<Matrix> dE_dt ;*/

	bool isPurelyElastic ;
    bool updated ;
	
	// constructor for pure elasticity
	LogarithmicCreep( const Matrix & rig) ;
	// standard constructor
	LogarithmicCreep( const Matrix & rig, const Matrix & v, double e) ;

	virtual ~LogarithmicCreep() { } ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual Form * getCopy() const ;

	virtual void print() const ;

	virtual void preProcess( double timeStep, ElementState & currentState ) ;
    virtual void step(double timestep, ElementState &s, double maxScore) ;

    virtual void setLogCreepAccumulator( LogCreepStressAccumulator acc) ;

    virtual void accumulateStress(double timeStep, ElementState & currentState) ;
    virtual void makeEquivalentViscosity(double timeStep, ElementState & currentState) ;

	virtual bool hasInducedForces() const { return fixCreepVariable ; }

	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;

struct LogarithmicCreepWithImposedDeformation : public LogarithmicCreep
{
	// size = 0 indicates no imposed deformation
	Vector imposed ;
    Vector prevImposed ;

	LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp ) ;
	LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, double e, const Vector & imp ) ;

	virtual ~LogarithmicCreepWithImposedDeformation() { } ;

	virtual Form * getCopy() const ;

	virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

	virtual bool hasInducedForces() const { return fixCreepVariable || imposed.size()>0 ; }

    virtual void step(double timestep, ElementState &s, double maxScore) ;

	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
} ;

struct LogarithmicCreepWithImposedDeformationAndFracture : public LogarithmicCreepWithImposedDeformation
{

	DamageModel * dfunc ;
	FractureCriterion * criterion ;
	bool noFracture ;

	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double e, const Vector & imp) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double e, const Vector & imp, FractureCriterion * c, DamageModel * d) ;

	virtual ~LogarithmicCreepWithImposedDeformationAndFracture() { if(dfunc){delete dfunc ;} if(criterion){delete criterion ;} } 

	virtual Form * getCopy() const ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

    virtual bool changed() const { return (noFracture ? false : dfunc->changed()) ; }
    virtual bool fractured() const { return (noFracture ? false : dfunc->fractured()) ; }

    virtual void step(double timestep, ElementState &s, double maxScore) ;

    virtual FractureCriterion * getFractureCriterion() const { return criterion ; }
    virtual DamageModel * getDamageModel() const { return dfunc ;}
    virtual void setFractureCriterion(FractureCriterion *frac) ;

    virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
    virtual Matrix getViscousTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;


} ;


#endif
