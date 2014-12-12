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

	virtual ~LogarithmicCreep() { } ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	virtual void applyViscous( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual Form * getCopy() const ;

	virtual void print() const ;

	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual void step(double timestep, ElementState &s, double maxScore) ;

	virtual bool hasInducedForces() const { return fixCreepVariable ; }

	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

} ;

struct LogarithmicCreepWithImposedDeformation : public LogarithmicCreep
{
	// size = 0 indicates no imposed deformation
	Vector imposed ;
	Vector prevImposed ;

	LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc = new LogCreepAccumulator() ) ;
	LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator() ) ;

	virtual ~LogarithmicCreepWithImposedDeformation() { } ;

	virtual Form * getCopy() const ;

	virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

	virtual bool hasInducedForces() const { return fixCreepVariable || imposed.size()>0 ; }

	virtual void preProcess( double timeStep, ElementState & currentState ) ;

	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
} ;

struct LogarithmicCreepWithImposedDeformationAndFracture : public LogarithmicCreepWithImposedDeformation
{

	DamageModel * dfunc ;
	FractureCriterion * criterion ;
	bool noFracture ;

	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc = new LogCreepAccumulator()) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator()) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc = new LogCreepAccumulator()) ;
	LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc = new TimeUnderLoadLogCreepAccumulator()) ;

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
