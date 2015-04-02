// C++ Interface: generalized_spacetime_viscoelasticity
//
// Description: Generalized visco-elastic behaviour 
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __FINITE_DIFFERENCE_VISCOELASTICTY_H_
#define __FINITE_DIFFERENCE_VISCOELASTICTY_H_

#include "viscoelasticity.h"

namespace Amie
{

typedef enum
{
	FORWARD_EULER,
	BACKWARD_EULER,
	CENTRAL_DIFFERENCE,
	NEWMARK,
	ZIENKIEWICZ,
} ViscoelasticFiniteDifferenceIntegration ;

struct FiniteDifferenceViscoelasticity : public LinearForm
{
	// viscoelastic model represented
	ViscoelasticModel model ;

	// list of tensors for easy access
	std::vector<Matrix> tensors ;

        Vector viscoelasticInternalForces ;

        ViscoelasticFiniteDifferenceIntegration scheme ;

        double theta ;
	
	std::vector<Variable> v ;
	
	// constructor for pure elasticity or pure viscosity
	FiniteDifferenceViscoelasticity( ViscoelasticModel model, const Matrix & rig, ViscoelasticFiniteDifferenceIntegration s = BACKWARD_EULER, double t = 0) ; 
	// constructor for elementary Kelvin-Voigt or Maxwell
	FiniteDifferenceViscoelasticity( ViscoelasticModel model, const Matrix & rig, const Matrix & eta, ViscoelasticFiniteDifferenceIntegration s = BACKWARD_EULER, double t = 0) ; 
	// constructor for Burger (KV + Maxwell in serial)
	FiniteDifferenceViscoelasticity( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx, ViscoelasticFiniteDifferenceIntegration s = BACKWARD_EULER, double t = 0) ; 
	// constructor for generalized KelvinVoigt or Maxwell
	FiniteDifferenceViscoelasticity( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches, ViscoelasticFiniteDifferenceIntegration s = BACKWARD_EULER, double t = 0) ;
	// constructor for generalized KelvinVoigt or Maxwell with 1 module only
	FiniteDifferenceViscoelasticity( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, ViscoelasticFiniteDifferenceIntegration s = BACKWARD_EULER, double t = 0) ;

	virtual ~FiniteDifferenceViscoelasticity() ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;

	virtual bool fractured() const { return false ; }
	virtual bool isViscous() const { return false ; }
        virtual bool hasInducedForces() const { return true ; }

	virtual ElementState * createElementState( IntegrableEntity * e) ;

	virtual void preProcess( double timeStep, ElementState & currentState ) ;

	virtual void updateElementState(double timestep, ElementState & currentState) const ;

	virtual Form * getCopy() const ;

	virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
	virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

	virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const { return param ; }	
	
        virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

	virtual void print() const ;

protected:
	void initializeThetaCoefficient() ;
  
} ;


} 


#endif
