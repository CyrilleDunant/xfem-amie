//
// C++ Interface: generalized maxwell model with finite difference time-stepping
//
// Description: 
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2012
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MAXWELL_H_
#define __MAXWELL_H_

#include "physics_base.h"
#include "stiffness.h"



namespace Mu
{
	struct NewmarkNumeroffMaxwell : public LinearForm
	{
		Matrix stiffness ;
		Vector decay ;
		int p ;
		double gamma ;
		Function affine, constant ;
		std::vector<Variable> v ;
		
		std::vector<Matrix> reducedStiffnessAtGaussPoints ;
		std::vector<Vector> imposedStressAtGaussPoints ;
		std::vector<Vector> fi, gi, li, pi ;
		
		NewmarkNumeroffMaxwell(const Matrix & rig, Vector d, int p = 10, double g = 0.5) ;
		NewmarkNumeroffMaxwell(const Matrix & rig, Vector d, Function affine, Function constant, int p = 10, double g = 0.5) ;
		virtual ~NewmarkNumeroffMaxwell() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const { return false ; } 
		
		virtual Form * getCopy() const ;

		virtual bool hasInducedForces() const { return true ; }

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = NULL, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = NULL, int g = -1) const ;

		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
		virtual void step(double timestep, ElementState & currentState) ;
		
		virtual ElementState * createElementState( IntegrableEntity * e) ;

		virtual void updateElementState(double timestep, ElementState & currentState) const ;

		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = NULL, int g = -1) const ;
		
		virtual void preProcess( double timeStep, ElementState & currentState ) ;

		
	protected:
		void preProcessAtGaussPoint(double timestep, ElementState & currentState, int g) ;	  
		void nextDelta(int k, double timestep, Vector & d0, Vector & d1, Vector & d2, Vector & d3, Vector & prev) ;
		Vector updateInternalStrain( size_t g, const Vector & eps) const ;
		Vector updateInternalStrainRate( size_t g, const Vector & eps) const ;
		void setNumberOfGaussPoints(size_t n) ;
	} ;
      
  

} ;

#endif