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

#ifndef __GENERALIZED_FD_MAXWELL_H_
#define __GENERALIZED_FD_MAXWELL_H_

#include "physics_base.h"
#include "stiffness.h"



namespace Mu
{
	struct MaxwellBranch
	{
		Matrix stiffness ;
		Vector decay ;
		int index ;
		int p ;
		double gamma ;
		Function affine, constant ;
		
		std::valarray<Matrix> reducedStiffnessAtGaussPoints ;
		std::valarray<Vector> imposedStressAtGaussPoints ;
		std::valarray<Vector> fi, gi, li, pi ;
		
		MaxwellBranch(const Matrix & rig, Vector d, int i, Function affine, Function constant, int p = 10, double g = 0.5) ;
		
		void step(double timestep, ElementState & currentState) ;
		void stepAtGaussPoint(double timestep, ElementState & currentState, int j) ;
		
		void nextDelta(int k, double timestep, Vector & d0, Vector & d1, Vector & d2, Vector & d3, Vector & prev) ;
		
		Vector updateInternalStrain( size_t g, const Vector & eps) const ;
		Vector updateInternalStrainRate( size_t g, const Vector & eps) const ;
	} ;
      
  
	struct GeneralizedFDMaxwell : public Stiffness
	{
		std::vector<MaxwellBranch> branches ;
		bool change ;
		
		GeneralizedFDMaxwell(const Matrix & rig, const std::vector<std::pair<Matrix, Vector> > & br , int p = 10, double g = 0.5) ;
		GeneralizedFDMaxwell(const Matrix & rig, const std::pair<Matrix, Vector> & br , int p = 10, double g = 0.5) ;

		virtual ~GeneralizedFDMaxwell() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const { return false ; } 
		
		virtual Form * getCopy() const ;

		virtual bool hasInducedForces() const { return true ; }

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = NULL) const ;

		std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
		virtual void step(double timestep, ElementState & currentState) ;
		
		std::vector<Matrix> makeStiffnessMatrixAtGaussPoints(size_t p) const ;
		std::vector<Vector> makeImposedStressAtGaussPoints(size_t p) const ;
		
		virtual Vector updateInternalStrain( size_t g, size_t i, Vector & eps) ;
		virtual Vector updateInternalStrainRate( size_t g, size_t i, Vector & eps) ;

		virtual ElementState * createElementState( IntegrableEntity * e) const ;

		virtual void updateElementState(double timestep, ElementState & currentState) const ;

		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = NULL) const ;
		
		virtual void preProcess( double timeStep, ElementState & currentState ) ;
	} ;

} ;

#endif