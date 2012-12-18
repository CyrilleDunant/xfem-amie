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

#ifndef __MAXWELL_AND_FRACTURE_H_
#define __MAXWELL_AND_FRACTURE_H_

#include "physics_base.h"
#include "stiffness.h"
#include "maxwell.h"



namespace Mu
{	
	struct GeneralizedIterativeMaxwellAndFracture : public LinearForm
	{
		std::vector<Variable> v ;
		std::vector<Vector> imposedStressAtGaussPoints ;
		std::vector<IterativeMaxwell *> branches ;
		Matrix r0 ;
		DamageModel * dfunc ;
		FractureCriterion * criterion ;		
		
		GeneralizedIterativeMaxwellAndFracture(const Matrix & rig, const std::vector<Matrix> & param, const std::vector<double> & chartime, FractureCriterion * crit, DamageModel * dmg = nullptr) ;
		virtual ~GeneralizedIterativeMaxwellAndFracture() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		virtual void step(double timestep, ElementState & currentState, double maxscore) ;
		virtual void updateElementState(double timestep, ElementState & currentState) const ;
		virtual void preProcess( double timeStep, ElementState & currentState ) ;
 		void preProcessAtGaussPoint(double timestep, ElementState & currentState, int g) ;	  

		
		virtual FractureCriterion * getFractureCriterion() const ;
		virtual DamageModel * getDamageModel() const ;
		virtual void setFractureCriterion(FractureCriterion * frac) ;
		virtual bool changed() const ;
		virtual bool fractured() const ;		
		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		
		virtual Form * getCopy() const ;
		virtual ElementState * createElementState( IntegrableEntity * e) ;
		virtual bool hasInducedForces() const { return true ; }

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

		void setNumberOfGaussPoints(size_t n) ;
		void syncNumberOfGaussPoints(ElementState & state) ;
		virtual void getCoefficients(double timestep) ;
		void getInstantaneousCoefficients() ;
	} ;
	

} ;

#endif