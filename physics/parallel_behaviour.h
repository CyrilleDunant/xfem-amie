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

#ifndef __PARALLEL_BEHAVIOUR_H_
#define __PARALLEL_BEHAVIOUR_H_

#include "physics_base.h"



namespace Amie
{
	struct ParallelBehaviour : public LinearForm
	{
		std::vector<Form * > branches ;
		std::vector<Variable> v ;

		ParallelBehaviour(  std::vector<Form *> b ) ;
		ParallelBehaviour(  Form * b1, Form * b2 ) ;

		virtual ~ParallelBehaviour() ;

		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const ;
		
		virtual Form * getCopy() const ;

		virtual bool hasInducedForces() const ;

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
		virtual void step(double timestep, ElementState & currentState, double maxscore) ;
		
		virtual ElementState * createElementState( IntegrableEntity * e) ;

		virtual void updateElementState(double timestep, ElementState & currentState) const ;

		virtual Matrix getTensor(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		
		virtual void preProcess( double timeStep, ElementState & currentState ) ;
		
	} ;

} ;

#endif