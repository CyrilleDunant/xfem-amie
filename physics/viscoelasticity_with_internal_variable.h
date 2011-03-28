// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VISCOELASTICITY_WITH_INTERNAL_VARIABLE_H_
#define __VISCOELASTICITY_WITH_INTERNAL_VARIABLE_H_

#include "physics_base.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"

namespace Mu
{


struct ViscoElasticity: public LinearForm
{
	Vector tau_g ;
	Vector tau_k ;
	Vector g ;
	Vector k ;
	Matrix a_g ;
	Vector a_k ;
	Vector average_delta_sigma ;

	ViscoElasticity( double _tau_k, double _tau_g, Vector g, Vector k);
	
	virtual ~ViscoElasticity() ;
		
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const;
	/** \todo remove usage of previousState. complement state instead*/
	virtual void step(double timestep, ElementState & currentState);
	virtual bool changed() {return true ;}
	virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	virtual bool hasInducedForces();
	
	virtual Form * getCopy() const ;
} ;

} ;

#endif // __ PHYSICS_H_



