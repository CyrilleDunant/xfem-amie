// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STRESS_DEFINED_STIFF_H_
#define __STRESS_DEFINED_STIFF_H_

#include "physics_base.h"
#include "non_linear_stiffness.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"

namespace Mu
{



struct StressDefinedStiffness : public NonLinearStiffness
{
	StressDefinedStiffness(Function f, double n, IntegrableEntity * parent) ;
	StressDefinedStiffness(Function f, double n, SpaceDimensionality dim) ;

	virtual ~StressDefinedStiffness() ;
	
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const;
	
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
	virtual bool isActive() const ;
	std::vector<BoundaryCondition * > StressDefinedStiffness::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	virtual Form * getCopy() const ;
	
	virtual bool fractured() const ;

} ;




} ;

#endif // __ PHYSICS_H_



