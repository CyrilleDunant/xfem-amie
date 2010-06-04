// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STRESS_DEFINED_STIFF_DEF_H_
#define __STRESS_DEFINED_STIFF_DEF_H_

#include "physics_base.h"
#include "non_linear_stiffness.h"
#include "stress_defined_stiffness.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"
#include "../polynomial/vm_base.h" 


namespace Mu
{



struct StressDefinedStiffnessWithImposedDeformation : public StressDefinedStiffness
{
	Function imposedX ;
	Function imposedY ;
	Function imposedZ ;

	StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defIso, IntegrableEntity * parent) ;
	StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defIso, SpaceDimensionality dim) ;
	StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defX, Function defY, Function defZ, IntegrableEntity * parent) ;
	StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defX, Function defY, Function defZ, SpaceDimensionality dim) ;

	virtual ~StressDefinedStiffnessWithImposedDeformation() ;
	
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const;

	virtual bool hasInducedForces() const ;
		
	virtual Vector getImposedStress(const Point & p) const ;
	
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
	std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	virtual bool isActive() const ;
	
	virtual Form * getCopy() const ;
	

} ;


} ;

#endif // __ PHYSICS_H_



