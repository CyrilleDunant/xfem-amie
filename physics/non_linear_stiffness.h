// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __NON_LINEAR_STIFFNESS_H_
#define __NON_LINEAR_STIFFNESS_H_

#include "physics_base.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"

namespace Mu
{




struct NonLinearStiffness : public NonLinearForm
{
	NonLinearStiffness(Function f, double n, IntegrableEntity * parent) ;
	NonLinearStiffness(Function f, double n, SpaceDimensionality dim) ;

	Function E ;
	double nu ;
	IntegrableEntity * parent ;
	
	virtual ~NonLinearStiffness() ;
	
	/** Apply the law.
	 * 
	 * @param p_i first basis polynomial.
	 * @param p_j second basis polynomial.
	 * @param e IntegrableEntity on which to perform the integration
	 * @return matrix resulting of the integration of \f$ \nabla H^T K \nabla H \f$.
	 */
	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const;
	
	virtual bool hasInducedForces() const;
	
	virtual bool hasInducedMatrix() const ;
	
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const;
	
	
	virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;

	virtual void setParent(IntegrableEntity * p) ; 
	
	virtual bool isActive() const ;
	
	virtual Form * getCopy() const ;
	
} ;

} ;

#endif // __ PHYSICS_H_



