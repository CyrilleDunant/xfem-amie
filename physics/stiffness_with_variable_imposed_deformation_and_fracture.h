//
// C++ Interface: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_WITH_VARIABLE_DEF_AND_FRAC
#define __STIFFNESS_WITH_VARIABLE_DEF_AND_FRAC

#include "physics_base.h"
#include "damagemodels/isotropiclineardamage.h"
#include "fracturecriteria/fracturecriterion.h"

namespace Mu
{

	
	/** \brief A linear Elastic Law
	*
	* The field param is the Cauchy-Green Strain Tensor
	* The imposed deformation are given as a vector and grow with time
	* A fracture criterion can be set
	*/
	struct StiffnessWithVariableImposedDeformationAndFracture : public LinearForm
	{
		Vector imposed ;
		Function x ;
		Function y ;
		IsotropicLinearDamage dfunc ;
		double previousDamage ;
		FractureCriterion * criterion ;
		bool frac ;
		bool change ; 
		double sigmaRupt ;
		double init ;
		double damage ;
		std::vector<Variable> v ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param imposedDef Vector of the original imposed deformations
		* @param criterion FractureCriterion to use.
		*/
		StiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * criterion) ;
		
		virtual ~StiffnessWithVariableImposedDeformationAndFracture() ;
		
		/** \brief  Apply the law.
		* 
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @return symbolic matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
		
		/** \brief  Apply the law.
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @param gp integration points
		* @param Jinv inverse Jacobian matrix at the integration points
		* @param ret, matrix to store the result
		* @param vm VirtualMachine to use to perform the computation
		* @return matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		/** \brief return true is the damage model has reached the breaking point.*/
		virtual bool fractured() const ;

		/** \brief return true if the damage state has changed*/
		virtual bool changed() const ;
		
		/** \brief return a copy of the bahaviour, including a copy of the fracture criterion and damage state.*/
		virtual Form * getCopy() const ;
		
		/** \brief return true */
		virtual bool hasInducedForces() const ;

		/** \brief return the imposed strain vector * the CG tensor. */
		virtual Vector getImposedStress(const Point & p) const ;
		
		/** \brief Increment the loading. The loading is increased by the timestep*a random value between 0 and 1 .*/
		virtual void step(double timestep, ElementState & currentState) ;

		virtual void artificialDamageStep(double d) ;
		
		/** \brief Return the virtual force resulting of the imposed stress
		* 
		* @param ElementState to consider 
		* @param p_i shape function to consider
		* @param gp GaussPointArray to consider for the quadrature
		* @param Jinv  inverse Jacobian matrix at the integration points
		* @param Vector v Vector in which to store the results
		* @return \f$ \int\nabla^T p_i K \varepsilon_{imposed} \f$ 
		*/
		virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
		
	} ;


} ;

#endif
