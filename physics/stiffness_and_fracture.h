//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_AND_FRACTURE_H
#define __STIFFNESS_AND_FRACTURE_H

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/lineardamage.h"
#include "damagemodels/anisotropicdamage.h"
#include "damagemodels/isotropiclineardamage.h"
#include "homogenization/homogenization_base.h"

namespace Mu
{

	/** \brief A linear Elastic Law with a pluggable fracture criterion
	*
	* The field param is the Cauchy-Green Strain Tensor
	*/
	struct StiffnessAndFracture : public LinearForm
	{
// 		AnisotropicLinearDamage * dfunc ;
		DamageModel * dfunc ;
// 		IsotropicLinearDamage dfunc ;
		double eps ;
		Vector previousPreviousDamage ;
		Vector intermediateDamage ;
		Vector previousDamage ;
		int count ;
		FractureCriterion * criterion ;
		bool change ; 
		bool previouschange ; 
		double sigmaRupt ;
		double init ;
		Vector damage ;
		std::vector<Variable> v ;

		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param c  FractureCriterion to use. The behaviour is responsible for deleting the criterion upon cleanup.
		*/
		StiffnessAndFracture(const Matrix & rig, FractureCriterion * c, double eps = 0.007)  ;
		void setNeighbourhoodRadius(double d ) ;

		virtual ~StiffnessAndFracture();
		
		/** Apply the law.
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @param gp integration points
		* @param Jinv inverse Jacobian matrix at the integration points
		* @param ret, matrix to store the result
		* @param vm VirtualMachine to use to perform the computation
		* @return matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		/** \brief Check for fracture
		*
		* @param timestep elapsed time
		* @param currentState state of the element
		* 
		* if the yield criterion is true, se fractured state to true
		*/
		virtual void step(double timestep, ElementState & currentState) ;
		
		/** Artificial damage step */
		virtual void artificialDamageStep(double d) ;

		/** \brief returns the previous damage */
		virtual Vector getPreviousDamage() {return previousDamage ; } ;

		/** \brief returns the previous previous damage */
		virtual Vector getPreviousPreviousDamage() {return previousPreviousDamage ; } ;

		/** \brief Impose previous and previous previous damage */
		virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) ;

		/** \brief Check for fracture state
		*
		* @return true if the element is fractured
		*/
		virtual bool fractured() const ;
		
		
		/** \brief get Copy of the behaviour
		*
		* The copy also copies the damage state.
		* @return pointer to the copy. Caller is responsible for cleaning memory
		*/
		virtual Form * getCopy() const ;

		/** \brief return true if the damage state has been modfied*/
		virtual bool changed() const ;
		
		/** \brief return to the previous damage state*/
		virtual void stepBack() ;

		/** \brief Return the (damaged) Stifness tensor*/
		virtual Matrix getTensor(const Point & p) const ;

		/** \brief Acessor, return the stifness criterion in use*/
		virtual FractureCriterion * getFractureCriterion() const ;

		virtual DamageModel * getDamageModel() const ;
				
		virtual Material toMaterial() ;
		
		virtual void scale(double f) ;
		
	} ;
	


} ;


#endif 
