//
// C++ Interface: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __WEI_STIFFNESS_WITH_VARIABLE_DEF_AND_FRAC
#define __WEI_STIFFNESS_WITH_VARIABLE_DEF_AND_FRAC

#include "physics_base.h"
#include "damagemodels/isotropiclineardamage.h"
#include "fracturecriteria/fracturecriterion.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "stiffness_with_variable_imposed_deformation_and_fracture.h"

namespace Mu
{

	
	/** \brief A linear Elastic Law with a fracture criterion and imposed strain. When applied to elements, properties are randomised.
	*
	* The field param is the Cauchy-Green Strain Tensor
	* The imposed deformation are given as a vector
	*/
	struct WeibullStiffnessWithVariableImposedDeformationAndFracture : public LinearForm
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
		double variability ;
		double materialRadius ;
		double neighbourhoodRadius ;
		std::vector<Variable> v ;
		/** Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param imposedDef original imposed deformation vector
		* @param criterion FractureCriterion to use
		*/
		WeibullStiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * criterion) ;
		
		virtual ~WeibullStiffnessWithVariableImposedDeformationAndFracture() ;
		
	/** \brief Apply the behaviour
		* This overloaded apply() is more efficient and is designed to minimise allocating and dealocating memory.
		* This method  is never used, however, because the behaviour only serves as a template for copying
		* 
		* @param p_i first shape function.
		* @param p_j second shape function.
		* @param gp Set of gauss points for numerical integration
		* @param Jinv inverse jacobian matrices at the gauss points
		* @param ret matrix in which to sore the results
		* @param vm pointer to the virtual machine dedicated for the computation
		*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		virtual bool fractured() const ;

		virtual bool changed() const ;
		
		/** \brief return a StiffnessWithVariableImposedDeformationAndFracture, which parameters are randomised using a Weibull distribution*/
		virtual Form * getCopy() const ;

		virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;
		
		virtual bool hasInducedForces() const ;
		
		virtual void step(double timestep, ElementState & currentState, double maxscore) ;

		/** \brief returns 0 */
		virtual Vector getPreviousDamage() ;
	
		/** \brief returns 0 */
		virtual Vector getPreviousPreviousDamage() {return Vector(0) ; } ;
		
		virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
	} ;


} ;

#endif
