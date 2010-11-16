//
// C++ Interface: weibull_distributed_stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __WEIBULL_STIFFNESS_H_
#define __WEIBULL_STIFFNESS_H_

#include "physics_base.h"
#include "damagemodels/damagemodel.h"
#include "fracturecriteria/fracturecriterion.h"

namespace Mu
{

	/** \brief A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	* Actual tensors are distributed around the prescribed value following a Weibull distribution
	*/
	struct WeibullDistributedStiffness : public LinearForm
	{
		double materialRadius ;
		double neighbourhoodRadius ;
		std::vector<Variable> v ;
		double variability ;
		double down ;
		double up ;
		MirrorState mirroring; 
		double dx ;
		double dy ;
		double dz ;
		DamageModel * damageModel ;
		/** Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param cri stress limit for the Mohr - Coulomb criterion to use
		*/
		WeibullDistributedStiffness(const Matrix & rig, double down, double up, MirrorState mirroring = NO_MIRROR, double dx = 0, double dy = 0, double dz = 0)  ;
		
		void setMaterialCharacteristicRadius(double r) {materialRadius = r ;} ;
		void setNeighbourhoodRadius(double r) {neighbourhoodRadius = r ;}
		
		virtual ~WeibullDistributedStiffness();
		
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

		/** \brief return a StifnessAndFracture Behaviour
		*
		* The parameters for the copy are determined randomly using a Weibull distribution. As this behaviour is spaceDependant, only the copies are used in elements
		*/
		virtual Form * getCopy() const ;
		
		void setDamageModel(DamageModel * newmodel) ;

		
	} ;

} ;

#endif
