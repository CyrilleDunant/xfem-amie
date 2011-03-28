//
// C++ Interface: weibull_distributed_stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __SPATIALLY_DISTRIBUTED_STIFFNESS_H_
#define __SPATIALLY_DISTRIBUTED_STIFFNESS_H_

#include "physics_base.h"
#include "../features/features.h"
#include "../geometry/geometry_base.h"

namespace Mu
{

	/** \brief A linear Elastic Law (without any fracture model)
	* The field param is the Cauchy-Green Strain Tensor
	* Actual tensors are distributed around the prescribed value following a Weibull distribution
	* Stiffness is decreased around the particles
	*/
	struct SpatiallyDistributedStiffness : public LinearForm
	{
		std::vector<Variable> v ;
		double variability ;
		Matrix pore ;
		double length ;
		double distance ;
		double criteriona ;
		double criterionb ;

		/** Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param inc inclusion vector
		*/
		SpatiallyDistributedStiffness(const Matrix & rig, const Matrix & pore, double l, double ca, double cb)  ;
		
		virtual ~SpatiallyDistributedStiffness();
		
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

		/** \brief return a Stifness Behaviour
		*
		* The parameters for the copy are determined randomly using a Weibull distribution. As this behaviour is spaceDependant, only the copies are used in elements
		*/
		virtual Form * getCopy() const ;

		void setDistance(double d) ;
	} ;

} ;

#endif
