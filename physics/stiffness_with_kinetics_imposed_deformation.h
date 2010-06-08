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

#ifndef __STIFFNESS_WITH_TIMED_DEF
#define __STIFFNESS_WITH_TIMED_DEF

#include "physics_base.h"
#include "stiffness_with_imposed_deformation.h"

namespace Mu
{

	
	/** \brief A linear elastic law with an imposed strain
	*
	* The field param is the Cauchy-Green Strain Tensor
	* The imposed deformation are given as a vector
	*/
	struct StiffnessWithKineticsImposedDeformation : public StiffnessWithImposedDeformation
	{
		std::vector<Variable> v ;
		Vector imposed ;
		Function kinetics ;
		Function x ;
		Function y ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		* @param imposedDef Imposed deformation
		*/
		StiffnessWithKineticsImposedDeformation(const Matrix & rig, Vector imposedDef, Function k) ;
		
		virtual ~StiffnessWithKineticsImposedDeformation() ;
		
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
		
		/** \brief return false */
		virtual bool fractured() const ;
		
		/** \brief Return a copy of this Behaviour*/
		virtual Form * getCopy() const ;
	
		
		/** \brief Return the Vector of imposed Stress at the considered point. As the imposed stress is uniform, the point is ignored*/
		virtual Vector getImposedStress(const Point & p) const ;

		/** \brief Increment the according to the function kinetics. The loading is increased by the timestep*a random value between 0 and 1 .*/
		virtual void step(double timestep, ElementState & currentState) ;
		
		std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
	} ;


} ;

#endif
