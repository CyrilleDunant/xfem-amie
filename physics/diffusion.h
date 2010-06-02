//
// C++ Interface: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __DIFFUSION_H_
#define __DIFFUSION_H_

#include "physics_base.h"
#include "homogenization/properties_base.h"

namespace Mu
{

	/** \brief A linear Diffusion Law
	* The field param is the diffusion matrix
	*/
	struct Diffusion : public LinearForm
	{
		std::vector<Variable> v ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the diffusion Tensor
		*/
		Diffusion(const Matrix & rig) ;
		
		virtual ~Diffusion() ;
		
		/** \brief Apply the law.
		* 
		* @param p_i first basis function.
		* @param p_j second basis function.
		* @param e Element in which to apply the law.
		* @return matrix resulting of the integration of \f$ \nabla \cdot p_i^T K \nabla \cdot p_j + \dot{p_i}p_j \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const; 
		
	 /** \brief Apply the law.
		 * 
		 * @param p_i first basis function.
		 * @param p_j second basis function.
		 * @param gp Gauss points to consider
		 * @param Jinv Inverse Jacobian Matrices at the Gauss points
		 * @param ret store the result
		 * @param vm VirtualMachine to use
		 */
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		/** \brief return false
		 * 
		 * @return false
		 */
		virtual bool fractured() const ;
		
		/** \brief Return a Copy of the behaviour
		 * 
		 * @return a new Diffusion
		 */
		virtual Form * getCopy() const ;
		
		/** \brief Return an empty Vector
		 * 
		 * @param s ElementState
		 * @param p_i basis function.
		 * @param gp Gauss points to consider
		 * @param Jinv Inverse Jacobian Matrices at the Gauss points
		 * @param v Vector to store the result
		 */
		virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
		
		virtual Material toMaterial() ;


	} ;

} ;

#endif
