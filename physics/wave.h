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

#ifndef __WAVE_H_
#define __WAVE_H_

#include "physics_base.h"

namespace Mu
{

	/** \brief A linear Elastic Law
	* The field param is the stiffness Tensor
	*/
	struct Wave : public LinearForm
	{
		std::vector<Variable> v ;
		/** \brief Constructor
		* 
		* @param rig Complete expression of the stiffness Tensor
		*/
		Wave(const Matrix & rig) ;
		
		virtual ~Wave() ;
		
/** \brief Apply the behaviour
	* This overloaded apply() is more efficient and is designed to minimise allocating and dealocating memory.
	* return the matrix resulting of \f$ \dot{H}\cdot \dot{H} - \nabla\cdot H^T K \nabla\cdot H \f$.
	* 
	* @param p_i first shape function.
	* @param p_j second shape function.
	* @param gp Set of gauss points for numerical integration
	* @param Jinv inverse jacobian matrices at the gauss points
	* @param ret matrix in which to sore the results
	* @param vm pointer to the virtual machine dedicated for the computation
	*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
	/** \brief return false*/
		virtual bool fractured() const ;
		
		/** \brief return a copy of the behaviour */
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif
