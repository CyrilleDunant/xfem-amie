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

#ifndef __STIFFNESS_WITH_DEF
#define __STIFFNESS_WITH_DEF

#include "physics_base.h"

namespace Mu
{

	
	/** A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	* The imposed deformation are given as a vector
	* This is used for the simulation of periodic conditions
	*/
	struct StiffnessWithImposedDeformation : public LinearForm
	{
		Vector imposed ;
		Function x ;
		Function y ;
		/** Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		*/
		StiffnessWithImposedDeformation(const Matrix & rig, Vector imposedDef) ;
		
		virtual ~StiffnessWithImposedDeformation() ;
		
		/** Apply the law.
		* 
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @return symbolic matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
		
		virtual Matrix apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
		virtual bool fractured() const ;
		
		virtual Form * getCopy() const ;
		
		virtual bool hasInducedForces() const ;
		
		
		virtual Vector getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;
		
	} ;


} ;

#endif
