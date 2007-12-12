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

namespace Mu
{

	/** A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	*/
	struct Diffusion : public LinearForm
	{
		
		/** Constructor
		* 
		* @param rig Complete expression of the Cauchy-Green Strain Tensor
		*/
		Diffusion(const Matrix & rig) ;
		
		virtual ~Diffusion() ;
		
		/** Apply the law.
		* 
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @return symbolic matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const; 
		
		virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const ;
		
		virtual bool fractured() const ;
		
		virtual Form * getCopy() const ;
		
		virtual Vector getForces(const ElementState & s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const ;
		
	} ;

} ;

#endif
