//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __STIFFNESS_AND_FRACTURE_H
#define __STIFFNESS_AND_FRACTURE_H

#include "physics_base.h"
#include "fracturecriterion.h"

namespace Mu
{

	struct StiffnessAndFracture : public LinearForm
	{
		Matrix previousParam ;
		FractureCriterion * criterion ;
		bool frac ;
		bool change ; 
		double sigmaRupt ;
		double init ;
		StiffnessAndFracture(const Matrix & rig, FractureCriterion *)  ;
		
		virtual ~StiffnessAndFracture();
		
		/** Perform the integration
		* 
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @param e element in which the integration is done.
		* @return matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const ;
		
		/** Apply the law.
		* @param p_i first basis polynomial.
		* @param p_j second basis polynomial.
		* @param gp integration points
		* @param Jinv inverse Jacobian matrix at the integration points
		* @return matrix resulting of \f$ \nabla H^T K \nabla H \f$.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const ;
		
		/** Check for fracture
		* @param timestep elapsed time
		* @param currentState state of the element
		* 
		* if the Von Mises yield criterion is true, se fractured state to true
		*/
		virtual void step(double timestep, ElementState & currentState) ;
		
		/** Check for fracture state
		*
		* @return true if the element is fractured
		*/
		virtual bool fractured() const ;
		
		
		/** get Copy of the behaviour
		*
		* @return pointer to the copy. Caller is responsible for cleaning memory
		*/
		virtual Form * getCopy() const ;
		
		virtual Vector getForces(const ElementState & s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const ;

		virtual bool changed() const ;
		
		virtual void stepBack() ;
		
	} ;
	


} ;


#endif 
