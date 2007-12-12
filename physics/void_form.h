//
// C++ Interface: void_form
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VOID_FORM_H_
#define __VOID_FORM_H_

#include "physics_base.h"

namespace Mu
{

	/** A void law
	*/
	class VoidForm : public LinearForm
	{
	public:
		VoidForm() ;
		
		Matrix constant ;
		
		virtual Matrix apply(const Function & p_i, const Function & p_j,const IntegrableEntity *e)  const ;
		
		virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const ;
		
		/** Step through time
		* 
		* @param timestep length of the timestep
		* @param currentState current state of the element -- behaviour can be dependant on it
		*/
		virtual void step(double timestep, const ElementState & currentState);
		
		virtual void updateElementState(double timestep, ElementState & s) const ;
		
		virtual bool fractured() const ;
		
		virtual ~VoidForm()  ;
		
		virtual Form * getCopy() const ;
		
		virtual Vector getForces(const ElementState & s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const ;
		
	} ;



} ;


#endif
