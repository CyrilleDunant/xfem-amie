//
// C++ Interface: void_form
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VOID_FORM_H_
#define __VOID_FORM_H_

#include "physics_base.h"

namespace Amie
{

	/** \brief A void law
	*/
	class VoidForm : public LinearForm
	{
	public:
		VoidForm() ;
		
		Matrix constant ;

		/** \brief Return a null Matrix*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const ;
		
		/** \brief Step through time, do nothing
		* 
		* @param timestep length of the timestep
		* @param currentState current state of the element -- behaviour can be dependant on it
		*/
		virtual void step(double timestep, ElementState & currentState, double maxScore);

		/** \brief do nothing*/
		virtual void updateElementState(double timestep, ElementState & s) const ;
		
		/** \brief return false*/
		virtual bool fractured() const ;
		
		virtual ~VoidForm()  ;
		
		/** \brief return another void form*/
		virtual Form * getCopy() const ;

	} ;



} ;


#endif
