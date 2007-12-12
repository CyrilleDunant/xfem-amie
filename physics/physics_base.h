//
// C++ Interface: physics_base
//
// Description: template behaviours : abstract linear and non-linear behaviours
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PHYSICS_BASE__H_
#define __PHYSICS_BASE__H_

#include "../matrixops.h"
#include "../polynomial/vm_base.h"
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_2D.h"

namespace Mu
{

	/** A Linear Form 
	*/
	class LinearForm : public Form
	{
	public:
		LinearForm(const Matrix & p, bool t = false, bool s = false, size_t numdof = 0 ) ;
		
		virtual ~LinearForm()  ;
		
		virtual bool fractured() const ; 
		
		virtual Form * getCopy() const = 0 ;
		
		virtual void step(double timestep, ElementState & s) ;
		
		virtual void updateElementState(double timestep, ElementState & currentState) const ;
		
	} ;
	
	/** A Non-Linear Form 
	*/
	class NonLinearForm : public Form
	{
	protected:
		std::vector<Point> hints ;

	public:
		NonLinearForm() ;
		
		virtual ~NonLinearForm() ;
		
		virtual std::vector<Point> getIntegrationHints() ;
		
		virtual Point getIntegrationHint(size_t i) ;
		
		virtual void setIntegrationHints(std::vector<Point> h) ;
		
		virtual void step(double timestep, ElementState & s) ;
		
		virtual void updateElementState(double timestep, ElementState & currentState) const ;
		
		virtual bool isActive() const = 0 ;
		
		virtual bool fractured() const ;
		
		virtual bool hasInducedMatrix() const = 0;
		
		virtual Form * getCopy() const =0 ;
		
	} ;
	
	/** A Linear Form and a constant
	*/
	class LinearFormAndConstant : public Form
	{
	public:
		LinearFormAndConstant(const Matrix & p, Matrix c) ;
		
		Matrix constant ;
		virtual ~LinearFormAndConstant() ;
	} ;


} ;


#endif
