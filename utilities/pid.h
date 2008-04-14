//
// C++ Interface: pid
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PID_H__
#define __PID_H__

namespace Mu
{

/**PID Controller.
 * This is a utility class provided for easey control of
 * typically load-controled experiments. It can be otherwise
 * used wherever robust control is useful. It is assumed that 
 * the controller has been dimensionned using the Ziegler-Nichols
 * method. Thus it has a single constant, for the proportional 
 * control, from which are deduced the integration and derivation
 * constants. This controller additionally does a priori control 
 * and is adaptative
 */

class PidController
{
protected:
	double k ;           ///< the proportional constant
	double p_d ;         ///< the derivative constant
	double p_i ;         ///< the integral constant
	double &control ;    ///< variable to control
	double &state ;      ///< current state
	double lastState ;   ///< previous state
	double target ;      ///< the value to reach
	double lastCommand ; ///< previous state

	double preverror ;
	double integralerror ;
	double derivativeerror ;
	
	double originalratio ;
public:
	PidController(const double & constant, double & var, double &state) ;
	void step() ;
	
	void setControlVariable(double & var) ;
	void setCommand(double c) ;
	void setStateVariable(double & var) ;
	
	void setProportionalConstant( double c ) ;
	void setDerivativeConstant( double c ) ;
	void setIntegralConstant( double c ) ;
	
	double getControlVariable() const ;
	double getCommand() const;
	double getStateVariable() const;
	
	double getProportionalConstant() const;
	double getDerivativeConstant() const;
	double getIntegralConstant() const;
} ;

} ;

#endif


