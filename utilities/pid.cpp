//
// C++ Implementation: pid
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "pid.h"
#include <cmath>

using namespace Mu ;


PidController::PidController(const double & constant, double & var, double &s) : k(constant), p_d(0.3*constant), p_i(constant), control(var), state(s), lastState(state), target(0.), lastCommand(0.), preverror(0.), integralerror(0.), derivativeerror(0.), originalratio(0)
{
	
}
	
void PidController::step()
{
	double apriori_command = 0 ;
	if(!(std::abs(lastState) < 1e-12))
		apriori_command = lastCommand / lastState * target;
	
	double error = target-state ;
	derivativeerror =  error-preverror ;
	integralerror += (error+preverror)*.5 ;
	
	double softeningFactor = 1. ;
	
	if(originalratio == 0 && std::abs(state) > 1e-12 && std::abs(control) > 1e-12)
		originalratio = control/state ;
	
	if(originalratio > 0 && std::abs(lastState) > 1e-12)
	{
		softeningFactor = (lastCommand / lastState ) / originalratio ;
	}
	double K_p = k*softeningFactor ;
	double K_i = p_i*softeningFactor ;
	double K_d = p_d*softeningFactor ;
	control = apriori_command + K_p*error + K_i*integralerror+ K_d*derivativeerror;
	preverror = error ;
	lastState = state ;
	lastCommand = control ;
	
}
	
void PidController::setControlVariable(double & var)
{
	control = var ;
}
void PidController::setCommand(double c)
{
	target = c ;
}
	
void PidController::setStateVariable(double & var)
{
	state = var ;
}
	
void PidController::setProportionalConstant( double c )
{
	k = c ;
}
	
void PidController::setDerivativeConstant( double c )
{
	p_d = c ;
}
	
void PidController::setIntegralConstant( double c )
{
	p_i = c ;
}
	
double PidController::getControlVariable() const 
{
	return control ;
}
	
double PidController::getCommand() const
{
	return target ;
}
	
double PidController::getStateVariable() const
{
	return state ;
}
	
double PidController::getProportionalConstant() const
{
	return k ;
}
double PidController::getDerivativeConstant() const
{
	return p_d ;
}
double PidController::getIntegralConstant() const
{
	return p_i ;
}
