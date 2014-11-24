#include "logcreep_accumulator.h"
#include "../logarithmic_creep.h"
#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

void LogCreepAccumulator::preProcess( double timeStep, ElementState & currentState ) 
{
	if(error.size() == 0)
	{
		error.push_back(1.) ;
		return ;
	}
	error.push_back( error[error.size() -1] + 1./( timeStep * log(1+currentState.getNodalCentralTime() ) ) ) ;
} ;


void RealTimeLogCreepAccumulator::preProcess( double timeStep, ElementState & currentState ) 
{
	LogCreepAccumulator::preProcess(timeStep, currentState) ;

	t = currentState.getNodalCentralTime()+timeStep*0.5 ;
	tau = dynamic_cast<LogarithmicCreep *>(currentState.getParent()->getBehaviour())->tau ;
}

double RealTimeLogCreepAccumulator::getKelvinVoigtSpringReduction() const 
{
//	double x = 1.+t/tau ;
	return fakeSpring ;//*x ;//1./(1.+std::log(x)) ;
}

double RealTimeLogCreepAccumulator::getKelvinVoigtDashpotReduction() const 
{
	return (1.+t/tau)*(1.-fakeSpring)/(1.+fakeSpring) ;
}

void TimeUnderLoadLogCreepAccumulator::preProcess( double timeStep, ElementState & currentState ) 
{
	LogCreepAccumulator::preProcess(timeStep, currentState) ;

	if(previousTimeStep < POINT_TOLERANCE_2D)
		previousTimeStep = timeStep ;

        Vector stress(2+(currentState.getParent()->spaceDimensions()==SPACE_THREE_DIMENSIONAL)) ; stress = 0 ;
	if(!currentState.getParent()->getBehaviour()->fractured())
	{
		dynamic_cast<GeneralizedSpaceTimeViscoElasticElementState &>(currentState).getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress, nullptr, -1, 1.) ;
//		double bulk = 0.5*(stress[0]+stress[1]) ;
		double dt = (timeStep+previousTimeStep)/2. ;
//		stress -= bulk ;
		stress = std::abs(stress) ;
		double prevStress = currentStress ;
		currentStress = stress.max() ;
		accumulatedStress += (currentStress+prevStress)*0.5*dt ;
		realtime = currentState.getNodalCentralTime()-dt ;
        }
        previousTimeStep = timeStep ;
}

double TimeUnderLoadLogCreepAccumulator::getKelvinVoigtSpringReduction() const 
{
	double t = 0. ;
	if(accumulatedStress > POINT_TOLERANCE_2D && currentStress > POINT_TOLERANCE_2D)
		t = accumulatedStress/currentStress ;
	double x = 1.+t/tau ;
	return 1./(1.+std::log(x)) ;
}

double TimeUnderLoadLogCreepAccumulator::getKelvinVoigtDashpotReduction() const 
{
	double t = 0. ;
	if(accumulatedStress > POINT_TOLERANCE_2D && currentStress > POINT_TOLERANCE_2D)
		t = accumulatedStress/currentStress ;
	return 1.+t/tau ;
}

