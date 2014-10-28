#include "logcreep_accumulator.h"
#include "../logarithmic_creep.h"
#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;


void RealTimeLogCreepAccumulator::preProcess( double timeStep, ElementState & currentState ) 
{
	t = currentState.getNodalCentralTime()+timeStep*0.5 ;
	tau = dynamic_cast<LogarithmicCreep *>(currentState.getParent()->getBehaviour())->tau ;
}

double RealTimeLogCreepAccumulator::getKelvinVoigtReduction() const 
{
	double x = 1.+t/tau ;
	return 1./(1./x+std::log(x)) ;
}

void TimeUnderLoadLogCreepAccumulator::preProcess( double timeStep, ElementState & currentState ) 
{
	if(previousTimeStep < POINT_TOLERANCE_2D)
		previousTimeStep = timeStep ;

        Vector stress(2+(currentState.getParent()->spaceDimensions()==SPACE_THREE_DIMENSIONAL)) ; stress = 0 ;
	if(!currentState.getParent()->getBehaviour()->fractured())
	{
		dynamic_cast<GeneralizedSpaceTimeViscoElasticElementState &>(currentState).getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress, nullptr, -1, 1.) ;
		double bulk = 0.5*(stress[0]+stress[1]) ;
		double dt = (timeStep+previousTimeStep)/2. ;
		stress -= bulk ;
		stress = std::abs(stress) ;
		double prevStress = currentStress ;
		currentStress = stress.max() ;
		accumulatedStress += (currentStress+prevStress)*0.5*dt ;
        }
        previousTimeStep = timeStep ;
}

