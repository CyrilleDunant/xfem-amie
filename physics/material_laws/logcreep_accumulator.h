#ifndef __LOGARITHMIC_CREEP_ACCUMULATOR_H_
#define __LOGARITHMIC_CREEP_ACCUMULATOR_H_

#include "../../elements/integrable_entity.h"

namespace Amie
{

struct LogCreepAccumulator
{
	LogCreepAccumulator() { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) { } ;
	virtual double getKelvinVoigtReduction() const { return 1. ; }
	virtual Function getKelvinVoigtPreviousFunction() const { return Function("0.5 1 t - *") ; }
	virtual Function getKelvinVoigtNextFunction() const { return Function("0.5 1 t + *") ; }
} ;

struct RealTimeLogCreepAccumulator : public LogCreepAccumulator
{
	double t ;
	double tau ;
	RealTimeLogCreepAccumulator() : LogCreepAccumulator(), t(0.), tau(1.) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtReduction() const ;
} ;

struct TimeUnderLoadLogCreepAccumulator : public LogCreepAccumulator
{
	double accumulatedStress ;
	double currentStress ;
	double tau ;
	double previousTimeStep ;
	TimeUnderLoadLogCreepAccumulator() : LogCreepAccumulator(), accumulatedStress(0.), currentStress(0.), previousTimeStep(0.), tau(1.) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtReduction() const ;
} ;

}

#endif
