#ifndef __LOGARITHMIC_CREEP_ACCUMULATOR_H_
#define __LOGARITHMIC_CREEP_ACCUMULATOR_H_

#include "../../elements/integrable_entity.h"

namespace Amie
{

struct LogCreepAccumulator
{
	std::vector<double> error ;

	LogCreepAccumulator() { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const { return 1. ; }
	virtual double getKelvinVoigtDashpotReduction() const { return 1. ; }
	virtual Function getKelvinVoigtPreviousFunction() const { return Function("0.5 1 t - *") ; }
	virtual Function getKelvinVoigtNextFunction() const { return Function("0.5 1 t + *") ; }
	virtual LogCreepAccumulator * getCopy() const { return new LogCreepAccumulator() ; }
} ;

struct RealTimeLogCreepAccumulator : public LogCreepAccumulator
{
	double t ;
	double tau ;
	double fakeSpring ;
	RealTimeLogCreepAccumulator(double f = 0.) : LogCreepAccumulator(), t(0.), tau(1.), fakeSpring(f) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const ;
	virtual double getKelvinVoigtDashpotReduction() const ;
	virtual LogCreepAccumulator * getCopy() const { return new RealTimeLogCreepAccumulator() ; }
} ;

struct TimeUnderLoadLogCreepAccumulator : public LogCreepAccumulator
{
	double accumulatedStress ;
	double currentStress ;
	double tau ;
	double previousTimeStep ;
	double realtime ;
	TimeUnderLoadLogCreepAccumulator() : LogCreepAccumulator(), accumulatedStress(0.), currentStress(0.), previousTimeStep(0.), tau(1.), realtime(0.) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const ;
	virtual double getKelvinVoigtDashpotReduction() const ;
	virtual LogCreepAccumulator * getCopy() const { return new TimeUnderLoadLogCreepAccumulator() ; }
} ;

}

#endif