#ifndef __LOGARITHMIC_CREEP_ACCUMULATOR_H_
#define __LOGARITHMIC_CREEP_ACCUMULATOR_H_

#include "../../elements/integrable_entity.h"

namespace Amie
{

struct LogCreepAccumulator
{
	std::vector<double> error ;

	LogCreepAccumulator() { } ;
        virtual ~LogCreepAccumulator() { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const { return 1. ; }
	virtual double getKelvinVoigtDashpotReduction() const { return 1. ; }
	virtual Function getKelvinVoigtPreviousFunction() const { return Function("0.5 1 t - *") ; }
	virtual Function getKelvinVoigtNextFunction() const { return Function("0.5 1 t + *") ; }
	virtual LogCreepAccumulator * getCopy() const { return new LogCreepAccumulator() ; }
} ;

/*PARSE No LogCreepAccumulator */
struct NoLogCreepAccumulator : public LogCreepAccumulator
{
	NoLogCreepAccumulator() : LogCreepAccumulator() { } ;
        virtual ~NoLogCreepAccumulator() { } ;
	virtual LogCreepAccumulator * getCopy() const { return new NoLogCreepAccumulator() ; }
} ;

/*PARSE StrainConsolidation LogCreepAccumulator */
struct StrainConsolidationLogCreepAccumulator : public LogCreepAccumulator
{
	double currentStrain = 0 ;
	double currentMaxwellStrain = 0 ;
	double fakeSpring = 0 ;
	double k = 0 ;
	double tau = 0 ;

	StrainConsolidationLogCreepAccumulator() : LogCreepAccumulator() { } ;
        virtual ~StrainConsolidationLogCreepAccumulator() { } ;
	virtual LogCreepAccumulator * getCopy() const { return new StrainConsolidationLogCreepAccumulator() ; }
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const ;
	virtual double getKelvinVoigtDashpotReduction() const ;
} ;


/*PARSE RealTime LogCreepAccumulator
	@value[viscous_flow] 0 // initial value of the viscous strain
*/
struct RealTimeLogCreepAccumulator : public LogCreepAccumulator
{
	double t ;
	double tau ;
	double fakeSpring = 0.;
        double flow ;
	RealTimeLogCreepAccumulator(double f = 0.) : LogCreepAccumulator(), t(0.), tau(1.), flow(f) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const ;
	virtual double getKelvinVoigtDashpotReduction() const ;
	virtual LogCreepAccumulator * getCopy() const { return new RealTimeLogCreepAccumulator() ; }
} ;

/*PARSE TimeUnderLoad LogCreepAccumulator */
struct TimeUnderLoadLogCreepAccumulator : public LogCreepAccumulator
{
	double accumulatedStress ;
	double currentStress ;
	double tau ;
	double previousTimeStep ;
	double realtime ;
	TimeUnderLoadLogCreepAccumulator() : LogCreepAccumulator(), accumulatedStress(0.), currentStress(0.), tau(1.), previousTimeStep(0.), realtime(0.) { } ;
	virtual void preProcess( double timeStep, ElementState & currentState ) ;
	virtual double getKelvinVoigtSpringReduction() const ;
	virtual double getKelvinVoigtDashpotReduction() const ;
	virtual LogCreepAccumulator * getCopy() const { return new TimeUnderLoadLogCreepAccumulator() ; }
} ;

}

#endif
