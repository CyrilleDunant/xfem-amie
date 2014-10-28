#include "logarithmic_creep_with_external_parameters.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Amie ;

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), c,d, acc), external(parseDefaultValues(args, sep)), plane(pt)
{
	noFracture = (c == nullptr) ;
	makeProperties(external) ;
}

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), acc), external(parseDefaultValues(args, sep))
{
	noFracture = true ;
	makeProperties(external) ;
}

void LogarithmicCreepWithExternalParameters::makeProperties(std::map<std::string, double> & values, double kvReduction) 
{
	prevParam = param ;
	prevEta = eta ;
	prevImposed = imposed ;

	if(values.find("bulk_modulus") != values.end())
	{
		double k_inst = values["bulk_modulus"] ;
		double mu_inst = values["shear_modulus"] ;
		C = Material::cauchyGreen( k_inst, mu_inst, false, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
	}
	else
	{
		double E_inst = values["young_modulus"] ;
		double nu_inst = values["poisson_ratio"] ;
		C = Material::cauchyGreen( E_inst, nu_inst, true, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
	}
	param = 0. ;
	placeMatrixInBlock( C, 0,0, param) ;
	if(values.find("creep_characteristic_time") != values.end())
	{
		isPurelyElastic = false ;
		if(values.find("creep_bulk") != values.end())
		{
			double k_visc = values["creep_bulk"] ;
			double mu_visc = values["creep_shear"] ;
			E = Material::cauchyGreen( k_visc, mu_visc, false, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
		}
		else
		{
			double E_visc = values["creep_modulus"] ;
			double nu_visc = values["creep_poisson"] ;
			E = Material::cauchyGreen( E_visc, nu_visc, true, (param.numCols()==6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
		}
		tau = values["creep_characteristic_time"] ;
		placeMatrixInBlock( C*(-1), 1,0, param) ;
		placeMatrixInBlock( C*(-1), 0,1, param) ;
		placeMatrixInBlock( C+E*kvReduction, 1,1, param) ;
		placeMatrixInBlock( E*kvReduction*tau, 1,1, eta) ;

		model = GENERALIZED_KELVIN_VOIGT ;

//		(param-eta/tau).print() ;

	}
	else
	{
		eta = 0. ;
		model = PURE_ELASTICITY ;
		isPurelyElastic = true ;
	}
	if(values.find("imposed_deformation") != values.end())
	{
		imposed.resize(C.numCols()) ;
		double a = values["imposed_deformation"] ;
		for(size_t i = 0 ; i < 2+(imposed.size()==6) ; i++)
		    imposed[i] = a ;
		if(std::abs(a) < POINT_TOLERANCE_2D)
		    imposed.resize(0) ;
	}
	else
		imposed.resize(0) ;

	if(imposed.size() > 0 && prevImposed.size() == 0)
	{
		prevImposed.resize(imposed.size()) ;
		prevImposed = 0. ;
	}

}

ElementState * LogarithmicCreepWithExternalParameters::createElementState( IntegrableEntity * e)
{
	return new GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables(e, external) ;
}

Form * LogarithmicCreepWithExternalParameters::getCopy() const
{
	LogarithmicCreepWithExternalParameters * copy ;
	std::string args ;
	for(auto ext : external)
	{
		args.append(ext.first) ;
		args.append(" = ") ;
		args.append(std::to_string(ext.second)) ;
		args.append(",") ;
	}
	args = args.substr(0, args.size()-1) ;

	if(noFracture)
	{
		copy = new LogarithmicCreepWithExternalParameters( args,  accumulator, (param.numCols() == 6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane)  ;
	}
	else
	{
		copy = new LogarithmicCreepWithExternalParameters( args,criterion->getCopy(), dfunc->getCopy(),  accumulator, (param.numCols() == 6) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane ) ;
		copy->dfunc->getState(true).resize(dfunc->getState().size());
		copy->dfunc->getState(true) = dfunc->getState() ;
		copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
		copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
		copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	}

	for(size_t i = 0 ; i < relations.size() ; i++)
		copy->relations.push_back(relations[i]) ;

	return copy ;

}

void LogarithmicCreepWithExternalParameters::step(double timestep, ElementState &s, double maxScore)
{
    for(size_t i = 0 ; i < relations.size() ; i++)
        relations[i]->step( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s), timestep) ;

    LogarithmicCreepWithImposedDeformationAndFracture::step(timestep, s, maxScore) ;
}

void LogarithmicCreepWithExternalParameters::preProcess( double timeStep, ElementState & currentState )
{
	if(reducedTimeStep > POINT_TOLERANCE_2D)
	{
		return ;
	}

	accumulator->preProcess(timeStep, currentState) ;
	std::map<std::string, double> prop = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState).getVariables() ;
	makeProperties( prop, accumulator->getKelvinVoigtReduction() ) ;
	currentState.getParent()->behaviourUpdated = true ;
}
