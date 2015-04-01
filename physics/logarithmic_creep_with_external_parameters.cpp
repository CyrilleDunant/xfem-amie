#include "logarithmic_creep_with_external_parameters.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"

using namespace Amie ;

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), c,d, acc), plane(pt), external(parseDefaultValues(args, sep))
{
	noFracture = (c == nullptr) ;
	if(args.size() > 0)
		makeProperties(external) ;

}

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, std::string ftension, std::string fcompression, DamageModel * d, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), nullptr,d, acc), plane(pt), external(parseDefaultValues(args, sep))
{
	double E = 0. ;	
	if(external.find("bulk_modulus") != external.end())
	{
		double k = external["bulk_modulus"] ;
		double mu = external["shear_modulus"] ;
		E = (9.*k*mu)/(3.*k+mu) ;
	}
	else
		E = external["young_modulus"] ;
	criterion = new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( ftension, fcompression, E ) ;

	noFracture = false ;
	if(args.size() > 0)
		makeProperties(external) ;

}

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), acc), plane(pt), external(parseDefaultValues(args, sep))
{
	noFracture = true ;
	if(args.size() > 0)
		makeProperties(external) ;
}

void LogarithmicCreepWithExternalParameters::makeProperties(std::map<std::string, double> & values, double kvSpringReduction, double kvDashpotReduction) 
{
	prevParam = param ;
	prevEta = eta ;
	prevImposed = imposed ;

	if(values.find("bulk_modulus") != values.end())
	{
		double k_inst = values["bulk_modulus"] ;
		double mu_inst = values["shear_modulus"] ;
		C = Tensor::cauchyGreen( k_inst, mu_inst, false, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
	}
	else
	{
		double E_inst = values["young_modulus"] ;
		double nu_inst = values["poisson_ratio"] ;
		C = Tensor::cauchyGreen( E_inst, nu_inst, true, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
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
			E = Tensor::cauchyGreen( k_visc, mu_visc, false, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
		}
		else
		{
			double E_visc = values["creep_modulus"] ;
			double nu_visc = values["creep_poisson"] ;
			E = Tensor::cauchyGreen( E_visc, nu_visc, true, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
		}
		double a = 0. ;
		if(values.find("recoverable_bulk") != values.end())
		{
			double k_rec = values["recoverable_bulk"] ;
			double mu_rec = values["recoverable_shear"] ;
			R = Tensor::cauchyGreen( k_rec, mu_rec, false, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
		}
		else
		{
			if(values.find("recoverable_modulus") != values.end())
			{
				double E_rec = values["recoverable_modulus"] ;
				double nu_rec = 0.2 ;
				if(values.find("recoverable_poisson") != values.end())
				{
					nu_rec = values["recoverable_poisson"] ;
				}
				else if(values.find("creep_poisson") != values.end())
				{
					nu_rec = values["creep_poisson"] ;
				}
				else
				{
					double k_visc = values["creep_bulk"] ;
					double mu_visc = values["creep_shear"] ;
					nu_rec = (3.*k_visc-2.*mu_visc)/(2.*(3*k_visc+mu_visc)) ;
				}
				R = Tensor::cauchyGreen( E_rec, nu_rec, true, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane) ;
			}
			else
			{
				R = E ;
			}
		}
		a = R[1][1]/E[1][1] ;
		tau = values["creep_characteristic_time"]/exp(-a) ;
		placeMatrixInBlock( C, 1,2, param) ;
		placeMatrixInBlock( C, 2,1, param) ;
		placeMatrixInBlock( C*(-1), 1,0, param) ;
		placeMatrixInBlock( C*(-1), 0,1, param) ;
		placeMatrixInBlock( C*(-1), 2,0, param) ;
		placeMatrixInBlock( C*(-1), 0,2, param) ;
		placeMatrixInBlock( C+E*kvSpringReduction, 1,1, param) ;
		placeMatrixInBlock( E*tau*kvDashpotReduction, 1,1, eta) ;
		placeMatrixInBlock( C+R, 2,2, param) ;
		placeMatrixInBlock( R*tau, 2,2, eta) ;

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
//		if(std::abs(a) < POINT_TOLERANCE)
//		    imposed.resize(0) ;

	}
	else
		imposed.resize(0) ;

	if(imposed.size() > 0 && prevImposed.size() == 0)
	{
		prevImposed.resize(imposed.size()) ;
		prevImposed = 0. ;
	}

	if(!noFracture)
	{
		AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion* crit = dynamic_cast< AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion* >(criterion) ;
		if(crit)
		{
			double tensileStrength = -1 ;
			double tensileUltimateStrength = -1 ;
			if(values.find("tensile_strength")!= values.end())
				tensileStrength = values["tensile_strength"] ;
			if(values.find("tensile_ultimate_strength")!= values.end())
				tensileUltimateStrength = values["tensile_ultimate_strength"] ;
			if(values.find("tensile_strength_decrease_factor")!= values.end())
			{
				if(values.find("tensile_strength")!= values.end())
					tensileUltimateStrength = tensileStrength*values["tensile_strength_decrease_factor"] ;
				else
				{
					tensileStrength = crit->getTensileStrength() ;
					if(tensileStrength > 0)
						tensileUltimateStrength = tensileStrength*(1.-values["tensile_strength_decrease_factor"]) ;
				}
			}
			double compressiveStrength = 1 ;
			double compressiveUltimateStrength = 1 ;
			if(values.find("compressive_strength")!= values.end())
				compressiveStrength = values["compressive_strength"] ;
			if(values.find("compressive_ultimate_strength")!= values.end())
				compressiveUltimateStrength = values["compressive_ultimate_strength"] ;
			if(values.find("compressive_strength_decrease_factor")!= values.end())
			{
				if(values.find("compressive_strength")!= values.end())
					compressiveUltimateStrength = compressiveStrength*(1.-values["compressive_strength_decrease_factor"]) ;
				else
				{
					compressiveStrength = crit->getTensileStrength() ;
					if(compressiveStrength < 0)
						compressiveUltimateStrength = compressiveStrength*values["compressive_strength_decrease_factor"] ;
				}
			}
			double tensileStrain = -1 ;
			double tensileUltimateStrain = -1 ;
			if(values.find("tensile_strain")!= values.end())
				tensileStrain = values["tensile_strain"] ;
			if(values.find("tensile_ultimate_strain")!= values.end())
				tensileUltimateStrain = values["tensile_ultimate_strain"] ;
			if(values.find("tensile_strain_increase_factor")!= values.end())
			{
				if(values.find("tensile_strain")!= values.end())
					tensileUltimateStrain = tensileStrain*values["tensile_strain_increase_factor"] ;
				else
				{
					tensileStrain = crit->getTensileStrain() ;
					if(tensileStrain > 0)
						tensileUltimateStrain = tensileStrain*values["tensile_strain_increase_factor"] ;
				}
			}
			double compressiveStrain = 1 ;
			double compressiveUltimateStrain = 1 ;
			if(values.find("compressive_strain")!= values.end())
				compressiveStrain = values["compressive_strain"] ;
			if(values.find("compressive_ultimate_strain")!= values.end())
				compressiveUltimateStrain = values["compressive_ultimate_strain"] ;
			if(values.find("compressive_strain_increase_factor")!= values.end())
			{
				if(values.find("compressive_strain")!= values.end())
					compressiveUltimateStrain = compressiveStrain*values["compressive_strain_increase_factor"] ;
				else
				{
					compressiveStrain = crit->getCompressiveStrain() ;
					if(compressiveStrain < 0)
						compressiveUltimateStrain = compressiveStrain*values["compressive_strain_increase_factor"] ;
				}
			}

			crit->setMaximumTensileStress( tensileStrength, tensileUltimateStrength ) ;
			crit->setMaximumTensileStrain( tensileStrain, tensileUltimateStrain ) ;
			crit->setMaximumCompressiveStress( compressiveStrength, compressiveUltimateStrength ) ;
			crit->setMaximumCompressiveStrain( compressiveStrain, compressiveUltimateStrain ) ;
		}
			
		
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
		std::stringstream stream ;
		stream << std::fixed << std::setprecision(16) << ext.second ;
		args.append(stream.str()) ;
		args.append(",") ;
	}
	args = args.substr(0, args.size()-1) ;

	if(noFracture)
	{
		copy = new LogarithmicCreepWithExternalParameters( args,  accumulator->getCopy(), ((int) param.numCols() == 3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane)  ;
	}
	else
	{
		copy = new LogarithmicCreepWithExternalParameters( args, criterion->getCopy(), dfunc->getCopy(),  accumulator->getCopy(), ((int) param.numCols() == 3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane ) ;
		copy->dfunc->getState(true).resize(dfunc->getState().size());
		copy->dfunc->getState(true) = dfunc->getState() ;
		copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
		copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
		copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	}

        copy->setBlocks(blocks) ;
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
	if( reducedTimeStep > POINT_TOLERANCE)
	{
		return ;
	}

	accumulator->preProcess(timeStep, currentState) ;
	dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState).synchronize(external) ;
	for(size_t i = 0 ; i < relations.size() ; i++)
		relations[i]->preProcess( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState), timeStep ) ;
	std::map<std::string, double> prop = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState).getVariables() ;
	makeProperties( prop, accumulator->getKelvinVoigtSpringReduction(), accumulator->getKelvinVoigtDashpotReduction() ) ;
	currentState.getParent()->behaviourUpdated = true ;
}

LogarithmicCreepWithExternalParameters::~LogarithmicCreepWithExternalParameters() 
{ 
/*	if(dfunc) 
		{delete dfunc ; } 
	if(criterion) 
		{delete criterion ;} */
	for(size_t i = 0 ; i < relations.size() ; i++)
	{
		if(relations[i])
			delete relations[i] ;
	}
	relations.resize(0) ;
}

std::vector<BoundaryCondition * > LogarithmicCreepWithExternalParameters::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = LogarithmicCreep::getBoundaryConditions(s, id, p_i, gp, Jinv ) ;
    if(imposed.size() == 0 || ((imposed.max() < POINT_TOLERANCE) && (imposed.min() > -POINT_TOLERANCE)) || fractured())
        return ret ;
    Vector istress = C*imposed ;
    if(dfunc)
        istress = dfunc->apply(C) * imposed   ;
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, static_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, static_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
    }
    if(v.size() == 4)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, static_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, static_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, static_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
    }
    return ret ;

}
