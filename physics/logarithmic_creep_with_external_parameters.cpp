#include "logarithmic_creep_with_external_parameters.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "fracturecriteria/mazars.h"
#include "../utilities/object_translator.h"

using namespace Amie ;

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::string args, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt, char sep) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), c,d, acc), plane(pt), external(parseDefaultValues(args, sep))
{
	noFracture = (c == nullptr) ;
	if(args.size() > 0)
		makeProperties(external) ;

}

LogarithmicCreepWithExternalParameters::LogarithmicCreepWithExternalParameters(std::map<std::string, double> values, ExternalMaterialLawList * mat, FractureCriterion * c, DamageModel * d, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt) : LogarithmicCreepWithImposedDeformationAndFracture( Matrix( 3+3*(dim == SPACE_THREE_DIMENSIONAL), 3+3*(dim == SPACE_THREE_DIMENSIONAL)), Vector(), c,d, acc), plane(pt), external( values )
{
	noFracture = (c == nullptr) ;
	if(values.size() > 0)
		makeProperties(external) ;

	if(mat != nullptr)
	{
		for(auto law = mat->begin() ; law != mat->end() ; law++)
			addMaterialLaw( *law ) ;
	}

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
	bool modelChange = false ;

	if(values.find("young_modulus") != values.end() && values.find("poisson_ratio") != values.end())
	{
		double E_inst = values["young_modulus"] ;
		double nu_inst = values["poisson_ratio"] ;
                if(values.find("angle") != values.end() && values.find("young_modulus_anisotropy_coefficient") != values.end())
			C = Tensor::orthotropicCauchyGreen( E_inst, E_inst*values["young_modulus_anisotropy_coefficient"], E_inst, nu_inst, values["angle"], plane) ;
		else
			C = Tensor::cauchyGreen( E_inst, nu_inst, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
	}
	else if(values.find("young_modulus") != values.end() && values.find("shear_modulus") != values.end())
	{
		double E_inst = values["young_modulus"] ;
		double G_inst = values["shear_modulus"] ;
		C = Tensor::cauchyGreen( E_inst, G_inst, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
	}
	else 

	{
		double k_inst = values["bulk_modulus"] ;
		double mu_inst = values["shear_modulus"] ;
		C = Tensor::cauchyGreen( k_inst, mu_inst, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, BULK_SHEAR) ;
	}
        Point angle ;
        if(values.find("angle_x") != values.end())
            angle.setX( values["angle_x"] ) ;
        if(values.find("angle_y") != values.end())
            angle.setY( values["angle_y"] ) ;
        if(values.find("angle_z") != values.end())
            angle.setZ( values["angle_z"] ) ;
        for(size_t i = 0 ; i < relations.size() ; i++)
            relations[i]->preProcess( C, angle, plane ) ;
	param = 0. ;
	placeMatrixInBlock( C, 0,0, param) ;

	if(values.find("creep_characteristic_time") != values.end())
	{
		isPurelyElastic = false ;
		if(values.find("creep_modulus") != values.end())
		{
			double E_visc = values["creep_modulus"] ;
			if(values.find("creep_poisson") != values.end())
			{
				double nu_visc = values["creep_poisson"] ;
				E = Tensor::cauchyGreen( E_visc, nu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
			}
			else if(values.find("creep_shear") != values.end())
			{
				double G_visc = values["creep_shear"] ;
				E = Tensor::cauchyGreen( E_visc, G_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
			}
			else if(values.find("poisson_ratio") != values.end())
			{
				double nu_visc = values["poisson_ratio"] ;
				E = Tensor::cauchyGreen( E_visc, nu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
			}
			else if(values.find("shear_modulus") != values.end())
			{
				double G_visc = values["shear_modulus"] ;
				E = Tensor::cauchyGreen( E_visc, G_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
			}
		}
		else
		{
			double k_visc = values["creep_bulk"] ;
			double mu_visc = values["creep_shear"] ;
			E = Tensor::cauchyGreen( k_visc, mu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, BULK_SHEAR) ;
		}
		double a = 0. ;
		if(values.find("recoverable_modulus") != values.end())
		{
			double E_rec = values["recoverable_modulus"] ;
			if(values.find("recoverable_poisson") != values.end())
			{
				double nu_visc = values["recoverable_poisson"] ;
				R = Tensor::cauchyGreen( E_rec, nu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
			}
			else if(values.find("recoverable_shear") != values.end())
			{
				double G_visc = values["recoverable_shear"] ;
				R = Tensor::cauchyGreen( E_rec, G_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
			}
			else if(values.find("creep_poisson") != values.end())
			{
				double nu_visc = values["creep_poisson"] ;
				R = Tensor::cauchyGreen( E_rec, nu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
			}
			else if(values.find("creep_shear") != values.end())
			{
				double G_visc = values["creep_shear"] ;
				R = Tensor::cauchyGreen( E_rec, G_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
			}
			else if(values.find("poisson_ratio") != values.end())
			{
				double nu_visc = values["poisson_ratio"] ;
				R = Tensor::cauchyGreen( E_rec, nu_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_POISSON) ;
			}
			else if(values.find("shear_modulus") != values.end())
			{
				double G_visc = values["shear_modulus"] ;
				R = Tensor::cauchyGreen( E_rec, G_visc, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, YOUNG_SHEAR) ;
			}
                        else
				R = E ;
		}
		else
		{
			if(values.find("recoverable_bulk") != values.end())
			{
				double k_rec = values["recoverable_bulk"] ;
				double mu_rec = values["recoverable_shear"] ;
				R = Tensor::cauchyGreen( k_rec, mu_rec, ((int) param.numCols()==3*blocks) ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL, plane, BULK_SHEAR) ;
			}
			else
				R = E ;
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

		modelChange = (model != GENERALIZED_KELVIN_VOIGT) ;
		model = GENERALIZED_KELVIN_VOIGT ;

//		(param-eta/tau).print() ;

	}
	else
	{
		eta = 0. ;
		modelChange = (model != PURE_ELASTICITY) ;
		model = PURE_ELASTICITY ;
		isPurelyElastic = true ;
	}
	if(values.find("imposed_deformation") != values.end())
	{
		imposed.resize(C.numCols()) ;
		double a = values["imposed_deformation"] ;
		for(size_t i = 0 ; i < 2+(imposed.size()==6) ; i++)
			imposed[i] = a ;
		if(values.find("imposed_deformation_xx") != values.end())
			imposed[0] += values["imposed_deformation_xx"] ;
		if(values.find("imposed_deformation_yy") != values.end())
			imposed[1] += values["imposed_deformation_yy"] ;
		if(values.find("imposed_deformation_zz") != values.end())
			imposed[2] += values["imposed_deformation_zz"] ;

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
		NonLocalSpaceTimeMazars* mazar = dynamic_cast< NonLocalSpaceTimeMazars* > (criterion) ;
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
		else if(mazar)
		{
			double thr = mazar->threshold ;
			double E = values["young_modulus"] ;
			double nu = values["poisson_ratio"] ;
			if(values.find("tensile_strain") == values.end())
				thr = values["tensile_stress"]/E ;
			else
				thr = values["tensile_strain"] ;
			double Gf = values["fracture_energy"] ;
			if(values.find("compressive_stress") != values.end() || values.find("compressive_strain") != values.end())
			{
				double cstress = values["compressive_stress"] ;
				double cstrain = values["compressive_strain"] ;
				mazar->reset( thr, E, nu, Gf, cstress, cstrain ) ;
			}
			else
				mazar->reset( thr, E, nu, Gf ) ;
		}
		else
		{
			std::map<std::string, std::string> str ;
			Object::resetFractureCriterion( criterion, values, str ) ;
		}
			
		
	}

	if(prevParam.array().max() < POINT_TOLERANCE)
	{
		prevParam = param ;
		prevEta = eta ;
		prevImposed = imposed ;
	}

	if(modelChange)
		makeBlockConnectivity() ;


}

ElementState * LogarithmicCreepWithExternalParameters::createElementState( IntegrableEntity * e)
{
	return new GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables(e, external, blocks) ;
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
/*		copy->dfunc->getState(true).resize(dfunc->getState().size());
		copy->dfunc->getState(true) = dfunc->getState() ;
		copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
		copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
                copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());*/
	}

        copy->setBlocks(blocks) ;
	for(size_t i = 0 ; i < relations.size() ; i++)
	{
		copy->relations.push_back(relations[i]) ;
		WeibullDistributedMaterialLaw * weib = dynamic_cast<WeibullDistributedMaterialLaw *>(relations[i]) ;
		if(weib != nullptr)
		{
			std::pair< std::string, double> w = weib->getWeibullVariable() ;
			if(external.find(w.first) == external.end())
				copy->addMaterialParameter(w.first, w.second) ;
		}
	}

	return copy ;

}

void LogarithmicCreepWithExternalParameters::step(double timestep, ElementState &s, double maxScore)
{
    for(size_t i = 0 ; i < relations.size() ; i++)
        relations[i]->step( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s), timestep) ;

    LogarithmicCreepWithImposedDeformationAndFracture::step(timestep, s, maxScore) ;
}

void LogarithmicCreepWithExternalParameters::replaceMaterialLaw( ExternalMaterialLaw * law, int index ) 
{
    if(law == nullptr)
    {
       removeMaterialLaw(index) ;
       return ;
    }
    if(index < 0 || index >= (int) relations.size())
    {
        addMaterialLaw( law ) ;
        return ;
    }
    relations[index] = law ;
}

void LogarithmicCreepWithExternalParameters::removeMaterialLaw( int index ) 
{
    std::vector<ExternalMaterialLaw *> tmp ;
    for(size_t j = 0 ; j < relations.size() ; j++)
    {
        if((int) j != index)
            tmp.push_back( relations[j] ) ;
    }
    relations.clear() ;
    relations = tmp ;
}


void LogarithmicCreepWithExternalParameters::preProcess( double timeStep, ElementState & currentState )
{
	if( reducedTimeStep > POINT_TOLERANCE)
	{
		return ;
	}

        #pragma omp critical(logcreep)
        {

		dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState).synchronize(external) ;
		for(size_t i = 0 ; i < relations.size() ; i++)
		{
			relations[i]->preProcess( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState), timeStep ) ;
		}
		accumulator->preProcess(timeStep, currentState) ;
		std::map<std::string, double> prop = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(currentState).getVariables() ;
		makeProperties( prop, accumulator->getKelvinVoigtSpringReduction(), accumulator->getKelvinVoigtDashpotReduction() ) ;

		if(dfunc) { dfunc->prepare() ; }
        }
	currentState.getParent()->behaviourUpdated = ( abs( param.array()-prevParam.array() ).max() > 1 )  ;
	currentState.getParent()->behaviourViscoUpdated = ( abs( eta.array()-prevEta.array() ).max() > 1 ) ;
        if( imposed.size() == prevImposed.size() && abs( imposed-prevImposed ).max() > POINT_TOLERANCE ) { currentState.getParent()->behaviourForcesUpdated = true ; }
}

LogarithmicCreepWithExternalParameters::~LogarithmicCreepWithExternalParameters() 
{ 
/*	if(dfunc) 
		{delete dfunc ; } 
	if(criterion) 
		{delete criterion ;} 
	for(size_t i = 0 ; i < relations.size() ; i++)
	{
		if(relations[i])
			delete relations[i] ;
	}
	relations.resize(0) ;*/
}

std::vector<BoundaryCondition * > LogarithmicCreepWithExternalParameters::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(!noFracture && fractured())
        return ret ;

    return LogarithmicCreepWithImposedDeformationAndFracture::getBoundaryConditions( s, id, p_i, gp, Jinv ) ;



}
