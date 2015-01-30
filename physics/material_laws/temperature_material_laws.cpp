#include "temperature_material_laws.h"
#include <stdlib.h>

namespace Amie
{

void ThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
    s.add("imposed_deformation", alpha*(T-T0)) ;
}

void RadiationInducedExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double kappa = s.get("radiation_expansion_delay",defaultValues) ;
    double epsmax = s.get("maximum_radiation_expansion",defaultValues) ;
    double delta = s.get("neutron_fluence_correction",defaultValues) ;
    double N = s.get("neutron_fluence",defaultValues) ;
    s.add("imposed_deformation", (kappa*epsmax*(exp(delta*N)-1)/(epsmax+kappa*exp(delta*N)))) ;
}

void DryingShrinkageMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double h = s.get("relative_humidity",defaultValues) ;
    double h0 = defaultValues["relative_humidity"] ;
    if(h > h0 - POINT_TOLERANCE_2D)
        return ;
    double ks = s.get("drying_shrinkage_coefficient",defaultValues) ;
    s.add("imposed_deformation", -ks*(h0-h)) ;
}

void LoadNonLinearCreepMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_modulus"))
	return ;

    Vector stress(3) ; stress=0. ;
    s.getAverageField(REAL_STRESS_FIELD, stress, nullptr, -1, 1.) ;
    double pstress =  stress.max() ;
    double strength = s.get("tensile_strength", defaultValues) ;
    if(stress.min() < 0 && -stress.min() > stress.max())
    {
	pstress = stress.min() ;
	strength = s.get("compressive_strength", defaultValues) ;
    }
    double C = s.get("young_modulus",defaultValues) ; //s.get("creep_modulus",defaultValues) ;
    double alpha = s.get("imposed_deformation", defaultValues) ;
    pstress += C*alpha ;
    double f = pstress/strength ;
    double factor = (1.-std::pow(f, 10.))/(1.+f*f) ;
    s.multiply("creep_modulus", factor) ;
    s.multiply("recoverable_modulus", factor) ;
}

void StrainRateDependentStrengthMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("tensile_strength") && !s.has("compressive_strength") )
	return ;

    if(!s.has("previous_principal_strain_0") || !s.has("previous_principal_strain_1"))
    {
	s.set("previous_principal_strain_0", 0) ;
	s.set("previous_principal_strain_1", 0) ;
    }

    Vector pstrain(2) ; pstrain=0. ;
    s.getAverageField(PRINCIPAL_STRAIN_FIELD, pstrain, nullptr, -1, 1.) ;
    if(std::abs(pstrain.max()) < POINT_TOLERANCE_2D*dt && std::abs(pstrain.min()) < POINT_TOLERANCE_2D*dt)
	return ;

    double strainrate = 0. ;
    if( std::abs(pstrain[0]) > std::abs(pstrain[1])  )
	strainrate = pstrain[0]-s.get("previous_principal_strain_0",defaultValues) ;
    else
	strainrate = pstrain[1]-s.get("previous_principal_strain_1",defaultValues) ;

    if(strainrate < 0)
	strainrate *= -1. ;

    strainrate /= s.getNodalDeltaTime() ;

    if(strainrate < strainRateRef)
    {

	    double f = strainrate/strainRateRef ;

	    if(s.has("tensile_strength"))
	    {
		s.multiply("tensile_strength", 0.6+0.4*std::pow(f, p)) ;
		s.multiply("ultimate_tensile_strain", 1.+std::log10(1./f)) ;
	    }

	    if(s.has("compressive_strength"))
	    {
		s.multiply("compressive_strength", 0.6+0.4*std::pow(f, p)) ;
		s.multiply("ultimate_compressive_strain", 1.+std::log10(1./f)) ;
	    }
    }

    s.set("previous_principal_strain_0",pstrain[0]);
    s.set("previous_principal_strain_1",pstrain[1]);

}

ArrheniusMaterialLaw::ArrheniusMaterialLaw(std::string a, std::string args, char sep) : ExternalMaterialLaw(args, sep), affected(a)
{
    coefficient = affected ; coefficient.append("_activation_energy") ;
}

void ArrheniusMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(defaultValues.find(affected) == defaultValues.end())
	defaultValues[affected] = s.get(affected, defaultValues) ;


    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double Ea = s.get(coefficient, defaultValues) ;
    s.multiply( affected, exp( Ea*(1./T-1./T0) ) ) ;
}

void CreepArrheniusMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_characteristic_time"))
        return ;

    if(defaultValues.find("creep_characteristic_time") == defaultValues.end())
	defaultValues["creep_characteristic_time"] = s.get("creep_characteristic_time", defaultValues) ;

    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double Ea = s.get("creep_activation_energy", defaultValues) ;
    s.multiply("creep_characteristic_time", exp( Ea*(1./T-1./T0) )) ;
}

void CreepRelativeHumidityMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_characteristic_time"))
        return ;

    double h = s.get("relative_humidity", defaultValues) ;
    double hc = s.get("creep_humidity_coefficient", defaultValues) ;
    double factor =  hc*(1-h) + exp(-hc*(1-h)) ;
    if(s.has("creep_modulus"))
    {
	s.multiply("creep_modulus", factor) ;
    } else {
	s.multiply("creep_bulk", factor) ;
	s.multiply("creep_shear", factor) ;
    }
    if(s.has("recoverable_modulus"))
    {
	s.multiply("recoverable_modulus", factor) ;
    } else {
	s.multiply("recoverable_bulk", factor) ;
	s.multiply("recoverable_shear", factor) ;
    }

}


} ;

