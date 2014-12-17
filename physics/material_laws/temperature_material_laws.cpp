#include "temperature_material_laws.h"
#include <stdlib.h>

namespace Amie
{

void ThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
    double imp = s.get("imposed_deformation", defaultValues) ;
    imp += alpha*(T-T0) ;
    s.set("imposed_deformation", imp) ;
}

void RadiationInducedExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double kappa = s.get("radiation_expansion_delay",defaultValues) ;
    double epsmax = s.get("maximum_radiation_expansion",defaultValues) ;
    double delta = s.get("neutron_fluence_correction",defaultValues) ;
    double N = s.get("neutron_fluence",defaultValues) ;
    double imp = s.get("imposed_deformation", defaultValues) ;
    imp += (kappa*epsmax*(exp(delta*N)-1)/(epsmax+kappa*exp(delta*N))) ;
    s.set("imposed_deformation", imp) ;
}

void DryingShrinkageMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double h = s.get("relative_humidity",defaultValues) ;
    double h0 = defaultValues["relative_humidity"] ;
    if(h > h0 - POINT_TOLERANCE_2D)
        return ;
    double ks = s.get("drying_shrinkage_coefficient",defaultValues) ;
    double imp = s.get("imposed_deformation", defaultValues) ;
    if(effective)
	s.set("effective_imposed_deformation", -ks*(h0-h)) ;
    imp -= ks*(h0-h) ;
    s.set("imposed_deformation", imp) ;
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
    double prop = defaultValues[affected] ;
    prop *= exp( Ea*(1./T-1./T0) ) ;
    s.set(affected, prop) ;
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
    double factor = exp( Ea*(1./T-1./T0) ) ;
    double tau = defaultValues["creep_characteristic_time"] ;
    s.set("creep_characteristic_time", tau*factor) ;
/*    if(s.has("creep_modulus"))
    {
        double E = defaultValues["creep_modulus"] ;
        s.set("creep_modulus", E*factor) ;
    } else {
        double k = defaultValues["creep_bulk"] ;
        s.set("creep_bulk", k*factor) ;
        double mu = defaultValues["creep_shear"] ;
        s.set("creep_shear", mu*factor) ;
    }*/
}

void CreepRelativeHumidityMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_characteristic_time"))
        return ;

    if(defaultValues.find("creep_modulus") == defaultValues.end() && s.has("creep_modulus"))
	defaultValues["creep_modulus"] = s.get("creep_modulus", defaultValues) ;
    else
    {
	if(defaultValues.find("creep_bulk") == defaultValues.end())
		defaultValues["creep_bulk"] = s.get("creep_bulk", defaultValues) ;
	if(defaultValues.find("creep_shear") == defaultValues.end())
		defaultValues["creep_shear"] = s.get("creep_shear", defaultValues) ;
    }

    if(defaultValues.find("recoverable_modulus") == defaultValues.end() && s.has("recoverable_modulus"))
	defaultValues["recoverable_modulus"] = s.get("recoverable_modulus", defaultValues) ;
    else
    {
	if(defaultValues.find("recoverable_bulk") == defaultValues.end())
		defaultValues["recoverable_bulk"] = s.get("recoverable_bulk", defaultValues) ;
	if(defaultValues.find("recoverable_shear") == defaultValues.end())
		defaultValues["recoverable_shear"] = s.get("recoverable_shear", defaultValues) ;
    }

    double h = s.get("relative_humidity", defaultValues) ;
    double hc = s.get("creep_humidity_coefficient", defaultValues) ;
    double factor =  hc*(1-h) + exp(-hc*(1-h)) ;
    if(s.has("creep_modulus"))
    {
        double E = defaultValues["creep_modulus"] ;
        s.set("creep_modulus", E*factor) ;
    } else {
        double k = defaultValues["creep_bulk"] ;
        s.set("creep_bulk", k*factor) ;
        double mu = defaultValues["creep_shear"] ;
        s.set("creep_shear", mu*factor) ;
    }
    if(s.has("recoverable_modulus"))
    {
        double E = defaultValues["recoverable_modulus"] ;
        s.set("recoverable_modulus", E*factor) ;
    } else {
        double k = defaultValues["recoverable_bulk"] ;
        s.set("recoverable_bulk", k*factor) ;
        double mu = defaultValues["recoverable_shear"] ;
        s.set("recoverable_shear", mu*factor) ;
    }

}


} ;

