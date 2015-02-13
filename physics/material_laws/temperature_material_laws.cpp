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


} ;

