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

void AnisotropicThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    if(!s.has("thermal_expansion_coefficient_xx"))
    {
        if(s.has("thermal_expansion_coefficient_11") && s.has("angle"))
        {
            Vector alpha(3) ;
            alpha[0] = s.get("thermal_expansion_coefficient_11", defaultValues) ;
            alpha[1] = s.get("thermal_expansion_coefficient_22", defaultValues) ;
            alpha = Tensor::rotate2ndOrderTensor2D( alpha, s.get("angle", defaultValues) ) ;
            alpha *= (T-T0) ;
            s.add("imposed_deformation_xx", alpha[0]*(T-T0)) ;
            s.add("imposed_deformation_yy", alpha[1]*(T-T0)) ;
            return ;
        }
        else
        {
            double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
            s.add("imposed_deformation", alpha*(T-T0)) ;
            return ;
        }
    }
    double alpha_xx = s.get("thermal_expansion_coefficient_xx", defaultValues) ;
    s.add("imposed_deformation_xx", alpha_xx*(T-T0)) ;
    double alpha_yy = s.get("thermal_expansion_coefficient_yy", defaultValues) ;
    s.add("imposed_deformation_yy", alpha_yy*(T-T0)) ;
}

void IncrementalThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("thermal_deformation"))
        s.set("thermal_deformation", 0.) ;

    double T = s.get("temperature", defaultValues) ;
    if(!s.has("temperature_previous"))
        s.set("temperature_previous", T) ;
    double T0 = s.get("temperature_previous", defaultValues) ;
    double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
    s.add("thermal_deformation", alpha*(T-T0)) ;
    s.add("imposed_deformation", s.get("thermal_deformation", defaultValues) ) ;
}

void AnisotropicIncrementalThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("thermal_deformation_xx"))
    {
        s.set("thermal_deformation_xx", 0.) ;
        s.set("thermal_deformation_yy", 0.) ;
    }

    double T = s.get("temperature", defaultValues) ;
    if(!s.has("temperature_previous"))
        s.set("temperature_previous", T) ;
    double T0 = s.get("temperature_previous", defaultValues) ;

    if(!s.has("thermal_expansion_coefficient_xx"))
    {
        if(s.has("thermal_expansion_coefficient_11") && s.has("angle"))
        {
            Vector alpha(3) ;
            alpha[0] = s.get("thermal_expansion_coefficient_11", defaultValues) ;
            alpha[1] = s.get("thermal_expansion_coefficient_22", defaultValues) ;
            alpha = Tensor::rotate2ndOrderTensor2D( alpha, s.get("angle", defaultValues) ) ;
            alpha *= (T-T0) ;
            s.add("thermal_deformation_xx", alpha[0]*(T-T0)) ;
            s.add("thermal_deformation_yy", alpha[1]*(T-T0)) ;
        }
        else
        {
            double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
            s.add("thermal_deformation_xx", alpha*(T-T0)) ;
            s.add("thermal_deformation_yy", alpha*(T-T0)) ;
        }
    }
    else
    {
        double alpha_xx = s.get("thermal_expansion_coefficient_xx", defaultValues) ;
        s.add("thermal_deformation_xx", alpha_xx*(T-T0)) ;
        double alpha_yy = s.get("thermal_expansion_coefficient_yy", defaultValues) ;
        s.add("thermal_deformation_yy", alpha_yy*(T-T0)) ;
    }

    s.add("imposed_deformation_xx", s.get("thermal_deformation_xx", defaultValues) ) ;
    s.add("imposed_deformation_yy", s.get("thermal_deformation_yy", defaultValues) ) ;
}

void RadiationInducedExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double kappa = s.get("radiation_expansion_delay",defaultValues) ;
    double epsmax = s.get("maximum_radiation_expansion",defaultValues) ;
    double delta = s.get("neutron_fluence_correction",defaultValues) ;
    double N = s.get("neutron_fluence",defaultValues) ;
    s.add("imposed_deformation", (kappa*epsmax*(exp(delta*N)-1)/(epsmax+kappa*exp(delta*N)))) ;
}

void TemperatureDependentRadiationInducedExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("rive_previous"))
        s.set("rive_previous", 0.) ;
    if(s.has("neutron_flux") && !s.has("neutron_fluence"))
        s.set("neutron_fluence", 0.) ;

    double eps_max = s.get("maximum_radiation_expansion",defaultValues) ;
    double eps_prev = s.get("rive_previous", defaultValues) ;
    if(eps_prev > eps_max-POINT_TOLERANCE)
    {
        s.add("imposed_deformation", eps_max) ;
        return ;
    }

    double ac = s.get("characteristic_fluence_coefficient", defaultValues) ;
    double bc = s.get("characteristic_fluence_correction", defaultValues) ;
    double al = s.get("latency_fluence_coefficient", defaultValues) ;
    double bl = s.get("latency_fluence_correction", defaultValues) ;
    double alpha = s.get("rive_integration_coefficient",defaultValues) ;

    double t_next = s.get("temperature", defaultValues)-273 ;
    if(!s.has("temperature_previous"))
        s.set("temperature_previous", t_next+273) ;
    double t_prev = s.get("temperature_previous", defaultValues)-273 ;
    double phi_prev = s.get("neutron_fluence_previous", defaultValues) ;
    double phi_next = phi_prev ;
    if(s.has("neutron_flux"))
    {
        phi_next += s.get("neutron_flux", defaultValues)*dt ;
        s.set("neutron_fluence", phi_next) ;
    }
    else
        phi_next = s.get("neutron_fluence", defaultValues) ;

    double delta_prev = 1./(ac*t_prev+bc) ;
    double delta_next = 1./(ac*t_next+bc) ;
    double kappa_prev = eps_max/exp((std::max(0.,al*t_prev+bl))/(ac*t_prev+bc)) ;
    double kappa_next = eps_max/exp((std::max(0.,al*t_next+bl))/(ac*t_next+bc)) ;

    double phieq_prev = (1./delta_prev) * std::log(  eps_max*(eps_prev + kappa_prev)/(kappa_prev*(eps_max-eps_prev)) ) ;
    double phieq_next = (1./delta_next) * std::log(  eps_max*(eps_prev + kappa_next)/(kappa_next*(eps_max-eps_prev)) ) ;

    double delta_eps_prev = kappa_prev*eps_max*delta_prev*exp( delta_prev*phieq_prev )*(eps_max+kappa_prev)/ ((eps_max + kappa_prev*exp(delta_prev*phieq_prev ))*(eps_max + kappa_prev*exp(delta_prev*phieq_prev ))) ;
    double delta_eps_next = kappa_next*eps_max*delta_next*exp( delta_next*phieq_next )*(eps_max+kappa_next)/ ((eps_max + kappa_next*exp(delta_next*phieq_next ))*(eps_max + kappa_next*exp(delta_next*phieq_next ))) ;

    double eps_next = std::min( eps_max, eps_prev + ( ( alpha*delta_eps_prev + (1.-alpha)*delta_eps_next ) * (phi_next - phi_prev) ) ) ;

    s.add("imposed_deformation", eps_next) ;
    s.set("rive_previous", eps_next) ;
}

ArrheniusMaterialLaw::ArrheniusMaterialLaw(std::string a, std::string args, char sep) : ExternalMaterialLaw(args, sep), affected(a)
{
    coefficient = affected ;
    coefficient.append("_activation_energy") ;
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


}

