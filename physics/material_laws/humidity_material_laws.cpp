#include "humidity_material_laws.h"
#include <stdlib.h>

namespace Amie
{

void DryingShrinkageMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    double h = s.get(input,defaultValues) ;
    double hmin = h ;
    double h0 = s.get(input+"_reference", defaultValues) ;
    if(!s.has(input+"_minimum"))
	s.set(input+"_minimum", h) ;
    hmin = s.get(input+"_minimum", defaultValues) ;
    if(h < hmin)
    {
	s.set(input+"_minimum", h) ;
	hmin = h ;
    }
    double ks = s.get("drying_shrinkage_coefficient",defaultValues) ;
    double thr = s.get("drying_shrinkage_irreversibility_threshold", defaultValues) ;
    double irr = s.get("drying_shrinkage_irreversibility_coefficient", defaultValues) ;
    double dh = 0. ;
    if(hmin < thr)
	dh = hmin-thr ;
    if(h < thr)
	dh += thr-h0 ;
    else
	dh += (h-h0)*(1.-irr)+(thr-h0)*irr ;
    s.add("imposed_deformation", -ks*(std::pow(dh, order)) ) ;
}

void KelvinCapillaryPressureMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double rw = 1e6 ; // water density
	double R = 8.3144621 ;
	double mw = s.get("water_molar_volume", defaultValues) ;
	double T = s.get("temperature", defaultValues) ;
	double h = s.get("relative_humidity", defaultValues) ;
	s.set("capillary_pressure", -(rw*R*T/mw)*std::log(h) ) ;
}

void BaroghelBounyWaterSaturationMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double pc = s.get("capillary_pressure", defaultValues) ;
	double a = s.get("capillary_pressure_coefficient", defaultValues) ;
	double b = s.get("capillary_pressure_exponent", defaultValues) ;
	s.set("water_saturation", a*std::pow( 1. + std::pow(pc/a, b/(b-1)), -1./b)  ) ;
}

void CapillaryPressureDryingShrinkageMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double pc = s.get("capillary_pressure", defaultValues) ;
	double k = s.get("bulk_modulus", defaultValues) ;
	double b = s.get("biot_coefficient", defaultValues) ;
	double xhi = 1. ;
	if(s.get("relative_humidity", defaultValues) > 0)
		xhi = s.get("water_saturation", defaultValues)/s.get("relative_humidity", defaultValues) ;
	s.add("imposed_deformation", -b/(3.*k)*xhi*pc) ;
}

void DisjoiningPressureDryingShrinkageMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double mw = s.get("mass_water_content", defaultValues) ;
	double ssa = s.get("specific_surface_area", defaultValues) ;
	double ap = s.get("asymptotic_disjoining_pressure", defaultValues) ;
	if(!s.has("disjoining_pressure_reference"))
	{
		double mw0 = s.get("mass_water_content_reference", defaultValues) ;
		s.set("disjoining_pressure_reference", ap*exp(-2.*mw0*1e3/(0.95*ssa))) ;
	}
	double p0 = s.get("disjoining_pressure_reference", defaultValues) ;
	double p = ap*exp(-2.*mw*1e3/(0.95*ssa)) ;
	if(p < p0)
		return ;
	double k = s.get("bulk_modulus", defaultValues) ;
	double w = s.get("water_fraction_reference", defaultValues) ;
	s.add("imposed_deformation", w/(3.*k)*(p0-p) ) ;
}

void BETIsothermMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double C = s.get("bet_adsorption_enthalpy_constant"+suffix, defaultValues) ;
	double k = s.get("bet_humidity_correction_factor"+suffix, defaultValues) ;
	double w = s.get("bet_water_layer_mass"+suffix, defaultValues) ;
	double h = s.get("relative_humidity", defaultValues) ;

	s.set("mass_water_content"+suffix, C*w*k*h/( (1.-k*h) * (1.+(C-1.)*k*h) ) ) ;
}

void BiExponentialIsothermMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double dry = s.get("dry_exponential_constant"+suffix, defaultValues) ;
	double sat = s.get("saturated_exponential_constant"+suffix, defaultValues) ;
	double dh = s.get("humidity_exponential_constant"+suffix, defaultValues) ;
	double h = s.get("relative_humidity", defaultValues) ;

	s.set("mass_water_content"+suffix, dry*(1.-exp(-h/dh))+sat*(exp(h/dh)-1.)  ) ;
}

void SorptionIsothermHysteresisMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double h = s.get("relative_humidity", defaultValues) ;
	bool first = false ;
	if(!s.has("relative_humidity_previous"))
	{
		s.set("relative_humidity_previous", h) ;
		first = true ;
	}
	double hprev = s.get("relative_humidity_previous", defaultValues) ;

	// constant humidity, do nothing
	if(!first && std::abs(h-hprev) < POINT_TOLERANCE)
		return ;

	if(!s.has("bet_adsorption_desorption_shift"))
		s.set("bet_adsorption_desorption_shift", 0.) ;
	if(!s.has("bet_adsorption_desorption_mode"))
		s.set("bet_adsorption_desorption_mode", -1.) ;

	desorption->preProcess(s, dt) ;
	adsorption->preProcess(s, dt) ;

	double shift = s.get("bet_adsorption_desorption_shift", defaultValues) ;
	double mode = s.get("bet_adsorption_desorption_mode", defaultValues) ;

	double w = 0. ;
	double wcd = s.get("mass_water_content_desorption", defaultValues) ;
	double wca = s.get("mass_water_content_adsorption", defaultValues) ;
	if(first)
	{
		s.set("mass_water_content_desorption_previous", wcd) ;
		s.set("mass_water_content_adsorption_previous", wca) ;
	}
	double wcdprev = s.get("mass_water_content_desorption_previous", defaultValues) ;
	double wcaprev = s.get("mass_water_content_adsorption_previous", defaultValues)  ;
	if(h-hprev < 0) // desorption mode
	{
		if(mode > 0) // mode change
			shift = wcdprev-wcaprev-shift ;

		w = std::min(std::max(wca, wcd-shift), wcd) ;
		s.set("bet_adsorption_desorption_mode", -1.) ;
	} 
	else
	{
		if(mode < 0) // mode change
			shift = wcdprev-wcaprev-shift ;

		w = std::max(std::min(wcd, wca+shift),wca) ;
		s.set("bet_adsorption_desorption_mode", 1.) ;
	}

	s.set("mass_water_content", w) ;
	if(!s.has("mass_water_content_reference"))
		s.set("mass_water_content_reference", w) ;
	s.set("water_saturation", w/s.get("mass_water_content_reference", defaultValues)) ;
	s.set("bet_adsorption_desorption_shift", shift) ;
	s.set("relative_humidity_previous", h) ;
	s.set("mass_water_content_desorption_previous", wcd) ;
	s.set("mass_water_content_adsorption_previous", wca) ;
}

void ThermalExpansionHumidityMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
	double R = 8.3144621 ;
	double m = s.get("water_molar_volume", defaultValues) ;
	double k = s.get("bulk_modulus", defaultValues) ;
	double ks = s.get("bulk_modulus_solid_skeleton", defaultValues) ;
	double h = s.get("relative_humidity", defaultValues) ;

	s.add("thermal_expansion_coefficient", -(h*std::log(h)*(R/m)*(1./(3.*k) - 1./ks)) ) ;
}


void CreepRelativeHumidityMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_characteristic_time"))
        return ;

    double h = s.get("relative_humidity", defaultValues) ;
    double hc = s.get("creep_humidity_coefficient", defaultValues) ;
    double factor =  (1-h)/hc + exp(-(1-h)/hc) ;
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

