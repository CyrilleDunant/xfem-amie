// C++ Interface: logarithmic_creep
//
// Description: Logarithmic visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __HUMIDITY_MATERIAL_LAW_H_
#define __HUMIDITY_MATERIAL_LAW_H_

#include "material_laws.h"

namespace Amie
{

struct ThermalExpansionHumidityMaterialLaw : public ExternalMaterialLaw
{
    ThermalExpansionHumidityMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { defaultValues["water_molar_volume"] = 1.8016e-5 ; }
    virtual ~ThermalExpansionHumidityMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct DryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    double order ;
    std::string input ;

    DryingShrinkageMaterialLaw(std::string in = "relative_humidity", double o = 1, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), order(o) , input(in)
    { 
	defaultValues["drying_shrinkage_irreversibility_threshold"] = -1. ;
	defaultValues["drying_shrinkage_irreversibility_factor"] = 0. ;
    }
    virtual ~DryingShrinkageMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct CapillaryPressureDryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    CapillaryPressureDryingShrinkageMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) {  }
    virtual ~CapillaryPressureDryingShrinkageMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

} ;

struct KelvinCapillaryPressureMaterialLaw : public ExternalMaterialLaw
{
    KelvinCapillaryPressureMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { 
        defaultValues["water_molar_volume"] = 1.8016e-2 ; 
        defaultValues["water_density"] = 1e3 ; 
    }
    virtual ~KelvinCapillaryPressureMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct VanGenuchtenCapillaryPressureMaterialLaw : public ExternalMaterialLaw
{
    VanGenuchtenCapillaryPressureMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~VanGenuchtenCapillaryPressureMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct VanGenuchtenWaterSaturationMaterialLaw : public ExternalMaterialLaw
{
    VanGenuchtenWaterSaturationMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~VanGenuchtenWaterSaturationMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct DisjoiningPressureDryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    DisjoiningPressureDryingShrinkageMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~DisjoiningPressureDryingShrinkageMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct BETIsothermMaterialLaw : public ExternalMaterialLaw
{
	std::string suffix ;

	BETIsothermMaterialLaw(std::string suf = std::string(), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), suffix(suf) { } ;
	virtual ~BETIsothermMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct BiExponentialIsothermMaterialLaw : public BETIsothermMaterialLaw
{
	BiExponentialIsothermMaterialLaw(std::string suf = std::string(), std::string args = std::string(), char sep = ',') : BETIsothermMaterialLaw(suf, args, sep) { } ;
	virtual ~BiExponentialIsothermMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct SorptionIsothermHysteresisMaterialLaw : public ExternalMaterialLaw
{
	ExternalMaterialLaw * desorption ;
	ExternalMaterialLaw * adsorption ;

	SorptionIsothermHysteresisMaterialLaw(ExternalMaterialLaw * d, ExternalMaterialLaw * a, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), desorption(d), adsorption(a) { } ;
	virtual ~SorptionIsothermHysteresisMaterialLaw() { if(desorption){ delete desorption ; delete adsorption ; } } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct BenboudjemaDryingCreepMaterialLaw : public ExternalMaterialLaw
{
    BenboudjemaDryingCreepMaterialLaw(std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep) { }
    virtual ~BenboudjemaDryingCreepMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};

struct HavlasekDryingCreepMaterialLaw : public ExternalMaterialLaw
{
    HavlasekDryingCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~HavlasekDryingCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};


struct WittmannCreepRelativeHumidityMaterialLaw : public ExternalMaterialLaw
{
    WittmannCreepRelativeHumidityMaterialLaw(std::string args = std::string("creep_humidity_coefficient = 0.2"), char sep = 'c') : ExternalMaterialLaw(args, sep) { }
    virtual ~WittmannCreepRelativeHumidityMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};


} 


#endif // __TEMPERATURE_MATERIAL_LAW_H_
