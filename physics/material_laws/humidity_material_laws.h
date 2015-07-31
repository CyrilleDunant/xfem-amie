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

/*PARSE HumidityDependentThermalExpansionCoefficient ExternalMaterialLaw */
struct HumidityDependentThermalExpansionCoefficientMaterialLaw : public ExternalMaterialLaw
{
    HumidityDependentThermalExpansionCoefficientMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { defaultValues["water_molar_volume"] = 1.8016e-5 ; }
    virtual ~HumidityDependentThermalExpansionCoefficientMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE DryingShrinkage ExternalMaterialLaw
    @string[parameter] relative_humidity // name of the parameter to use to calculate the drying shrinkage
 */
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

/*PARSE CapillaryPressureDrivenDryingShrinkage ExternalMaterialLaw */
struct CapillaryPressureDrivenDryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    CapillaryPressureDrivenDryingShrinkageMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) {  }
    virtual ~CapillaryPressureDrivenDryingShrinkageMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

} ;

/*PARSE KelvinCapillaryPressure ExternalMaterialLaw */
struct KelvinCapillaryPressureMaterialLaw : public ExternalMaterialLaw
{
    KelvinCapillaryPressureMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { 
        defaultValues["water_molar_volume"] = 1.8016e-2 ; 
        defaultValues["water_density"] = 1e3 ; 
    }
    virtual ~KelvinCapillaryPressureMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE VanGenuchtenCapillaryPressure ExternalMaterialLaw */
struct VanGenuchtenCapillaryPressureMaterialLaw : public ExternalMaterialLaw
{
    VanGenuchtenCapillaryPressureMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~VanGenuchtenCapillaryPressureMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE VanGenuchtenWaterSaturation ExternalMaterialLaw */
struct VanGenuchtenWaterSaturationMaterialLaw : public ExternalMaterialLaw
{
    VanGenuchtenWaterSaturationMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~VanGenuchtenWaterSaturationMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE DisjoiningPressureDrivenDryingShrinkage ExternalMaterialLaw */
struct DisjoiningPressureDrivenDryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    DisjoiningPressureDrivenDryingShrinkageMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~DisjoiningPressureDrivenDryingShrinkageMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE BETIsotherm ExternalMaterialLaw */
struct BETIsothermMaterialLaw : public ExternalMaterialLaw
{
	std::string suffix ;

	BETIsothermMaterialLaw(std::string suf = std::string(), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), suffix(suf) { } ;
	virtual ~BETIsothermMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE BiExponentialIsotherm ExternalMaterialLaw */
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

/*PARSE WaterVapourSaturationPressure ExternalMaterialLaw */
struct WaterVapourSaturationPressureMaterialLaw : public ExternalMaterialLaw
{
    std::string suffix ;

    WaterVapourSaturationPressureMaterialLaw(std::string s = std::string(), std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), suffix(s) { }
    virtual ~WaterVapourSaturationPressureMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};

struct ClausiusClapeyronRelativeHumidityMaterialLaw : public ExternalMaterialLaw
{
    ExternalMaterialLaw * reference ;
    ExternalMaterialLaw * vapour ;

    ClausiusClapeyronRelativeHumidityMaterialLaw(ExternalMaterialLaw * bet, std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), reference(bet) { 
        vapour = new WaterVapourSaturationPressureMaterialLaw("_T1") ;
    }
    virtual ~ClausiusClapeyronRelativeHumidityMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};



/*PARSE BenboudjemaDryingCreep ExternalMaterialLaw */
struct BenboudjemaDryingCreepMaterialLaw : public ExternalMaterialLaw
{
    BenboudjemaDryingCreepMaterialLaw(std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep) { }
    virtual ~BenboudjemaDryingCreepMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};

/*PARSE BazantRelativeHumidityDependentCreep ExternalMaterialLaw */
struct BazantRelativeHumidityDependentCreepMaterialLaw : public ExternalMaterialLaw
{
    BazantRelativeHumidityDependentCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~BazantRelativeHumidityDependentCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE HavlasekDryingCreep ExternalMaterialLaw */
struct HavlasekDryingCreepMaterialLaw : public ExternalMaterialLaw
{
    HavlasekDryingCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~HavlasekDryingCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};


/*PARSE WittmannRelativeHumidityDependentCreep ExternalMaterialLaw */
struct WittmannRelativeHumidityDependentCreepMaterialLaw : public ExternalMaterialLaw
{
    WittmannRelativeHumidityDependentCreepMaterialLaw(std::string args = std::string("creep_humidity_coefficient = 0.2"), char sep = 'c') : ExternalMaterialLaw(args, sep) { }
    virtual ~WittmannRelativeHumidityDependentCreepMaterialLaw() { } ;

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};


} 


#endif // __TEMPERATURE_MATERIAL_LAW_H_
