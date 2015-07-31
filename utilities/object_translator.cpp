/* this is an auto-generated file created on 31/6/2015 at 17:27  */

#include "object_translator.h"
#include "enumeration_translator.h"
#include "../physics/material_laws/humidity_material_laws.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/mechanical_material_laws.h"
#include "../physics/material_laws/logcreep_accumulator.h"
#include "../physics/material_laws/temperature_material_laws.h"

namespace Amie
{

    ExternalMaterialLaw * Object::getExternalMaterialLaw(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, std::vector<std::string>> & stringlists, std::map<std::string, double> & values)
    {
        // parsed from header file: ../physics/material_laws/humidity_material_laws.h
        if( type == "HumidityDependentThermalExpansionCoefficient" ) { return new HumidityDependentThermalExpansionCoefficientMaterialLaw() ; }
        if( type == "DryingShrinkage" )
        { 
            if( strings.find("parameter") == strings.end() )
                strings["parameter"] = "relative_humidity" ;
            return new DryingShrinkageMaterialLaw(strings["parameter"]) ;
        }
        if( type == "CapillaryPressureDrivenDryingShrinkage" ) { return new CapillaryPressureDrivenDryingShrinkageMaterialLaw() ; }
        if( type == "KelvinCapillaryPressure" ) { return new KelvinCapillaryPressureMaterialLaw() ; }
        if( type == "VanGenuchtenCapillaryPressure" ) { return new VanGenuchtenCapillaryPressureMaterialLaw() ; }
        if( type == "VanGenuchtenWaterSaturation" ) { return new VanGenuchtenWaterSaturationMaterialLaw() ; }
        if( type == "DisjoiningPressureDrivenDryingShrinkage" ) { return new DisjoiningPressureDrivenDryingShrinkageMaterialLaw() ; }
        if( type == "BETIsotherm" ) { return new BETIsothermMaterialLaw() ; }
        if( type == "BiExponentialIsotherm" ) { return new BiExponentialIsothermMaterialLaw() ; }
        if( type == "WaterVapourSaturationPressure" ) { return new WaterVapourSaturationPressureMaterialLaw() ; }
        if( type == "BenboudjemaDryingCreep" ) { return new BenboudjemaDryingCreepMaterialLaw() ; }
        if( type == "BazantRelativeHumidityDependentCreep" ) { return new BazantRelativeHumidityDependentCreepMaterialLaw() ; }
        if( type == "HavlasekDryingCreep" ) { return new HavlasekDryingCreepMaterialLaw() ; }
        if( type == "WittmannRelativeHumidityDependentCreep" ) { return new WittmannRelativeHumidityDependentCreepMaterialLaw() ; }
   
        // parsed from header file: ../physics/material_laws/material_laws.h
        if( type == "TimeDerivative" ) { return new TimeDerivativeMaterialLaw(strings["parameter"]) ; }
        if( type == "TimeIntegral" ) { return new TimeIntegralMaterialLaw(strings["parameter"]) ; }
        if( type == "Minimum" )
        { 
            if( strings.find("operation") == strings.end() )
                strings["operation"] = "SET" ;
            return new MinimumMaterialLaw(strings["output"], stringlists["parameters"], Enum::getEMLOperation(strings["operation"])) ;
        }
        if( type == "StoreMaximumValue" )
        { 
            if( values.find("starting_value") == values.end() )
                values["starting_value"] = 0 ;
            return new StoreMaximumValueMaterialLaw(strings["parameter"], values["starting_value"]) ;
        }
        if( type == "GetParticleOrientation" )
        { 
            if( strings.find("output") == strings.end() )
                strings["output"] = "angle" ;
            return new GetParticleOrientationMaterialLaw(strings["output"]) ;
        }
        if( type == "StorePreviousValue" ) { return new StorePreviousValueMaterialLaw(strings["parameter"]) ; }
        if( type == "StoreMinimumValue" )
        { 
            if( values.find("starting_value") == values.end() )
                values["starting_value"] = 0 ;
            return new StoreMinimumValueMaterialLaw(strings["parameter"], values["starting_value"]) ;
        }
        if( type == "Maximum" )
        { 
            if( strings.find("operation") == strings.end() )
                strings["operation"] = "SET" ;
            return new MaximumMaterialLaw(strings["output"], stringlists["parameters"], Enum::getEMLOperation(strings["operation"])) ;
        }
        if( type == "ExponentialDecay" ) { return new ExponentialDecayMaterialLaw(strings["output"], strings["target"]) ; }
        if( type == "GetField" ) { return new GetFieldMaterialLaw(strings["field"]) ; }
        if( type == "WeibullDistributed" )
        { 
            if( strings.find("weibull_variable") == strings.end() )
                strings["weibull_variable"] = "weibull_variable" ;
            return new WeibullDistributedMaterialLaw(stringlists["parameters"], strings["weibull_variable"]) ;
        }
   
        // parsed from header file: ../physics/material_laws/mechanical_material_laws.h
        if( type == "BulkShearConversion" ) { return new BulkShearConversionMaterialLaw() ; }
        if( type == "AdjustStrainStressCurve" ) { return new AdjustStrainStressCurveMaterialLaw() ; }
        if( type == "BazantLoadNonLinearCreep" ) { return new BazantLoadNonLinearCreepMaterialLaw() ; }
        if( type == "TensionCompressionCreep" ) { return new TensionCompressionCreepMaterialLaw() ; }
   
        // parsed from header file: ../physics/material_laws/temperature_material_laws.h
        if( type == "ThermalExpansion" ) { return new ThermalExpansionMaterialLaw() ; }
        if( type == "RadiationDependentThermalExpansionCoefficient" ) { return new RadiationDependentThermalExpansionCoefficientMaterialLaw() ; }
        if( type == "RadiationDependentPoissonRatio" ) { return new RadiationDependentPoissonRatioMaterialLaw() ; }
        if( type == "IncrementalThermalExpansion" ) { return new IncrementalThermalExpansionMaterialLaw() ; }
        if( type == "RadiationInducedVolumetricExpansion" ) { return new RadiationInducedVolumetricExpansionMaterialLaw() ; }
        if( type == "TemperatureDependentRadiationInducedVolumetricExpansion" ) { return new TemperatureDependentRadiationInducedVolumetricExpansionMaterialLaw() ; }
        if( type == "Arrhenius" ) { return new ArrheniusMaterialLaw(strings["parameter"]) ; }
        if( type == "CreepArrhenius" ) { return new CreepArrheniusMaterialLaw() ; }
   
        return nullptr ;
    }

    void Object::resetExternalMaterialLaw(ExternalMaterialLaw * target)
    {
        // parsed from header file: ../physics/material_laws/humidity_material_laws.h
   
        // parsed from header file: ../physics/material_laws/material_laws.h
   
        // parsed from header file: ../physics/material_laws/mechanical_material_laws.h
   
        // parsed from header file: ../physics/material_laws/temperature_material_laws.h
   
    }

}

