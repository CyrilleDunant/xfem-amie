/* this is an auto-generated file created on 21/0/2016 at 12:25  */

#include "object_translator.h"
#include "enumeration_translator.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/material_laws/mechanical_material_laws.h"
#include "../physics/material_laws/humidity_material_laws.h"
#include "../physics/void_form.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness_with_imposed_stress.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/damagemodels/isotropiclineardamage.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/confinedmohrcoulombwithstrain.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/fracturecriteria/confinedvonmises.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/mazars.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/fracturecriteria/spacetimeflowrule.h"
#include "../physics/fracturecriteria/confinedmohrcoulomb.h"
#include "../physics/fracturecriteria/limitstrains.h"
#include "../physics/material_laws/logcreep_accumulator.h"
#include "../features/microstructuregenerator.h"
#include "../utilities/granulo.h"
#include "../utilities/inclusion_family.h"
#include "../features/enrichmentmanagers/gelmanager.h"
#include "../geometry/sampler/regular_sampler.h"
#include "../geometry/sampler/gradient_sampler.h"
#include "../geometry/sampler/sampler.h"

namespace Amie
{

    ExternalMaterialLaw * Object::getExternalMaterialLaw(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, std::vector<std::string>> & stringlists, std::map<std::string, double> & values)
    {
        // parsed from header file: ../physics/material_laws/material_laws.h
        if( type == "Eval" )
        { 
            if( strings.find("operation") == strings.end() ) { strings["operation"] = "SET" ; } ; 
            return new EvalMaterialLaw(strings["output"], strings["function"], Enum::getEMLOperation(strings["operation"])) ;
        }
        if( type == "LinearInterpolated" ) { return new LinearInterpolatedMaterialLaw(strings["output"], strings["input"], strings["file_name"], Enum::getEMLOperation(strings["operation"])) ; }
        if( type == "TimeDerivative" ) { return new TimeDerivativeMaterialLaw(strings["parameter"]) ; }
        if( type == "TimeIntegral" ) { return new TimeIntegralMaterialLaw(strings["parameter"]) ; }
        if( type == "Minimum" )
        { 
            if( strings.find("operation") == strings.end() ) { strings["operation"] = "SET" ; } ; 
            return new MinimumMaterialLaw(strings["output"], stringlists["parameters"], Enum::getEMLOperation(strings["operation"])) ;
        }
        if( type == "StoreMaximumValue" )
        { 
            if( values.find("starting_value") == values.end() ) { values["starting_value"] = 0 ; } ; 
            return new StoreMaximumValueMaterialLaw(strings["parameter"], values["starting_value"]) ;
        }
        if( type == "GetParticleOrientation" )
        { 
            if( strings.find("output") == strings.end() ) { strings["output"] = "angle" ; } ; 
            return new GetParticleOrientationMaterialLaw(strings["output"]) ;
        }
        if( type == "StorePreviousValue" ) { return new StorePreviousValueMaterialLaw(strings["parameter"]) ; }
        if( type == "StoreMinimumValue" )
        { 
            if( values.find("starting_value") == values.end() ) { values["starting_value"] = 0 ; } ; 
            return new StoreMinimumValueMaterialLaw(strings["parameter"], values["starting_value"]) ;
        }
        if( type == "Maximum" )
        { 
            if( strings.find("operation") == strings.end() ) { strings["operation"] = "SET" ; } ; 
            return new MaximumMaterialLaw(strings["output"], stringlists["parameters"], Enum::getEMLOperation(strings["operation"])) ;
        }
        if( type == "ExponentialDecay" ) { return new ExponentialDecayMaterialLaw(strings["output"], strings["parameter"]) ; }
        if( type == "GetField" ) { return new GetFieldMaterialLaw(strings["parameter"]) ; }
        if( type == "WeibullDistributed" )
        { 
            if( strings.find("weibull_variable") == strings.end() ) { strings["weibull_variable"] = "weibull_variable" ; } ; 
            return new WeibullDistributedMaterialLaw(stringlists["parameters"], strings["weibull_variable"]) ;
        }
        if( type == "UniformDistributedPerParticle" )
        { 
            if( values.find("minimum") == values.end() ) { values["minimum"] = 0 ; } ; 
            if( values.find("maximum") == values.end() ) { values["maximum"] = 1 ; } ; 
            if( strings.find("operation") == strings.end() ) { strings["operation"] = "SET" ; } ; 
            return new UniformDistributedPerParticleMaterialLaw(strings["output"], values["minimum"], values["maximum"], Enum::getEMLOperation(strings["operation"])) ;
        }
   
        // parsed from header file: ../physics/material_laws/temperature_material_laws.h
        if( type == "ThermalExpansion" ) { return new ThermalExpansionMaterialLaw() ; }
        if( type == "RadiationDependentThermalExpansionCoefficient" ) { return new RadiationDependentThermalExpansionCoefficientMaterialLaw() ; }
        if( type == "RadiationDependentPoissonRatio" ) { return new RadiationDependentPoissonRatioMaterialLaw() ; }
        if( type == "IncrementalThermalExpansion" ) { return new IncrementalThermalExpansionMaterialLaw() ; }
        if( type == "RadiationInducedVolumetricExpansion" ) { return new RadiationInducedVolumetricExpansionMaterialLaw() ; }
        if( type == "TemperatureDependentRadiationInducedVolumetricExpansion" ) { return new TemperatureDependentRadiationInducedVolumetricExpansionMaterialLaw() ; }
        if( type == "Arrhenius" ) { return new ArrheniusMaterialLaw(strings["parameter"]) ; }
        if( type == "CreepArrhenius" ) { return new CreepArrheniusMaterialLaw() ; }
   
        // parsed from header file: ../physics/material_laws/mechanical_material_laws.h
        if( type == "BulkShearConversion" ) { return new BulkShearConversionMaterialLaw() ; }
        if( type == "AdjustStrainStressCurve" ) { return new AdjustStrainStressCurveMaterialLaw() ; }
        if( type == "BazantLoadNonLinearCreep" ) { return new BazantLoadNonLinearCreepMaterialLaw() ; }
        if( type == "TensionCompressionCreep" ) { return new TensionCompressionCreepMaterialLaw() ; }
        if( type == "Mineral" )
        { 
            if( strings.find("separators") == strings.end() ) { strings["separators"] = ".-" ; } ; 
            if( values.find("index") == values.end() ) { values["index"] = -1 ; } ; 
            if( values.find("factor") == values.end() ) { values["factor"] = 1 ; } ; 
            if( strings.find("force") == strings.end() ) { strings["force"] = "false" ; } ; 
            if( strings.find("cutting_plane") == strings.end() ) { strings["cutting_plane"] = "ZETA" ; } ; 
            return new MineralMaterialLaw(strings["file_name"], strings["separators"], values["index"], values["factor"], Enum::getbool(strings["force"]), Enum::getVariable(strings["cutting_plane"])) ;
        }
   
        // parsed from header file: ../physics/material_laws/humidity_material_laws.h
        if( type == "HumidityDependentThermalExpansionCoefficient" ) { return new HumidityDependentThermalExpansionCoefficientMaterialLaw() ; }
        if( type == "DryingShrinkage" )
        { 
            if( strings.find("parameter") == strings.end() ) { strings["parameter"] = "relative_humidity" ; } ; 
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
   
        return nullptr ;
    }

    bool Object::isExternalMaterialLaw(std::string type)
    {
        // parsed from header file: ../physics/material_laws/material_laws.h
        if( type == "Eval" ) { return true ; }
        if( type == "LinearInterpolated" ) { return true ; }
        if( type == "TimeDerivative" ) { return true ; }
        if( type == "TimeIntegral" ) { return true ; }
        if( type == "Minimum" ) { return true ; }
        if( type == "StoreMaximumValue" ) { return true ; }
        if( type == "GetParticleOrientation" ) { return true ; }
        if( type == "StorePreviousValue" ) { return true ; }
        if( type == "StoreMinimumValue" ) { return true ; }
        if( type == "Maximum" ) { return true ; }
        if( type == "ExponentialDecay" ) { return true ; }
        if( type == "GetField" ) { return true ; }
        if( type == "WeibullDistributed" ) { return true ; }
        if( type == "UniformDistributedPerParticle" ) { return true ; }
   
        // parsed from header file: ../physics/material_laws/temperature_material_laws.h
        if( type == "ThermalExpansion" ) { return true ; }
        if( type == "RadiationDependentThermalExpansionCoefficient" ) { return true ; }
        if( type == "RadiationDependentPoissonRatio" ) { return true ; }
        if( type == "IncrementalThermalExpansion" ) { return true ; }
        if( type == "RadiationInducedVolumetricExpansion" ) { return true ; }
        if( type == "TemperatureDependentRadiationInducedVolumetricExpansion" ) { return true ; }
        if( type == "Arrhenius" ) { return true ; }
        if( type == "CreepArrhenius" ) { return true ; }
   
        // parsed from header file: ../physics/material_laws/mechanical_material_laws.h
        if( type == "BulkShearConversion" ) { return true ; }
        if( type == "AdjustStrainStressCurve" ) { return true ; }
        if( type == "BazantLoadNonLinearCreep" ) { return true ; }
        if( type == "TensionCompressionCreep" ) { return true ; }
        if( type == "Mineral" ) { return true ; }
   
        // parsed from header file: ../physics/material_laws/humidity_material_laws.h
        if( type == "HumidityDependentThermalExpansionCoefficient" ) { return true ; }
        if( type == "DryingShrinkage" ) { return true ; }
        if( type == "CapillaryPressureDrivenDryingShrinkage" ) { return true ; }
        if( type == "KelvinCapillaryPressure" ) { return true ; }
        if( type == "VanGenuchtenCapillaryPressure" ) { return true ; }
        if( type == "VanGenuchtenWaterSaturation" ) { return true ; }
        if( type == "DisjoiningPressureDrivenDryingShrinkage" ) { return true ; }
        if( type == "BETIsotherm" ) { return true ; }
        if( type == "BiExponentialIsotherm" ) { return true ; }
        if( type == "WaterVapourSaturationPressure" ) { return true ; }
        if( type == "BenboudjemaDryingCreep" ) { return true ; }
        if( type == "BazantRelativeHumidityDependentCreep" ) { return true ; }
        if( type == "HavlasekDryingCreep" ) { return true ; }
        if( type == "WittmannRelativeHumidityDependentCreep" ) { return true ; }
   
        return false ;
    }

    void Object::resetExternalMaterialLaw(ExternalMaterialLaw * target)
    {
        // parsed from header file: ../physics/material_laws/material_laws.h
   
        // parsed from header file: ../physics/material_laws/temperature_material_laws.h
   
        // parsed from header file: ../physics/material_laws/mechanical_material_laws.h
   
        // parsed from header file: ../physics/material_laws/humidity_material_laws.h
   
    }

    Form * Object::getForm(std::string type, std::map<std::string, double> & values, std::map<std::string, std::string> & strings, std::map<std::string, ExternalMaterialLawList*> & externalmateriallawlists, std::map<std::string, FractureCriterion*> & fracturecriterions, std::map<std::string, DamageModel*> & damagemodels, std::map<std::string, LogCreepAccumulator*> & logcreepaccumulators)
    {
        // parsed from header file: ../physics/void_form.h
        if( type == "VoidForm" ) { return new VoidForm() ; }
   
        // parsed from header file: ../physics/stiffness.h
        if( type == "Stiffness" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new Stiffness(values["young_modulus"], values["poisson_ratio"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
        if( type == "WeibullDistributedElasticStiffness" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new WeibullDistributedElasticStiffness(values["young_modulus"], values["poisson_ratio"], values["variability"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/stiffness_with_imposed_deformation.h
        if( type == "StiffnessWithImposedDeformation" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new StiffnessWithImposedDeformation(values["young_modulus"], values["poisson_ratio"], values["imposed_deformation"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/logarithmic_creep_with_external_parameters.h
        if( type == "LogarithmicCreepWithExternalParameters" )
        { 
            if( externalmateriallawlists.find("relations") == externalmateriallawlists.end() ) { externalmateriallawlists["relations"] = nullptr ; } ; 
            if( fracturecriterions.find("fracture_criterion") == fracturecriterions.end() ) { fracturecriterions["fracture_criterion"] = nullptr ; } ; 
            if( damagemodels.find("damage_model") == damagemodels.end() ) { damagemodels["damage_model"] = nullptr ; } ; 
            if( logcreepaccumulators.find("accumulator") == logcreepaccumulators.end() ) { logcreepaccumulators["accumulator"] = new RealTimeLogCreepAccumulator() ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new LogarithmicCreepWithExternalParameters(values, externalmateriallawlists["relations"], fracturecriterions["fracture_criterion"], damagemodels["damage_model"], logcreepaccumulators["accumulator"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/weibull_distributed_stiffness.h
        if( type == "WeibullDistributedStiffness" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( values.find("variability") == values.end() ) { values["variability"] = 0.2 ; } ; 
            if( values.find("material_characteristic_radius") == values.end() ) { values["material_characteristic_radius"] = 0.001 ; } ; 
            return new WeibullDistributedStiffness(values["young_modulus"], values["poisson_ratio"], Enum::getSpaceDimensionality(strings["dimension"]), values["compressive_strength"], values["tensile_strength"], Enum::getplaneType(strings["plane_type"]), values["variability"], values["material_characteristic_radius"]) ;
        }
   
        // parsed from header file: ../physics/logarithmic_creep.h
        if( type == "LogarithmicCreep" )
        { 
            if( values.find("creep_modulus") == values.end() ) { values["creep_modulus"] = -1 ; } ; 
            if( values.find("creep_poisson") == values.end() ) { values["creep_poisson"] = -1 ; } ; 
            if( values.find("creep_characteristic_time") == values.end() ) { values["creep_characteristic_time"] = -1 ; } ; 
            if( values.find("recoverable_modulus") == values.end() ) { values["recoverable_modulus"] = -1 ; } ; 
            if( values.find("recoverable_poisson") == values.end() ) { values["recoverable_poisson"] = -1 ; } ; 
            if( logcreepaccumulators.find("accumulator") == logcreepaccumulators.end() ) { logcreepaccumulators["accumulator"] = new RealTimeLogCreepAccumulator() ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new LogarithmicCreep(values["young_modulus"], values["poisson_ratio"], values["creep_modulus"], values["creep_poisson"], values["creep_characteristic_time"], values["recoverable_modulus"], values["recoverable_poisson"], logcreepaccumulators["accumulator"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/viscoelasticity.h
        if( type == "Viscoelasticity" )
        { 
            if( values.find("creep_characteristic_time") == values.end() ) { values["creep_characteristic_time"] = 1 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new Viscoelasticity(Enum::getViscoelasticModel(strings["model"]), values["young_modulus"], values["poisson_ratio"], values["creep_characteristic_time"], strings["file_name"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
        if( type == "StandardViscoelasticity" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new StandardViscoelasticity(Enum::getCreepComplianceModel(strings["model"]), values["young_modulus"], values["poisson_ratio"], values["creep_modulus"], values["creep_poisson"], values["tau"], values["branches"], values, Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/viscoelasticity_and_fracture.h
        if( type == "ViscoelasticityAndFracture" )
        { 
            if( fracturecriterions.find("fracture_criterion") == fracturecriterions.end() ) { fracturecriterions["fracture_criterion"] = nullptr ; } ; 
            if( damagemodels.find("damage_model") == damagemodels.end() ) { damagemodels["damage_model"] = nullptr ; } ; 
            if( values.find("creep_characteristic_time") == values.end() ) { values["creep_characteristic_time"] = 1 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new ViscoelasticityAndFracture(Enum::getViscoelasticModel(strings["model"]), values["young_modulus"], values["poisson_ratio"], fracturecriterions["fracture_criterion"], damagemodels["damage_model"], values["creep_characteristic_time"], strings["file_name"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/viscoelasticity_and_imposed_deformation.h
        if( type == "ViscoelasticityAndImposedDeformation" )
        { 
            if( values.find("creep_characteristic_time") == values.end() ) { values["creep_characteristic_time"] = 1 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new ViscoelasticityAndImposedDeformation(Enum::getViscoelasticModel(strings["model"]), values["young_modulus"], values["poisson_ratio"], values["imposed_deformation"], values["creep_characteristic_time"], strings["file_name"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/stiffness_and_fracture.h
        if( type == "StiffnessAndFracture" )
        { 
            if( damagemodels.find("damage_model") == damagemodels.end() ) { damagemodels["damage_model"] = nullptr ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new StiffnessAndFracture(values["young_modulus"], values["poisson_ratio"], fracturecriterions["fracture_criterion"], damagemodels["damage_model"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/stiffness_with_imposed_stress.h
        if( type == "StiffnessWithImposedStress" )
        { 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new StiffnessWithImposedStress(values["young_modulus"], values["poisson_ratio"], values["imposed_stress"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/materials/aggregate_behaviour.h
        if( type == "AggregateBehaviour" )
        { 
            if( strings.find("elastic") == strings.end() ) { strings["elastic"] = "FALSE" ; } ; 
            if( strings.find("space_time") == strings.end() ) { strings["space_time"] = "FALSE" ; } ; 
            if( values.find("young_modulus") == values.end() ) { values["young_modulus"] = 59e9 ; } ; 
            if( values.find("poisson_ratio") == values.end() ) { values["poisson_ratio"] = 0.3 ; } ; 
            if( values.find("tensile_strength") == values.end() ) { values["tensile_strength"] = 10e6 ; } ; 
            if( values.find("material_characteristic_radius") == values.end() ) { values["material_characteristic_radius"] = 0.00025 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( values.find("variability") == values.end() ) { values["variability"] = 0.2 ; } ; 
            if( values.find("blocks") == values.end() ) { values["blocks"] = 0 ; } ; 
            return new AggregateBehaviour(Enum::getbool(strings["elastic"]), Enum::getbool(strings["space_time"]), values["young_modulus"], values["poisson_ratio"], values["tensile_strength"], values["material_characteristic_radius"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"]), values["variability"], values["blocks"]) ;
        }
   
        // parsed from header file: ../physics/materials/concrete_behaviour.h
        if( type == "ConcreteBehaviour" )
        { 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( strings.find("redistribution") == strings.end() ) { strings["redistribution"] = "UPPER_BOUND" ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            return new ConcreteBehaviour(values["young_modulus"], values["poisson_ratio"], values["compressive_strength"], Enum::getplaneType(strings["plane_type"]), Enum::getRedistributionType(strings["redistribution"]), Enum::getSpaceDimensionality(strings["dimension"])) ;
        }
   
        // parsed from header file: ../physics/materials/paste_behaviour.h
        if( type == "PasteBehaviour" )
        { 
            if( strings.find("elastic") == strings.end() ) { strings["elastic"] = "FALSE" ; } ; 
            if( strings.find("space_time") == strings.end() ) { strings["space_time"] = "FALSE" ; } ; 
            if( values.find("young_modulus") == values.end() ) { values["young_modulus"] = 12e9 ; } ; 
            if( values.find("poisson_ratio") == values.end() ) { values["poisson_ratio"] = 0.3 ; } ; 
            if( values.find("tensile_strength") == values.end() ) { values["tensile_strength"] = 3e6 ; } ; 
            if( values.find("short_term_creep_modulus") == values.end() ) { values["short_term_creep_modulus"] = 3.6e9 ; } ; 
            if( values.find("long_term_creep_modulus") == values.end() ) { values["long_term_creep_modulus"] = 4e9 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( values.find("material_characteristic_radius") == values.end() ) { values["material_characteristic_radius"] = 0.000175 ; } ; 
            if( values.find("variability") == values.end() ) { values["variability"] = 0.2 ; } ; 
            if( values.find("blocks") == values.end() ) { values["blocks"] = 0 ; } ; 
            return new PasteBehaviour(Enum::getbool(strings["elastic"]), Enum::getbool(strings["space_time"]), values["young_modulus"], values["poisson_ratio"], values["tensile_strength"], values["short_term_creep_modulus"], values["long_term_creep_modulus"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"]), values["material_characteristic_radius"], values["variability"], values["blocks"]) ;
        }
        if( type == "LogCreepPasteBehaviour" )
        { 
            if( strings.find("elastic") == strings.end() ) { strings["elastic"] = "FALSE" ; } ; 
            if( values.find("young_modulus") == values.end() ) { values["young_modulus"] = 12e9 ; } ; 
            if( values.find("poisson_ratio") == values.end() ) { values["poisson_ratio"] = 0.3 ; } ; 
            if( values.find("tensile_strength") == values.end() ) { values["tensile_strength"] = 3e6 ; } ; 
            if( values.find("creep_modulus") == values.end() ) { values["creep_modulus"] = 40e9 ; } ; 
            if( values.find("creep_characteristic_time") == values.end() ) { values["creep_characteristic_time"] = 2 ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( values.find("material_characteristic_radius") == values.end() ) { values["material_characteristic_radius"] = 0.000175 ; } ; 
            if( values.find("variability") == values.end() ) { values["variability"] = 0.2 ; } ; 
            return new LogCreepPasteBehaviour(Enum::getbool(strings["elastic"]), values["young_modulus"], values["poisson_ratio"], values["tensile_strength"], values["creep_modulus"], values["creep_characteristic_time"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"]), values["material_characteristic_radius"], values["variability"]) ;
        }
   
        // parsed from header file: ../physics/materials/gel_behaviour.h
        if( type == "GelBehaviour" )
        { 
            if( strings.find("space_time") == strings.end() ) { strings["space_time"] = "FALSE" ; } ; 
            if( strings.find("dimension") == strings.end() ) { strings["dimension"] = "SPACE_TWO_DIMENSIONAL" ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            if( values.find("blocks") == values.end() ) { values["blocks"] = 0 ; } ; 
            return new GelBehaviour(Enum::getbool(strings["space_time"]), values["young_modulus"], values["poisson_ratio"], values["imposed_deformation"], Enum::getSpaceDimensionality(strings["dimension"]), Enum::getplaneType(strings["plane_type"]), values["blocks"]) ;
        }
   
        return nullptr ;
    }

    bool Object::isForm(std::string type)
    {
        // parsed from header file: ../physics/void_form.h
        if( type == "VoidForm" ) { return true ; }
   
        // parsed from header file: ../physics/stiffness.h
        if( type == "Stiffness" ) { return true ; }
        if( type == "WeibullDistributedElasticStiffness" ) { return true ; }
   
        // parsed from header file: ../physics/stiffness_with_imposed_deformation.h
        if( type == "StiffnessWithImposedDeformation" ) { return true ; }
   
        // parsed from header file: ../physics/logarithmic_creep_with_external_parameters.h
        if( type == "LogarithmicCreepWithExternalParameters" ) { return true ; }
   
        // parsed from header file: ../physics/weibull_distributed_stiffness.h
        if( type == "WeibullDistributedStiffness" ) { return true ; }
   
        // parsed from header file: ../physics/logarithmic_creep.h
        if( type == "LogarithmicCreep" ) { return true ; }
   
        // parsed from header file: ../physics/viscoelasticity.h
        if( type == "Viscoelasticity" ) { return true ; }
        if( type == "StandardViscoelasticity" ) { return true ; }
   
        // parsed from header file: ../physics/viscoelasticity_and_fracture.h
        if( type == "ViscoelasticityAndFracture" ) { return true ; }
   
        // parsed from header file: ../physics/viscoelasticity_and_imposed_deformation.h
        if( type == "ViscoelasticityAndImposedDeformation" ) { return true ; }
   
        // parsed from header file: ../physics/stiffness_and_fracture.h
        if( type == "StiffnessAndFracture" ) { return true ; }
   
        // parsed from header file: ../physics/stiffness_with_imposed_stress.h
        if( type == "StiffnessWithImposedStress" ) { return true ; }
   
        // parsed from header file: ../physics/materials/aggregate_behaviour.h
        if( type == "AggregateBehaviour" ) { return true ; }
   
        // parsed from header file: ../physics/materials/concrete_behaviour.h
        if( type == "ConcreteBehaviour" ) { return true ; }
   
        // parsed from header file: ../physics/materials/paste_behaviour.h
        if( type == "PasteBehaviour" ) { return true ; }
        if( type == "LogCreepPasteBehaviour" ) { return true ; }
   
        // parsed from header file: ../physics/materials/gel_behaviour.h
        if( type == "GelBehaviour" ) { return true ; }
   
        return false ;
    }

    void Object::resetForm(Form * target)
    {
        // parsed from header file: ../physics/void_form.h
   
        // parsed from header file: ../physics/stiffness.h
   
        // parsed from header file: ../physics/stiffness_with_imposed_deformation.h
   
        // parsed from header file: ../physics/logarithmic_creep_with_external_parameters.h
   
        // parsed from header file: ../physics/weibull_distributed_stiffness.h
   
        // parsed from header file: ../physics/logarithmic_creep.h
   
        // parsed from header file: ../physics/viscoelasticity.h
   
        // parsed from header file: ../physics/viscoelasticity_and_fracture.h
   
        // parsed from header file: ../physics/viscoelasticity_and_imposed_deformation.h
   
        // parsed from header file: ../physics/stiffness_and_fracture.h
   
        // parsed from header file: ../physics/stiffness_with_imposed_stress.h
   
        // parsed from header file: ../physics/materials/aggregate_behaviour.h
   
        // parsed from header file: ../physics/materials/concrete_behaviour.h
   
        // parsed from header file: ../physics/materials/paste_behaviour.h
   
        // parsed from header file: ../physics/materials/gel_behaviour.h
   
    }

    DamageModel * Object::getDamageModel(std::string type, std::map<std::string, double> & values)
    {
        // parsed from header file: ../physics/damagemodels/isotropiclineardamage.h
        if( type == "Isotropic" ) { return new IsotropicLinearDamage() ; }
   
        // parsed from header file: ../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h
        if( type == "SpaceTimeFiberBasedIsotropic" )
        { 
            if( values.find("damage_increment") == values.end() ) { values["damage_increment"] = 0.1 ; } ; 
            if( values.find("time_tolerance") == values.end() ) { values["time_tolerance"] = 0.001 ; } ; 
            if( values.find("maximum_damage") == values.end() ) { values["maximum_damage"] = 0.6 ; } ; 
            return new SpaceTimeFiberBasedIsotropicLinearDamage(values["damage_increment"], values["time_tolerance"], values["maximum_damage"]) ;
        }
   
        return nullptr ;
    }

    bool Object::isDamageModel(std::string type)
    {
        // parsed from header file: ../physics/damagemodels/isotropiclineardamage.h
        if( type == "Isotropic" ) { return true ; }
   
        // parsed from header file: ../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h
        if( type == "SpaceTimeFiberBasedIsotropic" ) { return true ; }
   
        return false ;
    }

    void Object::resetDamageModel(DamageModel * target)
    {
        // parsed from header file: ../physics/damagemodels/isotropiclineardamage.h
   
        // parsed from header file: ../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h
   
    }

    FractureCriterion * Object::getFractureCriterion(std::string type, std::map<std::string, double> & values, std::map<std::string, std::string> & strings)
    {
        // parsed from header file: ../physics/fracturecriteria/mcft.h
        if( type == "NonLocalMCFT" ) { return new NonLocalMCFT(values["compressive_strength"], values["young_modulus"], values["material_characteristic_radius"], Enum::getRedistributionType(strings["redistribution_type"])) ; }
        if( type == "NonLocalSpaceTimeMCFT" ) { return new NonLocalSpaceTimeMCFT(values["compressive_strength"], values["young_modulus"], values["material_characteristic_radius"], Enum::getRedistributionType(strings["redistribution_type"])) ; }
   
        // parsed from header file: ../physics/fracturecriteria/maxstrain.h
        if( type == "MaximumStrain" ) { return new MaximumStrain(values["tensile_strain"]) ; }
        if( type == "SpaceTimeNonLocalMaximumStrain" ) { return new SpaceTimeNonLocalMaximumStrain(values["tensile_strain"]) ; }
        if( type == "SpaceTimeNonLocalLinearSofteningMaximumStrain" ) { return new SpaceTimeNonLocalLinearSofteningMaximumStrain(values["tensile_strain"], values["tensile_strength"], values["tensile_ultimate_strain"]) ; }
        if( type == "SpaceTimeNonLocalMaximumStress" ) { return new SpaceTimeNonLocalMaximumStress(values["tensile_strength"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulombwithstrain.h
        if( type == "ConfinedMohrCoulombWithStrainLimit" ) { return new ConfinedMohrCoulombWithStrainLimit(values["tensile_strength"], values["compressive_strength"], values["tensile_strain"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/vonmises.h
        if( type == "VonMises" ) { return new VonMises(values["tensile_strength"]) ; }
        if( type == "VonMisesStrain" ) { return new VonMisesStrain(values["tensile_strain"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedvonmises.h
        if( type == "ConfinedVonMises" ) { return new ConfinedVonMises(values["tensile_strength"], values["compressive_strength"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/mohrcoulomb.h
        if( type == "MohrCoulomb" ) { return new MohrCoulomb(values["tensile_strength"], values["compressive_strength"]) ; }
        if( type == "NonLocalMohrCoulomb" ) { return new NonLocalMohrCoulomb(values["tensile_strength"], values["compressive_strength"], values["young_modulus"]) ; }
        if( type == "SpaceTimeNonLocalMohrCoulomb" ) { return new SpaceTimeNonLocalMohrCoulomb(values["tensile_strength"], values["compressive_strength"], values["young_modulus"]) ; }
        if( type == "NonLocalLinearlyDecreasingMohrCoulomb" ) { return new NonLocalLinearlyDecreasingMohrCoulomb(values["tensile_strength"], values["compressive_strength"], values["tensile_ultimate_strain"], values["compressive_ultimate_strain"], values["young_modulus"]) ; }
        if( type == "NonLocalExponentiallyDecreasingMohrCoulomb" ) { return new NonLocalExponentiallyDecreasingMohrCoulomb(values["tensile_strength"], values["compressive_strength"], values["tensile_ultimate_strain"], values["compressive_ultimate_strain"], values["young_modulus"]) ; }
        if( type == "NonLocalInverseRootMohrCoulomb" ) { return new NonLocalInverseRootMohrCoulomb(values["tensile_strain"], values["tensile_ultimate_strain"], values["young_modulus"], values["tensile_strain_coefficient"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/mazars.h
        if( type == "NonLocalMazars" )
        { 
            if( values.find("compressive_strain") == values.end() ) { values["compressive_strain"] = 1 ; } ; 
            if( values.find("compressive_stress") == values.end() ) { values["compressive_stress"] = 1 ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new NonLocalMazars(values["tensile_strain"], values["young_modulus"], values["poisson_ratio"], values["tensile_fracture_energy"], values["compressive_strain"], values["compressive_stress"], values["material_characteristic_radius"], Enum::getplaneType(strings["plane_type"])) ;
        }
        if( type == "NonLocalSpaceTimeMazars" )
        { 
            if( values.find("compressive_strain") == values.end() ) { values["compressive_strain"] = 1 ; } ; 
            if( values.find("compressive_stress") == values.end() ) { values["compressive_stress"] = 1 ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            return new NonLocalSpaceTimeMazars(values["tensile_strain"], values["young_modulus"], values["poisson_ratio"], values["tensile_fracture_energy"], values["compressive_strain"], values["compressive_stress"], values["material_characteristic_radius"], Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h
        if( type == "SpaceTimeNonLocalMultiLinearSofteningFractureCriterion" )
        { 
            if( values.find("strain_renormalization_factor") == values.end() ) { values["strain_renormalization_factor"] = 1e4 ; } ; 
            if( values.find("stress_renormalization_factor") == values.end() ) { values["stress_renormalization_factor"] = 1e-6 ; } ; 
            return new SpaceTimeNonLocalMultiLinearSofteningFractureCriterion(strings["tension_file_name"], values["young_modulus"], values["strain_renormalization_factor"], values["stress_renormalization_factor"]) ;
        }
        if( type == "AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion" )
        { 
            if( values.find("strain_renormalization_factor") == values.end() ) { values["strain_renormalization_factor"] = 1e4 ; } ; 
            if( values.find("stress_renormalization_factor") == values.end() ) { values["stress_renormalization_factor"] = 1e-6 ; } ; 
            return new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion(strings["tension_file_name"], strings["compression_file_name"], values["young_modulus"], values["strain_renormalization_factor"], values["stress_renormalization_factor"]) ;
        }
   
        // parsed from header file: ../physics/fracturecriteria/boundedvonmises.h
        if( type == "BoundedVonMises" ) { return new BoundedVonMises(values["tensile_strength"], values["maximum_damage"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/spacetimeflowrule.h
        if( type == "SpaceTimeNonLocalDamageFlowRule" ) { return new SpaceTimeNonLocalDamageFlowRule(strings["file_name"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulomb.h
        if( type == "ConfinedMohrCoulomb" ) { return new ConfinedMohrCoulomb(values["tensile_strength"], values["compressive_strength"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/limitstrains.h
        if( type == "LimitStrains" ) { return new LimitStrains(values["tensile_strain"], values["compressive_strain"]) ; }
   
        return nullptr ;
    }

    bool Object::isFractureCriterion(std::string type)
    {
        // parsed from header file: ../physics/fracturecriteria/mcft.h
        if( type == "NonLocalMCFT" ) { return true ; }
        if( type == "NonLocalSpaceTimeMCFT" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/maxstrain.h
        if( type == "MaximumStrain" ) { return true ; }
        if( type == "SpaceTimeNonLocalMaximumStrain" ) { return true ; }
        if( type == "SpaceTimeNonLocalLinearSofteningMaximumStrain" ) { return true ; }
        if( type == "SpaceTimeNonLocalMaximumStress" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulombwithstrain.h
        if( type == "ConfinedMohrCoulombWithStrainLimit" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/vonmises.h
        if( type == "VonMises" ) { return true ; }
        if( type == "VonMisesStrain" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedvonmises.h
        if( type == "ConfinedVonMises" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/mohrcoulomb.h
        if( type == "MohrCoulomb" ) { return true ; }
        if( type == "NonLocalMohrCoulomb" ) { return true ; }
        if( type == "SpaceTimeNonLocalMohrCoulomb" ) { return true ; }
        if( type == "NonLocalLinearlyDecreasingMohrCoulomb" ) { return true ; }
        if( type == "NonLocalExponentiallyDecreasingMohrCoulomb" ) { return true ; }
        if( type == "NonLocalInverseRootMohrCoulomb" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/mazars.h
        if( type == "NonLocalMazars" ) { return true ; }
        if( type == "NonLocalSpaceTimeMazars" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h
        if( type == "SpaceTimeNonLocalMultiLinearSofteningFractureCriterion" ) { return true ; }
        if( type == "AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/boundedvonmises.h
        if( type == "BoundedVonMises" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/spacetimeflowrule.h
        if( type == "SpaceTimeNonLocalDamageFlowRule" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulomb.h
        if( type == "ConfinedMohrCoulomb" ) { return true ; }
   
        // parsed from header file: ../physics/fracturecriteria/limitstrains.h
        if( type == "LimitStrains" ) { return true ; }
   
        return false ;
    }

    void Object::resetFractureCriterion(FractureCriterion * target, std::map<std::string, double> & values, std::map<std::string, std::string> & strings)
    {
        // parsed from header file: ../physics/fracturecriteria/mcft.h
   
        // parsed from header file: ../physics/fracturecriteria/maxstrain.h
        if(dynamic_cast<SpaceTimeNonLocalLinearSofteningMaximumStrain *>(target) != nullptr) { dynamic_cast<SpaceTimeNonLocalLinearSofteningMaximumStrain *>(target)->reset(values["tensile_strain"], values["tensile_strength"], values["tensile_ultimate_strain"]) ; }
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulombwithstrain.h
   
        // parsed from header file: ../physics/fracturecriteria/vonmises.h
   
        // parsed from header file: ../physics/fracturecriteria/confinedvonmises.h
   
        // parsed from header file: ../physics/fracturecriteria/mohrcoulomb.h
   
        // parsed from header file: ../physics/fracturecriteria/mazars.h
        if(dynamic_cast<NonLocalMazars *>(target) != nullptr)
        { 
            if( values.find("compressive_strain") == values.end() ) { values["compressive_strain"] = 1 ; } ; 
            if( values.find("compressive_stress") == values.end() ) { values["compressive_stress"] = 1 ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            dynamic_cast<NonLocalMazars *>(target)->reset(values["tensile_strain"], values["young_modulus"], values["poisson_ratio"], values["tensile_fracture_energy"], values["compressive_strain"], values["compressive_stress"], values["material_characteristic_radius"], Enum::getplaneType(strings["plane_type"])) ;
        }
        if(dynamic_cast<NonLocalSpaceTimeMazars *>(target) != nullptr)
        { 
            if( values.find("compressive_strain") == values.end() ) { values["compressive_strain"] = 1 ; } ; 
            if( values.find("compressive_stress") == values.end() ) { values["compressive_stress"] = 1 ; } ; 
            if( strings.find("plane_type") == strings.end() ) { strings["plane_type"] = "PLANE_STRESS" ; } ; 
            dynamic_cast<NonLocalSpaceTimeMazars *>(target)->reset(values["tensile_strain"], values["young_modulus"], values["poisson_ratio"], values["tensile_fracture_energy"], values["compressive_strain"], values["compressive_stress"], values["material_characteristic_radius"], Enum::getplaneType(strings["plane_type"])) ;
        }
   
        // parsed from header file: ../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h
   
        // parsed from header file: ../physics/fracturecriteria/boundedvonmises.h
   
        // parsed from header file: ../physics/fracturecriteria/spacetimeflowrule.h
   
        // parsed from header file: ../physics/fracturecriteria/confinedmohrcoulomb.h
   
        // parsed from header file: ../physics/fracturecriteria/limitstrains.h
   
    }

    LogCreepAccumulator * Object::getLogCreepAccumulator(std::string type, std::map<std::string, double> & values)
    {
        // parsed from header file: ../physics/material_laws/logcreep_accumulator.h
        if( type == "No" ) { return new NoLogCreepAccumulator() ; }
        if( type == "StrainConsolidation" ) { return new StrainConsolidationLogCreepAccumulator() ; }
        if( type == "RealTime" )
        { 
            if( values.find("viscous_flow") == values.end() ) { values["viscous_flow"] = 0 ; } ; 
            return new RealTimeLogCreepAccumulator(values["viscous_flow"]) ;
        }
        if( type == "TimeUnderLoad" ) { return new TimeUnderLoadLogCreepAccumulator() ; }
   
        return nullptr ;
    }

    bool Object::isLogCreepAccumulator(std::string type)
    {
        // parsed from header file: ../physics/material_laws/logcreep_accumulator.h
        if( type == "No" ) { return true ; }
        if( type == "StrainConsolidation" ) { return true ; }
        if( type == "RealTime" ) { return true ; }
        if( type == "TimeUnderLoad" ) { return true ; }
   
        return false ;
    }

    void Object::resetLogCreepAccumulator(LogCreepAccumulator * target)
    {
        // parsed from header file: ../physics/material_laws/logcreep_accumulator.h
   
    }

    InclusionGenerator * Object::getInclusionGenerator(std::string type, std::map<std::string, double> & values)
    {
        // parsed from header file: ../features/microstructuregenerator.h
        if( type == "Circular" )
        { 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new CircularInclusionGenerator(values["placement_rotation"]) ;
        }
        if( type == "XFEM" ) { return new XFEMInclusionGenerator() ; }
        if( type == "SpaceTimeXFEM" ) { return new SpaceTimeXFEMInclusionGenerator() ; }
        if( type == "Ellipsoidal" )
        { 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("shape_factor_variability") == values.end() ) { values["shape_factor_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new EllipsoidalInclusionGenerator(values["shape_factor"], values["orientation"], values["orientation_variability"], values["shape_factor_variability"], values["placement_rotation"]) ;
        }
        if( type == "Rectangular" )
        { 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("shape_factor_variability") == values.end() ) { values["shape_factor_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new RectangularInclusionGenerator(values["shape_factor"], values["orientation"], values["orientation_variability"], values["shape_factor_variability"], values["placement_rotation"]) ;
        }
        if( type == "Polygonal" )
        { 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("vertex_variability") == values.end() ) { values["vertex_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new PolygonalInclusionGenerator(values["vertex"], values["orientation"], values["orientation_variability"], values["vertex_variability"], values["placement_rotation"]) ;
        }
        if( type == "VoronoiPolygonal" )
        { 
            if( values.find("shape_factor") == values.end() ) { values["shape_factor"] = 1 ; } ; 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("shape_factor_variability") == values.end() ) { values["shape_factor_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new VoronoiPolygonalInclusionGenerator(values["box_width"], values["grains"], values["spacing"], values["shape_factor"], values["orientation"], values["orientation_variability"], values["shape_factor_variability"], values["placement_rotation"]) ;
        }
        if( type == "GravelPolygonal" )
        { 
            if( values.find("amplitude_factor") == values.end() ) { values["amplitude_factor"] = 0.9 ; } ; 
            if( values.find("amplitude_exponent") == values.end() ) { values["amplitude_exponent"] = 1.9 ; } ; 
            if( values.find("degree") == values.end() ) { values["degree"] = 2 ; } ; 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("vertex_variability") == values.end() ) { values["vertex_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new GravelPolygonalInclusionGenerator(values["amplitude_factor"], values["amplitude_exponent"], values["degree"], values["vertex"], values["orientation"], values["orientation_variability"], values["vertex_variability"], values["placement_rotation"]) ;
        }
        if( type == "CrushedPolygonal" )
        { 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("vertex_variability") == values.end() ) { values["vertex_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new CrushedPolygonalInclusionGenerator(values["shape_factor"], values["vertex"], values["orientation"], values["orientation_variability"], values["vertex_variability"], values["placement_rotation"]) ;
        }
        if( type == "CrushedSubtendedPolygonal" )
        { 
            if( values.find("orientation") == values.end() ) { values["orientation"] = 0 ; } ; 
            if( values.find("orientation_variability") == values.end() ) { values["orientation_variability"] = 3.141592 ; } ; 
            if( values.find("vertex_variability") == values.end() ) { values["vertex_variability"] = 0 ; } ; 
            if( values.find("placement_rotation") == values.end() ) { values["placement_rotation"] = 0 ; } ; 
            return new CrushedSubtendedPolygonalInclusionGenerator(values["shape_factor"], values["angle_variability"], values["vertex"], values["orientation"], values["orientation_variability"], values["vertex_variability"], values["placement_rotation"]) ;
        }
   
        return nullptr ;
    }

    bool Object::isInclusionGenerator(std::string type)
    {
        // parsed from header file: ../features/microstructuregenerator.h
        if( type == "Circular" ) { return true ; }
        if( type == "XFEM" ) { return true ; }
        if( type == "SpaceTimeXFEM" ) { return true ; }
        if( type == "Ellipsoidal" ) { return true ; }
        if( type == "Rectangular" ) { return true ; }
        if( type == "Polygonal" ) { return true ; }
        if( type == "VoronoiPolygonal" ) { return true ; }
        if( type == "GravelPolygonal" ) { return true ; }
        if( type == "CrushedPolygonal" ) { return true ; }
        if( type == "CrushedSubtendedPolygonal" ) { return true ; }
   
        return false ;
    }

    void Object::resetInclusionGenerator(InclusionGenerator * target)
    {
        // parsed from header file: ../features/microstructuregenerator.h
   
    }

    ParticleSizeDistribution * Object::getParticleSizeDistribution(std::string type, std::map<std::string, double> & values, std::map<std::string, std::string> & strings)
    {
        // parsed from header file: ../utilities/granulo.h
        if( type == "PSDBolomeA" ) { return new PSDBolomeA() ; }
        if( type == "PSDBolomeB" ) { return new PSDBolomeB() ; }
        if( type == "PSDBolomeC" ) { return new PSDBolomeC() ; }
        if( type == "PSDBolomeD" ) { return new PSDBolomeD() ; }
        if( type == "PSDFuller" )
        { 
            if( values.find("radius_minimum") == values.end() ) { values["radius_minimum"] = 0 ; } ; 
            if( values.find("exponent") == values.end() ) { values["exponent"] = 0.5 ; } ; 
            return new PSDFuller(values["radius_minimum"], values["exponent"]) ;
        }
        if( type == "ConstantSizeDistribution" ) { return new ConstantSizeDistribution() ; }
        if( type == "GranuloFromCumulativePSD" )
        { 
            if( values.find("factor") == values.end() ) { values["factor"] = 1 ; } ; 
            if( values.find("radius_maximum") == values.end() ) { values["radius_maximum"] = -1 ; } ; 
            if( values.find("radius_minimum") == values.end() ) { values["radius_minimum"] = -1 ; } ; 
            return new GranuloFromCumulativePSD(strings["file_name"], Enum::getPSDSpecificationType(strings["specification"]), values["factor"], values["radius_maximum"], values["radius_minimum"]) ;
        }
   
        return nullptr ;
    }

    bool Object::isParticleSizeDistribution(std::string type)
    {
        // parsed from header file: ../utilities/granulo.h
        if( type == "PSDBolomeA" ) { return true ; }
        if( type == "PSDBolomeB" ) { return true ; }
        if( type == "PSDBolomeC" ) { return true ; }
        if( type == "PSDBolomeD" ) { return true ; }
        if( type == "PSDFuller" ) { return true ; }
        if( type == "ConstantSizeDistribution" ) { return true ; }
        if( type == "GranuloFromCumulativePSD" ) { return true ; }
   
        return false ;
    }

    void Object::resetParticleSizeDistribution(ParticleSizeDistribution * target)
    {
        // parsed from header file: ../utilities/granulo.h
   
    }

    InclusionFamily * Object::getInclusionFamily(std::string type, std::map<std::string, double> & values, std::map<std::string, ParticleSizeDistribution*> & particlesizedistributions, std::map<std::string, InclusionGenerator*> & inclusiongenerators, std::map<std::string, std::string> & strings)
    {
        // parsed from header file: ../utilities/inclusion_family.h
        if( type == "." ) { return new InclusionFamily(values["number"], values["radius_maximum"], values["surface"], particlesizedistributions["particle_size_distribution"], inclusiongenerators["geometry"]) ; }
        if( type == "Embedded" ) { return new EmbeddedInclusionFamily(values["number"], values["radius_maximum"], values["surface"], particlesizedistributions["particle_size_distribution"], inclusiongenerators["geometry"]) ; }
        if( type == "Masked" ) { return new MaskedInclusionFamily(values["number"], values["radius_maximum"], values["surface"], particlesizedistributions["particle_size_distribution"], inclusiongenerators["geometry"]) ; }
        if( type == "Concentric" ) { return new ConcentricInclusionFamily(values["layer_width"]) ; }
        if( type == "Voronoi" )
        { 
            if( values.find("outside_layer") == values.end() ) { values["outside_layer"] = -1 ; } ; 
            if( values.find("interface") == values.end() ) { values["interface"] = 0 ; } ; 
            if( values.find("maximum_vertex") == values.end() ) { values["maximum_vertex"] = 20 ; } ; 
            return new VoronoiInclusionFamily(values["radius"], values["surface_fraction"], values["correction_factor"], values["number_of_grains"], values["outside_layer"], values["interface"], values["maximum_vertex"]) ;
        }
        if( type == "FileDefinedCircle" )
        { 
            if( strings.find("column_1") == strings.end() ) { strings["column_1"] = "radius" ; } ; 
            if( strings.find("column_2") == strings.end() ) { strings["column_2"] = "center_x" ; } ; 
            if( strings.find("column_3") == strings.end() ) { strings["column_3"] = "center_y" ; } ; 
            return new FileDefinedCircleInclusionFamily(values["number"], strings["file_name"], strings["column_1"], strings["column_2"], strings["column_3"]) ;
        }
        if( type == "FileDefinedPolygon" ) { return new FileDefinedPolygonInclusionFamily(values["number"], strings["file_name"]) ; }
   
        return nullptr ;
    }

    bool Object::isInclusionFamily(std::string type)
    {
        // parsed from header file: ../utilities/inclusion_family.h
        if( type == "." ) { return true ; }
        if( type == "Embedded" ) { return true ; }
        if( type == "Masked" ) { return true ; }
        if( type == "Concentric" ) { return true ; }
        if( type == "Voronoi" ) { return true ; }
        if( type == "FileDefinedCircle" ) { return true ; }
        if( type == "FileDefinedPolygon" ) { return true ; }
   
        return false ;
    }

    void Object::resetInclusionFamily(InclusionFamily * target)
    {
        // parsed from header file: ../utilities/inclusion_family.h
   
    }

    EnrichmentManager * Object::getEnrichmentManager(std::string type, std::map<std::string, FeatureTree*> & featuretrees, std::map<std::string, InclusionFamily*> & inclusionfamilys, std::map<std::string, double> & values, std::map<std::string, std::string> & strings)
    {
        // parsed from header file: ../features/enrichmentmanagers/gelmanager.h
        if( type == "Gel" )
        { 
            if( values.find("reactive_fraction") == values.end() ) { values["reactive_fraction"] = 0.1 ; } ; 
            if( values.find("radius_increment") == values.end() ) { values["radius_increment"] = -1 ; } ; 
            return new GelManager(featuretrees["feature_tree"], inclusionfamilys["zones"], values["reactive_fraction"], values["radius_increment"]) ;
        }
        if( type == "FunctionBasedGel" )
        { 
            if( values.find("reactive_fraction") == values.end() ) { values["reactive_fraction"] = 0.1 ; } ; 
            return new FunctionBasedGelManager(featuretrees["feature_tree"], inclusionfamilys["zones"], strings["radius"], values["reactive_fraction"]) ;
        }
        if( type == "SpaceTimeGel" )
        { 
            if( values.find("reactive_fraction") == values.end() ) { values["reactive_fraction"] = 0.1 ; } ; 
            return new SpaceTimeGelManager(featuretrees["feature_tree"], inclusionfamilys["zones"], strings["radius"], values["reactive_fraction"]) ;
        }
   
        return nullptr ;
    }

    bool Object::isEnrichmentManager(std::string type)
    {
        // parsed from header file: ../features/enrichmentmanagers/gelmanager.h
        if( type == "Gel" ) { return true ; }
        if( type == "FunctionBasedGel" ) { return true ; }
        if( type == "SpaceTimeGel" ) { return true ; }
   
        return false ;
    }

    void Object::resetEnrichmentManager(EnrichmentManager * target)
    {
        // parsed from header file: ../features/enrichmentmanagers/gelmanager.h
   
    }

    Sampler * Object::getSampler(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, Point> & points, std::map<std::string, double> & values)
    {
        // parsed from header file: ../geometry/sampler/regular_sampler.h
        if( type == "Regular" )
        { 
            if( strings.find("force") == strings.end() ) { strings["force"] = "false" ; } ; 
            return new RegularSampler(Enum::getbool(strings["force"])) ;
        }
   
        // parsed from header file: ../geometry/sampler/gradient_sampler.h
        if( type == "Gradient" )
        { 
            if( values.find("randomize") == values.end() ) { values["randomize"] = -1 ; } ; 
            return new GradientSampler(points["start_point"], points["end_point"], values["start_factor"], values["end_factor"], values["randomize"]) ;
        }
   
        // parsed from header file: ../geometry/sampler/sampler.h
        if( type == "." ) { return new Sampler() ; }
   
        return nullptr ;
    }

    bool Object::isSampler(std::string type)
    {
        // parsed from header file: ../geometry/sampler/regular_sampler.h
        if( type == "Regular" ) { return true ; }
   
        // parsed from header file: ../geometry/sampler/gradient_sampler.h
        if( type == "Gradient" ) { return true ; }
   
        // parsed from header file: ../geometry/sampler/sampler.h
        if( type == "." ) { return true ; }
   
        return false ;
    }

    void Object::resetSampler(Sampler * target)
    {
        // parsed from header file: ../geometry/sampler/regular_sampler.h
   
        // parsed from header file: ../geometry/sampler/gradient_sampler.h
   
        // parsed from header file: ../geometry/sampler/sampler.h
   
    }

}

