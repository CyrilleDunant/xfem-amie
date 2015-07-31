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

#ifndef __TEMPERATURE_MATERIAL_LAW_H_
#define __TEMPERATURE_MATERIAL_LAW_H_

#include "material_laws.h"

namespace Amie
{

/* Material law for thermal expansion
 * The material parameters must include "temperature", "thermal_expansion_coefficient" and "imposed_deformation"
 * The reference temperature must be given in the defaults values
 * Temperatures must be set in Kelvin for compatibility with the Arrhenius law
 */
/*PARSE ThermalExpansion 0 VOID */
struct ThermalExpansionMaterialLaw : public ExternalMaterialLaw
{
    ThermalExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~ThermalExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE RadiationDependentThermalExpansionCoefficient 0 VOID */
struct RadiationDependentThermalExpansionCoefficientMaterialLaw : public ExternalMaterialLaw
{
    RadiationDependentThermalExpansionCoefficientMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~RadiationDependentThermalExpansionCoefficientMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE RadiationDependentPoissonRatio 0 VOID */
struct RadiationDependentPoissonRatioMaterialLaw : public ExternalMaterialLaw
{
    RadiationDependentPoissonRatioMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~RadiationDependentPoissonRatioMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct AnisotropicThermalExpansionMaterialLaw : public ThermalExpansionMaterialLaw
{
    AnisotropicThermalExpansionMaterialLaw(std::string args, char sep = ',') : ThermalExpansionMaterialLaw(args, sep) { }
    virtual ~AnisotropicThermalExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE IncrementalThermalExpansion 0 VOID */
struct IncrementalThermalExpansionMaterialLaw : public ExternalMaterialLaw
{
    IncrementalThermalExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~IncrementalThermalExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct AnisotropicIncrementalThermalExpansionMaterialLaw : public ExternalMaterialLaw
{
    AnisotropicIncrementalThermalExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~AnisotropicIncrementalThermalExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE RadiationInducedVolumetricExpansion 0 VOID */
struct RadiationInducedVolumetricExpansionMaterialLaw : public ExternalMaterialLaw
{
    RadiationInducedVolumetricExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~RadiationInducedVolumetricExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE TemperatureDependentRadiationInducedVolumetricExpansion 0 VOID */
struct TemperatureDependentRadiationInducedVolumetricExpansionMaterialLaw : public ExternalMaterialLaw
{
    TemperatureDependentRadiationInducedVolumetricExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~TemperatureDependentRadiationInducedVolumetricExpansionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/* Material law for a generic Arrhenius equation affecting a single parameter
 * The material parameters must include the affected parameter, the activation energy for that parameter (in Kelvin^{-1}, and written as parameter"_activation_energy"), "temperature" and "temperature_reference" (both in Kelvin)
 * The reference value of the affected parameter and "temperature" must be given in the defaults values
 */
/*PARSE Arrhenius 1 VOID */
struct ArrheniusMaterialLaw : public ExternalMaterialLaw
{
    std::string affected ;
    std::string coefficient ;
    ArrheniusMaterialLaw(std::string a, std::string args = std::string(), char sep = ',') ;
    virtual ~ArrheniusMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/* Material law for the Arrhenius equation affecting all creep parameters
 * The material parameters must include "creep_characteristic_time", "creep_modulus" (or the pair "creep_bulk" and "creep_shear"), "creep_activation_energy" (in Kelvin^{-1}), "temperature" and "temperature_reference"
 * The reference value of "creep_characteristic_time","creep_modulus" (or the pair "creep_bulk" and "creep_shear"), and "temperature" must be given in the default values
 */
/*PARSE CreepArrhenius 0 VOID */
struct CreepArrheniusMaterialLaw : public ExternalMaterialLaw
{
    CreepArrheniusMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }
    virtual ~CreepArrheniusMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};


} 


#endif // __TEMPERATURE_MATERIAL_LAW_H_
