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
struct ThermalExpansionMaterialLaw : public ExternalMaterialLaw
{
    ThermalExpansionMaterialLaw(std::string args, char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct RadiationInducedExpansionMaterialLaw : public ExternalMaterialLaw
{
    RadiationInducedExpansionMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct DryingShrinkageMaterialLaw : public ExternalMaterialLaw
{
    bool effective ;

    DryingShrinkageMaterialLaw(bool ef = false, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/* Material law for a generic Arrhenius equation affecting a single parameter
 * The material parameters must include the affected parameter, the activation energy for that parameter (in Kelvin^{-1}, and written as parameter"_activation_energy"), "temperature" and "temperature_reference" (both in Kelvin)
 * The reference value of the affected parameter and "temperature" must be given in the defaults values
 */
struct ArrheniusMaterialLaw : public ExternalMaterialLaw
{
    std::string affected ;
    std::string coefficient ;
    ArrheniusMaterialLaw(std::string a, std::string args, char sep = ',') ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/* Material law for the Arrhenius equation affecting all creep parameters
 * The material parameters must include "creep_characteristic_time", "creep_modulus" (or the pair "creep_bulk" and "creep_shear"), "creep_activation_energy" (in Kelvin^{-1}), "temperature" and "temperature_reference"
 * The reference value of "creep_characteristic_time","creep_modulus" (or the pair "creep_bulk" and "creep_shear"), and "temperature" must be given in the default values
 */
struct CreepArrheniusMaterialLaw : public ExternalMaterialLaw
{
    CreepArrheniusMaterialLaw(std::string args, char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct CreepRelativeHumidityMaterialLaw : public ExternalMaterialLaw
{
    CreepRelativeHumidityMaterialLaw(std::string args = std::string("creep_humidity_coefficient = 5"), char sep = 'c') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt);
};


} ;


#endif // __TEMPERATURE_MATERIAL_LAW_H_
