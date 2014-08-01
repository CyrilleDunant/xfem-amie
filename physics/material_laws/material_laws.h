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

#ifndef __MATERIAL_LAW_H_
#define __MATERIAL_LAW_H_

#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"

namespace Amie
{

/* Parser for the default values for material laws
 * The string must be ordered as a series of "parameter = value" separated with the appropriate character.
 * Spaces are removed during parsing.
 */
std::map<std::string, double> parseDefaultValues(std::string args, char sep = ',') ;

/* Basic structure for material law. This structure is empty ; the preProcess() and/or step() methods must be surcharged for the law to do something */
struct ExternalMaterialLaw
{
    std::map<std::string, double> defaultValues ;

    ExternalMaterialLaw(std::string args = std::string(), char sep = ',') : defaultValues(parseDefaultValues(args, sep)) { }
    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) { }
    virtual void step( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s) { }
};

/* Generic material law to set one or several variables to a constant
 * The variables affected and their values are the variables specified in the arguments
 */
struct ConstantExternalMaterialLaw : public ExternalMaterialLaw
{
    ConstantExternalMaterialLaw( std::string args, char sep = ',' ) : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;

};

/* Generic material law to set an external variable from an analytical function of the space-time coordinates */
struct SpaceTimeDependentExternalMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    Function f ;
    SpaceTimeDependentExternalMaterialLaw( std::string e, const char *f_, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_) { }
    SpaceTimeDependentExternalMaterialLaw( std::string e, Function & f_, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;


};

/* Generic material law to set an external variable as a function of another one variable */
struct SimpleDependentExternalMaterialLaw: public SpaceTimeDependentExternalMaterialLaw
{
    std::string coordinate ;

    SimpleDependentExternalMaterialLaw( std::string e, std::string c, const char *f_, std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, args, sep), coordinate(c) { }
    SimpleDependentExternalMaterialLaw( std::string e, std::string c,  Function & f_, std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, args, sep), coordinate(c) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;

};

/* Generic material law to set an external variable from an analytical function of other variables and/or the space-time coordinates
 * "coordinates" maps the coordinate axis (valid values are: x y z t u v w) to the name of the external variable to consider on each axis
 * If useSpaceTimeCoordinates is set to true, any x,y,z,t coordinates not specified as external variables will be set as the space-time coordinates of the element's center
 *
 * For functions of a single parameter, use the SimpleDependentExternalMaterialLaw instead.
 */
struct VariableDependentExternalMaterialLaw: public SpaceTimeDependentExternalMaterialLaw
{
    std::map<char, std::string> coordinates ;
    bool useSpaceTimeCoordinates ;
    VariableDependentExternalMaterialLaw( std::string e, const char *f_,  std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_,  std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, const char *f_,  std::map<char, std::string> c, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_,  std::map<char, std::string> c, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_,args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;

    bool has(char c) const { return !(coordinates.empty() || coordinates.find(c) == coordinates.end()) ; }
    void setAsX(std::string s) { coordinates['x'] = s ; }
    void setAsY(std::string s) { coordinates['y'] = s ; }
    void setAsZ(std::string s) { coordinates['z'] = s ; }
    void setAsT(std::string s) { coordinates['t'] = s ; }
    void setAsU(std::string s) { coordinates['u'] = s ; }
    void setAsV(std::string s) { coordinates['v'] = s ; }
    void setAsW(std::string s) { coordinates['w'] = s ; }

};

/* Material law for thermal expansion
 * The material parameters must include "temperature", "thermal_expansion_coefficient" and "imposed_deformation"
 * The reference temperature must be given in the defaults values
 * Temperatures must be set in Kelvin for compatibility with the Arrhenius law
 */
struct ThermalExpansionMaterialLaw : public ExternalMaterialLaw
{
    ThermalExpansionMaterialLaw(std::string args, char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;
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

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;
};

/* Material law for the Arrhenius equation affecting all creep parameters
 * The material parameters must include "creep_characteristic_time", "creep_modulus" (or the pair "creep_bulk" and "creep_shear"), "creep_activation_energy" (in Kelvin^{-1}), "temperature" and "temperature_reference"
 * The reference value of "creep_characteristic_time","creep_modulus" (or the pair "creep_bulk" and "creep_shear"), and "temperature" must be given in the default values
 */
struct CreepArrheniusMaterialLaw : public ExternalMaterialLaw
{
    CreepArrheniusMaterialLaw(std::string args, char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s ) ;
};



} ;


#endif // __MATERIAL_LAW_H_
