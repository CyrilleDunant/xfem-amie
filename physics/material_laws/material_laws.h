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
#include "../../features/features.h"

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
    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s , double dt) { }
    virtual void step( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt) { }
};

/* Generic material law to set one or several variables to a constant
 * The variables affected and their values are the variables specified in the arguments
 */
struct ConstantExternalMaterialLaw : public ExternalMaterialLaw
{
    ConstantExternalMaterialLaw( std::string args, char sep = ',' ) : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

};

/* Generic material law to set an external variable from an analytical function of the space-time coordinates */
struct SpaceTimeDependentExternalMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    Function f ;
    bool add ;

    SpaceTimeDependentExternalMaterialLaw( std::string e, const char *f_, bool a = false, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_), add(a) { }
    SpaceTimeDependentExternalMaterialLaw( std::string e, Function & f_, bool a = false, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_), add(a) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;


};

/* Generic material law to set an external variable as a function of another one variable */
struct SimpleDependentExternalMaterialLaw: public SpaceTimeDependentExternalMaterialLaw
{
    std::string coordinate ;

    SimpleDependentExternalMaterialLaw( std::string e, std::string c, const char *f_, bool a = false,  std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, a,args, sep), coordinate(c) { }
    SimpleDependentExternalMaterialLaw( std::string e, std::string c,  Function & f_,bool a = false, std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinate(c) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

} ;

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
    VariableDependentExternalMaterialLaw( std::string e, const char *f_, bool a = false, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_,a, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_, bool a = false, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, const char *f_, std::map<char, std::string> c, bool a = false, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_, std::map<char, std::string> c, bool a = false, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

    bool has(char c) const { return !(coordinates.empty() || coordinates.find(c) == coordinates.end()) ; }
    void setAsX(std::string s) { coordinates['x'] = s ; }
    void setAsY(std::string s) { coordinates['y'] = s ; }
    void setAsZ(std::string s) { coordinates['z'] = s ; }
    void setAsT(std::string s) { coordinates['t'] = s ; }
    void setAsU(std::string s) { coordinates['u'] = s ; }
    void setAsV(std::string s) { coordinates['v'] = s ; }
    void setAsW(std::string s) { coordinates['w'] = s ; }

};

/* Generic material law to set an external variable as a function of another one variable, using a linear interpolation.
 * "external" contains the pair of variable (first the input variable, second the variable to be set
 * "values" contains the list of values for the linear interpolation (first the input variable, second the output variable) [this vector can be read from a two-columns file, without delimiters]
 * The interpolation does not extrapolate outside the given bounds: the bounding values are used instead
 */
struct LinearInterpolatedExternalMaterialLaw : public ExternalMaterialLaw
{
    std::pair<std::string, std::string> external ;
    std::pair<Vector, Vector> values ;

    LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e,std::pair<Vector, Vector> v, std::string args = std::string(), char sep = ',' ) : ExternalMaterialLaw(args, sep), external(e), values(v) { }
    LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string file, std::string args = std::string(), char sep = ',' ) ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    double get(double x) const ;

};

struct CopyFromFeatureTreeExternalMaterialLaw : public ExternalMaterialLaw
{
    FeatureTree  * source ;
    std::vector<std::string> external ;

    CopyFromFeatureTreeExternalMaterialLaw(FeatureTree * f, std::string e, std::string args = std::string(), char s = ',') ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

};

struct TimeDerivativeMaterialLaw : public ExternalMaterialLaw
{
    std::string rate ;
    std::string base ;
    double previous ;

    TimeDerivativeMaterialLaw(std::string b, std::string r, double init = 0., std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), base(b), rate(r), previous(init) { }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct TimeIntegralMaterialLaw : public ExternalMaterialLaw
{
    std::string integral ;
    std::string base ;

    TimeIntegralMaterialLaw(std::string b, std::string i, double init = 0., std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), base(b), integral(i) { defaultValues[integral] = init ; }

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};


} ;


#endif // __MATERIAL_LAW_H_
