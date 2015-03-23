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
#include "../../utilities/random.h"

namespace Amie
{

/* Parser for the default values for material laws
 * The string must be ordered as a series of "parameter = value" separated with the appropriate character.
 * Spaces are removed during parsing.
 */
std::map<std::string, double> parseDefaultValues(std::string args, char sep = ',') ;

typedef enum : char
{ 
	SET,
	ADD,
	MULTIPLY,
	SUBSTRACT,
	DIVIDE,
}  EMLOperation ;

/* Basic structure for material law. This structure is empty ; the preProcess() and/or step() methods must be surcharged for the law to do something */
struct ExternalMaterialLaw
{
    std::map<std::string, double> defaultValues ;

    ExternalMaterialLaw(std::string args = std::string(), char sep = ',') : defaultValues(parseDefaultValues(args, sep)) { }
    ExternalMaterialLaw(const ExternalMaterialLaw & law) : defaultValues(law.defaultValues) { } 
    virtual ~ExternalMaterialLaw() { } ;
    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s , double dt) { }
    virtual void step( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt) { }
    void setDefaultValue(std::string str, double v) ;

    bool has(std::string test) const { return defaultValues.find(test) != defaultValues.end() ; }
};

/* Generic material law to set one or several variables to a constant
 * The variables affected and their values are the variables specified in the arguments
 */
struct ConstantExternalMaterialLaw : public ExternalMaterialLaw
{
    ConstantExternalMaterialLaw( std::string args, char sep = ',' ) : ExternalMaterialLaw(args, sep) { }
    virtual ~ConstantExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

};

/* Generic material law to set an external variable from an analytical function of the space-time coordinates */
struct SpaceTimeDependentExternalMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    Function f ;
    EMLOperation op ;

    SpaceTimeDependentExternalMaterialLaw( std::string e, const char *f_, EMLOperation o = SET, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_), op(o) { }
    SpaceTimeDependentExternalMaterialLaw( std::string e, Function & f_, EMLOperation o = SET, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(e), f(f_), op(o) { }
    virtual ~SpaceTimeDependentExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;


};

/* Generic material law to set an external variable as a function of another one variable */
struct SimpleDependentExternalMaterialLaw: public SpaceTimeDependentExternalMaterialLaw
{
    std::string coordinate ;

    SimpleDependentExternalMaterialLaw( std::string e, std::string c, const char *f_,  EMLOperation a = SET,  std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, a,args, sep), coordinate(c) { }
    SimpleDependentExternalMaterialLaw( std::string e, std::string c,  Function & f_, EMLOperation a = SET, std::string args = std::string(), char sep = ',') : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinate(c) { }
    virtual ~SimpleDependentExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

} ;

struct AssignExternalMaterialLaw: public ExternalMaterialLaw
{
    std::string input ;
    std::string target ;

    AssignExternalMaterialLaw( std::string e, std::string c, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), input(e), target(c) { }
    virtual ~AssignExternalMaterialLaw() { } ;

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
    VariableDependentExternalMaterialLaw( std::string e, const char *f_, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_,a, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, const char *f_, std::map<char, std::string> c, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    VariableDependentExternalMaterialLaw( std::string e, Function & f_, std::map<char, std::string> c, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    virtual ~VariableDependentExternalMaterialLaw() { } ;

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
    EMLOperation op ;

    LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e,std::pair<Vector, Vector> v, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) : ExternalMaterialLaw(args, sep), external(e), values(v), op(o) { }
    LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string file, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) ;
    virtual ~LinearInterpolatedExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    double get(double x) const ;

};

struct CopyFromFeatureTreeExternalMaterialLaw : public ExternalMaterialLaw
{
    FeatureTree  * source ;
    std::vector<std::string> external ;

    CopyFromFeatureTreeExternalMaterialLaw(FeatureTree * f, std::string e, std::string args = std::string(), char s = ',') ;
    virtual ~CopyFromFeatureTreeExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

};

struct TimeDerivativeMaterialLaw : public ExternalMaterialLaw
{
    std::string rate ;
    std::string base ;
    double previous ;

    TimeDerivativeMaterialLaw(std::string b, std::string r, double init = 0., std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), rate(r), base(b), previous(init) { }
    virtual ~TimeDerivativeMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct TimeIntegralMaterialLaw : public ExternalMaterialLaw
{
    std::string integral ;
    std::string base ;

    TimeIntegralMaterialLaw(std::string b, std::string i, double init = 0., std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), integral(i),  base(b){ defaultValues[integral] = init ; }
    virtual ~TimeIntegralMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct MinimumMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::vector<std::string> coordinates ;
    EMLOperation op ;

    MinimumMaterialLaw(std::string out, std::vector<std::string> coord, EMLOperation o = SET, std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), external(out), coordinates(coord), op(o) { } ;
    virtual ~MinimumMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct MaximumMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::vector<std::string> coordinates ;
    EMLOperation op ;

    MaximumMaterialLaw(std::string out, std::vector<std::string> coord, EMLOperation o = SET, std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), external(out), coordinates(coord), op(o) { } ;
    virtual ~MaximumMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct ExponentiallyDecreasingMaterialLaw : public ExternalMaterialLaw
{
    std::string output ;
    std::string target ;
    std::string coefficient ;

    ExponentiallyDecreasingMaterialLaw(std::string out, std::string tar, std::string coef, std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), output(out), target(tar), coefficient(coef) { } ;
    virtual ~ExponentiallyDecreasingMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct GetFieldMaterialLaw : public ExternalMaterialLaw
{
    FieldType field ;
    std::string base ;

    GetFieldMaterialLaw(FieldType f, std::string b, std::string args = std::string(), char sep = 'c') : ExternalMaterialLaw(args, sep), field(f), base(b) { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    virtual ~GetFieldMaterialLaw() { } ;

    std::string getParameterName(size_t i) const ;
} ;

struct WeibullDistributedMaterialLaw : public ExternalMaterialLaw
{
    std::vector<std::string> affected ;
    std::string weib ;
    double shape ;
    double scale ;
    double variability ;

    WeibullDistributedMaterialLaw( std::string a, std::string w, double sh = 5, double sc = 1, double v = 0.2, std::string args = std::string(), char sep = 'c') ;
    WeibullDistributedMaterialLaw( std::vector<std::string> aff, std::string w, double sh = 5, double sc = 1, double v = 0.2, std::string args = std::string(), char sep = 'c') ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    virtual ~WeibullDistributedMaterialLaw() { } ;
    

} ;

} 


#endif // __MATERIAL_LAW_H_
