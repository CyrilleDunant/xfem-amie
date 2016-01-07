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

typedef enum : char
{ 
	SET,
	ADD,
	MULTIPLY,
	SUBSTRACT,
	DIVIDE,
}  EMLOperation ;

/* Basic structure for material law. This structure is empty ; the preProcess() and/or step() methods must be surcharged for the law to do something */
/*SOURCE ExternalMaterialLaw */
struct ExternalMaterialLaw
{
    std::map<std::string, double> defaultValues ;

    ExternalMaterialLaw(std::string args = std::string(), char sep = ',') : defaultValues(parseDefaultValues(args, sep)) { }
    ExternalMaterialLaw(const ExternalMaterialLaw & law) : defaultValues(law.defaultValues) { } 
    virtual ~ExternalMaterialLaw() { } ;
    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s , double dt) { }
    virtual void preProcess( Matrix & stiffness, Point angle, planeType pt ) { } 
    virtual void step( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt) { }
    void setDefaultValue(std::string str, double v) ;

    bool has(std::string test) const { return defaultValues.find(test) != defaultValues.end() ; }
};

typedef std::vector<ExternalMaterialLaw *> ExternalMaterialLawList ;

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
/*PARSE Eval ExternalMaterialLaw
    @string[output] // name of the output
    @string[function] // natural form of the function
    @string<EMLOperation>[operation] SET // operation to apply
*/
struct EvalMaterialLaw: public SpaceTimeDependentExternalMaterialLaw
{
    std::map<char, std::string> coordinates ;
    bool useSpaceTimeCoordinates ;
    EvalMaterialLaw( std::string e, const char *f_, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_,a, args, sep), useSpaceTimeCoordinates(u) { }
    EvalMaterialLaw( std::string e, Function & f_, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), useSpaceTimeCoordinates(u) { }
    EvalMaterialLaw( std::string e, const char *f_, std::map<char, std::string> c, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    EvalMaterialLaw( std::string e, Function & f_, std::map<char, std::string> c, EMLOperation a = SET, std::string args = std::string(),bool u = true, char sep = ',' ) : SpaceTimeDependentExternalMaterialLaw(e,f_, a, args, sep), coordinates(c),useSpaceTimeCoordinates(u) { }
    EvalMaterialLaw( std::string out, std::string expr, EMLOperation a = SET) ;
    virtual ~EvalMaterialLaw() { } ;

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
/*PARSE LinearInterpolated ExternalMaterialLaw
    @string[output] // output parameter
    @string[input] // x coordinates for the linear interpolation
    @string[file_name] // address written as "file_name(input)"
    @string<EMLOperation>[operation] // operation to apply
*/
struct LinearInterpolatedMaterialLaw : public ExternalMaterialLaw
{
    std::pair<std::string, std::string> external ;
    std::pair<Vector, Vector> values ;
    EMLOperation op ;

    LinearInterpolatedMaterialLaw(std::pair<std::string, std::string> e,std::pair<Vector, Vector> v, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) : ExternalMaterialLaw(args, sep), external(e), values(v), op(o) { }
    LinearInterpolatedMaterialLaw(std::pair<std::string, std::string> e, std::string file, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) ;
    LinearInterpolatedMaterialLaw(std::string out, std::string in, std::string file, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) ;
    virtual ~LinearInterpolatedMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    double get(double x) const ;

};

/* Generic material law to set an external variable as a function of two other variables, using linear interpolations.
 * "input" contains the pair of input variables
 * "values" contains the list of values for the two input variables
 * "output" is the name of the output variable
 * "data" is the table of output values at each pair of input variable
 * The interpolation does not extrapolate outside the given bounds: the bounding values are used instead
 */
struct LinearBiInterpolatedExternalMaterialLaw : public ExternalMaterialLaw
{
    std::pair<std::string, std::string> input ;
    std::string output ;
    std::pair<Vector, Vector> values ;
    Matrix data ;
    EMLOperation op ;

    LinearBiInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string out, std::pair<Vector, Vector> v, Matrix & dat, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) : ExternalMaterialLaw(args, sep), input(e), output(out), values(v), data(dat), op(o) { }
    LinearBiInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string out, std::pair<Vector, Vector> v, std::string datafile, EMLOperation o = SET, std::string args = std::string(), char sep = ',' ) ;
    virtual ~LinearBiInterpolatedExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    double get(double x, double y) const ;
};

struct CopyFromFeatureTreeExternalMaterialLaw : public ExternalMaterialLaw
{
    FeatureTree  * source ;
    std::vector<std::string> external ;

    CopyFromFeatureTreeExternalMaterialLaw(FeatureTree * f, std::string e, std::string args = std::string(), char s = ',') ;
    virtual ~CopyFromFeatureTreeExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;

};


/*PARSE TimeDerivative ExternalMaterialLaw
    @string[parameter] // name of the parameter to derive
 */
struct TimeDerivativeMaterialLaw : public ExternalMaterialLaw
{
    std::string rate ;
    std::string base ;
    std::string previous ;

    TimeDerivativeMaterialLaw(std::string b, std::string r, std::string p, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), rate(r), base(b), previous(p) { }
    TimeDerivativeMaterialLaw(std::string b) : ExternalMaterialLaw(), rate(b+"_rate"), base(b), previous(b+"_previous") { }
    virtual ~TimeDerivativeMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE TimeIntegral ExternalMaterialLaw 
    @string[parameter] // name of the parameter to integrate
 */
struct TimeIntegralMaterialLaw : public ExternalMaterialLaw
{
    std::string integral ;
    std::string base ;

    TimeIntegralMaterialLaw(std::string b, std::string i, double init = 0., std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), integral(i),  base(b){ defaultValues[integral] = init ; }
    TimeIntegralMaterialLaw(std::string b) : ExternalMaterialLaw(), integral(b+"_integral"),  base(b){ defaultValues[integral] = 0 ; }
    virtual ~TimeIntegralMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE Minimum ExternalMaterialLaw
    @string[output] // parameter in which the output is stored
    @stringlist[parameters] // list of parameters which minima are taken
    @string<EMLOperation>[operation] SET // operation to apply
*/
struct MinimumMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::vector<std::string> coordinates ;
    EMLOperation op ;

    MinimumMaterialLaw(std::string out, std::vector<std::string> coord, EMLOperation o = SET, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(out), coordinates(coord), op(o) { } ;
    virtual ~MinimumMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE StoreMaximumValue ExternalMaterialLaw 
    @string[parameter] // parameter which maximum value is recorded
    @value[starting_value] 0 // initial value of the minimum 
*/
struct StoreMaximumValueMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::string max ;
    double start ;
    
    StoreMaximumValueMaterialLaw( std::string in, double s = 0., std::string append = std::string("_maximum"), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(in), max(in+append), start(s) { } ;

    virtual ~StoreMaximumValueMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE GetParticleOrientation ExternalMaterialLaw 
    @string[output] angle // parameter in which the orientation is stored
*/
struct GetParticleOrientationMaterialLaw : public ExternalMaterialLaw
{
    std::string variable ;
    bool orthogonal ;

    GetParticleOrientationMaterialLaw( std::string var = std::string("angle"), bool ortho = false, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), variable(var), orthogonal(ortho) { } ;
    GetParticleOrientationMaterialLaw( std::string var, EMLOperation op) : ExternalMaterialLaw(), variable(var), orthogonal(false) { } ;

    virtual ~GetParticleOrientationMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE StorePreviousValue ExternalMaterialLaw 
    @string[parameter] // parameter which previous value is recorded
*/
struct StorePreviousValueMaterialLaw : public ExternalMaterialLaw
{
    std::vector<std::string>  external ;
    std::string append ;
    
    StorePreviousValueMaterialLaw( std::vector<std::string> ext, std::string app = std::string("_previous"), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(ext), append(app) { } ;
    StorePreviousValueMaterialLaw( std::string ext, std::string app = std::string("_previous"), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), append(app) { external.push_back(ext) ; } ;
    virtual ~StorePreviousValueMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;


/*PARSE StoreMinimumValue ExternalMaterialLaw 
    @string[parameter] // parameter which minimum value is recorded
    @value[starting_value] 0 // initial value of the minimum 
*/
struct StoreMinimumValueMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::string max ;
    double start ;
    
    StoreMinimumValueMaterialLaw( std::string in, double s = 0., std::string append = std::string("_minimum"), std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(in), max(in+append), start(s) { } ;

    virtual ~StoreMinimumValueMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE Maximum ExternalMaterialLaw 
    @string[output] // parameter in which the output is stored
    @stringlist[parameters] // list of parameters which maxima are taken
    @string<EMLOperation>[operation] SET // operation to apply
*/
struct MaximumMaterialLaw : public ExternalMaterialLaw
{
    std::string external ;
    std::vector<std::string> coordinates ;
    EMLOperation op ;

    MaximumMaterialLaw(std::string out, std::vector<std::string> coord, EMLOperation o = SET, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), external(out), coordinates(coord), op(o) { } ;
    virtual ~MaximumMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE ExponentialDecay ExternalMaterialLaw 
    @string[output] // parameter in which the output is stored
    @string[parameter] // target towards which the output tends to with an exponential decay
*/
struct ExponentialDecayMaterialLaw : public ExternalMaterialLaw
{
    std::string output ;
    std::string target ;

    ExponentialDecayMaterialLaw(std::string out, std::string tar, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), output(out), target(tar) { } ;
    virtual ~ExponentialDecayMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE GetField ExternalMaterialLaw 
    @string[parameter] // name of the field to extract
*/
struct GetFieldMaterialLaw : public ExternalMaterialLaw
{
    FieldType field ;
    std::string base ;

    GetFieldMaterialLaw(std::string field, std::string args = std::string(), char sep = ',')  ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    virtual ~GetFieldMaterialLaw() { } ;

    std::string getParameterName(size_t i) const ;
} ;

/*PARSE WeibullDistributed ExternalMaterialLaw 
    @stringlist[parameters] // list of parameters that are affected by the Weibull distribution
    @string[weibull_variable] weibull_variable // name of the parameter in which the weibull variable is stored
*/
struct WeibullDistributedMaterialLaw : public ExternalMaterialLaw
{
    std::vector<std::string> affected ;
    std::string weib ;
    std::string shape ;
    std::string scale ;
    std::string variability ;

    WeibullDistributedMaterialLaw( std::string a, std::string w, std::string args = std::string(), char sep = ',') ;
    WeibullDistributedMaterialLaw( std::vector<std::string> aff, std::string w, std::string args = std::string(), char sep = ',') ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    virtual ~WeibullDistributedMaterialLaw() { } ;

    std::pair<std::string, double> getWeibullVariable() ;
    

} ;

/*PARSE UniformDistributedPerParticle ExternalMaterialLaw
    @string[output] // name of the parameter set to a random value
    @value[minimum] 0 // minimum value of the distribution
    @value[maximum] 1 // minimum value of the distribution
    @string<EMLOperation>[operation] SET // operation to apply
*/
struct UniformDistributedPerParticleMaterialLaw : public ExternalMaterialLaw
{
    std::string variable ;
    std::map<const Geometry *, double> values ;
    double minimum ;
    double maximum ;
    EMLOperation op ;

    UniformDistributedPerParticleMaterialLaw( std::string aff, double min, double max, EMLOperation o = SET, std::string args = std::string(), char sep = ',') ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
    virtual ~UniformDistributedPerParticleMaterialLaw() { } ;
} ;


} 


#endif // __MATERIAL_LAW_H_
