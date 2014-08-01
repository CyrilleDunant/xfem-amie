#include "material_laws.h"
#include <stdlib.h>

namespace Amie
{

std::map<std::string, double> parseDefaultValues(std::string args, char sep)
{
    std::map<std::string, double> ret ;
    std::string options = args ;
    if(options.size() == 0 || options.find('=') == std::string::npos)
        return ret ;


    size_t foundSpace = options.find(' ') ;
    while(foundSpace != std::string::npos)
    {
        options.erase( foundSpace, 1) ;
        foundSpace = options.find(' ') ;
    }

    size_t foundSep = options.find(sep) ;
    while(foundSep != std::string::npos)
    {
        std::string current = options.substr(0, foundSep) ;
        std::string external = current.substr(0, current.find('=') ) ;
        double value = atof(current.substr( current.find('=')+1 ).c_str()) ;
        ret[external] = value ;
        options.erase(0, foundSep+1) ;
        foundSep = options.find(sep) ;
    }
    std::string external = options.substr(0, options.find('=') ) ;
    double value = atof(options.substr( options.find('=')+1 ).c_str()) ;
    ret[external] = value ;

    return ret ;
}

void ConstantExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s )
{
    for(auto iter = defaultValues.begin() ; iter != defaultValues.end() ; iter++)
        s.set(iter->first, iter->second) ;
}

void SpaceTimeDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s )
{
    Point p = s.getParent()->getCenter() ;
    p.setT( s.getTime() ) ;
    s.set(external, VirtualMachine().eval(f, p)) ;
}

void SimpleDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s )
{
    double x = s.get(coordinate, defaultValues) ;
    s.set(external, VirtualMachine().eval(f, x)) ;
}

void VariableDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s )
{
    Point p ;
    if(useSpaceTimeCoordinates)
    {
        p = s.getParent()->getCenter() ;
        p.setT( s.getTime() ) ;
    }
    if( has('x')) { p.setX( s.get(coordinates['x'], defaultValues)) ; }
    if( has('y')) { p.setY( s.get(coordinates['y'], defaultValues)) ; }
    if( has('z')) { p.setZ( s.get(coordinates['z'], defaultValues)) ; }
    if( has('t')) { p.setT( s.get(coordinates['t'], defaultValues)) ; }
    Point q ;
    if( has('u')) { p.setX( s.get(coordinates['u'], defaultValues)) ; }
    if( has('v')) { p.setY( s.get(coordinates['v'], defaultValues)) ; }
    if( has('w')) { p.setZ( s.get(coordinates['w'], defaultValues)) ; }
    s.set(external, VirtualMachine().eval(f, p, q)) ;
}

void ThermalExpansionMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s)
{
    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double alpha = s.get("thermal_expansion_coefficient", defaultValues) ;
    double imp = s.get("imposed_deformation", defaultValues) ;
    imp += alpha*(T-T0) ;
    s.set("imposed_deformation", imp) ;
}

ArrheniusMaterialLaw::ArrheniusMaterialLaw(std::string a, std::string args, char sep) : ExternalMaterialLaw(args, sep), affected(a)
{
    coefficient = affected ; coefficient.append("_activation_energy") ;
}

void ArrheniusMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s)
{
    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double Ea = s.get(coefficient, defaultValues) ;
    double prop = defaultValues[affected] ;
    prop *= exp( Ea*(1./T-1./T0) ) ;
    s.set(affected, prop) ;
}

void CreepArrheniusMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s)
{
    if(!s.has("creep_characteristic_time"))
        return ;

    double T = s.get("temperature", defaultValues) ;
    double T0 = defaultValues["temperature"] ;
    double Ea = s.get("creep_activation_energy", defaultValues) ;
    double factor = exp( Ea*(1./T-1./T0) ) ;
    double tau = defaultValues["creep_characteristic_time"] ;
    s.set("creep_characteristic_time", tau*factor) ;
    if(s.has("creep_modulus"))
    {
        double E = defaultValues["creep_modulus"] ;
        s.set("creep_modulus", E*factor) ;
    } else {
        double k = defaultValues["creep_bulk"] ;
        s.set("creep_bulk", k*factor) ;
        double mu = defaultValues["creep_shear"] ;
        s.set("creep_shear", mu*factor) ;
    }
}


} ;

