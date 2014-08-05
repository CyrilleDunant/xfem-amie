#include "material_laws.h"
#include <stdlib.h>
#include <fstream>

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

void ConstantExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    for(auto iter = defaultValues.begin() ; iter != defaultValues.end() ; iter++)
        s.set(iter->first, iter->second) ;
}

void SpaceTimeDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    Point p = s.getParent()->getCenter() ;
    p.setT( s.getTime() ) ;
    s.set(external, VirtualMachine().eval(f, p)) ;
}

void SimpleDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double x = s.get(coordinate, defaultValues) ;
    s.set(external, VirtualMachine().eval(f, x)) ;
}

void VariableDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
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

LinearInterpolatedExternalMaterialLaw::LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string file, std::string args, char sep) : ExternalMaterialLaw(args, sep), external(e)
{
    std::vector<double> x ;
    std::vector<double> y ;

    std::fstream input ;
    input.open(file.c_str(), std::ios::in) ;
    double buffer ;

    while(!input.eof())
    {
        input >> buffer ;
        x.push_back(buffer) ;
        input >> buffer ;
        y.push_back(buffer) ;
    }

    Vector first(x.size()) ;
    Vector second(y.size()) ;
    for(size_t i = 0 ; i < x.size() ; i++)
    {
        first[i] = x[i] ;
        second[i] = y[i] ;
    }
    values = std::make_pair(first,second) ;

}


void LinearInterpolatedExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    if(external.first == std::string("x"))
        s.set(external.second, get(s.getParent()->getCenter().getX())) ;
    else if(external.first == std::string("y"))
        s.set(external.second, get(s.getParent()->getCenter().getY())) ;
    else if(external.first == std::string("z"))
        s.set(external.second, get(s.getParent()->getCenter().getZ())) ;
    else if(external.first == std::string("t"))
        s.set(external.second, get(s.getTime())) ;
    else
        s.set(external.second, get(s.get(external.second, defaultValues))) ;
}

double LinearInterpolatedExternalMaterialLaw::get(double x) const
{
    if(x < values.first[0])
        return values.second[0] ;
    if(x > values.first[values.first.size()-1])
        return values.second[values.second.size()-1] ;

    int i = 0 ;
    while(values.first[i]<x)
        i++ ;

    return values.second[i-1]+(values.second[i]-values.second[i-1])*(x-values.first[i-1])/(values.first[i]-values.first[i-1]) ;
}



} ;

