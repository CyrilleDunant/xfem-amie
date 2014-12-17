#include "material_laws.h"
#include "../../utilities/itoa.h"
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

void ExternalMaterialLaw::setDefaultValue(std::string str, double d)
{
	defaultValues[str] = d ;
}

void ConstantExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    for(auto iter = defaultValues.begin() ; iter != defaultValues.end() ; iter++)
        s.set(iter->first, iter->second) ;
}

void SpaceTimeDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    Point p = s.getParent()->getCenter() ;
    p.setT( s.getNodalCentralTime() ) ;
    double ret = VirtualMachine().eval(f, p) ;
    if(add)
        ret += s.get(external, defaultValues) ;
    s.set(external, ret) ;
}

void SimpleDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double x = s.get(coordinate, defaultValues) ;
    double ret = VirtualMachine().eval(f, x) ;
    if(add)
        ret += s.get(external, defaultValues) ;
    s.set(external, ret) ;
}

void AssignExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    s.set(target, s.get(input, defaultValues)) ;
}


void VariableDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    Point p ;
    if(useSpaceTimeCoordinates)
    {
        p = s.getParent()->getCenter() ;
        p.setT( s.getNodalCentralTime() ) ;
    }
    if( has('x')) { p.setX( s.get(coordinates['x'], defaultValues)) ; }
    if( has('y')) { p.setY( s.get(coordinates['y'], defaultValues)) ; }
    if( has('z')) { p.setZ( s.get(coordinates['z'], defaultValues)) ; }
    if( has('t')) { p.setT( s.get(coordinates['t'], defaultValues)) ; }
    Point q ;
    if( has('u')) { p.setX( s.get(coordinates['u'], defaultValues)) ; }
    if( has('v')) { p.setY( s.get(coordinates['v'], defaultValues)) ; }
    if( has('w')) { p.setZ( s.get(coordinates['w'], defaultValues)) ; }
    double ret = VirtualMachine().eval(f, p, q) ;
    if(add)
        ret += s.get(external, defaultValues) ;
    s.set(external, ret) ;
}

LinearInterpolatedExternalMaterialLaw::LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string file, std::string args, char sep) : ExternalMaterialLaw(args, sep), external(e)
{
    std::vector<double> x ;
    std::vector<double> y ;

    std::fstream input(file) ;
    if(!input.is_open())
    {
        std::cout << "file " << file << " doesn't exist!" << std::endl ;
        exit(0) ;
    }

    double buffer ;

    do {
        input >> buffer ;
        x.push_back(buffer) ;
        input >> buffer ;
        y.push_back(buffer) ;
    } while(!input.eof()) ;
    x.pop_back() ;
    y.pop_back() ;

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
        s.set(external.second, get(s.getNodalCentralTime())) ;
    else
        s.set(external.second, get(s.get(external.first, defaultValues))) ;
}

double LinearInterpolatedExternalMaterialLaw::get(double x) const
{
    if(x < values.first[0])
        return values.second[0] ;
    if(x > values.first[values.first.size()-1])
        return values.second[values.second.size()-1] ;

    int i = 0 ;
    while(values.first[i]<x)
    {
        if(std::abs(values.first[i]-x) < POINT_TOLERANCE_2D)
            return values.second[i] ;
        i++ ;
    }

    return values.second[i-1]+(values.second[i]-values.second[i-1])*(x-values.first[i-1])/(values.first[i]-values.first[i-1]) ;
}

CopyFromFeatureTreeExternalMaterialLaw::CopyFromFeatureTreeExternalMaterialLaw(FeatureTree * f, std::string e, std::string args, char s) : ExternalMaterialLaw(args, s), source(f)
{
    std::string options = e ;
    if(options.size() == 0)
        return ;

    size_t foundSpace = options.find(' ') ;
    while(foundSpace != std::string::npos)
    {
        options.erase( foundSpace, 1) ;
        foundSpace = options.find(' ') ;
    }

    size_t foundSep = options.find(s) ;
    while(foundSep != std::string::npos)
    {
        external.push_back(options.substr(0, foundSep)) ;
        options.erase(0, foundSep+1) ;
        foundSep = options.find(s) ;
    }
    external.push_back(options) ;

}

void CopyFromFeatureTreeExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables corresponding = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(dynamic_cast<DelaunayTriangle *>(source->get2DMesh()->getInTree(dynamic_cast<DelaunayTriangle *>(s.getParent())->index))->getState()) ;
    for(size_t i = 0 ; i < external.size() ; i++)
        s.set(external[i], corresponding.get(external[i], defaultValues)) ;
}

void TimeDerivativeMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double current = s.get(base, defaultValues) ;
    s.set(rate, (current-previous)/dt) ;
    previous = current ;
}

void TimeIntegralMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double current = s.get(integral, defaultValues) ;
    double rate = s.get(base, defaultValues) ;
    s.set(integral, current + rate*dt) ;
}

void MaximumMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
	Vector values(coordinates.size()) ;
	values = 0. ;
	for(size_t i = 0 ; i < coordinates.size() ; i++)
		values[i] = s.get( coordinates[i], defaultValues) ;
	double ret = values.max() ;
	if(add)
		ret += s.get(external, defaultValues) ;
	s.set(external, ret) ;
}

void MinimumMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
	Vector values(coordinates.size()) ;
	values = 0. ;
	for(size_t i = 0 ; i < coordinates.size() ; i++)
		values[i] = s.get( coordinates[i], defaultValues) ;
	double ret = values.min() ;
	if(add)
		ret += s.get(external, defaultValues) ;
	s.set(external, ret) ;
}

void ExponentiallyDecreasingMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
	double targetValue = s.get(target, defaultValues) ;
	double currentValue = s.get(output, defaultValues) ;
	double tau = s.get(coefficient, defaultValues) ;
	double nextValue = targetValue + (currentValue-targetValue)*(1-exp(-dt/tau)) ;
	s.set( output, nextValue ) ;

}

void GetFieldMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
	Vector v( fieldTypeElementarySize(field, SPACE_TWO_DIMENSIONAL, 1) ) ;
	s.getAverageField( field, v, nullptr, 0, 1.) ;
	for(size_t i = 0 ; i < v.size() ; i++)
	{
		s.set( getParameterName(i), v[i] ) ;
	}
}

std::string GetFieldMaterialLaw::getParameterName(size_t i) const 
{
	std::string b = base ;
	b.append("_") ;
	b.append(itoa(i)) ;
	return b ;
}



} ;

