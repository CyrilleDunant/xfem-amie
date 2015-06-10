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
    switch(op)
    {
    case SET:
        s.set(external, VirtualMachine().eval(f, p) ) ;
        break ;
    case ADD:
        s.add( external, VirtualMachine().eval(f, p) ) ;
        break ;
    case MULTIPLY:
        s.multiply( external, VirtualMachine().eval(f, p) ) ;
        break ;
    case SUBSTRACT:
        s.add( external, -1. * VirtualMachine().eval(f, p) ) ;
        break ;
    case DIVIDE:
        s.multiply( external, 1. / VirtualMachine().eval(f, p) ) ;
        break ;
    default:
        s.set(external, VirtualMachine().eval(f, p) ) ;
        break ;
    }
}

void SimpleDependentExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double x = s.get(coordinate, defaultValues) ;
    switch(op)
    {
    case SET:
        s.set(external, VirtualMachine().eval(f, x) ) ;
        break ;
    case ADD:
        s.add( external, VirtualMachine().eval(f, x) ) ;
        break ;
    case MULTIPLY:
        s.multiply( external, VirtualMachine().eval(f, x) ) ;
        break ;
    case SUBSTRACT:
        s.add( external, -1. * VirtualMachine().eval(f, x) ) ;
        break ;
    case DIVIDE:
        s.multiply( external, 1. / VirtualMachine().eval(f, x) ) ;
        break ;
    default:
        s.set(external, VirtualMachine().eval(f, x) ) ;
        break ;
    }
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
    if( has('x')) {
        p.setX( s.get(coordinates['x'], defaultValues)) ;
    }
    if( has('y')) {
        p.setY( s.get(coordinates['y'], defaultValues)) ;
    }
    if( has('z')) {
        p.setZ( s.get(coordinates['z'], defaultValues)) ;
    }
    if( has('t')) {
        p.setT( s.get(coordinates['t'], defaultValues)) ;
    }
    Point q ;
    if( has('u')) {
        p.setX( s.get(coordinates['u'], defaultValues)) ;
    }
    if( has('v')) {
        p.setY( s.get(coordinates['v'], defaultValues)) ;
    }
    if( has('w')) {
        p.setZ( s.get(coordinates['w'], defaultValues)) ;
    }
    switch(op)
    {
    case SET:
        s.set(external, VirtualMachine().eval(f, p, q) ) ;
        break ;
    case ADD:
        s.add( external, VirtualMachine().eval(f, p, q) ) ;
        break ;
    case MULTIPLY:
        s.multiply( external, VirtualMachine().eval(f, p, q) ) ;
        break ;
    case SUBSTRACT:
        s.add( external, -1. * VirtualMachine().eval(f, p, q) ) ;
        break ;
    case DIVIDE:
        s.multiply( external, 1. / VirtualMachine().eval(f, p, q) ) ;
        break ;
    default:
        s.set(external, VirtualMachine().eval(f, p, q) ) ;
        break ;
    }
}

LinearInterpolatedExternalMaterialLaw::LinearInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string file, EMLOperation o, std::string args, char sep) : ExternalMaterialLaw(args, sep), external(e), op(o)
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

LinearBiInterpolatedExternalMaterialLaw::LinearBiInterpolatedExternalMaterialLaw(std::pair<std::string, std::string> e, std::string out, std::pair<Vector, Vector> v, std::string file, EMLOperation o, std::string args, char sep) : ExternalMaterialLaw(args, sep), input(e), output(out), values(v), op(o)
{
    data.resize( values.first.size(), values.second.size() ) ;
    data = 0. ;

    std::fstream in(file) ;
    if(!in.is_open())
    {
        std::cout << "file " << file << " doesn't exist!" << std::endl ;
        exit(0) ;
    }

    double buffer ;

    size_t i = 0 ;
    size_t j = 0 ;
    do {
        in >> buffer ;
        data[i][j] = buffer ; 
        j++ ;
        if(j == data.numCols())
        {
            i++ ;
            j = 0 ;
        }
    } while(!in.eof() && i < data.numRows()) ;
}


void LinearInterpolatedExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double v = 0. ;
    if(external.first == std::string("x"))
        v = get(s.getParent()->getCenter().getX()) ;
    else if(external.first == std::string("y"))
        v = get(s.getParent()->getCenter().getY()) ;
    else if(external.first == std::string("z"))
        v = get(s.getParent()->getCenter().getZ()) ;
    else if(external.first == std::string("t"))
        v = get(s.getNodalCentralTime()) ;
    else
        v = get(s.get(external.first, defaultValues)) ;
    switch(op)
    {
    case SET:
        s.set(external.second, v) ;
        break ;
    case ADD:
        s.add(external.second, v) ;
        break ;
    case MULTIPLY:
        s.multiply(external.second, v) ;
        break ;
    case SUBSTRACT:
        s.add(external.second, -1.*v) ;
        break ;
    case DIVIDE:
        s.multiply(external.second, 1./v) ;
        break ;
    default:
        s.set(external.second, v) ;
        break ;

    }
}

void LinearBiInterpolatedExternalMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    double x = 0. ;
    if(input.first == std::string("x"))
        x = s.getParent()->getCenter().getX() ;
    else if(input.first == std::string("y"))
        x = s.getParent()->getCenter().getY() ;
    else if(input.first == std::string("z"))
        x = s.getParent()->getCenter().getZ() ;
    else if(input.first == std::string("t"))
        x = s.getNodalCentralTime() ;
    else
        x = s.get(input.first, defaultValues) ;

    double y = 0. ;
    if(input.second == std::string("x"))
        y = s.getParent()->getCenter().getX() ;
    else if(input.second == std::string("y"))
        y = s.getParent()->getCenter().getY() ;
    else if(input.second == std::string("z"))
        y = s.getParent()->getCenter().getZ() ;
    else if(input.second == std::string("t"))
        y = s.getNodalCentralTime() ;
    else
        y = s.get(input.second, defaultValues) ;

    double v = get(x,y) ;

    switch(op)
    {
    case SET:
        s.set(output, v) ;
        break ;
    case ADD:
        s.add(output, v) ;
        break ;
    case MULTIPLY:
        s.multiply(output, v) ;
        break ;
    case SUBSTRACT:
        s.add(output, -1.*v) ;
        break ;
    case DIVIDE:
        s.multiply(output, 1./v) ;
        break ;
    default:
        s.set(output, v) ;
        break ;

    }
}

double LinearInterpolatedExternalMaterialLaw::get(double x) const
{
    if(x < values.first[0]+POINT_TOLERANCE)
        return values.second[0] ;
    if(x > values.first[values.first.size()-1]-POINT_TOLERANCE)
        return values.second[values.second.size()-1] ;

    int i = 0 ;
    while(values.first[i]<x)
    {
        if(std::abs(values.first[i]-x) < POINT_TOLERANCE)
            return values.second[i] ;
        i++ ;
    }

    return values.second[i-1]+(values.second[i]-values.second[i-1])*(x-values.first[i-1])/(values.first[i]-values.first[i-1]) ;
}

double LinearBiInterpolatedExternalMaterialLaw::get(double x, double y) const
{
    Vector current( values.second ) ;
    int index = -1 ;
    double coef = 1 ;
    if(x < values.first[0]+POINT_TOLERANCE)
        index = 0 ;
    if(x > values.first[values.first.size()-1]-POINT_TOLERANCE)
        index = values.first.size()-1 ;

    if(index == -1)
    {
        int i = 0 ;
        bool at = false ;
        while(values.first[i]<x)
        {
            if(std::abs(values.first[i]-x) < POINT_TOLERANCE)
            {
                index = i ;
                at = true ;
                break ;
            }
            i++ ;
        }
        index = i ;
        if(!at)
            coef = (x-values.first[i-1])/(values.first[i]-values.first[i-1]) ;
    }

    for(size_t i = 0 ; i < current.size() ; i++)
    {
        if(coef == 1)
            current[i] = data[index][i] ;
        else
            current[i] = data[index-1][i] + coef*( data[index][i] - data[index-1][i] ) ;
    }

    if(y < values.second[0]+POINT_TOLERANCE)
        return current[0] ;
    if(y > values.second[current.size()-1]-POINT_TOLERANCE)
        return current[current.size()-1] ;

    int i = 0 ;
    while(values.second[i]<y)
    {
        if(std::abs(values.second[i]-y) < POINT_TOLERANCE)
            return current[i] ;
        i++ ;
    }

    return current[i-1]+(current[i]-current[i-1])*(y-values.second[i-1])/(values.second[i]-values.second[i-1]) ;

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
    if(!s.has(previous))
        s.set( previous, current ) ;
    double prev = s.get(previous, defaultValues) ;

    if(dt > POINT_TOLERANCE)
        s.set(rate, (current-prev)/dt) ;
    else
        s.set( rate, 0. ) ;
    s.set(previous, current) ;
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
    double v = values.max() ;
    switch(op)
    {
    case SET:
        s.set(external, v) ;
        break ;
    case ADD:
        s.add(external, v) ;
        break ;
    case MULTIPLY:
        s.multiply(external, v) ;
        break ;
    case SUBSTRACT:
        s.add(external, -1.*v) ;
        break ;
    case DIVIDE:
        s.multiply(external, 1./v) ;
        break ;
    default:
        s.set(external, v) ;
        break ;
    }
}


void MinimumMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    Vector values(coordinates.size()) ;
    values = 0. ;
    for(size_t i = 0 ; i < coordinates.size() ; i++)
        values[i] = s.get( coordinates[i], defaultValues) ;
    double v = values.min() ;
    switch(op)
    {
    case SET:
        s.set(external, v) ;
        break ;
    case ADD:
        s.add(external, v) ;
        break ;
    case MULTIPLY:
        s.multiply(external, v) ;
        break ;
    case SUBSTRACT:
        s.add(external, -1.*v) ;
        break ;
    case DIVIDE:
        s.multiply(external, 1./v) ;
        break ;
    default:
        s.set(external, v) ;
        break ;
    }
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

WeibullDistributedMaterialLaw::WeibullDistributedMaterialLaw( std::string a, std::string w, double sh, double sc, double v, std::string args, char sep) : ExternalMaterialLaw(args, sep), weib(w), shape(sh), scale(sc)
{
    affected.push_back(a) ;
    defaultValues[weib] = 1. ;
}

WeibullDistributedMaterialLaw::WeibullDistributedMaterialLaw( std::vector<std::string> aff, std::string w, double sh, double sc, double v, std::string args, char sep)  : ExternalMaterialLaw(args, sep), weib(w), shape(sh), scale(sc), variability(v)
{
    for(size_t i = 0 ; i < aff.size() ; i++)
        affected.push_back(aff[i]) ;
    defaultValues[weib] = 1. ;
}

void WeibullDistributedMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    if(!s.has(weib))
    {
        s.set(weib, std::max(0., 1. - variability + variability*RandomNumber().weibull(1.,5.)) ) ;
    }
    double w = s.get(weib, defaultValues) ;
    for(size_t i = 0 ; i < affected.size() ; i++)
        s.multiply(affected[i], w) ;
}

}

