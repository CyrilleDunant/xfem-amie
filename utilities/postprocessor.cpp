//
// C++ Implementation: post-processor
//
// Description: 
//
//
// Author: Alain Giorla <alain.giorla@gmail.com>, (C) 2008-2016
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#include <fstream>
#include "enumeration_translator.h"
#include "postprocessor.h"

using namespace Amie ;

void PostProcessor::write( std::string file, FeatureTree * F, std::vector<PostProcessor *> posts, bool console, bool erase ) 
{
    std::fstream out ;
    if(erase)
        out.open( file.c_str(), std::ios::out ) ;
    else
        out.open( file.c_str(), std::ios::out | std::ios::app ) ;

    out << F->getCurrentTime() << "\t" ;
    if(console) { std::cout << F->getCurrentTime() << "\t" ; } 

    for(size_t i = 0 ; i < posts.size() ; i++)
    {
        Vector v = posts[i]->postProcess(F) ;
        for(size_t j = 0 ; j < v.size() ; j++)
        {
            out << v[j] << "\t" ;
            if(console) { std::cout << v[j] << "\t" ; } 
        }
    }

    out << std::endl ;
    if(console) { std::cout << std::endl ; }

    out.close() ;
}

FieldPostProcessor::FieldPostProcessor( std::string var, int c, double t) : PostProcessor(c,t), field( DISPLACEMENT_FIELD ), variable(std::string())
{
     bool isFieldType = true ;
     FieldType ft = Enum::getFieldType( var, &isFieldType ) ;
     if( isFieldType )
         field = ft ;
     else
         variable = var ;
}

Vector AverageFieldPostProcessor::postProcess( FeatureTree * F ) 
{
    if(cacheIndex < 0) { cacheIndex = F->get2DMesh()->getAllElementsCacheID() ; }
    if( variable.length() == 0)
    {
        if(all)
            return F->getAverageField( field, instant ) ;
        return F->get2DMesh()->getField( field, cacheIndex, instant ) ;
    }
    Vector ret(1) ; ret[0] = 0 ;
    ret[0] = F->get2DMesh()->getField( variable, cacheIndex, correction ) ;
    return ret ;
}

Vector MinimumFieldPostProcessor::postProcess( FeatureTree * F ) 
{
    if(cacheIndex < 0) { cacheIndex = F->get2DMesh()->getAllElementsCacheID() ; }
    if(variable.length() == 0)
    {
        if(all)
            return F->getFieldMinMax( field, instant ).first ;
    
        Vector ret(0) ;
        VirtualMachine vm ;
        size_t cacheSize = F->get2DMesh()->getCache( cacheIndex ).size() ;
        for(size_t i = 0 ; i < cacheSize; i++)
        {
            DelaunayTriangle * elem = F->get2DMesh()->getElement( cacheIndex, i ) ;
            Vector f(0) ;
            elem->getState().getAverageField( field, f, &vm, instant ) ;
            if( ret.size() == 0 )
                ret.resize(f.size()) ;
            for( size_t j = 0 ; j < f.size() ; j++)
                ret[j] = std::min( ret[j], f[j] ) ;
        }
        return ret ;
    }

    Vector ret(1) ; ret[0] = 0 ;
    bool init = false ;
    size_t cacheSize = F->get2DMesh()->getCache( cacheIndex ).size() ;
    std::map<std::string, double> dummy ;
    for(size_t i = 0 ; i < cacheSize ; i++)
    {
        GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * state = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables *>(&( F->get2DMesh()->getElement( cacheIndex, i )->getState())) ;
        if(state == nullptr) { continue ; }
        if(!state->has( variable )) { continue ; }
        if(!init) { ret[0] = state->get( variable, dummy ) ; }
        else { ret[0] = std::min( ret[0], state->get( variable, dummy ) ) ; }
        init = true ;
    }
    return ret ;
}

Vector MaximumFieldPostProcessor::postProcess( FeatureTree * F ) 
{
    if(cacheIndex < 0) { cacheIndex = F->get2DMesh()->getAllElementsCacheID() ; }
    if(variable.length() == 0)
    {
        if(all)
            return F->getFieldMinMax( field, instant ).second ;
    
        Vector ret(0) ;
        VirtualMachine vm ;
        size_t cacheSize = F->get2DMesh()->getCache( cacheIndex ).size() ;
        for(size_t i = 0 ; i < cacheSize; i++)
        {
            DelaunayTriangle * elem = F->get2DMesh()->getElement( cacheIndex, i ) ;
            Vector f(0) ;
            elem->getState().getAverageField( field, f, &vm, instant ) ;
            if( ret.size() == 0 )
                ret.resize(f.size()) ;
            for( size_t j = 0 ; j < f.size() ; j++)
                ret[j] = std::max( ret[j], f[j] ) ;
        }
        return ret ;
    }

    Vector ret(1) ; ret[0] = 0 ;
    bool init = false ;
    size_t cacheSize = F->get2DMesh()->getCache( cacheIndex ).size() ;
    std::map<std::string, double> dummy ;
    for(size_t i = 0 ; i < cacheSize ; i++)
    {
        GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * state = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables *>(&( F->get2DMesh()->getElement( cacheIndex, i )->getState())) ;
        if(state == nullptr) { continue ; }
        if(!state->has( variable )) { continue ; }
        if(!init) { ret[0] = state->get( variable, dummy ) ; }
        else { ret[0] = std::max( ret[0], state->get( variable, dummy ) ) ; }
        init = true ;
    }
    return ret ;
}


Vector LocalFieldPostProcessor::postProcess( FeatureTree * F ) 
{
    if( trg == nullptr )
    {
        std::vector<DelaunayTriangle *> candidates = F->get2DMesh()->getConflictingElements( &p ) ;
	if(candidates.size() > 0) { trg = candidates[0] ; }
        for(size_t i = 0 ; i < candidates.size() ; i++)
        {
            if(candidates[i]->in(p)) { trg = candidates[i] ; }
        }
    }
    
    Vector ret(0) ;
    if(trg == nullptr) { return ret ; }

    if( variable.length() == 0 )
    {
        VirtualMachine vm ;
        trg->getState().getAverageField( field, ret, &vm, instant ) ;
        return ret ;
    }
    ret.resize(1) ; ret[0] = 0 ;
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * state = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables *>(&( trg->getState())) ;
    std::map<std::string, double> dummy ;
    if(state == nullptr ) { return ret ; }
    if(! state->has( variable )) { return ret ; } 
    ret[0] = state->get( variable, dummy ) ; 
    return ret ;
}

PostProcessor * ProfileFieldPostProcessor::getPostProcessor(size_t i) 
{
    Point p = start ;
    if(i >= div) { p = end ; }
    else if(i > 0) { p = start + (end-start)*((double) i)/((double) div ) ; }
    if(variable.length() == 0)
        return new LocalFieldPostProcessor( field, p.x, p.y, instant ) ;
    else
        return new LocalFieldPostProcessor( variable, p.x, p.y, instant ) ;
    return new DoNothingPostProcessor() ;
}

Vector ProfileFieldPostProcessor::postProcess( FeatureTree * F )
{
    if(gauges.size() == 0)
    {
        for(size_t i = 0 ; i < div+1 ; i++)
            gauges.push_back( getPostProcessor(i) ) ;
    }

    std::vector<double> tmp ;
    for(size_t i = 0 ; i < gauges.size() ; i++)
    {
        Vector v = gauges[i]->postProcess(F) ;
        for(size_t j = 0 ; j < v.size() ; j++)
            tmp.push_back(v[j]) ;
    }
    Vector ret(tmp.size()) ; ret = 0 ;
    for(size_t j = 0 ; j < tmp.size() ; j++)
        ret[j] = tmp[j] ;
    return ret ;    
}

Vector LinearStrainGaugePostProcessor::postProcess( FeatureTree * F ) 
{
    Vector dleft = left.postProcess(F) ;
    Vector dright = right.postProcess(F) ;

    double base = (left.p-right.p).norm() ; 
    Point pleft = left.p+Point( dleft[0], dleft[1] ) ;
    Point pright = left.p+Point( dright[0], dright[1] ) ;
    double stretch = (pleft-pright).norm() ;

    Vector ret(1) ;
    ret[0] = stretch/base - 1. ;
    return ret ;
}


std::pair<Point*, double> MacroscopicStrainPostProcessor::getClosestBoundingPoint( DelaunayTriangle * tri, Point p, double instant, Variable axis )
{
    double distance = tri->getRadius()*100. ;
    double t = tri->getState().getNodalCentralTime() ;
    double dt = tri->getState().getNodalDeltaTime() ;
    size_t j = tri->getBoundingPoints().size() + 1 ;
    for(size_t i = 0 ; i < tri->getBoundingPoints().size() ; i++)
    {
        if( (tri->timePlanes() == 1) || (std::abs( (t+ (instant*dt*0.5)) - tri->getBoundingPoint(i).getT() ) < POINT_TOLERANCE) )
        {
            if( (axis == XI) && (std::abs( p.getX() - tri->getBoundingPoint(i).getX() ) > POINT_TOLERANCE ) ) { continue ; }
            if( (axis == ETA) && (std::abs( p.getY() - tri->getBoundingPoint(i).getY() ) > POINT_TOLERANCE ) ) { continue ; }
            if( (axis == ZETA) && (std::abs( p.getZ() - tri->getBoundingPoint(i).getZ() ) > POINT_TOLERANCE ) ) { continue ; }

            double d = dist( p, tri->getBoundingPoint(i) ) ;            
            if( d < distance )
            {
                distance = d ;
                j = i ;
            }
        }
    }
    if( j < tri->getBoundingPoints().size() )
        return std::make_pair( &( tri->getBoundingPoint(j) ), distance ) ;
    return std::make_pair( nullptr, distance ) ;
}

int MacroscopicStrainPostProcessor::getTargetID( FeatureTree * F, Point p, double instant, Variable axis ) 
{
    std::vector<DelaunayTriangle *> candidates = F->get2DMesh()->getConflictingElements( &p ) ;
    if(candidates.size() == 0) { return -1 ; }
    std::pair<Point *, double> target = getClosestBoundingPoint( candidates[0], p, instant, axis ) ;
    for(size_t i = 1 ; i < candidates.size() ; i++)
    {
        std::pair<Point *, double> target_tmp = getClosestBoundingPoint( candidates[i], p, instant, axis ) ;
        if((target_tmp.second < target.second) && (target_tmp.first != nullptr))
            target = target_tmp ;
    }
    if( target.first == nullptr ) { return -1 ; }
    return target.first->getId() ;
}

Vector MacroscopicStrainPostProcessor::postProcess( FeatureTree * F )
{
    Vector ret(3) ; ret = 0 ;

    if(left < 0 || right < 0 || bottom < 0 || top < 0)
    {
        Feature * box = F->getFeature(0) ;
        std::vector<Point> bbox = box->getBoundingBox() ;
        dim = Point( std::abs(bbox[2].getX()-bbox[0].getX()), std::abs(bbox[0].getY()-bbox[2].getY()) ) ;

        if( left < 0 )
        {
            Point q( bbox[0].getX(), p.getY() ) ;
            left = MacroscopicStrainPostProcessor::getTargetID( F, q, instant, XI ) ;
            if( left < 0 ) { return ret ; }
        }

        if( right < 0 )
        {
            Point q( bbox[2].getX(), p.getY() ) ;
            right = MacroscopicStrainPostProcessor::getTargetID( F, q, instant, XI ) ;
            if( right < 0 ) { return ret ; }
        }

        if( bottom < 0 )
        {
            Point q( p.getX(), bbox[2].getY() ) ;
            bottom = MacroscopicStrainPostProcessor::getTargetID( F, q, instant, ETA ) ;
            if( bottom < 0 ) { return ret ; }
        }

        if( top < 0 )
        {
            Point q( p.getX(), bbox[0].getY() ) ;
            top = MacroscopicStrainPostProcessor::getTargetID( F, q, instant, ETA ) ;
            if( top < 0 ) { return ret ; }
        }
    }

    if(ndof == 0)
    {
        int id = F->get2DMesh()->getAllElementsCacheID() ;
        size_t i = 0 ;
        while(ndof == 0 && i < F->get2DMesh()->getCache(id).size() )
        {
            DelaunayTriangle * test = F->get2DMesh()->getElement(id, i) ;
            if(test == nullptr) { i++ ; continue ; }
            if(test->getBehaviour() == nullptr) { i++ ; continue ; }
            if(test->getBehaviour()->type == VOID_BEHAVIOUR) { i++ ; continue ; }
            ndof = test->getBehaviour()->getNumberOfDegreesOfFreedom() ;
            break ;
        }
        if( ndof == 0 ) { return ret ; }
    }

    if( dim.getX() < POINT_TOLERANCE || dim.getY() < POINT_TOLERANCE ) { return ret ; }

    Vector disp = F->getDisplacements() ;
    ret[0] = ( disp[ right*ndof ] - disp[ left*ndof ] ) / dim.getX() ;
    ret[1] = ( disp[ top*ndof+1 ] - disp[ bottom*ndof+1 ] ) / dim.getY() ;
    ret[2] = 0.5*( ( disp[ right*ndof+1 ] - disp[ left*ndof+1 ] ) / dim.getY() 
                 + ( disp[ top*ndof ] - disp[ bottom*ndof ] ) / dim.getX()  ) ;

    return ret ;

}

Vector MaximumMacroscopicStrainPostProcessor::postProcess( FeatureTree * F )
{
    if( gauges.size() < ndiv )
    {
        Feature * box = F->getFeature(0) ;
        std::vector<Point> bbox = box->getBoundingBox() ;
        double xmin = bbox[0].getX() ;
        double xmax = bbox[2].getX() ;
        double dx = (xmax-xmin) / (ndiv-1) ;
        double ymin = bbox[2].getY() ;
        double ymax = bbox[0].getY() ;
        double dy = (ymax-ymin) / (ndiv-1) ;

        for(size_t i = 0 ; i < ndiv ; i++)
            gauges.push_back( new MacroscopicStrainPostProcessor( xmin+dx*i, ymin+dy*i, instant ) ) ;
    }

    Vector ret(3) ; ret = 0 ;
    for(size_t i = 0 ; i < gauges.size() ; i++)
    {
        Vector str = gauges[i]->postProcess(F) ;
        if(std::abs(str[0]) > std::abs(ret[0])) { ret[0] = str[0] ; }
        if(std::abs(str[1]) > std::abs(ret[1])) { ret[1] = str[1] ; }
        if(std::abs(str[2]) > std::abs(ret[2])) { ret[2] = str[2] ; }
    }

    return ret ;
}




