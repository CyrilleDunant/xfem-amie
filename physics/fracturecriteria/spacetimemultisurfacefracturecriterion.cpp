//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "spacetimemultisurfacefracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include <fstream>

namespace Amie {



FractureCriterion * SpaceTimeMultiSurfaceFractureCriterion::getCopy() const 
{
    SpaceTimeMultiSurfaceFractureCriterion * ret = new SpaceTimeMultiSurfaceFractureCriterion( ) ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
        ret->add( surfaces[i]->getCopy() ) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

double SpaceTimeMultiSurfaceFractureCriterion::grade(ElementState &s)
{    
    double g = -1 ;
    double ms = -1 ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
    {
//        std::cout << g << " ";
        g = std::max( g, surfaces[i]->grade(s)) ;
        ms = std::max( ms, surfaces[i]->getScoreAtTimeStepEnd()  ) ;
    }
//    std::cout << g << std::endl ;
    scoreAtTimeStepEnd = ms ;
    return g ;
/*

    double atEnd = -1 ;
    double time = -1 ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
    {
        double localTime = -1 ;
        double localEnd = -1 ;
        if( !s.getParent()->getBehaviour()->getDamageModel()->fractured(i) )
        {
            localTime = surfaces[i]->grade(s) ;
            localEnd = surfaces[i]->getScoreAtTimeStepEnd() ;
        }

        instants[i] = localTime ;
        
        if( localTime > time )
            time = localTime ;
        if( localEnd > atEnd )
            atEnd = localEnd ;
    }

    scoreAtTimeStepEnd = atEnd ;

    return time ;*/
}

double SpaceTimeMultiSurfaceFractureCriterion::gradeAtTime(ElementState &s, double t)  
{
    double g = -1 ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
        g = std::max( g, surfaces[i]->gradeAtTime(s, t) ) ;
    return g ;
}

void SpaceTimeMultiSurfaceFractureCriterion::initialiseCache( ElementState& s ) 
{
    this->FractureCriterion::initialiseCache(s) ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
        surfaces[i]->initialiseCache(s) ;
}

void SpaceTimeMultiSurfaceFractureCriterion::updateCache( ElementState & s)
{
    this->FractureCriterion::updateCache(s) ;
    for(size_t i = 0 ; i < surfaces.size() ; i++)
        surfaces[i]->updateCache(s) ;
}

bool SpaceTimeMultiSurfaceFractureCriterion::directionMet(size_t direction, double t) 
{
//    std::cout << instants[0] << "\t" << instants[1] << "\t" << t << std::endl ;
    return t >= instants[direction] ;
}



}

