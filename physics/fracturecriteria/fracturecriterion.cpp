//
// C++ Implementation: fracturecriterion
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../../mesher/delaunay.h"
#include "../../physics/viscoelasticity.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
#include <omp.h>
using namespace Amie ;

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) :
    neighbourhoodvolume(-1),
    physicalCharacteristicRadius(.008),
    scoreAtState(-1),
    nonLocalScoreAtState(-1),
    metAtStep(false),
    mirroring(mirroring), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z),
    deltaScoreAtState(0),
    criterionDamageDifferential(0),
    mesh2d(nullptr), mesh3d(nullptr),
    stable(true), checkpoint(true), inset(false),inIteration(false),
    scoreTolerance(1e-3),
    initialScore(1),
    cachedInfluenceRatio(1),
    minDeltaInNeighbourhood(1),
    maxScoreInNeighbourhood(0),
    maxModeInNeighbourhood(-1),
    maxAngleShiftInNeighbourhood(0),
    smoothingType(GAUSSIAN_NONCOMPACT),
    cacheID(-1), cachecoreID(-1)
{
}

Vector FractureCriterion::getSmoothedField(FieldType f0,  ElementState &s ,double t)
{
    if(mesh2d)
        return mesh2d->getSmoothedField(f0, cachecoreID, s.getParent(), 0, t) ;
    return mesh3d->getSmoothedField(f0, cachecoreID, s.getParent(), 0, t) ;
    
}

std::pair<Vector, Vector> FractureCriterion::getSmoothedFields( FieldType f0, FieldType f1, ElementState &s ,double t)
{
    if(mesh2d)
        return mesh2d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), 0, t) ;
    return mesh3d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), 0, t) ;
  
}

void FractureCriterion::initialiseCache( ElementState & s)
{
    if(s.getMesh2D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
        
        double overlap = (smoothingType == QUARTIC_COMPACT)?6.:8. ;
        Circle epsilonAll( std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*2. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Circle epsilonReduced(physicalCharacteristicRadius*1.1+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        mesh2d = s.getMesh2D() ;
        cacheID = mesh2d->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;
        cachecoreID = mesh2d->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
       
    }
    if(s.getMesh3D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function z = Function("z")-s.getParent()->getCenter().getZ() ;
        Function rr = x*x+y*y+z*z ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
        
        double overlap = (smoothingType == QUARTIC_COMPACT)?6.:8. ;
        Sphere epsilonAll( std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*2. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Sphere epsilonReduced(physicalCharacteristicRadius*1.1+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        mesh3d = s.getMesh3D() ;
        cacheID = mesh3d->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;
        cachecoreID = mesh3d->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        
    }
}

double FractureCriterion::getMaxScoreInNeighbourhood() const
{
    double maxScore = nonLocalScoreAtState ;
    for(size_t i = 0 ; i< cache.size() ; i++)
    {
        DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i])) ;
        if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
        {
            double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
            maxScore = std::max(maxScore, renormScore) ;
        }
    }
    return maxScore ;
}

std::pair<double, double> FractureCriterion::setChange( ElementState &s, double thresholdScore)
{
    stable = true ;

    if( !s.getParent()->getBehaviour()->getDamageModel())
        return std::make_pair(0.,0.) ;
    if(!mesh2d && !mesh3d)
    {
        initialiseCache(s);
    }
    
    if(mesh2d)
    {
// 		thresholdScore = POINT_TOLERANCE_2D ;
        // outside of the checkpoints, we only care about the order of the elements in
        // term of their score. At the checkpoint, we consider the elements which
        // have met their criterion
        if(checkpoint) //new iteration
        {
            inset = false ;
            inIteration = false ;
            damagingSet.clear();
            proximitySet.clear() ;

            std::vector<unsigned int> newSet ;
            std::set<unsigned int> newProximity ;
            std::multimap<double, DelaunayTriangle *> sortedElements ;
            maxScoreInNeighbourhood = nonLocalScoreAtState ;

            for(size_t i = 0 ; i< mesh2d->getCache(cacheID).size() ; i++)
            {
                DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(mesh2d->getCache(cacheID)[i])) ;
                if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
                {
                    double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
                    sortedElements.insert( std::make_pair(-renormScore, ci)) ;
                    maxScoreInNeighbourhood = std::max(maxScoreInNeighbourhood, renormScore) ;
                }
            }

            initialScore = 1. ;
            if(thresholdScore > 0 && s.getParent()->getState().getDeltaTime() > POINT_TOLERANCE_2D)
                initialScore = 1. ;//std::max(1.+thresholdScore, POINT_TOLERANCE_2D) ;
            double minscore = thresholdScore ;
            double maxscore = 0 ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            if(!sortedElements.empty() && thresholdScore > 0 )
            {
                for(auto i = sortedElements.begin() ; i != sortedElements.end() ; i++ )
                {
                    if(std::abs(-i->first-thresholdScore) <= scoreTolerance*initialScore*.25 && -i->first > 0)
                    {
                        if(i->second == s.getParent() && met())
                            inset = true ;
                        if(i->second->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE_2D)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getDelta()) ;
                        if(!i->second->getBehaviour()->getDamageModel())
                            continue ;
                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(i->second->index);
                        minscore = std::min(-i->first, thresholdScore) ;
                    }
                    else
                    {
                        newProximity.insert(i->second->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = -i->first ;
                            foundmaxscore = true ;
                        }
                    }
                }
            }
// 			std::cout << newSet.size() << std::endl ;

            if(!inset)
            {
                damagingSet.clear();
                proximitySet.clear() ;
                if(std::abs(nonLocalScoreAtState-thresholdScore) < 4.*scoreTolerance*initialScore)
                    inIteration = true ;
                return std::make_pair(0.,0.) ;
            }
            inIteration = true ;
            if(!newSet.empty())
                std::stable_sort(newSet.begin(), newSet.end());

            damagingSet = newSet ;
            proximitySet.insert(proximitySet.end(), newProximity.begin(), newProximity.end()) ;
            for(size_t i = 0 ; i < proximitySet.size() ; i++)
                static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i]))->getBehaviour()->getFractureCriterion()->inIteration = true ;

            return std::make_pair(minscore - maxscore /*+ scoreTolerance*2.*initialScore*/, thresholdScore - minscore /*- scoreTolerance*initialScore*/) ;
        }
        else if (inset)
        {
            checkpoint = false ;
            DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(damagingSet[0])) ;
            double maxscore = 0 ;
            if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged )
            {
                maxscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
            }
            maxModeInNeighbourhood = ci->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = ci->getBehaviour()->getDamageModel()->getAngleShift() ;

            for(size_t i = 1 ; i < damagingSet.size() ; i++)
            {
                ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(damagingSet[i])) ;
                if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged)
                {
                    double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
                    maxscore = std::min(maxscore,nls) ;
                }
                maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;

            }


            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                    minscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
//					std::cout << i << " " << std::flush ;
                    ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i])) ;
                    if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged)
                    {
                        double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }

            return std::make_pair(maxscore - minscore /*+ scoreTolerance*2.*initialScore*/, thresholdScore - maxscore /*- scoreTolerance*initialScore*/) ;
        }
    }
    else
    {
        // outside of the checkpoints, we only care about the order of the elements in
        // term of their score. At the checkpoint, we consider the elements which
        // have met their criterion

        if(checkpoint) //new iteration
        {
            inset = false ;
            inIteration = false ;
            damagingSet.clear();
            proximitySet.clear() ;

            std::vector<unsigned int> newSet ;
            std::set<unsigned int> newProximity ;
            std::multimap<double, DelaunayTetrahedron *> sortedElements ;
            maxScoreInNeighbourhood = nonLocalScoreAtState ;

            for(size_t i = 0 ; i< mesh3d->getCache(cacheID).size() ; i++)
            {
                DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(mesh3d->getCache(cacheID)[i])) ;
                if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
                {
                    double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
                    sortedElements.insert( std::make_pair(-renormScore, ci)) ;
                    maxScoreInNeighbourhood = std::max(maxScoreInNeighbourhood, renormScore) ;
                }
            }

            if(thresholdScore > 0 && s.getParent()->getState().getDeltaTime() > POINT_TOLERANCE_2D)
                initialScore = std::max(1.+thresholdScore, POINT_TOLERANCE_2D) ;
            double minscore = thresholdScore ;
            double maxscore = 0 ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            if(!sortedElements.empty() && thresholdScore > 0 )
            {
                for(auto i = sortedElements.begin() ; i != sortedElements.end() ; i++ )
                {
                    if(std::abs(-i->first-thresholdScore) <= scoreTolerance*initialScore && -i->first > 0)
                    {
                        if(i->second == s.getParent() && met())
                            inset = true ;
                        if(i->second->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE_2D)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getDelta()) ;
                        if(!i->second->getBehaviour()->getDamageModel())
                            continue ;
                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(i->second->index);
                        minscore = std::min(-i->first, thresholdScore) ;
                    }
                    else
                    {
                        newProximity.insert(i->second->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = -i->first ;
                            foundmaxscore = true ;
                        }
                    }
                }
            }
// 			std::cout << newSet.size() << std::endl ;

            if(!inset)
            {
                damagingSet.clear();
                proximitySet.clear() ;
                if(std::abs(nonLocalScoreAtState-thresholdScore) < 4.*scoreTolerance*initialScore)
                    inIteration = true ;
                return std::make_pair(0.,0.) ;
            }
            inIteration = true ;
            if(!newSet.empty())
                std::stable_sort(newSet.begin(), newSet.end());

            damagingSet = newSet ;
            proximitySet.insert(proximitySet.end(), newProximity.begin(), newProximity.end()) ;
            for(size_t i = 0 ; i < proximitySet.size() ; i++)
                static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i]))->getBehaviour()->getFractureCriterion()->inIteration = true ;

            return std::make_pair(minscore - maxscore /*+ scoreTolerance*2.*initialScore*/, thresholdScore - minscore/* - scoreTolerance*initialScore*/) ;
        }
        else if (inset && !s.getParent()->getBehaviour()->getDamageModel()->converged)
        {
//			std::cout << "a" << std::flush ;
            checkpoint = false ;
            DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;
            double maxscore = 0 ;
            if(ci->getBehaviour()->getFractureCriterion())
            {
                maxscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
            }
            maxModeInNeighbourhood = ci->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = ci->getBehaviour()->getDamageModel()->getAngleShift() ;
//			std::cout << "c" << std::flush ;

            for(size_t i = 1 ; i < damagingSet.size() ; i++)
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                {
                    double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
                    maxscore = std::min(maxscore,nls) ;
                }
                maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;

            }
//			std::cout << "d" << std::flush ;

            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                    minscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
//					std::cout << i << " " << std::flush ;
                    ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                    if(ci->getBehaviour()->getFractureCriterion())
                    {
                        double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }
//			std::cout << "b" << std::endl ;
        }

    }

    //shut up the compiler
    return std::make_pair(0., 0.) ;
}

void FractureCriterion::step(ElementState &s)
{
    if(cache.empty() )
        initialiseCache(s) ;


    if(s.getParent()->getBehaviour()->fractured())
    {
        scoreAtState = -1 ;
        return ;
    }
    if(checkpoint || inIteration)
    {
        scoreAtState = grade(s) ;
    }

}

void FractureCriterion::computeNonLocalState(ElementState &s)
{
    metAtStep = scoreAtState > 0 ;
    nonLocalScoreAtState = scoreAtState ;

}

FractureCriterion::~FractureCriterion()
{
}


void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
    physicalCharacteristicRadius = r ;
    cache.resize(0) ;
    factors.resize(0);
}

bool FractureCriterion::met() const
{
    return metAtStep ;
}

Material FractureCriterion::toMaterial()
{
    Material mat ;
    return mat ;
}



