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
#include "../../polynomial/vm_function_extra.h"
#include <omp.h>
using namespace Amie ;

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) :
    physicalCharacteristicRadius(.008),
    scoreAtState(-1),
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
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;
    
    if(mesh2d)
        return mesh2d->getSmoothedField(f0, cachecoreID, s.getParent(), -1., t) ;
    return mesh3d->getSmoothedField(f0, cachecoreID, s.getParent(), -1., t) ;
    
}

std::pair<Vector, Vector> FractureCriterion::getSmoothedFields( FieldType f0, FieldType f1, ElementState &s ,double t)
{
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;
    
    if(mesh2d)
        return mesh2d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), -1., t) ;
    return mesh3d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), -1., t) ;
  
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
        Circle epsilonAll( std::max(physicalCharacteristicRadius*1.1, s.getParent()->getRadius()*3. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
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

double FractureCriterion::getMaxScoreInNeighbourhood(ElementState & s)
{
    double maxScore = scoreAtState ;
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;
    
    if(mesh2d)
    {
        for(size_t i = 0 ; i< mesh2d->getCache(cacheID).size() ; i++)
        {
            DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(mesh2d->getCache(cacheID)[i])) ;
            if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
            {
                double renormScore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
        }
    }
    else
    {
        for(size_t i = 0 ; i< mesh3d->getCache(cacheID).size() ; i++)
        {
            DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(mesh3d->getCache(cacheID)[i])) ;
            if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
            {
                double renormScore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
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
        std::cout << "cache init !" << std::endl ;
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

            initialScore = 1. ;
            double minscore = thresholdScore ;
            double maxscore = 0 ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            if(thresholdScore > 0 )
            {
                for(size_t i = 0 ; i< mesh2d->getCache(cacheID).size() ; i++)
                {
                    DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(mesh2d->getCache(cacheID)[i])) ;
                    if(!ci->getBehaviour()->getDamageModel())
                         continue ;
                        
                     if(ci->getBehaviour()->fractured())
                         continue ;
                     
                    if(std::abs(ci->getBehaviour()->getFractureCriterion()->scoreAtState-thresholdScore) <= scoreTolerance*initialScore*.25 && 
                        ci->getBehaviour()->getFractureCriterion()->scoreAtState > 0)
                    {
 
                        if(ci == s.getParent() && met())
                            inset = true ;
                        
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE_2D)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getDamageModel()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                        minscore = std::min(ci->getBehaviour()->getFractureCriterion()->scoreAtState, thresholdScore) ;
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                            foundmaxscore = true ;
                        }
                    }
                }
            };

            if(!inset)
            {
                damagingSet.clear();
                proximitySet.clear() ;
                if(std::abs(scoreAtState-thresholdScore) < 4.*scoreTolerance*initialScore)
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
                maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
            }
            maxModeInNeighbourhood = ci->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = ci->getBehaviour()->getDamageModel()->getAngleShift() ;

            for(size_t i = 1 ; i < damagingSet.size() ; i++)
            {
                ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(damagingSet[i])) ;
                if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged)
                {
                    double nls = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
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
                    minscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
                    ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i])) ;
                    if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged)
                    {
                        double nls = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

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

            initialScore = 1. ;
            double minscore = thresholdScore ;
            double maxscore = 0 ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            if(thresholdScore > 0 )
            {
                for(size_t i = 0 ; i< mesh3d->getCache(cacheID).size() ; i++)
                {
                    DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(mesh3d->getCache(cacheID)[i])) ;
                    if(!ci->getBehaviour()->getDamageModel())
                         continue ;
                        
                     if(ci->getBehaviour()->fractured())
                         continue ;
                     
                    if(std::abs(ci->getBehaviour()->getFractureCriterion()->scoreAtState-thresholdScore) <= scoreTolerance*initialScore*.25 && 
                        ci->getBehaviour()->getFractureCriterion()->scoreAtState > 0)
                    {
                        if(ci == s.getParent() && met())
                            inset = true ;
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE_2D)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getDamageModel()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                        minscore = std::min(ci->getBehaviour()->getFractureCriterion()->scoreAtState, thresholdScore) ;
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                            foundmaxscore = true ;
                        }
                    }
                }
            };

            if(!inset)
            {
                damagingSet.clear();
                proximitySet.clear() ;
                if(std::abs(scoreAtState-thresholdScore) < 4.*scoreTolerance*initialScore)
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

            return std::make_pair(minscore - maxscore /*+ scoreTolerance*2.*initialScore*/, thresholdScore - minscore /*- scoreTolerance*initialScore*/) ;
        }
        else if (inset && !s.getParent()->getBehaviour()->getDamageModel()->converged)
        {
            checkpoint = false ;
            DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;
            double maxscore = 0 ;
            if(ci->getBehaviour()->getFractureCriterion())
            {
                maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
            }
            maxModeInNeighbourhood = ci->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = ci->getBehaviour()->getDamageModel()->getAngleShift() ;

            for(size_t i = 1 ; i < damagingSet.size() ; i++)
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                {
                    double nls = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                    maxscore = std::min(maxscore,nls) ;
                }
                maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;

            }

            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                    minscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
                    ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                    if(ci->getBehaviour()->getFractureCriterion())
                    {
                        double nls = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }
        }

    }

    //shut up the compiler
    return std::make_pair(0., 0.) ;
}

void FractureCriterion::step(ElementState &s)
{
    if(!mesh2d && !mesh3d )
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
    metAtStep = scoreAtState > 0 ;

}

    double FractureCriterion::getScoreAtState() const {
        return scoreAtState ;
    }

FractureCriterion::~FractureCriterion()
{
}


void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
    physicalCharacteristicRadius = r ;
    mesh2d = nullptr ;
    mesh3d = nullptr ;
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



