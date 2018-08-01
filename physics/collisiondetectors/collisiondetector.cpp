//
// C++ Implementation: collision detector
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "collisiondetector.h"
#include "../geometryBasedEffects/geometryBasedEffect.h"
#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../../mesher/delaunay.h"
#include "../../physics/viscoelasticity.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
#include "../../polynomial/vm_function_extra.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
namespace Amie
{


CollisionDetector::CollisionDetector() : FractureCriterion()
{
}



double CollisionDetector::getMaxScoreInNeighbourhood(ElementState & s)
{
    double maxScore = scoreAtState ;
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(mesh2d)
    {
        for(auto i =  mesh2d->begin(cacheID) ; i !=  mesh2d->end(cacheID) ; i++)
        {
            if(i->getBehaviour()->getCollisionDetection())
            {
                double renormScore = i->getBehaviour()->getCollisionDetection()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
        }
    }
    else
    {
        for(auto i =  mesh3d->begin(cacheID) ; i !=  mesh3d->end(cacheID) ; i++)
        {
            if(i->getBehaviour()->getCollisionDetection() )
            {
                double renormScore = i->getBehaviour()->getCollisionDetection()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
        }
    }
    return maxScore ;
}

std::pair<double, double> CollisionDetector::setChange( ElementState &s, double thresholdScore)
{
    stable = true ;

    if( !s.getParent()->getBehaviour()->getGeometryBasedEffect())
        return std::make_pair(0.,0.) ;
    if(!mesh2d && !mesh3d)
    {
        initialiseCache(s);
    }

    if(mesh2d)
    {
// 		thresholdScore = POINT_TOLERANCE ;
        // outside of the checkpoints, we only care about the order of the elements in
        // term of their score. At the checkpoint, we consider the elements which
        // have met their criterion

        if( checkpoint ) //new iteration
        {
            
//             std::cout << "ccheckpoint " << thresholdScore << "  " << getScoreAtState() <<  std::endl ;
            inset = false ;
            inIteration = false ;
            damagingSet.clear();
            proximitySet.clear() ;
// 	    initialScore = std::max(scoreAtState, scoreTolerance*scoreTolerance) ;

            std::vector<unsigned int> newSet ;
            std::set<unsigned int> newProximity ;

            double maxscore = std::max(0.,thresholdScore) ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            if(s.getParent()->getBehaviour()->getFractureCriterion() && s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > getScoreAtState())
                return std::make_pair(0.,0.) ;
            if(thresholdScore > 0 )
            {
                for(auto ci = mesh2d->begin(cacheID) ; ci != mesh2d->end(cacheID) ; ci++)
                {
                    if(!ci->getBehaviour()->getGeometryBasedEffect())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(thresholdScore-ci->getBehaviour()->getCollisionDetection()->getScoreAtState() <= 4.*getScoreTolerance()*initialScore &&
                            ci->getBehaviour()->getCollisionDetection()->met())
                    {
  

                        if(ci == s.getParent() && met())
                        {
                            inset = true ;
                        }
                        if(ci->getBehaviour()->getGeometryBasedEffect()->getDelta() > POINT_TOLERANCE)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
                            foundmaxscore = true ;
                        }
                        else
                            maxscore = std::max(ci->getBehaviour()->getCollisionDetection()->getScoreAtState(), maxscore) ;

                    }
                }
            }
            
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
            {
                static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection()->inIteration = true ;
            }

            return std::make_pair(scoreAtState - maxscore, thresholdScore - scoreAtState - POINT_TOLERANCE) ;
        }
        else if (inset)
        {
            checkpoint = false ;
            DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(damagingSet[0])) ;

            maxModeInNeighbourhood = s.getParent()->getBehaviour()->getGeometryBasedEffect()->getMode() ;
            maxAngleShiftInNeighbourhood = s.getParent()->getBehaviour()->getGeometryBasedEffect()->getAngleShift() ;


            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getCollisionDetection())
                    minscore = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
                    ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i])) ;
                    if(ci->getBehaviour()->getCollisionDetection() && !ci->getBehaviour()->getGeometryBasedEffect()->converged)
                    {
                        double nls = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }

            return std::make_pair(scoreAtState - minscore, thresholdScore - scoreAtState - POINT_TOLERANCE) ;
        }
    }
    else
    {
        // outside of the checkpoints, we only care about the order of the elements in
        // term of their score. At the checkpoint, we consider the elements which
        // have met their criterion

         if( checkpoint ) //new iteration
        {
            inset = false ;
            inIteration = false ;
            damagingSet.clear();
            proximitySet.clear() ;
            
// 	    initialScore = std::max(scoreAtState, scoreTolerance*scoreTolerance) ;

            std::vector<unsigned int> newSet ;
            std::set<unsigned int> newProximity ;

            double maxscore = std::max(0.,thresholdScore) ;
            bool foundmaxscore = false ;
            minDeltaInNeighbourhood = 1 ;
            maxModeInNeighbourhood = -1 ;
            maxAngleShiftInNeighbourhood = 0 ;
            
            if(s.getParent()->getBehaviour()->getFractureCriterion() && s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > getScoreAtState())
                return std::make_pair(0.,0.) ;
            if(thresholdScore > 0 )
            {
                for(auto ci = mesh3d->begin(cacheID) ; ci != mesh3d->end(cacheID) ; ci++)
                {
                    if(!ci->getBehaviour()->getGeometryBasedEffect())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(thresholdScore-ci->getBehaviour()->getCollisionDetection()->getScoreAtState() <= 4.*getScoreTolerance()*initialScore &&
                            ci->getBehaviour()->getCollisionDetection()->met())
                    {
  

                        if(ci == s.getParent() && met())
                        {
                            inset = true ;
                        }
                        if(ci->getBehaviour()->getGeometryBasedEffect()->getDelta() > POINT_TOLERANCE)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getGeometryBasedEffect()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
                            foundmaxscore = true ;
                        }
                        else
                            maxscore = std::max(ci->getBehaviour()->getCollisionDetection()->getScoreAtState(), maxscore) ;
                    }
                }
            }
            
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
            {
                static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection()->inIteration = true ;
            }

            return std::make_pair(scoreAtState - maxscore, thresholdScore - scoreAtState) ;
        }
        else if (inset)
        {
            checkpoint = false ;
            DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;

            maxModeInNeighbourhood = s.getParent()->getBehaviour()->getGeometryBasedEffect()->getMode() ;
            maxAngleShiftInNeighbourhood = s.getParent()->getBehaviour()->getGeometryBasedEffect()->getAngleShift() ;


            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getCollisionDetection())
                    minscore = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
                    ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i])) ;
                    if(ci->getBehaviour()->getCollisionDetection() && !ci->getBehaviour()->getGeometryBasedEffect()->converged)
                    {
                        double nls = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }

            return std::make_pair(scoreAtState - minscore, thresholdScore - scoreAtState - scoreTolerance) ;
        }

    }

    //shut up the compiler
    return std::make_pair(0., 0.) ;
}

void CollisionDetector::step(ElementState &s)
{
    if(!mesh2d && !mesh3d )
        initialiseCache(s) ;

    if(mesh2d)
        dynamic_cast<DelaunayTriangle *>(s.getParent())->getSubTriangulatedGaussPoints() ;
    if(mesh3d)
        dynamic_cast<DelaunayTetrahedron *>(s.getParent())->getSubTriangulatedGaussPoints() ;

    if(s.getParent()->getBehaviour()->fractured())
    {
        scoreAtState = -1 ;
        scoreAtTimeStepEnd = -1 ;
        metAtStep = false ;
        return ;
    }
//     if(checkpoint || inIteration)
//     {
        scoreAtState = grade(s) ;
//     }

    metAtStep = scoreAtState > 0 ;    

}


CollisionDetector::~CollisionDetector()
{
}



}

