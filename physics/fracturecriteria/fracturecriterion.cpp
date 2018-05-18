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
#include "../collisiondetectors/collisiondetector.h"
#include "../contactmodels/contactmodel.h"
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

double getSmoothingKernelSize( SmoothingFunctionType type )
{
    return 2.5 ;
    switch(type)
    {
    case QUARTIC_COMPACT:
        return 2.5 ;
    case GAUSSIAN_NONCOMPACT:
        return 2.5 ;
    case LINEAR_COMPACT:
        return 2.5 ;
    }
    return 1. ;
}

Function getSmoothingKernelFunction( SmoothingFunctionType type, Function & rrn )
{
    switch(type)
    {
    case QUARTIC_COMPACT:
    {
        rrn/=2. ;
        return (rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
    }
    case GAUSSIAN_NONCOMPACT:
    {
        return f_exp(rrn*-1) ;
    }
    case LINEAR_COMPACT:
    {
        return f_sqrt(1.-rrn)*f_positivity(1.-rrn) ;
    }
    }
    return Function("1") ;
}

FractureCriterion::FractureCriterion() :
    restrictionSource(nullptr),
    initialScore(1),
    physicalCharacteristicRadius(.008),
    scoreAtState(0),
    deltaScoreAtState(0),
    metAtStep(false),
    stable(true),
    minDeltaInNeighbourhood(1),
    maxModeInNeighbourhood(-1),
    maxScoreInNeighbourhood(0),
    maxAngleShiftInNeighbourhood(0),
    scoreTolerance(2e-2),
    checkpoint(true),
    smoothingType(QUARTIC_COMPACT),
    cachedInfluenceRatio(1),
    cacheID(-1),
    cachecoreID(-1),
    needRestrictionUpdate(false),
    inIteration(false),
    inset(false),
    mesh2d(nullptr), mesh3d(nullptr)
{
    overlap = getSmoothingKernelSize( smoothingType ) ;
}

void FractureCriterion::setSmoothingFunctionType( SmoothingFunctionType type, bool over)
{
    smoothingType = type ;
    if(over)
        overlap = getSmoothingKernelSize( smoothingType ) ;
}

void FractureCriterion::copyEssentialParameters( const FractureCriterion * frac )
{
    setScoreTolerance( frac->getScoreTolerance() ) ;
    setMaterialCharacteristicRadius( frac->getMaterialCharacteristicRadius() ) ;
    setSmoothingFunctionType( frac->getSmoothingFunctionType(), false ) ;
    setSmoothingFunctionOverlap( frac->getSmoothingFunctionOverlap() ) ;
}

Vector FractureCriterion::getSmoothedField(FieldType f0,  ElementState &s ,double t)
{
//   std::cout << "a" << std::endl ;
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

//   std::cout << "b" << std::endl ;
    if(mesh2d)
        return mesh2d->getSmoothedField(f0, cachecoreID, s.getParent(), t, restriction) ;
    return mesh3d->getSmoothedField(f0, cachecoreID, s.getParent(), t, restriction) ;

}

std::pair<Vector, Vector> FractureCriterion::getSmoothedFields( FieldType f0, FieldType f1, ElementState &s ,double t)
{
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(mesh2d)
        return mesh2d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(),  t, restriction) ;
    return mesh3d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), t, restriction) ;

}

void FractureCriterion::updateRestriction(ElementState &s)
{
    if(!restrictionSource || !needRestrictionUpdate)
        return ;

    needRestrictionUpdate = false ;
    updateCache(s) ;

    Function xtransform = s.getParent()->getXTransform() ;
    Function ytransform = s.getParent()->getYTransform() ;
    Function ztransform = Function("0") ;
    Function ttransform = Function("0") ;
    if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
        ztransform = s.getParent()->getZTransform() ;
    if(s.getParent()->timePlanes() > 1)
        ttransform = s.getParent()->getTTransform() ;
    else
        ttransform = Function("0") ;

    VirtualMachine vm ;
    restriction.clear();

    for(size_t i= 0 ; i < s.getParent()->getGaussPoints().gaussPoints.size() ; i++)
    {
        Point p = s.getParent()->getGaussPoints().gaussPoints[i].first ;
        Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()),
                           vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()),
                           vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()),
                           vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

        restriction.push_back(restrictionSource->in(test));
    }
// }
}

void FractureCriterion::setRestriction(const Geometry * g,ElementState &s)
{
    restrictionSource = g ;
    restriction.clear();
    needRestrictionUpdate = true ;
}

void FractureCriterion::initialiseCache( ElementState & s)
{
    physicalCharacteristicRadius = std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*1.1) ;
    if(s.getMesh2D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        //it is desirable to always have a meaningful non-local law
//         physicalCharacteristicRadius = std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*1.5) ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth = getSmoothingKernelFunction( smoothingType, rrn ) ;

        Circle epsilonAll( std::max(physicalCharacteristicRadius*6., s.getParent()->getRadius()*4. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Circle epsilonReduced(physicalCharacteristicRadius*overlap+s.getParent()->getRadius()*1.5,s.getParent()->getCenter()) ;
        mesh2d = s.getMesh2D() ;
        cachecoreID = mesh2d->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        if(s.getParent()->timePlanes() > 1)
            cacheID = cachecoreID ;
        else
            cacheID = mesh2d->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;
    }
    if(s.getMesh3D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function z = Function("z")-s.getParent()->getCenter().getZ() ;
        Function rr = x*x+y*y+z*z ;
//         physicalCharacteristicRadius = std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*1.5) ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth = getSmoothingKernelFunction( smoothingType, rrn ) ;
        Sphere epsilonAll( std::max(physicalCharacteristicRadius*6., s.getParent()->getRadius()*4. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Sphere epsilonReduced(physicalCharacteristicRadius*overlap+s.getParent()->getRadius()*1.5,s.getParent()->getCenter()) ;
        mesh3d = s.getMesh3D() ;
        cachecoreID = mesh3d->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        if(s.getParent()->timePlanes() > 1)
            cacheID = cachecoreID ;
        else
            cacheID = mesh3d->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;

    }
    
//     std::cout << mesh2d->begin(cacheID)->size() << std::endl ;
}

void FractureCriterion::updateCache( ElementState & s)
{
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(s.getMesh2D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  getSmoothingKernelFunction( smoothingType, rrn ) ;

        mesh2d = s.getMesh2D() ;
        mesh2d->updateCache(cachecoreID, s.getParent()->getBehaviour()->getSource(), smooth) ;
        if(s.getParent()->timePlanes() <= 1)
            mesh2d->updateCache(cacheID) ;
    }
    if(s.getMesh3D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function z = Function("z")-s.getParent()->getCenter().getZ() ;
        Function rr = x*x+y*y+z*z ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  getSmoothingKernelFunction( smoothingType, rrn ) ;

        mesh3d = s.getMesh3D() ;
        mesh3d->updateCache(cachecoreID, s.getParent()->getBehaviour()->getSource(), smooth) ;
        if(s.getParent()->timePlanes() > 1)
            mesh3d->updateCache(cacheID) ;

    }
}

double FractureCriterion::getMaxScoreInNeighbourhood(ElementState & s)
{
    double maxScore = scoreAtState ;
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(mesh2d)
    {
        for(auto i =  mesh2d->begin(cacheID) ; i !=  mesh2d->end(cacheID) ; i++)
        {
            if(i->getBehaviour()->getFractureCriterion() && !i->getBehaviour()->fractured())
            {
                double renormScore = i->getBehaviour()->getFractureCriterion()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
            if(i->getBehaviour()->getCollisionDetection())
            {
                double renormScore = i->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
                maxScore = std::max(maxScore, renormScore) ;
            }
        }
    }
    else
    {
        for(auto i =  mesh3d->begin(cacheID) ; i !=  mesh3d->end(cacheID) ; i++)
        {
            if(i->getBehaviour()->getFractureCriterion() && !i->getBehaviour()->fractured())
            {
                double renormScore = i->getBehaviour()->getFractureCriterion()->scoreAtState ;
                maxScore = std::max(maxScore, renormScore) ;
            }
            if(i->getBehaviour()->getCollisionDetection())
            {
                double renormScore = i->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
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
    }

    if(mesh2d)
    {
// 		thresholdScore = POINT_TOLERANCE ;
        // outside of the checkpoints, we only care about the order of the elements in
        // term of their score. At the checkpoint, we consider the elements which
        // have met their criterion

        if( checkpoint ) //new iteration
        {
//             std::cout << "checkpoint " << thresholdScore  << "  "<< scoreAtState << "  "<<grade(s) << std::endl ;
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
            if(thresholdScore > 0 )
            {
                
                for(auto ci = mesh2d->begin(cacheID) ; ci != mesh2d->end(cacheID) ; ci++)
                {
                    
                    if(!ci->getBehaviour()->getDamageModel())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(thresholdScore-ci->getBehaviour()->getFractureCriterion()->scoreAtState <= 4.*scoreTolerance*initialScore &&
                            ci->getBehaviour()->getFractureCriterion()->met())
                    {
//                         std::cout << "met " << scoreAtState <<std::endl ;

                        if(ci == s.getParent() && met())
                        {
                            inset = true ;
                            if(ci->getBehaviour()->getCollisionDetection())
                                ci->getBehaviour()->getCollisionDetection()->inset = true ;
                        }
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getDamageModel()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                            foundmaxscore = true ;
                        }
                        else
                            maxscore = std::max(ci->getBehaviour()->getFractureCriterion()->scoreAtState, maxscore) ;
                        
                         if(ci->getBehaviour()->getCollisionDetection())
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
                static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i]))->getBehaviour()->getFractureCriterion()->inIteration = true ;
                if(static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection())
                    static_cast<DelaunayTriangle *>( mesh2d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection()->inIteration = true ;
            }


            return std::make_pair(scoreAtState - maxscore, thresholdScore - scoreAtState) ;
        }
        else if (inset)
        {
            checkpoint = false ;
            DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(damagingSet[0])) ;

            maxModeInNeighbourhood = s.getParent()->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = s.getParent()->getBehaviour()->getDamageModel()->getAngleShift() ;


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
                    if(ci->getBehaviour()->getContactModel() && !ci->getBehaviour()->getContactModel()->converged)
                    {
                        double nls = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }

            return std::make_pair(scoreAtState - minscore , thresholdScore - scoreAtState - scoreTolerance) ;
        }
    }
    else
    {
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
            if(thresholdScore > 0 )
            {
                for(auto ci = mesh3d->begin(cacheID) ; ci != mesh3d->end(cacheID) ; ci++)
                {
                    if(!ci->getBehaviour()->getDamageModel())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(thresholdScore-ci->getBehaviour()->getFractureCriterion()->scoreAtState <= 4.*scoreTolerance*initialScore &&
                            ci->getBehaviour()->getFractureCriterion()->met())
                    {


                        if(ci == s.getParent() && met())
                        {
                            inset = true ;
                            if(ci->getBehaviour()->getCollisionDetection())
                                ci->getBehaviour()->getCollisionDetection()->inset = true ;
                        }
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE)
                            minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, ci->getBehaviour()->getDamageModel()->getDelta()) ;

                        maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
                        maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;
                        newSet.push_back(ci->index);
                    }
                    else
                    {
                        newProximity.insert(ci->index) ;
                        if(!foundmaxscore)
                        {
                            maxscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;
                            foundmaxscore = true ;
                        }
                        else
                            maxscore = std::max(ci->getBehaviour()->getFractureCriterion()->scoreAtState, maxscore) ;
                        
                         if(ci->getBehaviour()->getCollisionDetection())
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
                static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i]))->getBehaviour()->getFractureCriterion()->inIteration = true ;
                if(static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection())
                    static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i]))->getBehaviour()->getCollisionDetection()->inIteration = true ;
            }


            return std::make_pair(scoreAtState - maxscore, thresholdScore - scoreAtState) ;
        }
        else if (inset)
        {
            checkpoint = false ;
            DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(damagingSet[0])) ;

            maxModeInNeighbourhood = s.getParent()->getBehaviour()->getDamageModel()->getMode() ;
            maxAngleShiftInNeighbourhood = s.getParent()->getBehaviour()->getDamageModel()->getAngleShift() ;


            double minscore = 0 ;
            if(!proximitySet.empty())
            {
                ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[0])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                    minscore = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

                for(size_t i = 1 ; i < proximitySet.size() ; i++)
                {
                    ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(proximitySet[i])) ;
                    if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->getDamageModel()->converged)
                    {
                        double nls = ci->getBehaviour()->getFractureCriterion()->scoreAtState ;

                        minscore = std::max(nls, minscore) ;
                    }
                    if(ci->getBehaviour()->getContactModel() && !ci->getBehaviour()->getContactModel()->converged)
                    {
                        double nls = ci->getBehaviour()->getCollisionDetection()->getScoreAtState() ;

                        minscore = std::max(nls, minscore) ;
                    }
                }
            }

            return std::make_pair(scoreAtState - minscore , thresholdScore - scoreAtState - scoreTolerance) ;
        }

    }

    //shut up the compiler
    return std::make_pair(0., 0.) ;
}

void FractureCriterion::step(ElementState &s)
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
    if(checkpoint || inIteration)
    {
        scoreAtState = grade(s) ;
//         std::cout << "score "<< scoreAtState << "  "<< getScoreAtState() <<std::endl ;
    }

//     
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

bool FractureCriterion::met(double threshold) const
{
    return scoreAtState > threshold ;
}

}

