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
    restrictionSource(nullptr),
    initialScore(1),
    physicalCharacteristicRadius(.008),
    scoreAtState(-1),
    deltaScoreAtState(0),
    criterionDamageDifferential(0),
    mirroring(mirroring),
    delta_x(delta_x), delta_y(delta_y), delta_z(delta_z),
    metAtStep(false),
    stable(true),
    minDeltaInNeighbourhood(1),
    maxModeInNeighbourhood(-1),
    maxScoreInNeighbourhood(0),
    maxAngleShiftInNeighbourhood(0),
    scoreTolerance(1e-5),
    checkpoint(true),
    inset(false),
    smoothingType(QUARTIC_COMPACT),
    cachedInfluenceRatio(1),
    cacheID(-1),
    cachecoreID(-1),
    inIteration(false),
    mesh2d(nullptr), mesh3d(nullptr)
{
}

void FractureCriterion::copyEssentialParameters( const FractureCriterion * frac ) 
{
    setScoreTolerance( frac->getScoreTolerance() ) ;
    setMaterialCharacteristicRadius( frac->getMaterialCharacteristicRadius() ) ;
    setSmoothingFunctionType( frac->getSmoothingFunctionType() ) ;
}

Vector FractureCriterion::getSmoothedField(FieldType f0,  ElementState &s ,double t)
{
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(mesh2d)
        return mesh2d->getSmoothedField(f0, cachecoreID, s.getParent(), -1., t, restriction) ;
    return mesh3d->getSmoothedField(f0, cachecoreID, s.getParent(), -1., t, restriction) ;

}

std::pair<Vector, Vector> FractureCriterion::getSmoothedFields( FieldType f0, FieldType f1, ElementState &s ,double t)
{
    if(!mesh2d && !mesh3d)
        initialiseCache(s) ;

    if(mesh2d)
        return mesh2d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), -1., t, restriction) ;
    return mesh3d->getSmoothedFields(f0, f1, cachecoreID, s.getParent(), -1., t, restriction) ;

}

void FractureCriterion::updateRestriction(ElementState &s)
{
    if(!restrictionSource)
        return ;
    
    if(mesh2d)
    {
        mesh2d->deleteCache(cachecoreID);
        mesh2d->deleteCache(cacheID);
    }
    if(mesh3d)
    {
        mesh3d->deleteCache(cachecoreID);
        mesh3d->deleteCache(cacheID);
    }
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
}

void FractureCriterion::setRestriction(const Geometry * g,ElementState &s)
{
    restrictionSource = g ;
    restriction.clear();
    if(mesh2d)
    {
        mesh2d->deleteCache(cachecoreID);
        mesh2d->deleteCache(cacheID);
    }
    if(mesh3d)
    {
        mesh3d->deleteCache(cachecoreID);
        mesh3d->deleteCache(cacheID);
    }
    
    if(!g)
        return ;
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
    
    for(size_t i= 0 ; i < s.getParent()->getGaussPoints().gaussPoints.size() ; i++)
    {
        Point p = s.getParent()->getGaussPoints().gaussPoints[i].first ;
        Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                           vm.eval(ytransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                           vm.eval(ztransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                           vm.eval(ttransform, p.getX(), p.getY(), p.getZ(), p.getT())) ;

        restriction.push_back(g->in(test));
    }
    initialiseCache(s) ;
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
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;

        double overlap = (smoothingType == QUARTIC_COMPACT)?6.:8. ;
        Sphere epsilonAll( std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*2. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Sphere epsilonReduced(physicalCharacteristicRadius*1.1+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        mesh3d = s.getMesh3D() ;
        cachecoreID = mesh3d->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        if(s.getParent()->timePlanes() > 1)
            cacheID = cachecoreID ;
        else
            cacheID = mesh3d->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;

    }
}
void FractureCriterion::updateCache( ElementState & s)
{
    if(s.getMesh2D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;

        mesh2d = s.getMesh2D() ;
        mesh2d->updateCache(cachecoreID, smooth) ;
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
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;

        mesh3d = s.getMesh3D() ;
        mesh3d->updateCache(cachecoreID, smooth) ;
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
                for(auto ci = mesh2d->begin(cacheID) ; ci != mesh2d->end(cacheID) ; ci++)
                {
                    if(!ci->getBehaviour()->getDamageModel())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(std::abs(ci->getBehaviour()->getFractureCriterion()->scoreAtState-thresholdScore) <= scoreTolerance*initialScore*.25 &&
                            ci->getBehaviour()->getFractureCriterion()->met())
                    {


                        if(ci == s.getParent() && met())
                                inset = true ;
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE)
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
                for(auto ci = mesh3d->begin(cacheID) ; ci != mesh3d->end() ; ci++)
                {
//                     DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(mesh3d->getCache(cacheID)[i])) ;
                    if(!ci->getBehaviour()->getDamageModel())
                        continue ;

                    if(ci->getBehaviour()->fractured())
                        continue ;

                    if(std::abs(ci->getBehaviour()->getFractureCriterion()->scoreAtState-thresholdScore) <= scoreTolerance*initialScore*.25 &&
                            ci->getBehaviour()->getFractureCriterion()->met())
                    {
                        if(ci == s.getParent() && met())
                            inset = true ;
                        if(ci->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE)
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
    
    if(mesh2d)
        dynamic_cast<DelaunayTriangle *>(s.getParent())->getSubTriangulatedGaussPoints() ;
    if(mesh3d)
        dynamic_cast<DelaunayTetrahedron *>(s.getParent())->getSubTriangulatedGaussPoints() ;
    updateRestriction(s) ;

    if(s.getParent()->getBehaviour()->fractured())
    {
        scoreAtState = -1 ;
        metAtStep = false ;
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



