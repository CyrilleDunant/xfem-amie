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
    deltaEnergyAtState(0),
    energyDamageDifferential(0),
    criterionDamageDifferential(0),
    energyIndexed(false),
    noEnergyUpdate(true),
    mesh2d(nullptr), mesh3d(nullptr),
    stable(true), checkpoint(true), inset(false),inIteration(false),
    scoreTolerance(1e-3),
    initialScore(1),
    cachedInfluenceRatio(1),
    currentAngle(0),
    minDeltaInNeighbourhood(1),
    maxScoreInNeighbourhood(0),
    maxModeInNeighbourhood(-1),
    maxAngleShiftInNeighbourhood(0),
    smoothingType(GAUSSIAN_NONCOMPACT),
    cacheID(-1), cachecoreID(-1)
{
}

double FractureCriterion::smoothedCrackAngle( ElementState &s) const
{
    double angle = 0 ;
    if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        double area = s.getParent()->area() ;
        double fact =  0 ; //area*s.getParent()->getBehaviour()->getDamageModel()->getState().max() ;

        for( size_t i = 0 ; i < cache.size() ; i++ )
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i]) ) ;
            double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
            if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()
                    || !ci->getBehaviour()->getFractureCriterion()
                    || ci->getBehaviour()->fractured()
                    || ci->getBehaviour()->type == VOID_BEHAVIOUR
                    || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource()
                    || dc >  4.*physicalCharacteristicRadius * physicalCharacteristicRadius)
            {
                continue ;
            }

            area = ci->area() ;
            //this is to eliminate scaling effects ;
            double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;

            double d = ci->getBehaviour()->getDamageModel()->getState().max() ;
            angle += atan2(ci->getCenter().getY()-s.getParent()->getCenter().getY(),ci->getCenter().getX()-s.getParent()->getCenter().getX())*area*d ;
            fact += area*d ;

            if( mirroring == MIRROR_X && std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
            {
                angle += atan2(ci->getCenter().getY()-s.getParent()->getCenter().getY(),ci->getCenter().getX()-s.getParent()->getCenter().getX())*area*d ;
                fact += area*d ;
            }

            if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
            {
                angle += atan2(ci->getCenter().getY()-s.getParent()->getCenter().getY(),ci->getCenter().getX()-s.getParent()->getCenter().getX())*area*d ;
                fact += area*d ;
            }

            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                angle += atan2(ci->getCenter().getY()-s.getParent()->getCenter().getY(),ci->getCenter().getX()-s.getParent()->getCenter().getX())*area*d ;
                fact += area*d ;
            }

            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                angle += atan2(ci->getCenter().getY()-s.getParent()->getCenter().getY(),ci->getCenter().getX()-s.getParent()->getCenter().getX())*area*d ;
                fact += area*d ;
            }
        }

        if(fact > POINT_TOLERANCE_2D)
            angle /= (fact+ s.getParent()->area()*s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
        else
            return 0. ;
    }
    return angle ;
}

double FractureCriterion::smoothedPrincipalStressAngle( ElementState &s, StressCalculationMethod m )
{
    double angle = 0 ;
    Vector ss = getSmoothedField(s, REAL_STRESS_FIELD) ;
    if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        return 0.5 * atan2( 0.5*ss[2], ss[0] - ss[1] ) ;
    }
    else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
    {
        return  0.5 * atan2(0.5*ss[3] , ss[0] - ss[1] ) ;
    }

    return angle ;
}
/*
Vector FractureCriterion::smoothedPrincipalStrain(ElementState &s)
{
    return smoothedPrincipalStressAndStrain(s).second ;
}

std::pair< Vector, Vector > FractureCriterion::smoothedPrincipalStressAndStrain(ElementState& s, StressCalculationMethod m , double t )
{
        std::pair< Vector, Vector > stressStrain = smoothedStressAndStrain(s,m, t) ;
        return std::make_pair(toPrincipal(stressStrain.first), toPrincipal(stressStrain.second)) ;
}*/

double FractureCriterion::smoothedScore(ElementState& s)
{

    double score = 0;
    if(factors.size()==0)
        initialiseFactors(s) ;
    double total = 0 ;
    auto fiterator = &factors[0] ;
    if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        score =scoreAtState*(*fiterator) ;
        total += (*fiterator);
        fiterator++ ;
        for( size_t i = 1 ; i < physicalcache.size() ; i++ )
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(physicalcache[i]) ) ;
            if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured())
            {
                score += ci->getBehaviour()->getFractureCriterion()->scoreAtState*(*fiterator) ;
                total += (*fiterator);
                fiterator++ ;
            }
        }
        score /= total ;

        return score ;
    }
    else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
    {
        score =scoreAtState*(*fiterator) ;
        total += (*fiterator);
        fiterator++ ;
        for( size_t i = 0 ; i < cache.size() ; i++ )
        {
            DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i]) ) ;
            if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured() && *fiterator >= POINT_TOLERANCE_2D)
            {
                score += ci->getBehaviour()->getFractureCriterion()->scoreAtState*(*fiterator) ;
                total += (*fiterator);
                fiterator++ ;
            }
        }
        score /= total ;
        return score ;
    }

    return 0. ;
}

void FractureCriterion::initialiseFactors(const ElementState & s)
{

    bool compact = (smoothingType == QUARTIC_COMPACT) ;

    if(cache.empty())
    {
        physicalcache.resize(0);
        initialiseCache(s);
    }
    if(!factors.size() == 0)
        factors.resize(0) ;

    std::vector<double> tmpfactors ;
    std::vector<unsigned int> tmpphysicalcache ;
    VirtualMachine vm ;
    if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        std::valarray< std::pair<Point, double> > fin(7) ;
        fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456), 0.062969590272413) ;
        fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456), 0.062969590272413) ;
        fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087), 0.062969590272413) ;
        fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770), 0.066197076394253) ;
        fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115), 0.066197076394253) ;
        fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115), 0.066197076394253) ;
        fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.1125) ;

        std::valarray< std::pair<Point, double> > fintmp = fin ;
        double j = s.getParent()->area()*2. ;//1./det(J) ;
        for(size_t i = 0 ; i < fintmp.size() ; i++)
        {
            fintmp[i].second *= j;
        }
        GaussPointArray gp(fintmp, QUINTIC) ;

        Function x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
        Function y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;

// 		for(double i = -2.*physicalCharacteristicRadius ; i < 2.*physicalCharacteristicRadius ; i += 0.01*physicalCharacteristicRadius)
// 			std::cout << vm.eval(smooth, i, s.getParent()->getCenter().getY()) << "  " <<   vm.eval(x, i, s.getParent()->getCenter().getY()) << "  " <<   vm.eval(y, i, s.getParent()->getCenter().getY())<< "  " <<   vm.eval(rrn, i, s.getParent()->getCenter().getY())<< std::endl ;
//
// 		exit(0) ;
        double weight = vm.ieval(smooth, gp) ;
        double fact = weight ;


        tmpfactors.push_back(weight);
        tmpphysicalcache.push_back(dynamic_cast<DelaunayTriangle *>(s.getParent())->index);

        if( mirroring == MIRROR_X && std::abs( s.getParent()->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
        {
            x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            rr = x*x+y*y ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, gp) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_Y &&  std::abs( s.getParent()->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
        {
            x = s.getParent()->getXTransform() ;
            y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getX() ;
            rr = x*x+y*y ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, gp) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            rr = x*x+y*y ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, gp) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getY() ;
            rr = x*x+y*y ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, gp) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }

        for( size_t i = 0 ; i < cache.size() ; i++ )
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i]) ) ;
            if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()
                    || ci->getBehaviour()->type == VOID_BEHAVIOUR
                    || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() )
            {
                continue ;
            }

            //this is to eliminate scaling effects ;
// 			double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ; ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;
// 			double d = dist(ci->getCenter(), s.getParent()->getCenter()) ;
// 			if(d > farthest)
// 				farthest = d ;
            x = ci->getXTransform()-s.getParent()->getCenter().getX() ;
            y = ci->getYTransform()-s.getParent()->getCenter().getY() ;
            rr = x*x+y*y ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius ) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;

            fintmp = fin ;
            double j = ci->area()*2. ;//1./det(J) ;
            for(size_t i = 0 ; i < fintmp.size() ; i++)
            {
                fintmp[i].second *= j;
            }
            GaussPointArray gp(fintmp, QUINTIC) ;

            weight = vm.ieval(smooth, gp) ;
            if(weight < POINT_TOLERANCE_2D)
                continue ;

            tmpphysicalcache.push_back(ci->index) ;
            tmpfactors.push_back(weight);

            fact += weight ;

            if( mirroring == MIRROR_X && std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
            {
                x = ci->getXTransform()*-1 -std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX();
                y = ci->getYTransform()-s.getParent()->getCenter().getY() ;
                rr = x*x+y*y ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth =!compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, gp) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }

            if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
            {
                x = ci->getXTransform()-s.getParent()->getCenter().getX() ;
                y = ci->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getY() ;
                rr = x*x+y*y ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, gp) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }

            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()*-1 -std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX();
                y = ci->getYTransform()-s.getParent()->getCenter().getY() ;
                rr = x*x+y*y ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth =!compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, gp) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }

            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()-s.getParent()->getCenter().getX() ;
                y = ci->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getY() ;
                rr = x*x+y*y ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, gp) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }

        }
        factors.resize(tmpfactors.size());
        physicalcache.resize(tmpphysicalcache.size());
        std::copy(tmpfactors.begin(), tmpfactors.end(), &factors[0]) ;
        std::copy(tmpphysicalcache.begin(), tmpphysicalcache.end(), &physicalcache[0]) ;

        return ;
    }
    else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
    {

        Function x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
        Function y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
        Function z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
        Function rr = x*x+y*y+z*z ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
        double weight = vm.ieval(smooth, s.getParent()) ;
        double fact = weight ;

        tmpfactors.push_back(weight);
        tmpphysicalcache.push_back(dynamic_cast<DelaunayTriangle *>(s.getParent())->index);

        double selfarea = s.getParent()->area() ;
        double farthest = 0 ;
        if( mirroring == MIRROR_X && std::abs( s.getParent()->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
        {
            x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_Y &&  std::abs( s.getParent()->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
        {
            x = s.getParent()->getXTransform() ;
            y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getX() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_Z &&  std::abs( s.getParent()->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
        {
            x = s.getParent()->getXTransform() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getX() ;
            z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().getZ()  - delta_z )-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            s.getParent()->setOrder(CUBIC) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =!compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XZ &&  std::abs( s.getParent()->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().getX()  - delta_x )-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_XZ &&  std::abs( s.getParent()->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().getZ()  - delta_z )-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_YZ &&  std::abs( s.getParent()->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().getY()  - delta_y )-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        if( mirroring == MIRROR_YZ &&  std::abs( s.getParent()->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
        {
            x = s.getParent()->getXTransform()-s.getParent()->getCenter().getX() ;
            y = s.getParent()->getYTransform()-s.getParent()->getCenter().getY() ;
            z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().getZ()  - delta_z )-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, s.getParent()) ;
            tmpfactors.back() += weight ;
            fact += weight ;
        }
        for( size_t i = 0 ; i < cache.size() ; i++ )
        {
            DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i]) ) ;
            if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()
                    || ci->getBehaviour()->type == VOID_BEHAVIOUR
                    || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() )
            {
                continue ;
            }
            double mindist = physicalCharacteristicRadius ;

            //this is to eliminate scaling effects ;
            double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ; ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;
            double d = dist(ci->getCenter(), s.getParent()->getCenter()) ;
            if(d > farthest)
                farthest = d ;
            x = ci->getXTransform()-s.getParent()->getCenter().getX() ;
            y = ci->getYTransform()-s.getParent()->getCenter().getY() ;
            z = ci->getZTransform()-s.getParent()->getCenter().getZ() ;
            rr = x*x+y*y+z*z ;
            rrn = rr/(mindist * mindist) ;
            smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
            weight = vm.ieval(smooth, ci) ;

            if(weight*factor < POINT_TOLERANCE_2D)
                continue ;

            tmpphysicalcache.push_back(ci->index) ;
            tmpfactors.push_back(weight*factor);

            fact += weight*factor ;

            if( mirroring == MIRROR_X && std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
            {
                x = ci->getXTransform()*-1-std::abs( ci->getCenter().getX()  - delta_x )-ci->getCenter().getX() ;
                y = ci->getYTransform()-ci->getCenter().getY() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
            {
                x = ci->getXTransform() ;
                y = ci->getYTransform()*-1-std::abs( ci->getCenter().getY()  - delta_y )-ci->getCenter().getX() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
            {
                x = ci->getXTransform() ;
                y = ci->getYTransform()-ci->getCenter().getX() ;
                z = ci->getZTransform()*-1-std::abs( ci->getCenter().getZ()  - delta_z )-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()*-1-std::abs( ci->getCenter().getX()  - delta_x )-ci->getCenter().getX() ;
                y = ci->getYTransform()-ci->getCenter().getY() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()-ci->getCenter().getX() ;
                y = ci->getYTransform()*-1-std::abs( ci->getCenter().getY()  - delta_y )-ci->getCenter().getY() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().getX()  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()*-1-std::abs( ci->getCenter().getX()  - delta_x )-ci->getCenter().getX() ;
                y = ci->getYTransform()-ci->getCenter().getY() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()-ci->getCenter().getX() ;
                y = ci->getYTransform()-ci->getCenter().getY() ;
                z = ci->getZTransform()*-1-std::abs( ci->getCenter().getZ()  - delta_z )-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().getY()  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()-ci->getCenter().getX() ;
                y = ci->getYTransform()*-1-std::abs( ci->getCenter().getY()  - delta_y )-ci->getCenter().getY() ;
                z = ci->getZTransform()-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth =  !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }
            if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().getZ()  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
            {
                x = ci->getXTransform()-ci->getCenter().getX() ;
                y = ci->getYTransform()-ci->getCenter().getY() ;
                z = ci->getZTransform()*-1-std::abs( ci->getCenter().getZ()  - delta_z )-ci->getCenter().getZ() ;
                rr = x*x+y*y+z*z ;
                rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
                smooth = !compact?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
                weight = vm.ieval(smooth, ci) ;
                tmpfactors.back() += weight ;
                fact += weight ;
            }

        }
        factors.resize(tmpfactors.size());
        physicalcache.resize(tmpphysicalcache.size());
        std::copy(tmpfactors.begin(), tmpfactors.end(), &factors[0]) ;
        std::copy(tmpphysicalcache.begin(), tmpphysicalcache.end(), &physicalcache[0]) ;
        return ;
    }

}

// Vector FractureCriterion::smoothedPrincipalStress( ElementState &s, StressCalculationMethod m)
// {
//     return smoothedPrincipalStressAndStrain(s, m).first ;
// }
// 
// 
// std::pair<Vector, Vector> FractureCriterion::smoothedStressAndStrain( ElementState &s , StressCalculationMethod m, double t)
// {
//     size_t vlength = 3 ;
//     if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
//         vlength = 6 ;
//     VirtualMachine vm ;
//     Vector str(0., vlength) ;
//     Vector stra(0., vlength) ;
//     Vector tmpstr(0.,vlength) ;
//     Vector tmpstra(0.,vlength) ;
// 
//     if(factors.size() == 0)
//         initialiseFactors(s) ;
//     double sumFactors = 0 ;
//     double stra0 = stra[0] ;
//     double stra1 = stra[1] ;
//     double stra2 = stra[2] ;
//     double str0 = str[0] ;
//     double str1 = str[1] ;
//     double str2 = str[2] ;
//     double stra3 = 0 ;
//     double stra4 = 0 ;
//     double stra5 = 0 ;
//     double str3 = 0 ;
//     double str4 = 0 ;
//     double str5 = 0 ;
// 
//     if(vlength == 6)
//     {
//         stra3 = stra[3] ;
//         stra4 = stra[4] ;
//         stra5 = stra[5] ;
//         str3 = str[3] ;
//         str4 = str[4] ;
//         str5 = str[5] ;
//     }
//     if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
//     {
//         if(s.getParent()->getOrder() < CONSTANT_TIME_LINEAR)
//         {
//             if(m == EFFECTIVE_STRESS)
//             {
//                 double iteratorValue = factors[0] ;
//                 s.getAverageField(STRAIN_FIELD,EFFECTIVE_STRESS_FIELD, tmpstra,tmpstr, &vm, 0, t);
// 
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
//                 #pragma omp parallel
//                 {
// #ifdef HAVE_OPENMP
//                     std::valarray<Vector> tmpstrap(tmpstra, omp_get_num_threads()) ;
//                     std::valarray<Vector> tmpstrp(tmpstr, omp_get_num_threads()) ;
//                     std::valarray<VirtualMachine> vmp(vm,omp_get_num_threads()) ;
// #else
//                     std::valarray<Vector> tmpstrap(tmpstra, 1) ;
//                     std::valarray<Vector> tmpstrp(tmpstr, 1) ;
//                     std::valarray<VirtualMachine> vmp(vm,1) ;
// #endif
//                     #pragma omp for reduction(+:stra0,stra1,stra2,str0,str1,str2,sumFactors)
//                     for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                     {
// #ifdef HAVE_OPENMP
//                         size_t thread = omp_get_thread_num() ;
// #else
//                         size_t thread = 0 ;
// #endif
//                         DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
//                         if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured())
//                         {
//                             ci->getState().getAverageField(STRAIN_FIELD,EFFECTIVE_STRESS_FIELD, tmpstrap[thread],tmpstrp[thread], &vmp[thread], 0, t);
//                             stra0 += tmpstrap[thread][0]*factors[i] ;
//                             stra1 += tmpstrap[thread][1]*factors[i] ;
//                             stra2 += tmpstrap[thread][2]*factors[i] ;
//                             str0 += tmpstrp[thread][0]*factors[i] ;
//                             str1 += tmpstrp[thread][1]*factors[i] ;
//                             str2 += tmpstrp[thread][2]*factors[i] ;
//                             sumFactors += factors[i] ;
//                         }
//                     }
//                 }
//             }
//             else
//             {
//                 double iteratorValue = factors[0] ;
//                 s.getAverageField(STRAIN_FIELD,REAL_STRESS_FIELD, tmpstra,tmpstr, &vm, 0, t);
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
//                 #pragma omp parallel
//                 {
// #ifdef HAVE_OPENMP
//                     std::valarray<Vector> tmpstrap(tmpstra, omp_get_num_threads()) ;
//                     std::valarray<Vector> tmpstrp(tmpstr, omp_get_num_threads()) ;
//                     std::valarray<VirtualMachine> vmp(vm,omp_get_num_threads()) ;
// #else
//                     std::valarray<Vector> tmpstrap(tmpstra, 1) ;
//                     std::valarray<Vector> tmpstrp(tmpstr, 1) ;
//                     std::valarray<VirtualMachine> vmp(vm,1) ;
// #endif
//                     #pragma omp for reduction(+:stra0,stra1,stra2,str0,str1,str2,sumFactors)
//                     for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                     {
// #ifdef HAVE_OPENMP
//                         size_t thread = omp_get_thread_num() ;
// #else
//                         size_t thread = 0 ;
// #endif
//                         DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
//                         if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured() )
//                         {
//                             ci->getState().getAverageField(STRAIN_FIELD,REAL_STRESS_FIELD, tmpstrap[thread],tmpstrp[thread], &vmp[thread], 0, t);
// 
//                             stra0 += tmpstrap[thread][0]*factors[i] ;
//                             stra1 += tmpstrap[thread][1]*factors[i] ;
//                             stra2 += tmpstrap[thread][2]*factors[i] ;
//                             str0 += tmpstrp[thread][0]*factors[i] ;
//                             str1 += tmpstrp[thread][1]*factors[i] ;
//                             str2 += tmpstrp[thread][2]*factors[i] ;
//                             sumFactors += factors[i] ;
//                         }
//                     }
//                 }
//             }
// 
//             str[0] += str0 ;
//             str[1] += str1 ;
//             str[2] += str2 ;
//             stra[0] += stra0 ;
//             stra[1] += stra1 ;
//             stra[2] += stra2 ;
//             currentAngle =  atan2( stra[0] - stra[1] , stra[2]) ;
// 
//             str /= sumFactors ;
//             stra /= sumFactors ;
// 
//             return std::make_pair(stra*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()), stra) ;
//         }
//         else
//         {
// //			std::cout << omp_get_thread_num() << std::endl ;
// 
// 
//             int blocks = 1. ;
//             Viscoelasticity * visco = static_cast<Viscoelasticity *>(s.getParent()->getBehaviour()) ;
//             if(visco)
//                 blocks = visco->blocks ;
// 
//             tmpstra.resize( stra.size()*blocks, 0.) ;
//             stra.resize( tmpstra.size(), 0. ) ;
// 
//             Vector strar(0., stra.size() ) ;
//             Vector tmpstrar(0., stra.size() ) ;
// // 			if(m == EFFECTIVE_STRESS)
// // 			{
// // 				#pragma omp parallel for
// // 				for( size_t i = 0 ; i < physicalcache.size() ; i++ )
// // 				{
// // 					Vector tmpstr(0.,vlength) ;
// // 					Vector tmpstra(0.,vlength*blocks) ;
// // 					Vector tmpstrar(0., vlength*blocks ) ;
// // 					DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[physicalcache[i]] ) ;
// // 					static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(ci->getState()).getEssentialAverageFields(EFFECTIVE_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,nullptr,t) ;
// // 				}
// // 			}
// // 			else
// // 			{
// // 				#pragma omp parallel for
// // 				for( size_t i = 0 ; i < physicalcache.size() ; i++ )
// // 				{
// // 					Vector tmpstr(0.,vlength) ;
// // 					Vector tmpstra(0.,vlength*blocks) ;
// // 					Vector tmpstrar(0., vlength*blocks ) ;
// // 					DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[physicalcache[i]] ) ;
// // 					static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(ci->getState()).getEssentialAverageFields(REAL_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,nullptr,t) ;
// // 				}
// // 			}
// 
//             if(m == EFFECTIVE_STRESS)
//             {
// 
//                 double iteratorValue = factors[0] ;
//                 static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(s).getEssentialAverageFields(EFFECTIVE_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,&vm,t) ;
// 
//                 strar = tmpstra*factors[0] ;
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
//                 for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                 {
//                     DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i])) ;
// 
//                     if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured())
//                     {
//                         Vector tmpstr(0.,vlength) ;
//                         Vector tmpstra(0.,vlength*blocks) ;
//                         Vector tmpstrar(0., vlength*blocks ) ;
//                         static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(ci->getState()).getEssentialAverageFields(EFFECTIVE_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,&vm,t) ;
// 
//                         #pragma omp critical
//                         {
//                             stra = stra+tmpstra*factors[i] ;
//                             strar = strar+tmpstrar*factors[i] ;
//                             str = str+tmpstr*factors[i] ;
//                             sumFactors += factors[i] ;
//                         }
//                     }
// 
//                 }
//             }
//             else
//             {
// 
//                 double iteratorValue = factors[0] ;
//                 static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(s).getEssentialAverageFields(REAL_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,&vm,t) ;
// 
// 
//                 strar = tmpstra*factors[0] ;
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
//                 for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                 {
//                     DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i])) ;
// 
//                     if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource() && !ci->getBehaviour()->fractured())
//                     {
// 
//                         Vector tmpstr(0.,vlength) ;
//                         Vector tmpstra(0.,vlength*blocks) ;
//                         Vector tmpstrar(0., vlength*blocks ) ;
//                         static_cast<GeneralizedSpaceTimeViscoElasticElementState &>(ci->getState()).getEssentialAverageFields(REAL_STRESS_FIELD,tmpstr,tmpstra,tmpstrar,&vm,t) ;
// 
//                         #pragma omp critical
//                         {
//                             stra = stra+tmpstra*factors[i] ;
//                             strar = strar+tmpstrar*factors[i] ;
//                             str = str+tmpstr*factors[i] ;
// 
//                             sumFactors += factors[i] ;
//                         }
//                     }
//                 }
// 
//             }
// 
// 
//             str /= sumFactors ;
//             stra /= sumFactors ;
//             strar /= sumFactors ;
// 
//             Vector strFromFullStrain = stra*s.getParent()->getBehaviour()->getTensor(Point()) ;
//             strFromFullStrain += (Vector) (strar*s.getParent()->getBehaviour()->getViscousTensor(Point())) ;
//             Vector strFromStrain(0., str.size()) ;
//             Vector straFromStrain(0., str.size()) ;
//             for(size_t i = 0 ; i < str.size() ; i++)
//             {
//                 strFromStrain[i] = strFromFullStrain[i] ;
//                 straFromStrain[i] = stra[i] ;
//             }
//             strFromStrain -= s.getParent()->getBehaviour()->getImposedStress(Point(0,0,0,t)) ;
//             currentAngle = -0.5*atan2( 2.*straFromStrain[2],  straFromStrain[0] -  straFromStrain[1] ) ;
//             return std::make_pair(strFromStrain, straFromStrain) ;
//         }
// 
//     }
//     else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
//     {
//         if(s.getParent()->getOrder() < CONSTANT_TIME_LINEAR)
//         {
// 
//             if(m == EFFECTIVE_STRESS)
//             {
//                 double iteratorValue = factors[0] ;
//                 s.getAverageField(STRAIN_FIELD,EFFECTIVE_STRESS_FIELD, tmpstra,tmpstr,&vm, 0, t);
//                 currentAngle = 0.5*atan2( tmpstra[3],  tmpstra[0] -  tmpstra[1] ) ;
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
//                 std::vector<Vector> tmpstrat ;
//                 std::vector<Vector> tmpstrt ;
// #ifdef HAVE_OPENMP
//                 for(size_t i = 0 ; i < omp_get_max_threads() ; i++)
//                 {
//                     tmpstrat.push_back(Vector(0., vlength));
//                     tmpstrt.push_back(Vector(0., vlength));
//                 }
// #else
//                 tmpstrat.push_back(Vector(0., vlength));
//                 tmpstrt.push_back(Vector(0., vlength));
// #endif
// 
//                 #pragma omp parallel for reduction(+:stra0,stra1,stra2,str0,str1,str2,sumFactors) firstprivate(tmpstrat,tmpstrt)
//                 for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                 {
//                     DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>(  mesh3d->getInTree(physicalcache[i])) ;
//                     if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
//                     {
//                         double correction = 1. ; //ci->getBehaviour()->getDamageModel()->getState().max() < .999 ;
// #ifdef HAVE_OPENMP
//                         int tnum = omp_get_thread_num() ;
// #else
//                         int tnum = 0 ;
// #endif
//                         ci->getState().getAverageField(STRAIN_FIELD,EFFECTIVE_STRESS_FIELD, tmpstrat[tnum],tmpstrt[tnum],&vm, 0, t);
// 
//                         stra0 += tmpstrat[tnum][0]*factors[i]*correction ;
//                         stra1 += tmpstrat[tnum][1]*factors[i]*correction ;
//                         stra2 += tmpstrat[tnum][2]*factors[i]*correction ;
//                         stra3 += tmpstrat[tnum][3]*factors[i]*correction ;
//                         stra4 += tmpstrat[tnum][4]*factors[i]*correction ;
//                         stra5 += tmpstrat[tnum][5]*factors[i]*correction ;
//                         str0 += tmpstrt[tnum][0]*factors[i]*correction ;
//                         str1 += tmpstrt[tnum][1]*factors[i]*correction ;
//                         str2 += tmpstrt[tnum][2]*factors[i]*correction ;
//                         str3 += tmpstrt[tnum][3]*factors[i]*correction ;
//                         str4 += tmpstrt[tnum][4]*factors[i]*correction ;
//                         str5 += tmpstrt[tnum][5]*factors[i]*correction ;
//                         sumFactors += factors[i]*correction ;
//                     }
//                 }
//             }
//             else
//             {
//                 double iteratorValue = factors[0] ;
//                 s.getAverageField(STRAIN_FIELD,REAL_STRESS_FIELD, tmpstra,tmpstr,&vm, 0, t);
//                 currentAngle = 0.5*atan2( tmpstra[3],  tmpstra[0] -  tmpstra[1] ) ;
//                 stra = tmpstra*factors[0] ;
//                 str = tmpstr*factors[0] ;
//                 sumFactors += factors[0] ;
// 
// 
//                 std::vector<Vector> tmpstrat ;
//                 std::vector<Vector> tmpstrt ;
// #ifdef HAVE_OPENMP
//                 for(size_t i = 0 ; i < omp_get_max_threads() ; i++)
//                 {
//                     tmpstrat.push_back(Vector(0., vlength));
//                     tmpstrt.push_back(Vector(0., vlength));
//                 }
// #else
//                 tmpstrat.push_back(Vector(0., vlength));
//                 tmpstrt.push_back(Vector(0., vlength));
// #endif
// 
//                 #pragma omp parallel for reduction(+:stra0,stra1,stra2,str0,str1,str2,sumFactors) firstprivate(tmpstrat,tmpstrt)
//                 for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//                 {
//                     DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
//                     if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
//                     {
//                         double correction = 1. ;
//                         // 					if(ci->getBehaviour()->getDamageModel())
//                         // 						correction = ci->getBehaviour()->getDamageModel()->getState().max() < .95 ;
// 
// #ifdef HAVE_OPENMP
//                         int tnum = omp_get_thread_num() ;
// #else
//                         int tnum = 0 ;
// #endif
//                         ci->getState().getAverageField(STRAIN_FIELD,REAL_STRESS_FIELD, tmpstrat[tnum],tmpstrt[tnum],&vm, 0, t);
// 
//                         stra0 += tmpstrat[tnum][0]*factors[i]*correction ;
//                         stra1 += tmpstrat[tnum][1]*factors[i]*correction ;
//                         stra2 += tmpstrat[tnum][2]*factors[i]*correction ;
//                         stra3 += tmpstrat[tnum][3]*factors[i]*correction ;
//                         stra4 += tmpstrat[tnum][4]*factors[i]*correction ;
//                         stra5 += tmpstrat[tnum][5]*factors[i]*correction ;
//                         str0 += tmpstrt[tnum][0]*factors[i]*correction ;
//                         str1 += tmpstrt[tnum][1]*factors[i]*correction ;
//                         str2 += tmpstrt[tnum][2]*factors[i]*correction ;
//                         str3 += tmpstrt[tnum][3]*factors[i]*correction ;
//                         str4 += tmpstrt[tnum][4]*factors[i]*correction ;
//                         str5 += tmpstrt[tnum][5]*factors[i]*correction ;
//                         sumFactors += factors[i]*correction ;
//                     }
//                 }
//             }
// 
//             str[0] += str0 ;
//             str[1] += str1 ;
//             str[2] += str2 ;
//             str[3] += str3 ;
//             str[4] += str4 ;
//             str[5] += str5 ;
//             stra[0] += stra0 ;
//             stra[1] += stra1 ;
//             stra[2] += stra2 ;
//             stra[3] += stra3 ;
//             stra[4] += stra4 ;
//             stra[5] += stra5 ;
// 
//             str /= sumFactors ;
//             stra /= sumFactors ;
// 
//             return std::make_pair(stra*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()), stra) ;
//         }
//         else
//         {
// 
//             int blocks = 1. ;
//             Viscoelasticity * visco = dynamic_cast<Viscoelasticity *>(s.getParent()->getBehaviour()) ;
//             if(visco)
//                 blocks = visco->blocks ;
// 
//             tmpstra.resize( stra.size()*blocks, 0.) ;
//             stra.resize( tmpstra.size(), 0. ) ;
// 
//             Vector strar(0., stra.size() ) ;
//             Vector tmpstrar(0., stra.size() ) ;
// 
//             if(m == EFFECTIVE_STRESS)
//                 s.getAverageField(EFFECTIVE_STRESS_FIELD, tmpstr,&vm, -1, t) ;
//             else
//                 s.getAverageField(REAL_STRESS_FIELD, tmpstr,&vm, -1, t) ;
// 
// 
//             s.getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, tmpstra,&vm, -1, t ) ;
//             s.getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, tmpstrar,&vm, -1, t ) ;
// 
//             currentAngle = 0.5*atan2( tmpstra[3],  tmpstra[0] -  tmpstra[1] ) ;
// 
//             strar = tmpstrar*factors[0] ;
//             stra = tmpstra*factors[0] ;
//             str = tmpstr*factors[0] ;
// 
//             sumFactors += factors[0] ;
// // 			cosangle = cos(0.5*atan2( tmpstra[2],  tmpstra[0] -  tmpstra[1] ))*iteratorValue ;
// 
//             //#pragma omp parallel for
//             for( size_t i = 1 ; i < physicalcache.size() ; i++ )
//             {
//                 DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
//                 if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
//                 {
//                     if(m == EFFECTIVE_STRESS)
//                         ci->getState().getAverageField(EFFECTIVE_STRESS_FIELD, tmpstr,&vm, -1, t) ;
//                     else
//                         ci->getState().getAverageField(REAL_STRESS_FIELD, tmpstr,&vm, -1, t) ;
// 
//                     ci->getState().getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, tmpstra,&vm, -1, t ) ;
//                     ci->getState().getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, tmpstrar,&vm, -1, t ) ;
// 
//                     //#pragma omp critical
//                     if(!ci->getBehaviour()->fractured() )
//                     {
//                         strar += tmpstrar*factors[i] ;
//                         stra += tmpstra*factors[i] ;
//                         str += tmpstr*factors[i] ;
//                         sumFactors += factors[i] ;
//                     }
//                 }
//             }
// 
//             str /= sumFactors ;
//             stra /= sumFactors ;
//             strar /= sumFactors ;
// 
// 
//             Vector strFromFullStrain = stra*s.getParent()->getBehaviour()->getTensor(Point()) ;
//             strFromFullStrain += (Vector) (strar*s.getParent()->getBehaviour()->getViscousTensor(Point())) ;
//             Vector strFromStrain(0., str.size()) ;
//             Vector straFromStrain(0., str.size()) ;
//             for(size_t i = 0 ; i < 6 ; i++)
//             {
//                 strFromStrain[i] = strFromFullStrain[i] ;
//                 straFromStrain[i] = stra[i] ;
//             }
// 
//             return std::make_pair(strFromStrain, straFromStrain) ;
//         }
// 
//     }
// 
//     return std::make_pair(smoothedPrincipalStress(s, m), smoothedPrincipalStrain(s)) ;
// }

Vector FractureCriterion::getSmoothedField( ElementState &s ,FieldType f0, double t)
{
    Vector first ;
    Vector strain ;
    Vector stress ;
    Vector strainrate ;
    Vector buffer ;
    int tsize = 3 ;
    int psize = 2 ;
    if(mesh3d)
    {
        tsize = 6 ;
        psize = 3 ;
    }
    bool spaceTime = s.getParent()->getOrder() >= CONSTANT_TIME_LINEAR ;
    VirtualMachine vm ;
    if(f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD )
    {
        //we first need to compute the strain field
        if(!spaceTime)
        {
            double sumFactors(0) ;
            
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                IntegrableEntity *ci ;
                if(mesh2d)
                    ci = static_cast<DelaunayTriangle *>(mesh2d->getInTree(physicalcache[i])) ;
                else
                    ci = static_cast<DelaunayTetrahedron *>(mesh3d->getInTree(physicalcache[i])) ;
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(STRAIN_FIELD, buffer, &vm, 0, t);
                    if(!strain.size())
                        strain.resize(0., buffer.size());
                    strain += buffer*factors[i] ;
                    sumFactors += factors[i] ;
                }
            }
            strain /= sumFactors ;
        }
        else
        {
            double sumFactors(0) ;
            Vector tmpstrain ;
            Vector tmpstrainrate ;
            
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                IntegrableEntity *ci ;
                if(mesh2d)
                    ci = static_cast<DelaunayTriangle *>(mesh2d->getInTree(physicalcache[i])) ;
                else
                    ci = static_cast<DelaunayTetrahedron *>(mesh3d->getInTree(physicalcache[i])) ;
                
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, &vm, 0, t);
                    if(!tmpstrain.size())
                        tmpstrain.resize(0., buffer.size());
                    tmpstrain += buffer*factors[i] ;
                    sumFactors += factors[i] ;
                }
            }
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, &vm, 0, t);
                    if(!tmpstrainrate.size())
                        tmpstrainrate.resize(0., buffer.size());
                    tmpstrainrate += buffer*factors[i] ;
                }
            }
            tmpstrain /= sumFactors ;
            tmpstrainrate /=sumFactors ;
            
            Vector tmpstress = tmpstrain*s.getParent()->getBehaviour()->getTensor(Point()) + (Vector) (tmpstrainrate*s.getParent()->getBehaviour()->getViscousTensor(Point())) ;
            stress.resize(tsize, 0.) ;
            strain.resize(tsize, 0.) ;
            for(size_t i = 0 ; i < tsize ; i++)
            {
                stress[i] = tmpstress[i] ;
                strain[i] = tmpstrain[i] ;
            }
        }
        
        if(f0 == PRINCIPAL_STRAIN_FIELD)
        {
            first.resize(psize);
            first = toPrincipal(strain) ;
        }
        if(f0 == REAL_STRESS_FIELD)
        {
            first.resize(tsize);
            if(!spaceTime)
                first = strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()) ;
            else
                first = stress ;
        }
        if(f0 == EFFECTIVE_STRESS_FIELD)
        {
            first.resize(tsize);
            first = strain*s.getParent()->getBehaviour()->param ;
        }
        if(f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
        {
            first.resize(psize);
            first = toPrincipal(strain*s.getParent()->getBehaviour()->param) ;
        }
        if(f0 == PRINCIPAL_REAL_STRESS_FIELD)
        {
            first.resize(psize);
            if(!spaceTime)
                first = toPrincipal(strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
            else
                first = toPrincipal(stress) ;
        }
        
    }
    else
    {
        double sumFactors(0) ;
        
        for(size_t i = 0 ; i < physicalcache.size() ; i++)
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
            if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
            {
                ci->getState().getAverageField(f0, buffer, &vm, 0, t);
                if(!first.size())
                    first.resize(0., buffer.size());
                first += buffer*factors[i] ;
                sumFactors += factors[i] ;
            }
        }

        first /= sumFactors ;
    }
    
    
    return first ;
}


std::pair<Vector, Vector> FractureCriterion::getSmoothedFields( ElementState &s ,FieldType f0, FieldType f1, double t)
{
    Vector first ;
    Vector second ;
    Vector strain ;
    Vector stress ;
    Vector strainrate ;
    Vector buffer ;
    int tsize = 3 ;
    int psize = 2 ;
    if(mesh3d)
    {
        tsize = 6 ;
        psize = 3 ;
    }
    bool spaceTime = s.getParent()->getOrder() >= CONSTANT_TIME_LINEAR ;
    VirtualMachine vm ;
    if(f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f0 == STRAIN_FIELD ||
        f1 == PRINCIPAL_STRAIN_FIELD || f1 == REAL_STRESS_FIELD || f1 == EFFECTIVE_STRESS_FIELD || f1 == PRINCIPAL_REAL_STRESS_FIELD || f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f1 == STRAIN_FIELD
    )
    {
        //we first need to compute the strain field
        if(!spaceTime)
        {
            buffer.resize(tsize, 0.);
            strain.resize(tsize, 0.);
            double sumFactors(0) ;
            
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                IntegrableEntity *ci ;
                if(mesh2d)
                    ci = static_cast<DelaunayTriangle *>(mesh2d->getInTree(physicalcache[i])) ;
                else
                    ci = static_cast<DelaunayTetrahedron *>(mesh3d->getInTree(physicalcache[i])) ;
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(STRAIN_FIELD, buffer, &vm, 0, t);
                    strain += buffer*factors[i] ;
                    sumFactors += factors[i] ;
                }
            }
            strain /= sumFactors ;
        }
        else
        {
            size_t blocks = s.getParent()->getBehaviour()->getNumberOfDegreesOfFreedom()/2 ;
            for(size_t i = 0 ; i < physicalcache.size() && !blocks; i++)
            {
                IntegrableEntity *ci ;
                if(mesh2d)
                    ci = static_cast<DelaunayTriangle *>(mesh2d->getInTree(physicalcache[i])) ;
                else
                    ci = static_cast<DelaunayTetrahedron *>(mesh3d->getInTree(physicalcache[i])) ;
                blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom()/(mesh2d?2:3) ;
            }
            Vector tmpstrain ;
            Vector tmpstrainrate ;
                
            tmpstrain.resize(fieldTypeElementarySize(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, mesh2d?SPACE_TWO_DIMENSIONAL:SPACE_THREE_DIMENSIONAL, blocks), 0.) ;
            buffer.resize(fieldTypeElementarySize(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, mesh2d?SPACE_TWO_DIMENSIONAL:SPACE_THREE_DIMENSIONAL, blocks), 0.);
            tmpstrainrate.resize(fieldTypeElementarySize(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, mesh2d?SPACE_TWO_DIMENSIONAL:SPACE_THREE_DIMENSIONAL, blocks), 0.) ;
            double sumFactors(0) ;
           
            
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                IntegrableEntity *ci ;
                if(mesh2d)
                    ci = static_cast<DelaunayTriangle *>(mesh2d->getInTree(physicalcache[i])) ;
                else
                    ci = static_cast<DelaunayTetrahedron *>(mesh3d->getInTree(physicalcache[i])) ;
                
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, &vm, 0, t);
                    if(!tmpstrain.size())
                        tmpstrain.resize(buffer.size(), 0.);
                    tmpstrain += buffer*factors[i] ;
                    sumFactors += factors[i] ;
                }
            }
            buffer.resize(fieldTypeElementarySize(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, mesh2d?SPACE_TWO_DIMENSIONAL:SPACE_THREE_DIMENSIONAL, blocks), 0.);
            for(size_t i = 0 ; i < physicalcache.size() ; i++)
            {
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
                if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
                {
                    ci->getState().getAverageField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, &vm, 0, t);
                    if(!tmpstrainrate.size())
                        tmpstrainrate.resize(buffer.size(), 0.);
                    tmpstrainrate += buffer*factors[i] ;
                }
            }
            tmpstrain /= sumFactors ;
            tmpstrainrate /=sumFactors ;
            
            Vector tmpstress = tmpstrain*s.getParent()->getBehaviour()->getTensor(Point()) + (Vector) (tmpstrainrate*s.getParent()->getBehaviour()->getViscousTensor(Point())) ;
            stress.resize(tsize, 0.) ;
            strain.resize(tsize, 0.) ;
            for(size_t i = 0 ; i < tsize ; i++)
            {
                stress[i] = tmpstress[i] ;
                strain[i] = tmpstrain[i] ;
            }
        }
        
        if(f0 == PRINCIPAL_STRAIN_FIELD)
        {
            first.resize(psize);
            first = toPrincipal(strain) ;
        }
        if(f1 == PRINCIPAL_STRAIN_FIELD)
        {
            second.resize(psize);
            second = toPrincipal(strain) ;
        }
        if(f0 == REAL_STRESS_FIELD)
        {
            first.resize(tsize);
            if(!spaceTime)
                first = strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()) ;
            else
                first = stress ;
        }
        if(f1 == REAL_STRESS_FIELD)
        {
            second.resize(tsize);
            if(!spaceTime)
                second = strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()) ;
            else
                second = stress ;
        }
        if(f0 == EFFECTIVE_STRESS_FIELD)
        {
            first.resize(tsize);
            first = strain*s.getParent()->getBehaviour()->param ;
        }
        if(f1 == EFFECTIVE_STRESS_FIELD)
        {
            second.resize(tsize);
            second = strain*s.getParent()->getBehaviour()->param ;
        }
        if(f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
        {
            first.resize(psize);
            first = toPrincipal(strain*s.getParent()->getBehaviour()->param) ;
        }
        if(f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD)
        {
            second.resize(psize);
            second = toPrincipal(strain*s.getParent()->getBehaviour()->param) ;
        }
        if(f0 == PRINCIPAL_REAL_STRESS_FIELD)
        {
            first.resize(psize);
            if(!spaceTime)
                first = toPrincipal(strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
            else
                first = toPrincipal(stress) ;
        }
        if(f1 == PRINCIPAL_REAL_STRESS_FIELD)
        {
            second.resize(psize);
            if(!spaceTime)
                second = toPrincipal(strain*s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())) ;
            else
                second = toPrincipal(stress) ;
        }
        
        if(f0 == STRAIN_FIELD)
        {
            first.resize(tsize) ;
            first = strain ;
        }
        if(f1 == STRAIN_FIELD)
        {
            second.resize(tsize) ;
            second = strain ;
        }
        
    }
    else
    {
        double sumFactors(0) ;
        
        for(size_t i = 0 ; i < physicalcache.size() ; i++)
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
            if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
            {
                ci->getState().getAverageField(f0, buffer, &vm, 0, t);
                if(!first.size())
                    first.resize(0., buffer.size());
                first += buffer*factors[i] ;
                sumFactors += factors[i] ;
            }
        }
        for(size_t i = 0 ; i < physicalcache.size() ; i++)
        {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *>(  mesh2d->getInTree(physicalcache[i]) ) ;
            if(ci->getBehaviour()->getSource() == s.getParent()->getBehaviour()->getSource())
            {
                ci->getState().getAverageField(f1, buffer, &vm, 0, t);
                if(!second.size())
                    second.resize(0., buffer.size());
                second += buffer*factors[i] ;
            }
        }
        first /= sumFactors ;
        second /= sumFactors ;
    }
    
    
    return std::make_pair(first, second) ;
}

double FractureCriterion::getDeltaEnergy(const ElementState & s, double delta_d)
{
    Assembly K ;
    double originalscore = 0 ;
    double originalenergy = 0 ;
    double score = 0 ;
    double energy = 0 ;
    double volume = 0 ;

    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTriangle *> elements ;
        std::vector<LinearForm *> behaviours ;
        TriElement father(LINEAR) ;

        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTriangle * tric = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTriangle(nullptr, nullptr,   tric->first,  tric->second,   tric->third,  nullptr) );
            elements.back()->setBehaviour(mesh2d, tric->getBehaviour()->getCopy()) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(mesh2d, false) ;
            K.add(elements.back());
            volume += elements.back()->area() ;
            for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
            {
                if(tric->getNeighbour(j)->isTriangle )
                {
                    DelaunayTriangle * tri = static_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*2];
                                    double ey = tri->getState().getDisplacements()[l*2+1];
                                    K.setPoint(ex, ey ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }
        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double a = elements[i]->area() ;
                originalenergy += elements[i]->getState().elasticEnergy()*a ;
            }
        }

        std::valarray<Point *> nularray(0) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            elements[i]->setBoundingPoints(nularray) ;
            delete elements[i] ;
        }

    }
    else
    {
        Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTetrahedron *> elements ;
        std::vector<LinearForm *> behaviours ;
        TetrahedralElement father ;
        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTetrahedron(nullptr, nullptr,   tet->first,  tet->second,   tet->third, tet->fourth,  nullptr) );
            elements.back()->setBehaviour(s.getMesh3D(),tet->getBehaviour()->getCopy()) ;

            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh3D(), false) ;
            K.add(elements.back());
            volume+= elements.back()->volume() ;

            for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
            {
                if(tet->getNeighbour(j)->isTetrahedron() )
                {
                    DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*3];
                                    double ey = tri->getState().getDisplacements()[l*3+1];
                                    double ez = tri->getState().getDisplacements()[l*3+2];
                                    K.setPoint(ex, ey, ez ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }


        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        double energy = 0 ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double v = elements[i]->volume() ;
                originalenergy += elements[i]->getState().elasticEnergy()*v ;
            }
        }
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            delete elements[i] ;
        }
    }



    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTriangle *> elements ;
        std::vector<LinearForm *> behaviours ;
        TriElement father(LINEAR) ;

        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTriangle * tric = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTriangle(nullptr, nullptr,   tric->first,  tric->second,   tric->third,  nullptr) );
            elements.back()->setBehaviour(s.getMesh2D(), tric->getBehaviour()->getCopy()) ;

            if(tric == s.getParent())
                elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

            elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh2D(),false) ;
            K.add(elements.back());
            volume += elements.back()->area() ;
            for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
            {
                if(tric->getNeighbour(j)->isTriangle )
                {
                    DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*2];
                                    double ey = tri->getState().getDisplacements()[l*2+1];
                                    K.setPoint(ex, ey ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double a = elements[i]->area() ;
                energy += elements[i]->getState().elasticEnergy()*a ;
            }
        }

        std::valarray<Point *> nularray(0) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            elements[i]->setBoundingPoints(nularray) ;
            delete elements[i] ;
        }
    }
    else
    {
        Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTetrahedron *> elements ;
        std::vector<LinearForm *> behaviours ;
        TetrahedralElement father ;
        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * >( mesh3d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTetrahedron(nullptr, nullptr,   tet->first,  tet->second,   tet->third, tet->fourth,  nullptr) );
            elements.back()->setBehaviour(s.getMesh3D(), tet->getBehaviour()->getCopy()) ;
            if(tet == s.getParent())
                elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
            elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh3D(),false) ;
            K.add(elements.back());
            volume+= elements.back()->volume() ;

            for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
            {
                if(tet->getNeighbour(j)->isTetrahedron() )
                {
                    DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                int id = tri->getBoundingPoint(k).getId() ;
                                double ex = tri->getState().getDisplacements()[k*3];
                                double ey = tri->getState().getDisplacements()[k*3+1];
                                double ez = tri->getState().getDisplacements()[k*3+2];
                                K.setPoint(ex, ey, ez ,id);
                            }
                        }
                    }
                }
            }


        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        double energy = 0 ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
            }
        }
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            delete elements[i] ;
        }
    }


    return (originalenergy-energy)/(delta_d) ;
}

std::pair<double, double> FractureCriterion::getCrackOpeningAndSlip(const ElementState & s)
{
    Matrix rotationMatrix(2,2) ;
    rotationMatrix[0][0] = cos(-currentAngle) ;
    rotationMatrix[0][1] = sin(-currentAngle) ;
    rotationMatrix[1][0] = -sin(-currentAngle) ;
    rotationMatrix[1][1] = cos(-currentAngle) ;

    Vector displacementLeft(0., 2) ;
    Vector displacementRight(0., 2) ;
    double countLeft = 0 ;
    double countRight = 0 ;

    for(size_t i = 0 ; i < physicalcache.size() ; i++)
    {
        DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(physicalcache[i])) ;
        for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
        {
            if((ci->getBoundingPoint(j) - s.getParent()->getCenter()).angle()-currentAngle > 0 &&
                    (ci->getBoundingPoint(j) - s.getParent()->getCenter()).angle()-currentAngle < M_PI
              )
            {
                displacementLeft[0] += ci->getState().getDisplacements()[2*j]*factors[i] ;
                displacementLeft[1] += ci->getState().getDisplacements()[2*j+1]*factors[i] ;
                countLeft += factors[i];
            }
            else
            {
                displacementRight[0] += ci->getState().getDisplacements()[2*j]*factors[i] ;
                displacementRight[1] += ci->getState().getDisplacements()[2*j+1]*factors[i] ;
                countRight += factors[i];
            }
        }
    }


    displacementLeft /= countLeft ;
    displacementRight /= countRight ;
    Vector delta = rotationMatrix*(displacementLeft-displacementRight) ;

    return std::make_pair(delta[0], delta[1]) ;
}

void FractureCriterion::initialiseCacheNew( FeatureTree * ft, const ElementState & s)
{
    if(ft->is2D())
    {
        Function x = Function("x")-s.getParent()->getCenter().getX() ;
        Function y = Function("y")-s.getParent()->getCenter().getY() ;
        Function rr = x*x+y*y ;
        Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
        Function smooth =  (smoothingType == GAUSSIAN_NONCOMPACT)?f_exp(rrn*-0.5):(rrn-1.)*(rrn-1.)*f_positivity(1.-rrn) ;
        
        double overlap = (smoothingType == QUARTIC_COMPACT)?6.:8. ;
        Circle epsilonAll( std::max(physicalCharacteristicRadius, s.getParent()->getRadius()*2. )*overlap+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        Circle epsilonReduced(physicalCharacteristicRadius*1.1+s.getParent()->getRadius(),s.getParent()->getCenter()) ;
        int cacheID = ft->get2DMesh()->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;
        int cachecoreID = ft->get2DMesh()->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        mesh2d = ft->get2DMesh() ;
    }
    if(ft->is3D())
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
        int cacheID = ft->get3DMesh()->generateCache(&epsilonAll, s.getParent()->getBehaviour()->getSource()) ;
        int cachecoreID = ft->get3DMesh()->generateCache(&epsilonReduced, s.getParent()->getBehaviour()->getSource(), smooth) ;
        mesh3d = ft->get3DMesh() ;
    }
}

void FractureCriterion::initialiseCache(const ElementState & s)
{
    DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
    DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
    if(testedTri)
    {
        if(!cache.empty())
        {
            cache.clear();
        }
        double overlap = (smoothingType == QUARTIC_COMPACT)?6.:8. ;
        Circle epsilon( std::max(physicalCharacteristicRadius, testedTri->getRadius()*2. )*overlap+testedTri->getRadius(),testedTri->getCenter()) ;
        if(!s.getMesh2D())
            return ;
        mesh2d = s.getMesh2D()  ;

        std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getNeighbouringElementsInGeometry(testedTri, &epsilon);
        std::vector<DelaunayTriangle *> neighbourhood ;
        std::vector<DelaunayTriangle *> neighbours = s.getMesh2D()->getNeighbourhood(testedTri) ;
        for(auto & n : neighbours)
        {
            if(n->getBehaviour() &&
                    n->getBehaviour()->type != VOID_BEHAVIOUR &&
                    n->getBehaviour()->getFractureCriterion())
            {
                neighbourhood.push_back(n);
                cache.push_back(n->index);
            }
        }
        for(size_t i = 0 ; i < tempcache.size() ; i++)
        {
            if(tempcache[i]->getBehaviour() &&
                    tempcache[i]->getBehaviour()->type != VOID_BEHAVIOUR &&
                    tempcache[i]->getBehaviour()->getFractureCriterion())
            {
                bool inNeighbourhood = false ;
                for(size_t j = 0 ; j < neighbourhood.size() ; j++)
                {
                    if(neighbourhood[j] == tempcache[i])
                    {
                        inNeighbourhood = true ;
                        break ;
                    }
                }
                if(!inNeighbourhood)
                {
                    cache.push_back(tempcache[i]->index);
                }
            }
        }


        if(cache.empty())
            cache.push_back(testedTri->index);

        initialiseFactors(s);

    }
    std::sort(cache.begin(), cache.end()) ;
}

std::pair<double, double> FractureCriterion::getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const
{
    Assembly K ;
    double originalscore = 0 ;
    double originalenergy = 0 ;
    double score = 0 ;
    double energy = 0 ;
    double volume = 0 ;

    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTriangle *> elements ;
        std::vector<LinearForm *> behaviours ;
        TriElement father(LINEAR) ;

        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ( mesh2d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTriangle(nullptr, nullptr,   tric->first,  tric->second,   tric->third,  nullptr) );
            elements.back()->setBehaviour(s.getMesh2D(), tric->getBehaviour()->getCopy()) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh2D(),false) ;
            K.add(elements.back());
            volume += elements.back()->area() ;
            for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
            {
                if(tric->getNeighbour(j)->isTriangle )
                {
                    DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*2];
                                    double ey = tri->getState().getDisplacements()[l*2+1];
                                    K.setPoint(ex, ey ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }
        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double a = elements[i]->area() ;
                double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
                if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
                    originalscore += 1./(1.-sc) ;
                originalenergy += elements[i]->getState().elasticEnergy()*a ;
            }
        }

        std::valarray<Point *> nularray(0) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            elements[i]->setBoundingPoints(nularray) ;
            delete elements[i] ;
        }

    }
    else
    {
        Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTetrahedron *> elements ;
        std::vector<LinearForm *> behaviours ;
        TetrahedralElement father ;
        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTetrahedron(nullptr, nullptr,   tet->first,  tet->second,   tet->third, tet->fourth,  nullptr) );
            elements.back()->setBehaviour(s.getMesh3D(), tet->getBehaviour()->getCopy()) ;

            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh3D(),false) ;
            K.add(elements.back());
            volume+= elements.back()->volume() ;

            for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
            {
                if(tet->getNeighbour(j)->isTetrahedron() )
                {
                    DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*3];
                                    double ey = tri->getState().getDisplacements()[l*3+1];
                                    double ez = tri->getState().getDisplacements()[l*3+2];
                                    K.setPoint(ex, ey, ez ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }


        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        double energy = 0 ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double v = elements[i]->volume() ;
                double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
                if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
                    originalscore += 1./(1.-sc) ;
                originalenergy += elements[i]->getState().elasticEnergy()*v ;
            }
        }
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            delete elements[i] ;
        }
    }



    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTriangle *> elements ;
        std::vector<LinearForm *> behaviours ;
        TriElement father(LINEAR) ;

        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ( mesh2d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTriangle(nullptr, nullptr,   tric->first,   tric->second,    tric->third,  nullptr) );
            elements.back()->setBehaviour(s.getMesh2D(),tric->getBehaviour()->getCopy()) ;

            if(tric == s.getParent())
                elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

            elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh2D(), false) ;
            K.add(elements.back());
            volume += elements.back()->area() ;
            for(size_t j = 0 ; j <  tric->neighbour.size() ; j++)
            {
                if( tric->getNeighbour(j)->isTriangle )
                {
                    DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>( tric->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
                                {
                                    int id = tri->getBoundingPoint(k).getId() ;
                                    double ex = tri->getState().getDisplacements()[l*2];
                                    double ey = tri->getState().getDisplacements()[l*2+1];
                                    K.setPoint(ex, ey ,id);
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double a = elements[i]->area() ;
                double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
                if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
                    score += 1./(1.-sc);
                energy += elements[i]->getState().elasticEnergy()*a ;
            }
        }

        std::valarray<Point *> nularray(0) ;

        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            elements[i]->setBoundingPoints(nularray) ;
            delete elements[i] ;
        }
    }
    else
    {
        Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

        std::vector<DelaunayTetrahedron *> elements ;
        std::vector<LinearForm *> behaviours ;
        TetrahedralElement father ;
        for(size_t i = 0 ; i < cache.size() ; i++)
        {
            DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * > ( mesh3d->getInTree(cache[i])) ;
            elements.push_back(new DelaunayTetrahedron(nullptr, nullptr,   tet->first,  tet->second,   tet->third, tet->fourth,  nullptr) );
            elements.back()->setBehaviour(s.getMesh3D(), tet->getBehaviour()->getCopy()) ;
            if(tet == s.getParent())
                elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
            elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
            elements.back()->refresh(&father);
            elements.back()->getState().initialize(s.getMesh3D(),false) ;
            K.add(elements.back());
            volume+= elements.back()->volume() ;

            for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
            {
                if(tet->getNeighbour(j)->isTetrahedron() )
                {
                    DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
                    if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
                        {
                            if(!c.in(tri->getBoundingPoint(k)))
                            {
                                int id = tri->getBoundingPoint(k).getId() ;
                                double ex = tri->getState().getDisplacements()[k*3];
                                double ey = tri->getState().getDisplacements()[k*3+1];
                                double ez = tri->getState().getDisplacements()[k*3+2];
                                K.setPoint(ex, ey, ez ,id);
                            }
                        }
                    }
                }
            }


        }

        K.cgsolve() ;
        for(size_t i = 0 ; i < elements.size() ; i++)
            elements[i]->step(0., &K.getDisplacements()) ;

        double energy = 0 ;
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            if(elements[i]->getBehaviour()->getFractureCriterion())
            {
                double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
                if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
                    score += 1./(1.-sc) ;
                energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
            }
        }
        for(size_t i = 0 ; i < elements.size() ; i++)
        {
            delete elements[i] ;
        }
    }


    return std::make_pair( (energy-originalenergy)/(delta_d),(score-originalscore)/(delta_d)) ;
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

std::pair<double, double> FractureCriterion::setChange(const ElementState &s, double thresholdScore)
{
    stable = true ;

    if( !s.getParent()->getBehaviour()->getDamageModel())
        return std::make_pair(0.,0.) ;
    if(cache.size() == 0)
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

            for(size_t i = 0 ; i< cache.size() ; i++)
            {
                DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i])) ;
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

            for(size_t i = 0 ; i< cache.size() ; i++)
            {
                DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;
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
    if(energyIndexed && s.getDeltaTime() > POINT_TOLERANCE_2D )
        noEnergyUpdate = true ;

    if(energyIndexed && noEnergyUpdate /*&& scoreAtState > -.5*/)
    {
        noEnergyUpdate = false ;
        currentEnergy = 0 ;
        bool buildvolume = false ;
        if(neighbourhoodvolume < 0)
        {
            neighbourhoodvolume = 0 ;
            buildvolume = true ;
        }

        if(mesh2d)
        {
            Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
            for(size_t i = 0 ; i < getCache().size() ; i++)
            {
                DelaunayTriangle *tri = static_cast<DelaunayTriangle *>( mesh2d->getInTree(getCache()[i]) ) ;
                if(c.in(*tri->first) || c.in(*tri->second) || c.in(*tri->third))
                {
                    currentEnergy += tri->getState().elasticEnergy()*tri->area() ;
                    if(buildvolume)
                        neighbourhoodvolume += tri->area() ;
                }
            }
        }
        else if(mesh3d)
        {
            Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

            for(size_t i = 0 ; i < getCache().size() ; i++)
            {
                DelaunayTetrahedron *tet = static_cast<DelaunayTetrahedron *>(  mesh3d->getInTree(getCache()[i]) ) ;
                if(c.in(*tet->first) || c.in(*tet->second) || c.in(*tet->third) || c.in(*tet->fourth))
                {
                    currentEnergy += tet->getState().elasticEnergy()*tet->volume() ;
                    if(buildvolume)
                        neighbourhoodvolume += tet->volume() ;
                }
            }
        }
        deltaEnergyAtState = (currentEnergy-previousEnergy) ;

        previousEnergy = currentEnergy ;

    }

    if(energyIndexed /*&& scoreAtState > -.5*/)
    {
        std::pair<double, double> dedc = getDeltaEnergyDeltaCriterion(s, 0.0001) ;
        energyDamageDifferential = dedc.first ;
        criterionDamageDifferential = dedc.second ;
    }

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

void FractureCriterion::computeNonLocalState(ElementState &s, NonLocalSmoothingType st)
{
    metAtStep = false ;

// 	if( s.getParent()->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel()->fractured())
// 	{
// 		metAtStep = 0 ; scoreAtState > 0 ;
// 		nonLocalScoreAtState = -1 ; scoreAtState ;
// 		return  ;
// 	}

    DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
    DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
    HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;

    switch (st)
    {
    case NULL_SMOOTH :
    {
        metAtStep = scoreAtState > 0 ;
        nonLocalScoreAtState = scoreAtState ;
        return ;
    }
    case MAX_PROXIMITY_SMOOTH :
    {
        if(testedTri)
        {
            if(scoreAtState <= 2.*scoreTolerance)
            {
                nonLocalScoreAtState = scoreAtState ;
                metAtStep = false ;
                return ;
            }

            if(!cache.empty())
            {
                nonLocalScoreAtState = scoreAtState ;
                for(size_t i = 0 ; i< cache.size() ; i++)
                {
                    DelaunayTriangle * ci = static_cast<DelaunayTriangle *>( mesh2d->getInTree(cache[i])) ;

                    if(squareDist2D(ci->getCenter(), s.getParent()->getCenter()) > physicalCharacteristicRadius*physicalCharacteristicRadius
                            || ci->getBehaviour()->getSource() !=  s.getParent()->getBehaviour()->getSource()
                            || ci->getBehaviour()->fractured()
                      )
                        continue ;

                    nonLocalScoreAtState = std::max(nonLocalScoreAtState, ci->getBehaviour()->getFractureCriterion()->scoreAtState) ;
                }
            }
            metAtStep = true ;
            return ;

        }

        if(testedTet)
        {
            nonLocalScoreAtState = scoreAtState ;

            if (scoreAtState <= 0)
            {
                metAtStep = false ;
                return  ;
            }

            double maxNeighbourhoodScore = 0 ;
// 				double matchedArea = 0 ;
// 				std::map<double, DelaunayTetrahedron *> scores ;
// 				std::map<DelaunayTetrahedron *, double> areatemp ;
// 				DelaunayTetrahedron * maxLocus = nullptr;

            if(!cache.empty())
            {
                for(size_t i = 0 ; i< cache.size() ; i++)
                {
                    DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;

                    if( !ci->getBehaviour()->fractured())
                    {
                        double s = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
// 							scores[-s] =  ci;
                        if(s > maxNeighbourhoodScore)
                        {
                            maxNeighbourhoodScore = s ;
// 								maxLocus = ci ;
                        }
                    }
// 						else if(ci->getBehaviour()->fractured())
// 						{
// 							double s = POINT_TOLERANCE_2D ;
// // 							scores[-s] =  ci;
// 						}
// 						areatemp[ci] = area[i] ;

                }
            }

            if(maxNeighbourhoodScore < 0)
            {
                metAtStep = false ;
                return  ;
            }

            std::vector<DelaunayTetrahedron *> maxloci ;

            for(size_t i = 0 ; i< cache.size() ; i++)
            {
                DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;
                if(ci->getBehaviour()->getFractureCriterion())
                    if(std::abs(ci->getBehaviour()->getFractureCriterion()->getScoreAtState()-maxNeighbourhoodScore) < scoreTolerance)
                        maxloci.push_back(ci) ;
            }


            for(size_t i = 0 ; i < maxloci.size() ; i++)
                if(squareDist3D(maxloci[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
                {
                    metAtStep = true ;
                    return  ;
                }
            metAtStep = false ;
            return  ;
        }
        else if(testedHex)
        {
            std::set<HexahedralElement *> neighbourhood ;
            std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
            for(size_t i = 0 ; i < neighbours.size() ; i++)
            {
                for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
                {
                    if(neighbours[i]->neighbourhood[j] != testedHex
                            && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
                        neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
                }
            }
            double score = grade(s) ;
            double maxNeighbourhoodScore = 0 ;
            if(!neighbourhood.empty())
            {
                for(auto i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
                {
                    if((*i)->getBehaviour()->getFractureCriterion()
                            && !(*i)->getBehaviour()->fractured())
                        maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
                                                         (*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
                    if((*i)->getBehaviour()->changed())
                    {
                        maxNeighbourhoodScore = 10.*score ;
                        break ;
                    }

                    if (maxNeighbourhoodScore > score)
                        break ;

                }
            }

            if( score > 0 )
            {
                if(score > maxNeighbourhoodScore)
                {
                    metAtStep = true ;
                    return  ;
                }
            }
            metAtStep = false ;
            return  ;
        }
        else
        {
            std::cout << " criterion not implemented for this kind of element" << std::endl ;
            metAtStep = false ;
            return  ;
        }
    }
    case GAUSSIAN_SMOOTH :
    {
        if(testedTri)
        {
            double smoothscore = smoothedScore(s) ;

            metAtStep =  (smoothscore > scoreTolerance) ;
            if(scoreAtState >= 0)
                std::cout << scoreAtState << "  " << smoothscore << std::endl ;
            nonLocalScoreAtState = smoothscore ;
            return ;

        }
        if(testedTet)
        {


            double str = 0 ;
            double fact = 0 ;

            for(size_t i = 0 ; i< cache.size() ; i++)
            {
                DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>( mesh3d->getInTree(cache[i])) ;
                double dc = squareDist3D(s.getParent()->getCenter(), ci->getCenter()) ;
                if(ci->getBehaviour()->getFractureCriterion())
                {
                    double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
                    double a = ci->volume() ;
                    str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                    fact+=d*a ;
                    if(mirroring == MIRROR_X && std::abs(ci->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_Y &&  std::abs(ci->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_Z &&  std::abs(ci->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                }
            }

            double smoothscore = str/fact ;
            metAtStep =  smoothscore > 0 ;
            nonLocalScoreAtState = smoothscore ;
            return ;
        }
        else if(testedHex)
        {
            std::set<HexahedralElement *> neighbourhood ;
            std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
            for(size_t i = 0 ; i < neighbours.size() ; i++)
            {
                for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
                {
                    neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
                }
            }
            double str = 0 ;
            double fact = 0 ;

            for(auto ci = neighbourhood.begin() ; ci != neighbourhood.end() ; ++ci)
            {
                double dc = squareDist3D(s.getParent()->getCenter(), (*ci)->getCenter()) ;
                if((*ci)->getBehaviour()->getFractureCriterion())
                {
                    double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
                    double a = (*ci)->volume() ;
                    str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                    fact+=d*a ;
                    if(mirroring == MIRROR_X && std::abs((*ci)->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_Y &&  std::abs((*ci)->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_Z &&  std::abs((*ci)->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().getX()  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().getY()  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                    if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().getZ()  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
                    {
                        str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
                        fact+=d*a ;
                    }
                }
            }

            double smoothscore = str/fact ;
            metAtStep =  smoothscore > 0 ;
            nonLocalScoreAtState = smoothscore ;
            return ;
        }
        else
        {
            std::cout << " criterion not implemented for this kind of element" << std::endl ;
            metAtStep = false ;
            return  ;
        }
    }
    }
    //shut up the compiler
    metAtStep = false ;
    return  ;
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



