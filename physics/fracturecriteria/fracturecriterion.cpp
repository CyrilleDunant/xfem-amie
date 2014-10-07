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
//     deltaEnergyAtState(0),
//     energyDamageDifferential(0),
    criterionDamageDifferential(0),
//     energyIndexed(false),
//     noEnergyUpdate(true),
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
    Vector ss = getSmoothedField( REAL_STRESS_FIELD, s) ;
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

void FractureCriterion::initialiseFactors( ElementState & s)
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

void FractureCriterion::initialiseCacheOld( ElementState & s)
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
//         std::pair<double, double> dedc = getDeltaEnergyDeltaCriterion(s, 0.0001) ;
//         energyDamageDifferential = dedc.first ;
//         criterionDamageDifferential = dedc.second ;
//     }

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



