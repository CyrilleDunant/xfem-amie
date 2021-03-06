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
#include "mohrcoulomb.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../damagemodels/damagemodel.h"

namespace Amie
{

MohrCoulomb::MohrCoulomb ( double up, double down): upVal(up), downVal(down)
{
    metInTension = false ;
    metInCompression = false ;
}


MohrCoulomb::~MohrCoulomb()
{
}

double MohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->getDamageModel()->fractured() )
        return -1 ;

    Vector pstress(0., s.getParent()->spaceDimensions()) ;
    s.getField( PRINCIPAL_REAL_STRESS_FIELD, s.getParent()->getCenter(), pstress, false) ;

    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;

    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = !metInCompression ;

    if( metInTension )
    {   
        return std::abs( maxStress/upVal )-1. ;
    }

    metInCompression = true ;
    return std::abs( minStress/downVal )-1. ;
}

FractureCriterion *MohrCoulomb::getCopy() const
{
    return new MohrCoulomb( *this ) ;
}

NonLocalMohrCoulomb::NonLocalMohrCoulomb( double up, double down) : upVal( up ), downVal( down )
{
    metInTension = false ;
    metInCompression = false ;
}


NonLocalMohrCoulomb::~NonLocalMohrCoulomb()
{

}

double NonLocalMohrCoulomb::grade( ElementState &s )
{
    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    Vector stress =  getSmoothedField(PRINCIPAL_REAL_STRESS_FIELD , s) ;
    
    double maxStress = stress.max() ;
    double minStress = stress.min() ;

    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = !metInCompression ;

//     if(std::abs(maxStress) > 1e-6 || std::abs(minStress) > 1e-6)
//     {
//         std::cout << maxStress << "  " << upVal <<" -> " << maxStress/upVal - 1. << " ; " <<minStress << "  " << downVal << " -> "<< minStress/downVal - 1. << std::endl ;
//     }
    return std::max(maxStress/upVal - 1., minStress/downVal - 1.) ;
}

double SpaceTimeNonLocalMohrCoulomb::grade( ElementState &s)
{
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;
    scoreAtTimeStepEnd = gradeAfter ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > 0)
    {
        return 1 ;
    }
    
    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 0 ;

    
    while(std::abs(upTime-downTime) > 1e-7)
    {
        double gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
            downTime = testTime ;
        else if(gradeTest > 0)
            upTime = testTime ;
        else
            return testTime ;
        
        testTime = 0.5*(downTime+upTime) ;
    }
    return 1.-(testTime*.5+.5) ;
}

double SpaceTimeNonLocalMohrCoulomb::gradeAtTime( ElementState &s, double t )
{

    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    Vector pstress( getSmoothedField(PRINCIPAL_REAL_STRESS_FIELD, s, t)) ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;

    metInTension = false ;
    metInCompression = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

    std::vector<double> scores ;
    scores.push_back(-1);
    if( maxStress >= upVal && maxStress > 0 )
    {
        metInTension = true;
        scores.push_back(1. - std::abs( upVal / maxStress ));
    }
    else if(maxStress > 0)
        scores.push_back(-1. + std::abs( maxStress / upVal ));

    if( minStress <= downVal && minStress < 0 )
    {
        metInCompression = true ;
        scores.push_back(1. - std::abs( downVal / minStress )) ;
    }
    else if(minStress < 0 )
    {
        scores.push_back(-1. + std::abs( minStress / downVal )) ;
    }
    std::sort(scores.begin(), scores.end()) ;

    return scores.back() ;

}

FractureCriterion *NonLocalMohrCoulomb::getCopy() const
{
    return new NonLocalMohrCoulomb( *this ) ;
}



NonLocalLinearlyDecreasingMohrCoulomb::NonLocalLinearlyDecreasingMohrCoulomb( double up, double down,double limittstrain, double limitcstrain, double E ) : upVal( up ), downVal( down ), stiffness(E),limittstrain(limittstrain),limitcstrain(limitcstrain)
{
    metInTension = false ;
    metInCompression = false ;
}


NonLocalLinearlyDecreasingMohrCoulomb::~NonLocalLinearlyDecreasingMohrCoulomb()
{

}

double NonLocalLinearlyDecreasingMohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->getDamageModel()->fractured() )
        return -1 ;

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_TOTAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;
    double minStrain = pstrain.min() ;
    
    metInTension = false ;
    metInCompression = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;
    double us = upVal/stiffness ;
    double ds = downVal/stiffness ;
    
    if(maxStrain < us)
        return std::abs(stiffness*maxStrain/upVal)-1. ;
    
    double tfactor = std::max(1.-std::abs((maxStrain-us)/(limittstrain-us)),0.) ;
    double cfactor = std::max(1.-std::abs((-minStrain+ds)/(-limitcstrain+ds)), 0.) ;

    double  upStress = tfactor*upVal +0.01e6;

    if(upStress < POINT_TOLERANCE)
    {
        if(maxStress < 1e-3*upVal)
            return -1 ;
        metInTension = true;
        metInCompression = false ;
        return  maxStress/upVal-1. ;
    }

    double  downStress = cfactor*downVal -0.01e6;
    if(downStress > -POINT_TOLERANCE)
    {
       if(minStress > 1e-3*downVal)
            return -1 ;
        metInTension = false;
        metInCompression = true ;
        return  std::abs(minStress/downVal)-1. ;
    }


    std::vector<double> scores ;
    scores.push_back(-1);
    if( maxStress > 0 )
    {
        metInTension = true;
        scores.push_back( maxStress / upStress-1.);
    }
    else if( minStress < 0 )
    {
        metInCompression = true ;
        scores.push_back( minStress / downStress-1.) ;
    }

    std::sort(scores.begin(), scores.end()) ;

    return scores.back() ;
}

FractureCriterion *NonLocalLinearlyDecreasingMohrCoulomb::getCopy() const
{
    return new NonLocalLinearlyDecreasingMohrCoulomb( *this ) ;
}

NonLocalExponentiallyDecreasingMohrCoulomb::NonLocalExponentiallyDecreasingMohrCoulomb( double up, double down,double limittstrain, double limitcstrain, double E ) :  upVal( up ), downVal( down ), stiffness(E),limittstrain(limittstrain),limitcstrain(limitcstrain)
{
    metInTension = false ;
    metInCompression = false ;
}


NonLocalExponentiallyDecreasingMohrCoulomb::~NonLocalExponentiallyDecreasingMohrCoulomb()
{

}

double NonLocalExponentiallyDecreasingMohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_TOTAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;

    metInTension = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;


    double tfactor = exp(-(maxStrain-upVal/stiffness)/(limittstrain-upVal/stiffness)) ;

    if(maxStrain <= upVal/stiffness)
        tfactor = 1 ;


    double  upStress = tfactor*upVal ;
    std::vector<double> scores ;
    scores.push_back(-1);

    if( maxStress >= upStress && maxStress > 0 )
    {
        metInTension = true;
        scores.push_back(std::abs( maxStress/upStress  )-1);
    }

    std::sort(scores.begin(), scores.end()) ;

    return scores.back() ;
}

FractureCriterion *NonLocalExponentiallyDecreasingMohrCoulomb::getCopy() const
{
    return new NonLocalExponentiallyDecreasingMohrCoulomb( *this ) ;
}



NonLocalInverseRootMohrCoulomb::NonLocalInverseRootMohrCoulomb(double limitstrain, double limitystrain, double E,double c) : stiffness(E)
    ,limitstrain(limitstrain),limitystrain(std::max(limitystrain,limitstrain)), c(c)
{
    metInTension = false ;
    upVal = E*limitstrain ;
}


NonLocalInverseRootMohrCoulomb::~NonLocalInverseRootMohrCoulomb()
{

}

double NonLocalInverseRootMohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_TOTAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
// 	double maxStress = pstress.max() ;
// 	double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;
//     double minStrain = pstrain.min() ;

    metInTension = false ;
    metInTension = std::abs( maxStrain / limitstrain ) > 1. ;

    double effectiveStiffness = stiffness ;
    if(s.getParent()->getBehaviour()->getDamageModel())
        effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

    double tfactor = (1+sqrt(c*(maxStrain-limitystrain))) ;
    if(maxStrain < limitystrain)
        tfactor = 1. ;

    double  upStrain = tfactor*upVal ;
    maxStrain *= effectiveStiffness ;
    std::vector<double> scores ;
    scores.push_back(-1);
    if( maxStrain <= upStrain )
    {
        metInTension = true;
        scores.push_back(1. - std::abs( upStrain / maxStrain ));
    }
    else if( upStrain > POINT_TOLERANCE)
        scores.push_back(-1. + std::abs( maxStrain / upStrain ));
    else if(maxStrain > 0)
        return POINT_TOLERANCE ;

    std::sort(scores.begin(), scores.end()) ;
    return scores.back() ;
}

FractureCriterion *NonLocalInverseRootMohrCoulomb::getCopy() const
{
    return new NonLocalInverseRootMohrCoulomb( *this ) ;
}

}
