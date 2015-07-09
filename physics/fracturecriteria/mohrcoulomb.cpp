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

MohrCoulomb::MohrCoulomb( double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
    , upVal( up ), downVal( down )
{
    metInTension = false ;
    metInCompression = false ;
}


MohrCoulomb::~MohrCoulomb()
{
}

double MohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->fractured() )
        return 0 ;

    Vector pstress(0., s.getParent()->spaceDimensions()) ;
    s.getField( PRINCIPAL_REAL_STRESS_FIELD, s.getParent()->getCenter(), pstress, false) ;

    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
    metInTension = false ;
    metInCompression = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

    if( maxStress >= upVal )
    {
        metInTension = true;

        if( minStress <= downVal )
            metInCompression = true ;

        return 1. - std::abs( upVal / maxStress ) ;
    }

    if( minStress <= downVal )
    {
        metInCompression = true ;
        return 1. - std::abs( downVal / minStress ) ;
    }

    double s0 = -1. + std::abs( maxStress / upVal );
    double s1 = -1. + std::abs( minStress / downVal ) ;

    if( minStress > 0 )
    {
        return s0 ;
    }

    if( maxStress < 0 )
    {
        return s1 ;
    }

    if( std::abs( s0 ) > std::abs( s1 ) )
        return s0 ;

    return s1;
}

FractureCriterion *MohrCoulomb::getCopy() const
{
    return new MohrCoulomb( *this ) ;
}

NonLocalMohrCoulomb::NonLocalMohrCoulomb( double up, double down, double E, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
    , upVal( up ), downVal( down ), stiffness(E)
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

    Vector stress =  getSmoothedField(PRINCIPAL_REAL_STRESS_FIELD, s) ;
    double score = stress.max() / upVal - 1.;
    if(std::abs(stress.max()) < std::abs(stress.min()))
        score = std::abs(stress.min() / downVal) - 1.; ;
    if(score > 0)
        metInTension = true ;
    else
        metInTension = false ;
    
    return score ;

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

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s, t)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;
//     double maxStrain = pstrain.max() ;
//     double minStrain = pstrain.min() ;

//  std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
    metInTension = false ;
    metInCompression = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

//  double effectiveStiffness = stiffness ;
//  if(s.getParent()->getBehaviour()->getDamageModel())
//      effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

    double  upStrain = upVal ;///effectiveStiffness ;
    double  downStrain = downVal ;///effectiveStiffness ;
    std::vector<double> scores ;
    scores.push_back(-1);
    if( maxStress >= upStrain*stiffness && maxStress > 0 )
    {
        metInTension = true;
        scores.push_back(1. - std::abs( upStrain*stiffness / maxStress ));
    }
    else if(maxStress > 0)
        scores.push_back(-1. + std::abs( maxStress / (upStrain*stiffness) ));

    if( minStress <= downStrain*stiffness && minStress < 0 )
    {
        metInCompression = true ;
        scores.push_back(1. - std::abs( downStrain*stiffness / minStress )) ;
    }
    else if(minStress < 0 )
    {
        scores.push_back(-1. + std::abs( minStress / (downStrain*stiffness) )) ;
    }
    std::sort(scores.begin(), scores.end()) ;

    return scores.back() ;

}

FractureCriterion *NonLocalMohrCoulomb::getCopy() const
{
    return new NonLocalMohrCoulomb( *this ) ;
}



NonLocalLinearlyDecreasingMohrCoulomb::NonLocalLinearlyDecreasingMohrCoulomb( double up, double down,double limittstrain, double limitcstrain, double E, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
    , upVal( up ), downVal( down ), stiffness(E),limittstrain(limittstrain),limitcstrain(limitcstrain)
{
    metInTension = false ;
    metInCompression = false ;
}


NonLocalLinearlyDecreasingMohrCoulomb::~NonLocalLinearlyDecreasingMohrCoulomb()
{

}

double NonLocalLinearlyDecreasingMohrCoulomb::grade( ElementState &s )
{

    if( s.getParent()->getBehaviour()->fractured() )
        return -1 ;

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;
    double minStrain = pstrain.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
    metInTension = false ;
    metInCompression = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

//     double d = s.getParent()->getBehaviour()->getDamageModel()->getState().max() ;
//     if(1.-d < POINT_TOLERANCE)
//     {
//         double c0 = std::abs(maxStress/upVal) ;
//         double c1 = std::abs(minStress/downVal) ;
//         return std::max(c0, c1) ;
//     }
// 
//     double effMaxStress = std::max(limittstrain/(limittstrain-upVal/stiffness-1./(stiffness*(1.-d))), 0.) ;
//     double effMinStress = std::max(std::abs(limitcstrain)/(std::abs(limitcstrain)-std::abs(downVal)/stiffness-1./(stiffness*(1.-d))), 0.) ;
// 
//     double c0 = std::abs(maxStress/upVal)-effMaxStress ;
//     double c1 = std::abs(minStress/downVal)-effMinStress ;
//     return std::max(c0, c1) ;

    double effectiveStiffness = stiffness ;
    if(s.getParent()->getBehaviour()->getDamageModel())
        effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

    double tfactor = 1.-(maxStrain-upVal/stiffness)/(limittstrain-upVal/stiffness) ;
// 	if(maxStrain > limittstrain)
// 		return POINT_TOLERANCE ;
// 	if(maxStrain <= upVal/stiffness)
// 		tfactor = 1 ;

    double cfactor = 1.-(-minStrain+downVal/stiffness)/(-limitcstrain+downVal/stiffness) ;
// 	if(minStrain < limitcstrain)
// 		return POINT_TOLERANCE ;
// 	if(minStrain > downVal/stiffness)
// 		cfactor = 1 ;
//
    double  upStrain = tfactor*upVal/effectiveStiffness ;
    double  downStrain = cfactor*downVal/effectiveStiffness ;
    std::vector<double> scores ;
    scores.push_back(-1);
    if( maxStrain >= upStrain && maxStrain > 0 )
    {
        metInTension = true;
        scores.push_back(std::abs( maxStrain / upStrain ));
    }
    else if(maxStrain > 0 && upStrain > POINT_TOLERANCE)
        scores.push_back(-1. + std::abs( maxStrain / upStrain ));
    else if(maxStrain > 0)
        return POINT_TOLERANCE ;

    if( minStrain <= downStrain && minStrain < 0 )
    {
        metInCompression = true ;
        scores.push_back(std::abs( minStrain / downStrain )) ;
    }
    else if(minStrain < 0  && downStrain < -POINT_TOLERANCE)
        scores.push_back(-1. + std::abs( minStrain / downStrain )) ;
    else if(minStrain < 0)
        return POINT_TOLERANCE ;

    std::sort(scores.begin(), scores.end()) ;
// 	if (scores.back() > .99)
// 	{
// 		std::cout << maxStrain << "  " << limittstrain << std::endl ;
// 		std::cout << scores.back() << std::endl ;
// 		exit(0) ;
// 	}
    return scores.back() ;
}

FractureCriterion *NonLocalLinearlyDecreasingMohrCoulomb::getCopy() const
{
    return new NonLocalLinearlyDecreasingMohrCoulomb( *this ) ;
}

NonLocalExponentiallyDecreasingMohrCoulomb::NonLocalExponentiallyDecreasingMohrCoulomb( double up, double down,double limittstrain, double limitcstrain, double E, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
    , upVal( up ), downVal( down ), stiffness(E),limittstrain(limittstrain),limitcstrain(limitcstrain)
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

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
    double maxStress = pstress.max() ;
    double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;
//     double minStrain = pstrain.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
    metInTension = false ;
    metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
    metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

    double effectiveStiffness = stiffness ;
    if(s.getParent()->getBehaviour()->getDamageModel())
        effectiveStiffness = stiffness*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

    double tfactor = exp(-(maxStrain-upVal/stiffness)/(limittstrain-upVal/stiffness)) ;
// 	if(maxStrain > limittstrain)
// 		return POINT_TOLERANCE ;
// 	if(tfactor < 1e-1)
// 		tfactor = 1e-1 ;
    if(maxStrain <= upVal/stiffness)
        tfactor = 1 ;

//     double cfactor = exp(-(-minStrain+downVal/stiffness)/(-limitcstrain+downVal/stiffness)) ;
// 	if(minStrain < limitcstrain)
// 		return POINT_TOLERANCE ;
//     if(minStrain > downVal/stiffness)
//         cfactor = 1 ;
//
    double  upStrain = tfactor*upVal/effectiveStiffness ;
//     double  downStrain = cfactor*downVal/effectiveStiffness ;
    std::vector<double> scores ;
    scores.push_back(-1);

    if( maxStrain >= upStrain && maxStrain > 0 )
    {
        metInTension = true;
        scores.push_back(1. - std::abs( upStrain / maxStrain ));
    }
    else if(maxStrain > 0 && upStrain > POINT_TOLERANCE)
        scores.push_back(-1. + std::abs( maxStrain / upStrain ));
    else if(maxStrain > 0)
        return -POINT_TOLERANCE ;

// 	if( minStrain <= downStrain && minStrain < 0 )
// 	{
// 		metInCompression = true ;
// 		scores.push_back(1. - std::abs( downStrain / minStrain )) ;
// 	}
// 	else if(minStrain < 0  && downStrain < -POINT_TOLERANCE)
// 		scores.push_back(-1. + std::abs( minStrain / downStrain )) ;
// 	else if(minStrain < 0)
// 		return POINT_TOLERANCE ;

    std::sort(scores.begin(), scores.end()) ;
// 	if (scores.back() > .99)
// 	{
// 		std::cout << maxStrain << "  " << limittstrain << std::endl ;
// 		std::cout << scores.back() << std::endl ;
// 		exit(0) ;
// 	}
    return scores.back() ;
}

FractureCriterion *NonLocalExponentiallyDecreasingMohrCoulomb::getCopy() const
{
    return new NonLocalExponentiallyDecreasingMohrCoulomb( *this ) ;
}



NonLocalInverseRootMohrCoulomb::NonLocalInverseRootMohrCoulomb(double limitstrain, double limitystrain, double E,double c, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z ), stiffness(E)
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

    std::pair<Vector, Vector> pstressStrain( getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s)) ;
    Vector pstress = pstressStrain.first ;
    Vector pstrain = pstressStrain.second ;
// 	double maxStress = pstress.max() ;
// 	double minStress = pstress.min() ;
    double maxStrain = pstrain.max() ;
//     double minStrain = pstrain.min() ;

// 	std::cout << maxStress << ", " << maxStrain << std::endl ;
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
