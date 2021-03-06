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
#include "maxstrain.h"
#include "../damagemodels/damagemodel.h"
#include "../damagemodels/spacetimebifurcation.h"
#include <fstream>

namespace Amie {

MaximumStrain::MaximumStrain(double up) :  upVal(up)
{
	metInCompression = false ;
	metInTension = false ;
}


MaximumStrain::~MaximumStrain()
{
}

double MaximumStrain::grade(ElementState &s)
{
	Vector pstrain(0., s.getParent()->getBoundingPoints().size()*(3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL))) ;
	s.getField( MECHANICAL_STRAIN_FIELD, s.getParent()->getBoundingPoints(), pstrain, false) ;
	double maxStrain = pstrain.max();
	metInCompression = false ;
	metInTension = false ;
	if(maxStrain > upVal)
	{
		metInTension = true ;
		return 1.-std::abs(upVal/maxStrain) ;
	}
	else
		return -1.+ std::abs(maxStrain/upVal);
	
}

FractureCriterion * MaximumStrain::getCopy() const
{
	return new MaximumStrain(*this) ;
}

double SpaceTimeNonLocalMaximumStrain::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;


	Vector stateBefore( getSmoothedField( PRINCIPAL_TOTAL_STRAIN_FIELD, s, -1) ) ;
	Vector stateAfter( getSmoothedField( PRINCIPAL_TOTAL_STRAIN_FIELD, s, 1) ) ;
	double maxStrainAfter = stateAfter.max() ;
	double maxStrainBefore = stateBefore.max() ;

	metInCompression = false ;
	metInTension = false ;

	if(maxStrainAfter > upVal)
	{
		metInTension = true ;
		scoreAtTimeStepEnd = 1.-std::abs(upVal/maxStrainAfter) ;
		return  std::min(1.,1. - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore)) ;
	}
	scoreAtTimeStepEnd = -1.+ std::abs(maxStrainAfter/upVal);
	return -1.+ maxStrainAfter/upVal;
	
	
}

double SpaceTimeLimitFirstStrainInvariant::grade(ElementState &s)
{
    if(s.getParent()->getBehaviour()->getDamageModel()->fractured(0))
    {
        scoreAtTimeStepEnd = -1 ;
        return -1 ;
    }

    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;

//    std::cout << gradeBefore << "\t" << gradeAfter << std::endl ;

    scoreAtTimeStepEnd = gradeAfter ;

    if(gradeBefore > 0)
    {
        return 1 ;
    }
    if(gradeAfter < 0)
    {
        return gradeAfter ;
    }

    double upTime = 1. ;
    double downTime = -1 ;
    double testTime = (upTime+downTime)*0.5 ;
    double gradeTest = 0 ;
    while(std::abs(upTime-downTime) > 1e-4)
    {
        gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
        {
            if(gradeBefore < 0)
                downTime = testTime ;
            else
                upTime = testTime ;
        }
        else if(gradeTest > 0)
        {
            if(gradeBefore < 0)
                upTime = testTime ;
            else
                downTime = testTime ;
        }
        else
        {
//            std::cout << gradeBefore << "\t" << gradeAfter << std::endl ;
            return 1.-(testTime*.5+.5) ;
        }
        testTime = (downTime+upTime)*0.5 ;
    }


//    std::cout << gradeBefore << "\t" << gradeAfter << std::endl ;
    return 1.-(testTime*.5+.5) ;
}



double SpaceTimeLimitFirstStrainInvariant::gradeAtTime(ElementState &s, double t)
{
        Vector strain = getSmoothedField( MECHANICAL_STRAIN_FIELD, s, t ) ;
        double trace = strain[0]+strain[1]+(strain.size() == 6 ? strain[2] : 0 ) ;
        if(upVal > 0)
        {
            return -(upVal-trace)/upVal ;
        }
        else
            return -(upVal-trace)/upVal ;
}

double SpaceTimeNonLocalLinearSofteningMaximumStrain::grade(ElementState &s)
{    
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;
    scoreAtTimeStepEnd = gradeAfter ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > -POINT_TOLERANCE)
        return 1 ;
    
    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 2.*(-gradeBefore)/(gradeAfter-gradeBefore)-1. ;
    double gradeDown = gradeBefore ;
    double gradeUp = gradeAfter ;
    double gradeTest = 0 ;

    while(std::abs(upTime-downTime) > 1e-4 )
    {
        gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
        {
            downTime = testTime ;
            gradeDown = gradeTest ;
        }
        else if(gradeTest > 0)
        {
            upTime = testTime ;
            gradeUp = gradeTest ;
        }
        else
            return 1.-(testTime*.5+.5) ;
        
	if(gradeDown > POINT_TOLERANCE)
	        testTime = downTime + (upTime-downTime)*(-gradeDown)/(gradeUp-gradeDown) ;
	else
		testTime = (downTime + upTime)*0.5 ;
	if(std::abs(testTime - upTime) < POINT_TOLERANCE && std::abs(gradeTest) < POINT_TOLERANCE) 
		return 1.-(testTime*.5+.5) ;
    }

    return 1.-(testTime*.5+.5) ;
}



double SpaceTimeNonLocalLinearSofteningMaximumStrain::gradeAtTime(ElementState &s, double t)
{
	if( s.getParent()->getBehaviour()->fractured() )
	{
		return -1 ;
	}

	std::pair<Vector, Vector> currentState = getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, t ) ;
	double Esoft = maxstress / ( yieldstrain - upVal) ;
	double Einst = maxstress/upVal * (1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;

	double epsMax = std::max(upVal, yieldstrain*Esoft/(Esoft+Einst) ) ;

// ;
// 	double Eprev = stateBefore.first.max() / stateBefore.second.max() ;

	double currentStrain = currentState.second.max() ;

	metInCompression = false ;
	metInTension = false ;
	if(currentStrain > epsMax)
	{
		metInTension = true ;
		return  std::min(1.,1. - epsMax / currentStrain ) ;
	}
	else if (currentStrain > 0)
		return -1.+ currentStrain/epsMax;
	return -1. ;

}

double SpaceTimeNonLocalMaximumStress::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, s, 1) ) ;
	double maxStressAfter = stateAfter.first.max() ;
	double maxStressBefore = stateBefore.first.max() ;

	metInCompression = false ;
	metInTension = false ;
	if(maxStressAfter > maxstress)
	{
		metInTension = true ;
		scoreAtTimeStepEnd = 1.-std::abs(maxstress/maxStressAfter) ;
		return std::min(1., 1. - (maxstress - maxStressBefore) / (maxStressAfter - maxStressBefore)) ;
	}
	scoreAtTimeStepEnd = -1.+ maxStressAfter/maxstress ;
	return -1.+ maxStressAfter/maxstress;
	
	
}



SpaceTimeNonLocalEllipsoidalMixedCriterion::SpaceTimeNonLocalEllipsoidalMixedCriterion(double up, double mstr, double E0, double Einf) : MaximumStrain(up),maxstress(mstr), E_inst(E0), E_relaxed(Einf)
{
	double strainrelaxed = maxstress/E_relaxed ;
	double stressinst = upVal*E_inst ;

	renormStrain=1 ;
	renormStress=1 ;

	Point center(strainrelaxed*renormStrain, stressinst*renormStress) ;
	Point mj(std::abs(upVal - strainrelaxed)*renormStrain, 0.) ;
	Point mn(0.,std::abs(stressinst - maxstress)*renormStress) ;

	surface = new Ellipse(center, mj, mn) ;
	surface->sampleBoundingSurface(512) ;
	base = std::abs(VirtualMachine().eval( getEllipseFormFunction(surface) )) ;
	upVal = surface->getCenter().getX()-surface->getMajorRadius() ;

}

double SpaceTimeNonLocalEllipsoidalMixedCriterion::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	metInCompression = false ;
	metInTension = false ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( REAL_STRESS_FIELD, MECHANICAL_STRAIN_FIELD, s, 1) ) ;

	Point before( stateBefore.second.max()*renormStrain/surface->getMinorRadius(), stateBefore.first.max()*renormStress/surface->getMajorRadius()) ;
	Point after( stateAfter.second.max()*renormStrain/surface->getMinorRadius() , stateAfter.first.max()*renormStress/surface->getMajorRadius()) ;
	
	
	if(stateAfter.second.max() <= POINT_TOLERANCE || stateAfter.first.max() <= POINT_TOLERANCE)
		return -1 ;
	
	if(after.getX() > surface->getCenter().getX()/surface->getMinorRadius()+1)
	{
		double maxStressAfter = stateAfter.first.max() ;
		double maxStressBefore = stateBefore.first.max() ;

		metInCompression = false ;
		metInTension = false ;
		if(maxStressAfter > maxstress)
		{
			metInTension = true ;
			
			return std::min(1., 1. - (maxstress - maxStressBefore) / (maxStressAfter - maxStressBefore)) ;
		}
		return -1.+ maxStressAfter/maxstress;
	}
	
	if(after.getX() < surface->getCenter().getX()/surface->getMinorRadius()-1  || after.getY() >= surface->getCenter().getY()/surface->getMajorRadius())
	{
		double maxStrainAfter = stateAfter.second.max() ;
		double maxStrainBefore = stateBefore.second.max() ;
		
		metInCompression = false ;
		metInTension = false ;
		if(maxStrainAfter > upVal)
		{
			metInTension = true ;
			double score = 1. - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore) ;
			if(score > 1)
			{
				double maxStressAfter = stateAfter.first.max() ;
				double maxStressBefore = stateBefore.first.max() ;

				metInCompression = false ;
				metInTension = false ;
				if(maxStressAfter > maxstress)
				{
					metInTension = true ;
					
					return std::min(1., 1. - (maxstress - maxStressBefore) / (maxStressAfter - maxStressBefore)) ;
				}
				return -1.+ maxStressAfter/maxstress;
			}
			return  1. - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore) ;
		}
		return -1.+ maxStrainAfter/upVal;
	}
	
	Segment evolution(before,after) ;
	Circle renormsurf(1.,surface->getCenter().getX()/surface->getMinorRadius(), surface->getCenter().getY()/surface->getMajorRadius()) ;
	
	metInTension = true ;

	std::vector<Point> inter = evolution.intersection( &renormsurf ) ;

	if(inter.size() == 0)
	{
		Point proj(after) ;
		renormsurf.project(&proj);
		return std::max(-1.,-squareDist2D(proj,after) );
	}

	metInTension = true ;
	int i = 0 ;
	if(inter.size() == 2 && squareDist2D( before, inter[1] ) < squareDist2D( before, inter[0] ) )
		i = 1 ;

	// should not happen, but just in case...
	if(!evolution.on(inter[i]))
	{
		Point proj(after) ;
		renormsurf.project(&proj);
		return std::max(-1.,-squareDist2D(proj,after) );
	}

	return std::min(1., std::sqrt(squareDist2D( inter[i], after )/squareDist2D( before, after ) )) ;

}

}

