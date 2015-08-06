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
#include "spacetimemultilinearsofteningfracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include <fstream>

namespace Amie {

SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string file, double E_, double e, double s_, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(E_, mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s_)
{
    std::fstream input(file) ;
    if(!input.is_open())
    {
        std::cout << "file " << file << " doesn't exist!" << std::endl ;
        exit(0) ;
    }

    double buffer ;
    std::vector<Point> p ;

    do {
        input >> buffer ;
        double x = buffer*renormStrain ;
        input >> buffer ;
        double y = buffer*renormStress ;
	p.push_back( Point(x,y) ) ;
    } while(!input.eof()) ;
    p.pop_back() ;


	double e1 = p[0].getX() ;
	upVal = e1 ;
//	double e2 = p[1].getX() ;
	double s1 = p[0].getY() ;
//	double s2 = p[1].getY() ;

	double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
	std::valarray<Point *> curve ;
	curve.resize( p.size()+1 ) ;
	Point * origin = new Point( 0., stressAtOrigin) ;
	curve[0] = origin ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		curve[i+1] = &p[i] ;
	}
//	double stressAtInfinity = curve[p.size()]->getY() ;
	Point horizontal(1.,0.) ;

	stressStrainCurve = new SegmentedLine( curve ) ;
	asymptote = new Line( *curve[p.size()], horizontal ) ;

}

SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & p, double E_, double e, double s_, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( p[0].getX(), mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s_)
{
	double e1 = p[0].getX()*renormStrain ;
	upVal = e1 ;
//	double e2 = p[1].getX()*renormStrain ;
	double s1 = p[0].getY()*renormStress ;
//	double s2 = p[1].getY()*renormStress ;

	double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
	std::valarray<Point *> curve ;
	curve.resize( p.size()+1 ) ;
	Point * origin = new Point( 0., stressAtOrigin) ;
	curve[0] = origin ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		curve[i+1] = new Point(p[i].getX()*renormStrain, p[i].getY()*renormStress) ;
	}
//	double stressAtInfinity = curve[p.size()]->getY()*renormStress ;
	Point horizontal(1.,0.) ;

	stressStrainCurve = new SegmentedLine( curve ) ;
	asymptote = new Line( *curve[p.size()], horizontal ) ;
}

SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( SegmentedLine * s, Line * l, double E_, double e, double s_, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(E_, mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s_), stressStrainCurve(s), asymptote(l)
{

}

FractureCriterion * SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getCopy() const 
{
	std::vector<Point> p ;
	for(size_t i = 1 ; i < stressStrainCurve->size() ; i++)
	{
		Point q( stressStrainCurve->getPoint(i).getX()/renormStrain, stressStrainCurve->getPoint(i).getY()/renormStress ) ;
		p.push_back(q) ;
	}
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion * ret = new SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( p, E, renormStrain, renormStress ) ;
	ret->copyEssentialParameters( this ) ;
	return ret ;
}

double SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::grade(ElementState &s)
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

    while(std::abs(upTime-downTime) > 1e-6)
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


double SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::gradeAtTime(ElementState &s, double t)  
{
	std::pair<Vector, Vector> currentState( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, t) ) ;
	Point current( currentState.second.max()*renormStrain, currentState.first.max()*renormStress ) ;
	Point inter ;
	bool found = false ;

	Line direction( Point(0,0), current ) ;
	if(direction.intersects( stressStrainCurve ) )
	{
		std::vector<Point> all = direction.intersection( stressStrainCurve ) ;
		if(all.size() > 0)
		{
			found = true ;
			inter = all[0] ;
			for(size_t i = 1 ; i < all.size() ; i++)
			{
				if( dist(all[i], current ) < dist(inter, current) )
					inter = all[i] ;
			}
		}
	}
	if(!found && direction.intersects( *asymptote ))
	{
		found = true ;
		inter = direction.intersection( *asymptote ) ;
	}

	Segment history( Point(0,0), current ) ;
	if( found )
	{
		if( history.on( inter ) )
			return std::min(1., 1.-inter.getX()/current.getX()) ;
		else
			return std::max(-1., -1.+current.getX()/inter.getX()) ;
	}
	
	if( current.getY() < stressStrainCurve->getPoint(0).getY() )
		return std::max( -1., -1.+current.getY()/stressStrainCurve->getPoint(0).getY() ) ;
	else
		return std::min( 1., 1.-stressStrainCurve->getPoint(0).getY()/current.getY() ) ;

}

void SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumStress(double s, double smax) 
{
	double targetStress = s ;
	double targetMaxStress = smax ;

	double currentStress = stressStrainCurve->getPoint(1).getY()/renormStress ;
	double currentMaxStress = stressStrainCurve->getPoint(stressStrainCurve->size()-1).getY()/renormStress ;

	if(s < 0)
		targetStress = currentStress ;
	if(smax < 0)
		targetMaxStress = currentMaxStress ;

	double factor = 0.; 
	if(std::abs(currentStress-currentMaxStress) > POINT_TOLERANCE)
		factor = (targetStress-targetMaxStress)/(currentStress-currentMaxStress) ;

	for(size_t i = 0 ; i < stressStrainCurve->size() ; i++)
	{
		double oldStress = stressStrainCurve->getPoint(i).getY()/renormStress ;
		double newStress = targetMaxStress + (oldStress-currentMaxStress)*factor ;
		stressStrainCurve->getPoint(i).setY( newStress*renormStress ) ;
	}
	Point p( asymptote->origin().getX(), stressStrainCurve->getPoint(stressStrainCurve->size()-1).getY()) ;
	asymptote->setOrigin(p) ;
}

void SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumStrain(double e, double emax) 
{
	double targetStrain = e ;
	double targetMaxStrain = emax ;

	double currentStrain = stressStrainCurve->getPoint(1).getX()/renormStrain ;
	double currentMaxStrain = stressStrainCurve->getPoint(stressStrainCurve->size()-1).getX()/renormStrain ;

	if(e < 0)
		targetStrain = currentStrain ;
	if(emax < 0)
		targetMaxStrain = currentMaxStrain ;

	double factor = 0.; 
	if(std::abs(currentStrain-currentMaxStrain) > POINT_TOLERANCE)
		factor = (targetStrain-targetMaxStrain)/(currentStrain-currentMaxStrain) ;

	for(size_t i = 0 ; i < stressStrainCurve->size() ; i++)
	{
		double oldStrain = stressStrainCurve->getPoint(i).getX()/renormStrain ;
		double newStrain = targetMaxStrain + (oldStrain-currentMaxStrain)*factor ;
		stressStrainCurve->getPoint(i).setX( newStrain*renormStrain ) ;
	}
}

AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string ftension, std::string fcompression, double E_, double e, double s, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(E_, mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s)
{
	tensileStressStrainCurve = nullptr ;
	compressiveStressStrainCurve = nullptr ;
	tensileAsymptote = nullptr ;
	compressiveAsymptote = nullptr ;

	if(ftension.size() > 0)
	{
	    std::fstream input(ftension) ;
	    if(!input.is_open())
	    {
		std::cout << "file " << ftension << " doesn't exist!" << std::endl ;
		exit(0) ;
	    }

	    double buffer ;
	    std::vector<Point> p ;

	    do {
		input >> buffer ;
		double x = buffer*renormStrain ;
		input >> buffer ;
		double y = buffer*renormStress ;
		p.push_back( Point(x,y) ) ;
	    } while(!input.eof()) ;
	    p.pop_back() ;


		double e1 = p[0].getX() ;
		upVal = e1 ;
//		double e2 = p[1].getX() ;
		double s1 = p[0].getY() ;
//		double s2 = p[1].getY() ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( p.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			curve[i+1] = &p[i] ;
		}
//		double stressAtInfinity = curve[p.size()]->getY() ;
		Point horizontal(1.,0.) ;

		tensileStressStrainCurve = new SegmentedLine( curve ) ;
		tensileAsymptote = new Line( *curve[p.size()], horizontal ) ;
	}

	if(fcompression.size() > 0)
	{
	    std::fstream input(fcompression) ;
	    if(!input.is_open())
	    {
		std::cout << "file " << fcompression << " doesn't exist!" << std::endl ;
		exit(0) ;
	    }

	    double buffer ;
	    std::vector<Point> p ;

	    do {
		input >> buffer ;
		double x = buffer*renormStrain ;
		input >> buffer ;
		double y = buffer*renormStress ;
		p.push_back( Point(x,y) ) ;
	    } while(!input.eof()) ;
	    p.pop_back() ;


//		double e1 = p[0].getX() ;
//		double e2 = p[1].getX() ;
		double s1 = p[0].getY() ;
//		double s2 = p[1].getY() ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( p.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			curve[i+1] = &p[i] ;
		}
//		double stressAtInfinity = curve[p.size()]->getY() ;
		Point horizontal(1.,0.) ;

		compressiveStressStrainCurve = new SegmentedLine( curve ) ;
		compressiveAsymptote = new Line( *curve[p.size()], horizontal ) ;
	}
}

AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & ptension, const std::vector<Point> & pcompression, double E_, double e, double s, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s)
{
	tensileStressStrainCurve = nullptr ;
	compressiveStressStrainCurve = nullptr ;
	tensileAsymptote = nullptr ;
	compressiveAsymptote = nullptr ;

	if(ptension.size() > 0)
	{
		double e1 = ptension[0].getX()*renormStrain ;
		upVal = e1 ;
//		double e2 = ptension[1].getX()*renormStrain ;
		double s1 = ptension[0].getY()*renormStress ;
//		double s2 = ptension[1].getY()*renormStress ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( ptension.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < ptension.size() ; i++)
		{
			curve[i+1] = new Point(ptension[i].getX()*renormStrain, ptension[i].getY()*renormStress) ;
		}
//		double stressAtInfinity = curve[ptension.size()]->getY()*renormStress ;
		Point horizontal(1.,0.) ;

		tensileStressStrainCurve = new SegmentedLine( curve ) ;
		tensileAsymptote = new Line( *curve[ptension.size()], horizontal ) ;
	}
	if(pcompression.size() > 0)
	{
//		double e1 = pcompression[0].getX()*renormStrain ;
//		double e2 = pcompression[1].getX()*renormStrain ;
		double s1 = pcompression[0].getY()*renormStress ;
//		double s2 = pcompression[1].getY()*renormStress ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( pcompression.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < pcompression.size() ; i++)
		{
			curve[i+1] = new Point(pcompression[i].getX()*renormStrain, pcompression[i].getY()*renormStress) ;
		}
//		double stressAtInfinity = curve[pcompression.size()]->getY()*renormStress ;
		Point horizontal(1.,0.) ;

		compressiveStressStrainCurve = new SegmentedLine( curve ) ;
		compressiveAsymptote = new Line( *curve[pcompression.size()], horizontal ) ;
	}
}

FractureCriterion * AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getCopy() const 
{
	std::vector<Point> ptension ;
	if(tensileStressStrainCurve)
	{
		for(size_t i = 1 ; i < tensileStressStrainCurve->size() ; i++)
		{
			Point q( tensileStressStrainCurve->getPoint(i).getX()/renormStrain, tensileStressStrainCurve->getPoint(i).getY()/renormStress ) ;
			ptension.push_back(q) ;
		}
	}
	std::vector<Point> pcompression ;
	if(compressiveStressStrainCurve)
	{
		for(size_t i = 1 ; i < compressiveStressStrainCurve->size() ; i++)
		{
			Point q( compressiveStressStrainCurve->getPoint(i).getX()/renormStrain, compressiveStressStrainCurve->getPoint(i).getY()/renormStress ) ;
			pcompression.push_back(q) ;
		}
	}
	AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion * ret = new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( ptension, pcompression, E, renormStrain, renormStress ) ;
	ret->copyEssentialParameters( this ) ;
	return ret ;
}

AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::~AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion() 
{ 
	if(tensileStressStrainCurve)
		delete tensileStressStrainCurve ; 
	if(tensileAsymptote)
		delete tensileAsymptote ; 
	if(compressiveStressStrainCurve)
		delete compressiveStressStrainCurve ; 
	if(compressiveAsymptote)
		delete compressiveAsymptote ; 
} 

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::grade(ElementState &s)
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
    double testTime = 2.*(-gradeBefore)/(gradeAfter-gradeBefore)-1. ;
    double gradeDown = gradeBefore ;
    double gradeUp = gradeAfter ;
    double gradeTest = 0 ;

    while(std::abs(upTime-downTime) > 1e-6 )
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
        
        if(gradeUp > POINT_TOLERANCE)
            testTime = downTime + (upTime-downTime)*(-gradeDown)/(gradeUp-gradeDown) ;
        else
            testTime = (downTime + upTime)*0.5 ;
    }

    return 1.-(testTime*.5+.5) ;
}

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::gradeAtTime(ElementState &s, double t)  
{
	std::pair<Vector, Vector> currentState( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, t) ) ;
	Point current( currentState.second.max()*renormStrain, currentState.first.max()*renormStress ) ;
	Point inter ;
	bool found = false ;

	double tension = -1. ;
	double compression = -1. ;

	Line direction( Point(0,0), current ) ;
	double cas = 0 ;

	if(tensileStressStrainCurve)
	{
		if(direction.intersects( tensileStressStrainCurve ) )
		{
			std::vector<Point> all = direction.intersection( tensileStressStrainCurve ) ;
			if(all.size() > 0)
			{
				found = true ;
				inter = all[0] ;
				for(size_t i = 1 ; i < all.size() ; i++)
				{
					if( dist(all[i], current ) < dist(inter, current) )
						inter = all[i] ;
				}
			}
		}
		if(!found && direction.intersects( *tensileAsymptote ))
		{
			inter = direction.intersection( *tensileAsymptote ) ;
			if( inter.getX() > tensileAsymptote->origin().getX() )
				found = true ;
		}
	
		Segment history( Point(0,0), current ) ;
		if( found )
		{
			if( history.on( inter ) )
			{
				tension = std::min(1., 1.-inter.getX()/current.getX()) ;
				cas = 1 ;
			}
			else
			{
				tension = std::max(-1., -1.+current.getX()/inter.getX()) ;
				cas = 2 ;
			}
		}
                else
                {
		
			if( current.getY() < tensileStressStrainCurve->getPoint(0).getY() )
			{
				tension = std::max( -1., -1.+current.getY()/tensileStressStrainCurve->getPoint(0).getY() ) ;
				cas = 3 ;
			}
			else
			{
				tension = std::min( 1., 1.-tensileStressStrainCurve->getPoint(0).getY()/current.getY() ) ;
				cas = 4 ;
			}
                }
	}

/*	if(t == -1 && tension > 0)
	{
		std::cout << currentState.first.min() << "\t" << currentState.first.max() << "\t"  << currentState.second.min() << "\t" << currentState.second.max() << std::endl ;
		current.print() ;
	}*/


	current = Point( currentState.second.min()*renormStrain, currentState.first.min()*renormStress ) ;
	direction.setVector( current ) ;
	found = false ;
	if(compressiveStressStrainCurve)
	{
		if(direction.intersects( compressiveStressStrainCurve ) )
		{
			std::vector<Point> all = direction.intersection( compressiveStressStrainCurve ) ;
			if(all.size() > 0)
			{
				found = true ;
				inter = all[0] ;
				for(size_t i = 1 ; i < all.size() ; i++)
				{
					if( dist(all[i], current ) < dist(inter, current) )
						inter = all[i] ;
				}
			}
		}
		if(!found && direction.intersects( *compressiveAsymptote ))
		{
			found = true ;
			inter = direction.intersection( *compressiveAsymptote ) ;
		}
	
		Segment history( Point(0,0), current ) ;
		if( found )
		{
			if( history.on( inter ) )
				compression = std::min(1., 1.-inter.getX()/current.getX()) ;
			
			else
				compression = std::max(-1., -1.+current.getX()/inter.getX()) ;
		}
		else
		{
		
			if( current.getY() > compressiveStressStrainCurve->getPoint(0).getY() )
				compression = std::max( -1., -1.+current.getY()/compressiveStressStrainCurve->getPoint(0).getY() ) ;
			else
				compression = std::min( 1., 1.-compressiveStressStrainCurve->getPoint(0).getY()/current.getY() ) ;
		}
	}
	

	return std::max( compression, tension ) ;

}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumTensileStress(double s, double smax) 
{
	if(!tensileStressStrainCurve || (s<0 && smax<0))
		return ;


	double targetStress = s ;
	double targetMaxStress = smax ;

	double currentStress = tensileStressStrainCurve->getPoint(1).getY()/renormStress ;
	double currentMaxStress = tensileStressStrainCurve->getPoint(tensileStressStrainCurve->size()-1).getY()/renormStress ;

	if(s < 0)
		targetStress = currentStress ;
	if(smax < 0)
		targetMaxStress = currentMaxStress ;

	double factor = 0.; 
	if(std::abs(currentStress-currentMaxStress) > POINT_TOLERANCE)
		factor = (targetStress-targetMaxStress)/(currentStress-currentMaxStress) ;

	for(size_t i = 0 ; i < tensileStressStrainCurve->size() ; i++)
	{
		double oldStress = tensileStressStrainCurve->getPoint(i).getY()/renormStress ;
		double newStress = targetMaxStress + (oldStress-currentMaxStress)*factor ;
		tensileStressStrainCurve->getPoint(i).setY( newStress*renormStress ) ;
	}
	Point p( tensileAsymptote->origin().getX(), tensileStressStrainCurve->getPoint(tensileStressStrainCurve->size()-1).getY()) ;
	tensileAsymptote->setOrigin(p) ;
}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumTensileStrain(double e, double emax) 
{
	if(!tensileStressStrainCurve || (e<0 && emax<0))
		return ;

	double targetStrain = e ;
	double targetMaxStrain = emax ;

	double currentStrain = tensileStressStrainCurve->getPoint(1).getX()/renormStrain ;
	double currentMaxStrain = tensileStressStrainCurve->getPoint(tensileStressStrainCurve->size()-1).getX()/renormStrain ;

	if(e < 0)
		targetStrain = currentStrain ;
	if(emax < 0)
		targetMaxStrain = currentMaxStrain ;

	double factor = 0.; 
	if(std::abs(currentStrain-currentMaxStrain) > POINT_TOLERANCE)
		factor = (targetStrain-targetMaxStrain)/(currentStrain-currentMaxStrain) ;

	for(size_t i = 0 ; i < tensileStressStrainCurve->size() ; i++)
	{
		double oldStrain = tensileStressStrainCurve->getPoint(i).getX()/renormStrain ;
		double newStrain = targetMaxStrain + (oldStrain-currentMaxStrain)*factor ;
		tensileStressStrainCurve->getPoint(i).setX( newStrain*renormStrain ) ;
	}

	Point p( tensileStressStrainCurve->getPoint(tensileStressStrainCurve->size()-1).getX(), tensileAsymptote->origin().getY() ) ;
	tensileAsymptote->setOrigin(p) ;

}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumCompressiveStress(double s, double smax) 
{
	if(!compressiveStressStrainCurve || (s>0 && smax>0))
		return ;

	double targetStress = s ;
	double targetMaxStress = smax ;

	double currentStress = compressiveStressStrainCurve->getPoint(1).getY()/renormStress ;
	double currentMaxStress = compressiveStressStrainCurve->getPoint(compressiveStressStrainCurve->size()-1).getY()/renormStress ;

	if(s > 0)
		targetStress = currentStress ;
	if(smax > 0)
		targetMaxStress = currentMaxStress ;

	double factor = 0.; 
	if(std::abs(currentStress-currentMaxStress) > POINT_TOLERANCE)
		factor = (targetStress-targetMaxStress)/(currentStress-currentMaxStress) ;

	for(size_t i = 0 ; i < compressiveStressStrainCurve->size() ; i++)
	{
		double oldStress = compressiveStressStrainCurve->getPoint(i).getY()/renormStress ;
		double newStress = targetMaxStress + (oldStress-currentMaxStress)*factor ;
		compressiveStressStrainCurve->getPoint(i).setY( newStress*renormStress ) ;
	}
	Point p( compressiveAsymptote->origin().getX(), compressiveStressStrainCurve->getPoint(compressiveStressStrainCurve->size()-1).getY()) ;
	compressiveAsymptote->setOrigin(p) ;
}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumCompressiveStrain(double e, double emax) 
{
	if(!compressiveStressStrainCurve || (e>0 && emax>0))
		return ;

	double targetStrain = e ;
	double targetMaxStrain = emax ;

	double currentStrain = compressiveStressStrainCurve->getPoint(1).getX()/renormStrain ;
	double currentMaxStrain = compressiveStressStrainCurve->getPoint(compressiveStressStrainCurve->size()-1).getX()/renormStrain ;

	if(e > 0)
		targetStrain = currentStrain ;
	if(emax > 0)
		targetMaxStrain = currentMaxStrain ;

	double factor = 0.; 
	if(std::abs(currentStrain-currentMaxStrain) > POINT_TOLERANCE)
		factor = (targetStrain-targetMaxStrain)/(currentStrain-currentMaxStrain) ;

	for(size_t i = 0 ; i < compressiveStressStrainCurve->size() ; i++)
	{
		double oldStrain = compressiveStressStrainCurve->getPoint(i).getX()/renormStrain ;
		double newStrain = targetMaxStrain + (oldStrain-currentMaxStrain)*factor ;
		compressiveStressStrainCurve->getPoint(i).setX( newStrain*renormStrain ) ;
	}
}

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getTensileStrength() const 
{
	if(!tensileStressStrainCurve)
		return -1 ;

	return tensileStressStrainCurve->getPoint(0).getY()/renormStress ;
}

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getTensileStrain() const 
{
	if(!tensileStressStrainCurve)
		return -1 ;

	return tensileStressStrainCurve->getPoint(0).getX()/renormStrain ;
}

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getCompressiveStrength() const 
{
	if(!compressiveStressStrainCurve)
		return -1 ;

	return compressiveStressStrainCurve->getPoint(0).getY()/renormStress ;
}

double AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::getCompressiveStrain() const 
{
	if(!compressiveStressStrainCurve)
		return -1 ;

	return compressiveStressStrainCurve->getPoint(0).getX()/renormStrain ;
}


}

