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
	double e2 = p[1].getX() ;
	double s1 = p[0].getY() ;
	double s2 = p[1].getY() ;

	double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
	std::valarray<Point *> curve ;
	curve.resize( p.size()+1 ) ;
	Point * origin = new Point( 0., stressAtOrigin) ;
	curve[0] = origin ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		curve[i+1] = &p[i] ;
	}
	double stressAtInfinity = curve[p.size()]->getY() ;
	Point horizontal(1.,0.) ;

	stressStrainCurve = new SegmentedLine( curve ) ;
	asymptote = new Line( *curve[p.size()], horizontal ) ;

}

SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & p, double E_, double e, double s_, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( p[0].getX(), mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s_)
{
	double e1 = p[0].getX()*renormStrain ;
	upVal = e1 ;
	double e2 = p[1].getX()*renormStrain ;
	double s1 = p[0].getY()*renormStress ;
	double s2 = p[1].getY()*renormStress ;

	double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
	std::valarray<Point *> curve ;
	curve.resize( p.size()+1 ) ;
	Point * origin = new Point( 0., stressAtOrigin) ;
	curve[0] = origin ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		curve[i+1] = new Point(p[i].getX()*renormStrain, p[i].getY()*renormStress) ;
	}
	double stressAtInfinity = curve[p.size()]->getY()*renormStress ;
	Point horizontal(1.,0.) ;

	stressStrainCurve = new SegmentedLine( curve ) ;
	asymptote = new Line( *curve[p.size()], horizontal ) ;
}

SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( SegmentedLine * s, Line * l, double E_, double e, double s_, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(E_, mirroring, delta_x, delta_y, delta_z), stressStrainCurve(s), E(E_),  asymptote(l), renormStrain(e), renormStress(s_)
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
	return new SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( p, E, renormStrain, renormStress ) ;
}

double SpaceTimeNonLocalMultiLinearSofteningFractureCriterion::grade(ElementState &s)  
{
	std::pair<Vector, Vector> stateBefore( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, 1) ) ;

/*	if(typeid(s) == typeid(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables))
	{
		std::map<std::string, double> dummy ;
		double imposed = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(s).get("imposed_deformation", dummy ) ;
		stateBefore.second -= imposed ;
		stateAfter.second -= imposed ;
	}*/

	Point before( stateBefore.second.max()*renormStrain, stateBefore.first.max()*renormStress ) ;
	Point after( stateAfter.second.max()*renormStrain, stateAfter.first.max()*renormStress ) ;

	bool compressive = stressStrainCurve->getPoint(0).getY() < 0 ;
	if( compressive )
	{
		before.setX( stateBefore.second.min()*renormStrain) ; before.setY(stateBefore.first.min()*renormStress ) ;
		after.setX( stateAfter.second.min()*renormStrain) ; after.setY( stateAfter.first.min()*renormStress ) ;
	}

	Segment history(before, after) ;

	std::vector<Point> intersection = history.intersection( stressStrainCurve ) ;
	Point inter ;

	if(!history.intersects(stressStrainCurve))
	{
		Line lhistory( history ) ;	
		std::vector<Point> lintersection = lhistory.intersection( stressStrainCurve ) ;
		if(lintersection.size() == 0)
		{
			if(!lhistory.intersects( *asymptote ))
				return -1. ;
			inter = lhistory.intersection( *asymptote ) ;
		}
		else
		{
			inter = lintersection[0] ;			
		}
	}
	else
	{
		inter = intersection[0] ;
	}

	if(compressive)
	{
		if(inter.getX() < before.getX() && inter.getX() > after.getX())
			return std::min(1., (after.getX()-inter.getX())/(after.getX()-before.getX()) ) ;
		if(inter.getX() < before.getX() && inter.getX() < after.getX())
			return std::min( -POINT_TOLERANCE, std::max(-1., -1.+(after.getX()-before.getX())/(inter.getX()-before.getX()) ) ) ;
	}
	else
	{
		if(inter.getX() > before.getX() && inter.getX() < after.getX())
			return std::min(1., (after.getX()-inter.getX())/(after.getX()-before.getX()) ) ;
		if(inter.getX() > before.getX() && inter.getX() > after.getX())
			return std::min( -POINT_TOLERANCE, std::max(-1., -1.+(after.getX()-before.getX())/(inter.getX()-before.getX()) ) ) ;
	}

	return -1. ;

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

AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string ftension, std::string fcompression, double E_, double f, double e, double s, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(E_, mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s), fmax(f), currentFraction(0.)
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
		double e2 = p[1].getX() ;
		double s1 = p[0].getY() ;
		double s2 = p[1].getY() ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( p.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			curve[i+1] = &p[i] ;
		}
		double stressAtInfinity = curve[p.size()]->getY() ;
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


		double e1 = p[0].getX() ;
		double e2 = p[1].getX() ;
		double s1 = p[0].getY() ;
		double s2 = p[1].getY() ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( p.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < p.size() ; i++)
		{
			curve[i+1] = &p[i] ;
		}
		double stressAtInfinity = curve[p.size()]->getY() ;
		Point horizontal(1.,0.) ;

		compressiveStressStrainCurve = new SegmentedLine( curve ) ;
		compressiveAsymptote = new Line( *curve[p.size()], horizontal ) ;
	}
}

AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & ptension, const std::vector<Point> & pcompression, double E_, double f, double e, double s, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain( 0., mirroring, delta_x, delta_y, delta_z), E(E_), renormStrain(e), renormStress(s), fmax(f), currentFraction(0.)
{
	tensileStressStrainCurve = nullptr ;
	compressiveStressStrainCurve = nullptr ;
	tensileAsymptote = nullptr ;
	compressiveAsymptote = nullptr ;

	if(ptension.size() > 0)
	{
		double e1 = ptension[0].getX()*renormStrain ;
		upVal = e1 ;
		double e2 = ptension[1].getX()*renormStrain ;
		double s1 = ptension[0].getY()*renormStress ;
		double s2 = ptension[1].getY()*renormStress ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( ptension.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < ptension.size() ; i++)
		{
			curve[i+1] = new Point(ptension[i].getX()*renormStrain, ptension[i].getY()*renormStress) ;
		}
		double stressAtInfinity = curve[ptension.size()]->getY()*renormStress ;
		Point horizontal(1.,0.) ;

		tensileStressStrainCurve = new SegmentedLine( curve ) ;
		tensileAsymptote = new Line( *curve[ptension.size()], horizontal ) ;
	}
	if(pcompression.size() > 0)
	{
		double e1 = pcompression[0].getX()*renormStrain ;
		double e2 = pcompression[1].getX()*renormStrain ;
		double s1 = pcompression[0].getY()*renormStress ;
		double s2 = pcompression[1].getY()*renormStress ;

		double stressAtOrigin = s1 ;//(s2 - (s1*e2/e1))/(1.-e2/e1) ;
		std::valarray<Point *> curve ;
		curve.resize( pcompression.size()+1 ) ;
		Point * origin = new Point( 0., stressAtOrigin) ;
		curve[0] = origin ;
		for(size_t i = 0 ; i < pcompression.size() ; i++)
		{
			curve[i+1] = new Point(pcompression[i].getX()*renormStrain, pcompression[i].getY()*renormStress) ;
		}
		double stressAtInfinity = curve[pcompression.size()]->getY()*renormStress ;
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
	return new AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( ptension, pcompression, E, fmax, renormStrain, renormStress ) ;
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
	if(currentFraction > fmax)
		return -1 ;

	double gtension = -1. ;
	double gcompression = -1. ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s, 1) ) ;

	if(tensileStressStrainCurve)
	{	
		bool found = false ;
		Point before( stateBefore.second.max()*renormStrain, stateBefore.first.max()*renormStress ) ;
		Point after( stateAfter.second.max()*renormStrain, stateAfter.first.max()*renormStress ) ;

		Segment history(before, after) ;

		std::vector<Point> intersection = history.intersection( tensileStressStrainCurve ) ;
		Point inter ;

		if(!history.intersects(tensileStressStrainCurve))
		{
			Line lhistory( history ) ;	
			std::vector<Point> lintersection = lhistory.intersection( tensileStressStrainCurve ) ;
			if(lintersection.size() == 0)
			{
				inter = lhistory.intersection( *tensileAsymptote ) ;
				if(lhistory.intersects( *tensileAsymptote ) && inter.getX() > tensileAsymptote->origin().getX())
					found = true ;
			}
			else
			{
				found = true ;
				inter = lintersection[0] ;			
			}
		}
		else
		{
			found = true ;
			inter = intersection[0] ;
		}

		if(found)
		{
			if(inter.getX() > before.getX() && inter.getX() < after.getX())
				gtension = std::min(1., (after.getX()-inter.getX())/(after.getX()-before.getX()) ) ;
			else if(inter.getX() > before.getX() && inter.getX() > after.getX())
				gtension = std::min( -POINT_TOLERANCE, std::max(-1., -1.+(after.getX()-before.getX())/(inter.getX()-before.getX()) ) ) ;
		}
	}

	if(compressiveStressStrainCurve)
	{	
		bool found = false ;
		Point before( stateBefore.second.min()*renormStrain, stateBefore.first.min()*renormStress ) ;
		Point after( stateAfter.second.min()*renormStrain, stateAfter.first.min()*renormStress ) ;

		Segment history(before, after) ;

		std::vector<Point> intersection = history.intersection( compressiveStressStrainCurve ) ;
		Point inter ;

		if(!history.intersects(compressiveStressStrainCurve))
		{
			Line lhistory( history ) ;	
			std::vector<Point> lintersection = lhistory.intersection( compressiveStressStrainCurve ) ;
			if(lintersection.size() == 0)
			{
				inter = lhistory.intersection( *compressiveAsymptote ) ;
				if(lhistory.intersects( *compressiveAsymptote ) && inter.getX() < compressiveAsymptote->origin().getX())
				{
					found = true ;
				}
			}
			else
			{
				found = true ;
				inter = lintersection[0] ;			
			}
		}
		else
		{
			found = true ;
			inter = intersection[0] ;
		}

		if(found)
		{
			if(inter.getX() < before.getX() && inter.getX() > after.getX())
			{
				gcompression = std::min(1., (after.getX()-inter.getX())/(after.getX()-before.getX()) ) ;				
			}
			else if(inter.getX() < before.getX() && inter.getX() < after.getX())
			{
				gcompression = std::min( -POINT_TOLERANCE, std::max(-1., -1.+(after.getX()-before.getX())/(inter.getX()-before.getX()) ) ) ;
			}
		}
	}
	
	return std::max(gcompression, gtension) ;

}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumTensileStress(double s, double smax) 
{
	if(!tensileStressStrainCurve)
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
	if(!tensileStressStrainCurve)
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
}

void AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion::setMaximumCompressiveStress(double s, double smax) 
{
	if(!compressiveStressStrainCurve)
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
	if(!compressiveStressStrainCurve)
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

