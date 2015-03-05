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
#include <fstream>

namespace Amie {

MaximumStrain::MaximumStrain(double up, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up)
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
	s.getField( STRAIN_FIELD, s.getParent()->getBoundingPoints(), pstrain, false) ;
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

Material MaximumStrain::toMaterial()
{
	Material mat ;
	return mat ;
}

double SpaceTimeNonLocalMaximumStrain::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;


	Vector stateBefore( getSmoothedField( PRINCIPAL_STRAIN_FIELD, s, -1) ) ;
	Vector stateAfter( getSmoothedField( PRINCIPAL_STRAIN_FIELD, s, 1) ) ;
	double maxStrainAfter = stateAfter.max() ;
	double maxStrainBefore = stateBefore.max() ;

	metInCompression = false ;
	metInTension = false ;
	if(maxStrainAfter > upVal)
	{
		metInTension = true ;
		return  std::min(1.,1. - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore)) ;
	}
	return -1.+ maxStrainAfter/upVal;
	
	
}

double SpaceTimeNonLocalLinearSofteningMaximumStrain::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_STRAIN_FIELD, s, 1) ) ;

	double Esoft = maxstress / ( yieldstrain - upVal) ;
	double Einst = stateAfter.first.max() / stateAfter.second.max() ; //maxstress/upVal * (1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
// ;
	double Eprev = stateBefore.first.max() / stateBefore.second.max() ;
	double epsMax = (Esoft / (Esoft+Einst))*yieldstrain ;

	double maxStrainBefore = stateBefore.second.max() ;
	double maxStrainAfter = stateAfter.second.max() ;
	double maxStressBefore = stateBefore.first.max() ;
	double maxStressAfter = stateAfter.first.max() ;

//		std::cout << maxStressBefore << "/" << epsMax << "/" << maxStressAfter << "\t" ;

		if( false ) //s.getParent()->getBoundingPoint(0).getId() == 9)
			std::cout << s.getNodalDeltaTime() << "\t" << s.getParent()->getBoundingPoint(0).getId() << "\t" << s.getParent()->getBoundingPoint(1).getId() << "\t" << s.getParent()->getBoundingPoint(2).getId() << "\t" <<  maxStrainBefore << "\t" << maxStrainAfter << "\t" << maxStressBefore << "\t" << maxStressAfter << "\t" << upVal << "\t" << yieldstrain << "\t" << epsMax << "\t" << epsMax*Einst << "\t" << Esoft << "\t" << Einst << "\t" <<  (Esoft / (Esoft+Einst)) << "\t" << std::min(1.,1. - (epsMax - maxStrainBefore) / (maxStrainAfter - maxStrainBefore)) <<  std::endl ;


	metInCompression = false ;
	metInTension = false ;
	if(maxStrainAfter > epsMax)
	{
		metInTension = true ;
/*		std::cout << (Esoft-Einst)/Esoft << std::endl ;
		Point before(maxStrainBefore, maxStressBefore) ;
		Point after(maxStrainAfter, maxStressAfter) ;
		Point intersect(epsMax, epsMax*Einst) ;
		before.print() ;
		intersect.print() ;
		after.print() ;*/
		return  std::min(1.,1. - (epsMax - maxStrainBefore) / (maxStrainAfter - maxStrainBefore)) ;
	}
//	std::cout << maxStrainBefore << "\t" << maxStrainAfter << "\t" <<  -1.+ maxStrainAfter/epsMax << std::endl ;
	return -1.+ maxStrainAfter/epsMax;

/*	Point before( (stateBefore.second).max(), (stateBefore.first).max()/1e6 ) ;
	Point after( (stateAfter.second).max(), (stateAfter.first).max()/1e6 ) ;

		SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string file, double E_, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::vector<Point> p, double E_, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;


	Segment history(before, after) ;
	Line current(history) ;

	Point x0( upVal, maxstress/1e6 ) ;
	Point x1( yieldstrain, 0. ) ;
	Line behaviour(Segment(x0, x1)) ;

	
	Point t = behaviour.intersection(current) ;


	if(t.getX() < after.getX())
	{
		if(history.on(t))
		{
			metInTension = true ;
			ret =  std::min(1.,1. - (t.getX() - before.getX()) / (after.getX() - before.getX())) ;
		}else
		ret = -1. + std::min(before.getX(),after.getX())/yieldstrain ;
	}
	else
		ret = -1.+ after.getX()/t.getX();
	
	if(false)//ret > 0)
	{
	before.print() ;
	after.print() ;
	x0.print() ;
	x1.print() ;
	t.print() ;
	std::cout << ret << std::endl ;
	}	

	return ret ;*/
}

double SpaceTimeNonLocalMaximumStress::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( REAL_STRESS_FIELD, STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( REAL_STRESS_FIELD, STRAIN_FIELD, s, 1) ) ;
	double maxStressAfter = stateAfter.first.max() ;
	double maxStressBefore = stateBefore.first.max() ;

	metInCompression = false ;
	metInTension = false ;
	if(maxStressAfter > maxstress)
	{
		metInTension = true ;
 //       std::cout << maxStressBefore << " " << maxStressAfter << std::endl ;
		return std::min(1., 1. - (maxstress - maxStressBefore) / (maxStressAfter - maxStressBefore)) ;
	}
	return -1.+ maxStressAfter/maxstress;
	
	
}



SpaceTimeNonLocalEllipsoidalMixedCriterion::SpaceTimeNonLocalEllipsoidalMixedCriterion(double up, double mstr, double E0, double Einf, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(up, mirroring, delta_x, delta_y, delta_z),maxstress(mstr), E_inst(E0), E_relaxed(Einf)
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
	base = std::abs(VirtualMachine().eval( surface->getEllipseFormFunction() )) ;
	upVal = surface->getCenter().getX()-surface->getMajorRadius() ;

}

double SpaceTimeNonLocalEllipsoidalMixedCriterion::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	metInCompression = false ;
	metInTension = false ;

	std::pair<Vector, Vector> stateBefore( getSmoothedFields( REAL_STRESS_FIELD, STRAIN_FIELD, s, -1) ) ;
	std::pair<Vector, Vector> stateAfter( getSmoothedFields( REAL_STRESS_FIELD, STRAIN_FIELD, s, 1) ) ;

	Point before( stateBefore.second.max()*renormStrain/surface->getMinorRadius(), stateBefore.first.max()*renormStress/surface->getMajorRadius()) ;
	Point after( stateAfter.second.max()*renormStrain/surface->getMinorRadius() , stateAfter.first.max()*renormStress/surface->getMajorRadius()) ;
	
	
	if(stateAfter.second.max() <= POINT_TOLERANCE_2D || stateAfter.first.max() <= POINT_TOLERANCE_2D)
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
// 			std::cout << upVal << "\t" << maxStrainAfter  << "\t" << maxStrainBefore << std::endl ;
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

