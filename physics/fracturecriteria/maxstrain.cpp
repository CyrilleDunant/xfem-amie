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

namespace Mu {

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


	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;
	double maxStrainAfter = stateAfter.second.max() ;
	double maxStrainBefore = stateBefore.second.max() ;

	metInCompression = false ;
	metInTension = false ;
	if(maxStrainAfter > upVal)
	{
		metInTension = true ;
//		std::cout << upVal << "\t" << maxStrainAfter  << "\t" << maxStrainBefore << std::endl ;
		return  1. - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore) ;
	}
	return -1.+ maxStrainAfter/upVal;
	
	
}

double SpaceTimeNonLocalMaximumStress::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;


	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;
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

SpaceTimeNonLocalEllipsoidalMixedCriterion::SpaceTimeNonLocalEllipsoidalMixedCriterion(double up, double mstr, double E0, double Einf, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : MaximumStrain(up, mirroring, delta_x, delta_y, delta_z),maxstress(mstr), E_inst(E0), E_relaxed(Einf)
{
	double strainrelaxed = maxstress/E_relaxed ;
	double stressinst = upVal*E_inst ;

	Point center(strainrelaxed, stressinst) ;
	Point mj(std::abs(upVal - strainrelaxed), 0.) ;
	Point mn(0.,std::abs(stressinst - maxstress)) ;

	Ellipse ell(center, mj, mn) ;
	surface = ell.getEllipseFormFunction() ;


}


double SpaceTimeNonLocalEllipsoidalMixedCriterion::grade(ElementState &s)
{
	if( s.getParent()->getBehaviour()->fractured() )
		return -1 ;

	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;

	Point after( stateAfter.second.max(), stateAfter.first.max() ) ;
	Point before( stateBefore.second.max(), stateBefore.first.max() ) ;

	double kafter = VirtualMachine().eval(surface, after) ;
	double kbefore = VirtualMachine().eval(surface, before) ;

	if(after.x > maxstress/E_relaxed)
		kafter = -2.-kafter ;
	if(after.y > upVal*E_inst)
		kafter = -2.-kafter ;

//	Segment transition(after, before) ;
	metInCompression = false ;
	metInTension = false ;
	if(kafter < 0)
	{
		metInTension = true ;
		return std::min(1., std::abs(kafter)/(std::abs(kbefore) + std::abs(kafter))) ;
	}

/*	if(transition.intersects(&surface))
	{
		metInTension = true ;
		Point k = transition.intersection(&surface)[0] ;
		double pos = (k-before).norm()/(after-before).norm() ;
		return std::min(1., 1.-pos) ;
	}*/

	if(kafter < kbefore)
		return (-1.)* kafter/kbefore;
	
	return -1. ;
}

}
