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
	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;
	double maxStrainAfter = stateAfter.second.max() ;
	double maxStrainBefore = stateBefore.second.max() ;
	metInCompression = false ;
	metInTension = false ;
	if(maxStrainAfter > upVal)
	{
		metInTension = true ;
		return std::min(1., 1 - (upVal - maxStrainBefore) / (maxStrainAfter - maxStrainBefore)) ;
	}
	else
		return -1.+ std::abs(maxStrainAfter/upVal);
	
}

}
