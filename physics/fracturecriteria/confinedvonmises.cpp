
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "confinedvonmises.h"
namespace Mu {

ConfinedVonMises::ConfinedVonMises(double threshdown, double threshup, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, thresholdup(threshup),thresholddown(threshdown)
{
}


ConfinedVonMises::~ConfinedVonMises()
{
}

double ConfinedVonMises::grade(ElementState &s)
{
	Vector v(0.,1) ;
	Vector pstress0(0., 3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
	s.getField( VON_MISES_REAL_STRESS_FIELD, PRINCIPAL_REAL_STRESS_FIELD, s.getParent()->getCenter(), v, pstress0, false) ;
	double maxStress = v[0] ;
	
	if(pstress0.max() < 0) //we are confined
	{
		if(maxStress > 1.25*thresholddown )
		{
			return 1. - std::abs(1.25*thresholddown/maxStress) ;
		}
		else 
		{
			return -1.+ std::abs(1.25*maxStress/thresholddown);
		}
	}
	if(std::abs(pstress0.min()) > std::abs(pstress0.max())) //we are in compression
	{
		if(maxStress > thresholddown )
		{
			return 1. - std::abs(thresholddown/maxStress) ;
		}
		else 
		{
			return -1.+ std::abs(maxStress/thresholddown);
		}
	}
	else // we are in tension
	{
		if(maxStress > thresholdup )
		{
			return 1. - std::abs(thresholdup/maxStress) ;
		}
		else 
		{
			return -1.+ std::abs(maxStress/thresholdup);
		}
	}
}

FractureCriterion * ConfinedVonMises::getCopy() const
{
	return new ConfinedVonMises(thresholdup, thresholddown) ;
}

Material ConfinedVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}

}
