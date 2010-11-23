
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "vonmises.h"
namespace Mu {

VonMises::VonMises(double thresh, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), threshold(thresh)
{
}


VonMises::~VonMises()
{
}

double VonMises::grade(const ElementState &s)
{
	metInCompression = true ;
	metInTension = true ;
	double maxStress = s.getMaximumVonMisesStress() ;

	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}
}

FractureCriterion * VonMises::getCopy() const
{
	return new VonMises(threshold) ;
}

Material VonMises::toMaterial()
{
	Material mat ;
	return mat ;
}


VonMisesStrain::VonMisesStrain(double thresh) : threshold(thresh)
{
}


VonMisesStrain::~VonMisesStrain()
{
}

double VonMisesStrain::grade(const ElementState &s)
{
	double maxStress = s.getVonMisesStrain(s.getParent()->getCenter()) ;

	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}
}

FractureCriterion * VonMisesStrain::getCopy() const
{
	return new VonMisesStrain(threshold) ;
}

Material VonMisesStrain::toMaterial()
{
	Material mat ;
	return mat ;
}

}
