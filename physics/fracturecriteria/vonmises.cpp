
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
#include "vonmises.h"
namespace Mu {

VonMises::VonMises(double thresh, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), threshold(thresh)
{
}


VonMises::~VonMises()
{
}

double VonMises::grade(ElementState &s)
{
	Vector maxStress(0.,1) ;
	s.getAverageField(VON_MISES_REAL_STRESS_FIELD, maxStress) ;

	if(maxStress[0] > threshold )
	{
		return 1. - std::abs(threshold/maxStress[0]) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress[0]/threshold);
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

double VonMisesStrain::grade(ElementState &s)
{
	Vector maxStrain(0.,1) ;
	s.getAverageField(VON_MISES_STRAIN_FIELD, maxStrain) ;

	if(maxStrain[0] > threshold )
	{
		return 1. - std::abs(threshold/maxStrain[0]) ;
	}
	else 
	{
		return -1.+ std::abs(maxStrain[0]/threshold);
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
