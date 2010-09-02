
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

VonMises::VonMises(double thresh) : threshold(thresh)
{
}


VonMises::~VonMises()
{
}

double VonMises::grade(const ElementState &s)
{
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
	mat(TAG_MAX_STRESS,threshold) ;
	return mat ;
}

}
