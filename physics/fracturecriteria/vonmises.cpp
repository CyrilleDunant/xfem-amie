
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
namespace Amie {

VonMises::VonMises(double thresh) : threshold(thresh)
{
}


VonMises::~VonMises()
{
}

double VonMises::grade(ElementState &s)
{
	Vector maxStress(0.,1) ;
	s.getAverageField(VON_MISES_REAL_STRESS_FIELD, maxStress, nullptr, 0) ;

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
	VonMises * ret = new VonMises(threshold) ;
	ret->copyEssentialParameters( this ) ;
	return ret ;
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
	s.getAverageField(VON_MISES_STRAIN_FIELD, maxStrain, nullptr, 0) ;

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
	VonMisesStrain * ret = new VonMisesStrain(threshold) ;
	ret->copyEssentialParameters( this ) ;
	return ret ;
}

}
