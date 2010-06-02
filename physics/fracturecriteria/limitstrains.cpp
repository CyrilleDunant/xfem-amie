//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "limitstrains.h"

namespace Mu {

LimitStrains::LimitStrains(double maxdown, double maxup)
	:  maxUpVal(maxup),maxDownVal(maxdown)
{
}


LimitStrains::~LimitStrains()
{
}

double LimitStrains::grade(const ElementState &s) const
{
	Vector pstrain = s.getPrincipalStrains(s.getParent()->getCenter()) ;
	double maxStrain = pstrain.max();
	if(maxUpVal < maxStrain)
		return 1.-std::abs(maxUpVal/maxStrain) ;
	else if(maxDownVal > maxStrain)
		return 1.-std::abs(maxDownVal/maxStrain) ;
	else if (maxStrain > 0)
		return -1.+ std::abs(maxStrain/maxUpVal);
	else
		return -1.+ std::abs(maxStrain/maxDownVal);
	
}

FractureCriterion * LimitStrains::getCopy() const
{
	return new LimitStrains(*this) ;
}

Material LimitStrains::toMaterial()
{
	Material mat ;
	mat(TAG_MAX_TENSILE_STRAIN,maxUpVal) ;
	mat(TAG_MAX_COMPRESSIVE_STRAIN,maxDownVal) ;
	return mat ;
}


}
