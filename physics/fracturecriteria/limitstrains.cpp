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
#include "limitstrains.h"

namespace Amie {

LimitStrains::LimitStrains(double maxdown, double maxup) :  maxUpVal(maxup),maxDownVal(maxdown)
{
	metInCompression = false ;
	metInTension = false ;
}


LimitStrains::~LimitStrains()
{
}

double LimitStrains::grade(ElementState &s)
{
	Vector pstrain(0., s.getParent()->spaceDimensions()) ;
	s.getField( PRINCIPAL_TOTAL_STRAIN_FIELD, s.getParent()->getCenter(), pstrain, false) ;
	double maxStrain = pstrain.max();
	metInCompression = false ;
	metInTension = false ;
	if(maxUpVal < maxStrain)
	{
		metInTension = true ;
		if(maxDownVal > maxStrain)
			metInCompression = true ;
		return 1.-std::abs(maxUpVal/maxStrain) ;
	}
	else if(maxDownVal > maxStrain)
	{
		metInCompression = true ;
		return 1.-std::abs(maxDownVal/maxStrain) ;
	}
	else if (maxStrain > 0)
		return -1.+ std::abs(maxStrain/maxUpVal);
	else
		return -1.+ std::abs(maxStrain/maxDownVal);
	
}

FractureCriterion * LimitStrains::getCopy() const
{
	return new LimitStrains(*this) ;
}


}
