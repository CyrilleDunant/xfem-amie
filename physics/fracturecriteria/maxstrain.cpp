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
#include "maxstrain.h"
#include "../mesher/delaunay.h"

namespace Mu {

MaximumStrain::MaximumStrain(double up)
	: upVal(up)
{
}


MaximumStrain::~MaximumStrain()
{
}

double MaximumStrain::grade(const ElementState &s) const
{
	Vector pstrain = s.getStrain(s.getParent()->getBoundingPoints()) ;
	double maxStrain = pstrain.max();
	if(maxStrain > upVal)
		return 1.-std::abs(upVal/maxStrain) ;
	else
		return 0 ;
	
}

FractureCriterion * MaximumStrain::getCopy() const
{
	return new MaximumStrain(*this) ;
}

}
