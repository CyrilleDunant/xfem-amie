//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MAX_STRAIN_H__
#define MAX_STRAIN_H__

#include "fracturecriterion.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The maximum (tensile) strain criterion is met when a strain limit is reached.
	
*/
class MaximumStrain : public FractureCriterion
{
	double upVal ;
public:
	MaximumStrain(double up);

	virtual ~MaximumStrain();

	virtual bool met(const ElementState * s) const ;
};

}

#endif
