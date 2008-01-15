//
// C++ Interface: ruptureenergy
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ENR_CRIT_H
#define ENR_CRIT_H

#include "fracturecriterion.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when one of the principal stresses is below or above the prescribed limits 
	
*/
class RuptureEnergy : public FractureCriterion
{
	

public:

	double energy ;

    RuptureEnergy(double energy);

    virtual ~RuptureEnergy();

	virtual bool met(const ElementState & s) const ;

	virtual FractureCriterion * getCopy() const;
};

}

#endif
