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
#ifndef MUMOHRCOULOMB_H
#define MUMOHRCOULOMB_H

#include "fracturecriterion.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when one of the principal stresses is below or above the prescribed limits 
	
*/
class MohrCoulomb : public FractureCriterion
{
	double upVal ;
	double downVal ;
public:
    MohrCoulomb(double up, double down);

    virtual ~MohrCoulomb();

	virtual bool met(const ElementState * s) const ;
};

}

#endif
