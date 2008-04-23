//
// C++ Interface: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MULINEARDAMAGE_H
#define MULINEARDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class LinearDamage : public DamageModel
{
protected:
	Vector state ;
	double thresholdDensity ;
public:
	LinearDamage(int numDof, double threshold) ;

	virtual ~LinearDamage();

	virtual const Vector & damageState() const ;
	virtual void step(ElementState & s) ;
	virtual Matrix apply(const Matrix & m) const;

};

}

#endif
