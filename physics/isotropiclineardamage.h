//
// C++ Interface: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_ISOTROPIC_LINEARDAMAGE_H
#define MU_ISOTROPIC_LINEARDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class IsotropicLinearDamage : public DamageModel
{
protected:
	Vector state ;
public:
	IsotropicLinearDamage(int numDof) ;

	virtual ~IsotropicLinearDamage();

	virtual const Vector & damageState() const ;
	virtual void step(ElementState & s) ;
	virtual Matrix apply(const Matrix & m) const;

};

}

#endif
