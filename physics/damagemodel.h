//
// C++ Interface: damagemodel
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __DAMAGE_MODEL_H__
#define __DAMAGE_MODEL_H__

#include "../elements/integrable_entity.h"
#include "../matrixops.h"

namespace Mu
{

class DamageModel
{
public:
	DamageModel() { isNull = true ; } ;
	virtual const Vector & damageState() const = 0 ;
	virtual void step(ElementState & s) = 0 ;
	virtual Matrix apply(const Matrix & m) const = 0 ;
	bool isNull ;
} ;

class NullDamage : public DamageModel
{
protected:
	Vector state ;
public:
	NullDamage() : state(0) { } ;
	virtual const Vector & damageState() const { return state ;} ;
	virtual void step(ElementState & s) { } ;
	virtual Matrix apply(const Matrix & m) const {return m ;} ;
} ;

}

#endif
