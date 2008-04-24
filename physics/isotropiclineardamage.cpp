//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "isotropiclineardamage.h"

namespace Mu {

IsotropicLinearDamage::IsotropicLinearDamage(int numDof)
 : state(1)
{
	isNull = false ;
}

const Vector & IsotropicLinearDamage::damageState() const
{
	return state ;
}

void IsotropicLinearDamage::step(ElementState & s)
{
	state[0] += std::max(state[0], 0.01) ;
// 	if(state[0] > .6)
// 		state[0] = .999 ;
	state[0] = std::min(.999, state[0]) ;

}

Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	return ret*(1.-state[0]) ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}


}
