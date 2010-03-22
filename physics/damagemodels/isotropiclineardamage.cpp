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
	state[0] = 0 ;
	isNull = false ;
}

const Vector & IsotropicLinearDamage::damageState() const
{
	return state ;
}

Vector & IsotropicLinearDamage::damageState()
{
	return state ;
}


void IsotropicLinearDamage::step(ElementState & s)
{
	double maxD = .999999 ; 

	state[0] += .1 ; 
	state[0] = std::min(maxD, state[0]) ;

}

Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	return ret*(1.-state[0]) ;
}


bool IsotropicLinearDamage::fractured() const 
{
	return state[0] >= .9 ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}


}
