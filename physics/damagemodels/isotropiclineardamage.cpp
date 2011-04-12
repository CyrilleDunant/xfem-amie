//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "isotropiclineardamage.h"

namespace Mu {

IsotropicLinearDamage::IsotropicLinearDamage() 
{
	getState().resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
}

Vector IsotropicLinearDamage::computeDamageIncrement(ElementState & s)
{
	Vector ret(1) ;
	ret[0] = 1 ;
// 	ret[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE-state[0], state[0]) ;
// 	ret[0] = std::min(.99999, state[0]) ;
// 	ret[0] = std::max(0., state[0]) ;
	return ret ;

}

void IsotropicLinearDamage::artificialDamageStep(double d)
{
	getState()[0] = std::min(getState()[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE_2D) ;
}


Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	return ret*(1.-getState()[0]) ;
}


Matrix IsotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	return ret*(1.-getPreviousState()[0]) ;
}

bool IsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity/fraction ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}


}
