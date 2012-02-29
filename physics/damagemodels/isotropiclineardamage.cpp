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
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > IsotropicLinearDamage::computeDamageIncrement( Mu::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void IsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getState()[0]) ;
}


Matrix IsotropicLinearDamage::applyPrevious(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getPreviousState()[0]) ;
}

bool IsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}

} ;

