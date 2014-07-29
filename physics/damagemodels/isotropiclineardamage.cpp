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

namespace Amie {

IsotropicLinearDamage::IsotropicLinearDamage() 
{
	getState(true).resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > IsotropicLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void IsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix IsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
	
	if(fractured())
		return m*0 ;
	
	return m*(1.-getState()[0]) ;
}

Matrix IsotropicLinearDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getState()[0]) ;
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

IsotropicLinearDamageRate::IsotropicLinearDamageRate() 
{
	getState(true).resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > IsotropicLinearDamageRate::computeDamageIncrement( Amie::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void IsotropicLinearDamageRate::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix IsotropicLinearDamageRate::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{
	
	if(fraction < 0)
		return m ;
	if(fractured())
		return m*0 ;
	
	double ratio = (p.getT()+1)*.5 ;
	
	return m*(1.-std::min(initalState[0] + ratio*getState()[0], 1.)) ;
}

Matrix IsotropicLinearDamageRate::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

	if(fractured())
		return m*0 ;
	if(fraction < 0)
		return m ;
	
	double ratio = (p.getT()+1)*.5 ;
	
	return m*(1.-std::min(initalState[0] + ratio*getState()[0], 1.)) ;
}

bool IsotropicLinearDamageRate::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity ;
}


IsotropicLinearDamageRate::~IsotropicLinearDamageRate()
{
}

} ;

