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

IsotropicLinearDamage::IsotropicLinearDamage(int numDof, double characteristicRadius) : DamageModel(characteristicRadius),
 state(1)
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
	if(fraction < 0)
	{
		double volume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			volume = s.getParent()->area() ;
		else
			volume = s.getParent()->volume() ;
		
		double charVolume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			charVolume = M_PI*characteristicRadius*characteristicRadius ;
		else
			charVolume = 4./3.*M_PI*characteristicRadius*characteristicRadius*characteristicRadius ;
		fraction = volume/charVolume ;
	}
	state[0] += damageDensityIncrement*fraction ; 
	state[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;

}

void IsotropicLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE) ;
}


Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0.00001 ;
	return ret*(1.-state[0]) ;
}


bool IsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state[0] >= thresholdDamageDensity/fraction ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}


}
