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
		if(s.getParent()->spaceDimensions() == 2)
			volume = s.getParent()->area() ;
		else
			volume = s.getParent()->volume() ;
		
		double charVolume ;
		if(s.getParent()->spaceDimensions() == 2)
			charVolume = M_PI*characteristicRadius*characteristicRadius ;
		else
			charVolume = 4./3*M_PI*characteristicRadius*characteristicRadius*characteristicRadius ;
		fraction = volume/charVolume ;
	}
	
	state[0] += damageDensityIncrement*fraction ; 
	state[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;

}

void IsotropicLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,0.9999) ;
}


Matrix IsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	return ret*(1.-state[0]) ;
}


bool IsotropicLinearDamage::fractured() const 
{
	return state[0] >= thresholdDamageDensity/fraction ;
}


IsotropicLinearDamage::~IsotropicLinearDamage()
{
}


}
