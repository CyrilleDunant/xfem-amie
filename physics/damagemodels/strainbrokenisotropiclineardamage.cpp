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
#include "strainbrokenisotropiclineardamage.h"

namespace Mu {

StrainBrokenIsotropicLinearDamage::StrainBrokenIsotropicLinearDamage(int numDof, double characteristicRadius, double limitStrain) : DamageModel(characteristicRadius), limitStrain(limitStrain)
{
	state.resize(1, 0.);
	state[0] = 0 ;
	isNull = false ;
}

const Vector & StrainBrokenIsotropicLinearDamage::damageState() const
{
	return state ;
}

Vector & StrainBrokenIsotropicLinearDamage::damageState()
{
	return state ;
}


void StrainBrokenIsotropicLinearDamage::step(ElementState & s)
{
	previousstate.resize(state.size());
	previousstate = state ;
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
		fraction = std::min(fraction, 1.) ;
	}
	if(s.getPrincipalStrains(s.getParent()->getCenter()).max() > limitStrain)
	{
		state[0] = thresholdDamageDensity/fraction+POINT_TOLERANCE ;
	}
	state[0] += damageDensityIncrement*fraction ; 
	state[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;

}

void StrainBrokenIsotropicLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE) ;
}


Matrix StrainBrokenIsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0.00001 ;
	return ret*(1.-state[0]) ;
}


Matrix StrainBrokenIsotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0.00001 ;
	return ret*(1.-previousstate[0]) ;
}

bool StrainBrokenIsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state[0] >= thresholdDamageDensity/fraction ;
}


StrainBrokenIsotropicLinearDamage::~StrainBrokenIsotropicLinearDamage()
{
}


}
