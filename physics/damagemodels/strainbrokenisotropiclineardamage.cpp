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
	previousstate.resize(1, 0.);
	state[0] = 0 ;
	isNull = false ;
}

Vector StrainBrokenIsotropicLinearDamage::computeDamageIncrement(ElementState & s)
{
	Vector ret(1) ; ret = 0 ;
	if(s.getPrincipalStrains(s.getParent()->getCenter()).max() > limitStrain)
	{
		ret[0] = thresholdDamageDensity/fraction+POINT_TOLERANCE ;
	}
	ret[0] += 0.5 ; 
// 	ret[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;
	return ret ;

}

void StrainBrokenIsotropicLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE) ;
}


Matrix StrainBrokenIsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*1e-6 ;
	return ret*(1.-state[0]) ;
}


Matrix StrainBrokenIsotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*1e-6 ;
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
