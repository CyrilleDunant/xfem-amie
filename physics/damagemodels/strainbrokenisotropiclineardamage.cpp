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
#include "strainbrokenisotropiclineardamage.h"

namespace Amie {

StrainBrokenIsotropicLinearDamage::StrainBrokenIsotropicLinearDamage(int numDof, double limitStrain) :  limitStrain(limitStrain)
{
	state.resize(1, 0.);
	state[0] = 0 ;
	isNull = false ;
}

std::pair<Vector, Vector> StrainBrokenIsotropicLinearDamage::computeDamageIncrement(ElementState & s)
{
	Vector ret(1) ; 
	Vector pstrain(0., s.getParent()->spaceDimensions()) ;
	s.getField( PRINCIPAL_TOTAL_STRAIN_FIELD, s.getParent()->getCenter(), pstrain, false) ;
	if(pstrain.max() > limitStrain)
	{
		ret[0] = (1.+damageDensityTolerance*64.)*thresholdDamageDensity/fraction - state[0] ; //thresholdDamageDensity/fraction+POINT_TOLERANCE ;
	}
// 	ret[0] += 0.5 ; 
// 	ret[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;
	return std::make_pair(state, ret) ;

}

void StrainBrokenIsotropicLinearDamage::computeDelta(ElementState & s)
{
	Vector ret(1) ; 
	Vector pstrain(0., s.getParent()->spaceDimensions()) ;
	s.getField( PRINCIPAL_TOTAL_STRAIN_FIELD, s.getParent()->getCenter(), pstrain, false) ;
	if(pstrain.max() > limitStrain)
	{
		ret[0] = (1.+damageDensityTolerance*64.)*thresholdDamageDensity/fraction - state[0] ; //thresholdDamageDensity/fraction+POINT_TOLERANCE ;
	}
	// 	ret[0] += 0.5 ; 
	// 	ret[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;
	delta = (ret-state).max() ;
}


Matrix StrainBrokenIsotropicLinearDamage::apply(const Matrix & m, const Point & p , const IntegrableEntity * e , int g) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*1e-6 ;
	return ret*(1.-state[0]) ;
}


bool StrainBrokenIsotropicLinearDamage::fractured(int direction) const 
{
	if(fraction < 0)
		return false ;
	return state[0] >= thresholdDamageDensity/fraction ;
}


StrainBrokenIsotropicLinearDamage::~StrainBrokenIsotropicLinearDamage()
{
}

DamageModel * StrainBrokenIsotropicLinearDamage::getCopy() const
{
    StrainBrokenIsotropicLinearDamage * ret = new StrainBrokenIsotropicLinearDamage(1,limitStrain) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}



}
