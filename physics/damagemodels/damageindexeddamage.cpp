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
#include "damageindexeddamage.h"

namespace Mu {

IndexedLinearDamage::IndexedLinearDamage(int numDof, FractureCriterion * e) : DamageModel(e->getMaterialCharacteristicRadius())
{
	state.resize(1, 0.);
	state[0] = 0 ;
	fixedDamage.resize(1,0.);
	isNull = false ;
}

const Vector & IndexedLinearDamage::damageState() const
{
	return state ;
}

Vector & IndexedLinearDamage::damageState()
{
	return state ;
}


void IndexedLinearDamage::step(ElementState & s)
{
	previousstate.resize(state.size());
	previousstate = state ;
	if(fraction < 0)
	{
		double volume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			volume = sqrt(s.getParent()->area()) ;
		else
			volume = pow(s.getParent()->volume(), 2./3.) ;
		
		double charVolume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			charVolume = sqrt(M_PI*characteristicRadius*characteristicRadius) ;
		else
			charVolume = pow(4./3.*M_PI*characteristicRadius*characteristicRadius*characteristicRadius, 2./3.) ;
		fraction = volume/charVolume ;
		if(fraction > 1)
			std::cout << "elements too large for damage characteristic radius!" << std::endl ;
		fraction = std::min(fraction, 1.) ;
	}
	double E_2 = s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())[0][0] ; E_2*=E_2 ;
	double l_2 = s.getParent()->area() ; 
	double maxincrement = std::abs((l_2*E_2-1.)/(l_2+l_2*E_2)) ;
	state[0] += std::min(damageDensityIncrement*fraction, maxincrement ) ; 
	state[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, state[0]) ;
	state[0] = std::min(.99999, state[0]) ;
	state[0] = std::max(0., state[0]) ;

}

void IndexedLinearDamage::artificialDamageStep(double d)
{
	state[0] = std::min(state[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE) ;
}


Matrix IndexedLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0. ;
	return ret*(1.-state[0]) ;
}


Matrix IndexedLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0. ;
	return ret*(1.-previousstate[0]) ;
}

bool IndexedLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return state[0] >= thresholdDamageDensity/fraction ;
}


IndexedLinearDamage::~IndexedLinearDamage()
{
}


}
