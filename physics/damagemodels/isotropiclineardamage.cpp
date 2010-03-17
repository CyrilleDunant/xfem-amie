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
	double maxD = .9999 ; //1- 0.0000005/s.getParent()->area() ;
// 	Vector pstrain = s.getPrincipalStresses(s.getParent()->getCenter()) ;
// 	Vector strain = s.getStrain(s.getParent()->getCenter()) ;
	
// 	bool inCompression = pstrain.min() < 0 && std::abs(pstrain.min()) > pstrain.max() ;
// 	
// 	if(inCompression)
// 		maxD = -.9999 ;
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		state[0] += 0.1 ; //5e-5*maxD/sqrt(s.getParent()->area()) ;
	else
		state[0] += 0.1 ; //5e-5*maxD/sqrt(s.getParent()->volume()) ;
// 	std::cout << 1e-5*maxD/sqrt(s.getParent()->area()) << std::endl ;
// 	if(!inCompression)
		state[0] = std::min(maxD, state[0]) ;
// 	else
// 		state[0] = std::max(maxD, state[0]) ;
// 	if(s.getStrain(s.getParent()->getBoundingPoints()).max() > .5)
// 		state[0] = maxD ;

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
