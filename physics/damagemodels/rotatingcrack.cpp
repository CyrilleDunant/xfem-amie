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
#include "rotatingcrack.h"

namespace Mu {

RotatingCrack::RotatingCrack(): damages(100, 0.)
{
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
	currentAngle = 0 ;
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( Mu::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

Matrix RotatingCrack::apply(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getState()[0]) ;
}


double RotatingCrack::getDamage()
{
	int index = round(currentAngle/M_PI*99);
	
}

Matrix RotatingCrack::applyPrevious(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getPreviousState()[0]) ;
}

bool RotatingCrack::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity ;
}


RotatingCrack::~RotatingCrack()
{
}



}
