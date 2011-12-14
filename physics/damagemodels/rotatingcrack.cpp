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
#include "damagemodel.h"
#include "../../elements/integrable_entity.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

RotatingCrack::RotatingCrack(): damages(0., 100)
{
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	isNull = false ;
	currentAngle = 0 ;
	currentDamage = 0 ;
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( ElementState &s)
{

	if(s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint())
	{
		currentAngle = s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle() ;
		if(currentAngle > M_PI*.5) 
			currentAngle -= M_PI ;
		if (currentAngle < -M_PI*.5)
			currentAngle += M_PI ;
		
		if(ceil(100.*(currentAngle+M_PI*.5)/M_PI - floor(100.*(currentAngle+M_PI*.5)/M_PI)-(currentAngle+M_PI*.5)/M_PI) > POINT_TOLERANCE_2D)
		{
			currentDamage = damages[std::max(std::min((int)ceil( 100.*(currentAngle+M_PI*.5)/M_PI)-1, (int)damages.size()-1), 0)]*((currentAngle+M_PI*.5)/M_PI-floor(100.*(currentAngle+M_PI*.5)/M_PI))
										+ damages[std::max(std::min((int)floor(100.*(currentAngle+M_PI*.5)/M_PI)-1, (int)damages.size()-1), 0)]*(ceil(100.*(currentAngle+M_PI*.5)/M_PI)-(currentAngle+M_PI*.5)/M_PI);
			currentDamage /= ceil(100.*(currentAngle+M_PI*.5)/M_PI - floor(100.*(currentAngle+M_PI*.5)/M_PI)-(currentAngle+M_PI*.5)/M_PI) ;
		}
		else
			currentDamage = damages[std::max(std::min((int)ceil( 100.*(currentAngle+M_PI*.5)/M_PI)-1, (int)damages.size()-1), 0)] ;
		
		state = currentDamage ;
	}
	
	return std::make_pair(state, Vector(1., 1)) ;
}

Matrix RotatingCrack::apply(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getState()[0]) ;
}


Matrix RotatingCrack::applyPrevious(const Matrix & m) const
{

	if(fractured())
		return m*0 ;
	
	return m*(1.-getPreviousState()[0]) ;
}

bool RotatingCrack::fractured() const 
{
// 	return false ;
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity ;
}

void RotatingCrack::postProcess()
{
	if(converged && getState()[0] > POINT_TOLERANCE_2D)
	{
		Vector damagesadded(0., damages.size()) ;
		
		for(double i = 25 ; i <  75 ; i++)
		{
			double cval = cos((M_PI*.5 - M_PI*(i/99.))*2.) ;
			cval*= getState()[0]-currentDamage ;

// 			double cval = (1.-std::abs(-1 + 2*(i/99.)))*getState()[0] ;
//			cval*= getState()[0]-currentDamage ;
			
// 				double cval = getState()[0]-currentDamage ;
			if(cval > 0)
				damagesadded[i] = cval ;
		}
		
		damagesadded.cshift((int)round((2.*currentAngle/M_PI)*100)) ;

		for(size_t i = 0 ; i <  damagesadded.size() ; i++)
			damages[i] = std::min(damagesadded[i]+damages[i], 1.) ;
		
// 		for(size_t i = 0 ; i <  damagesadded.size() ; i++)
// 			std::cout << damages[i] << "  "<< std::flush ;
// 		std::cout << std::endl ;
// 		exit(0) ;
	}
}

RotatingCrack::~RotatingCrack()
{
}



}
