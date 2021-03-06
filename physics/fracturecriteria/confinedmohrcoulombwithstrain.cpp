//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "confinedmohrcoulombwithstrain.h"

namespace Amie {

ConfinedMohrCoulombWithStrainLimit::ConfinedMohrCoulombWithStrainLimit(double up, double down, double strainLimit) : 
	upVal(up), downVal(down), strainLimit(strainLimit)
{
	metInCompression = false ;
	metInTension = false ;
}


ConfinedMohrCoulombWithStrainLimit::~ConfinedMohrCoulombWithStrainLimit()
{
}

double ConfinedMohrCoulombWithStrainLimit::grade(ElementState &s)
{

	if(s.getParent()->getBehaviour()->fractured())
		return 0 ;
	
	Point p_ = s.getParent()->getCenter() ;
	Vector pstress(0., s.getParent()->spaceDimensions()) ;
	Vector pstrain(0., s.getParent()->spaceDimensions()) ;
	s.getField(PRINCIPAL_TOTAL_STRAIN_FIELD, PRINCIPAL_REAL_STRESS_FIELD, p_, pstrain, pstress, false) ;

	double maxStress = pstress.max();
	double minStress = pstress.min();
	
	if(pstrain.max() > strainLimit)
		return 1.-strainLimit/pstrain.max() ;
	
	metInCompression = false ;
	metInTension = false ;
	
	if(maxStress < 0 && minStress < 0) //we are in a confined situation
	{
		metInCompression = true ;
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			if( std::abs(minStress) >= std::abs(5.*downVal) )
			{
				return 1. - std::abs((5.*downVal)/minStress) ;
			}
			return -1. + std::abs(minStress/(5.*downVal)) ;
		}
		else
		{
			if( std::abs(minStress) >= std::abs(1.25*downVal) )
			{
				return 1. - std::abs((1.25*downVal)/minStress) ;
			}
			return -1. + std::abs(minStress/(1.25*downVal)) ;
		}
	}
	
	if( maxStress >= upVal )
	{
		metInTension = true ;
		if(minStress <= downVal)
			metInCompression = true ;
		return 1. - std::abs(upVal/maxStress) ;
	}
		
	if( minStress <= downVal )
	{
		metInCompression = true ;
		return 1. - std::abs(downVal/minStress) ;
	}
	
	double s0 = -1. + std::abs(maxStress/upVal);
	double s1 = -1. + std::abs(minStress/downVal) ;
	
	if(minStress > 0)
	{
		return s0 ;
	}
	
	if(maxStress < 0)
	{
		return s1 ;
	}
	
	
	if(std::abs(s0) > std::abs(s1))
		return s0 ;

	return s1;
}

FractureCriterion * ConfinedMohrCoulombWithStrainLimit::getCopy() const
{
	return new ConfinedMohrCoulombWithStrainLimit(*this) ;
}

}
