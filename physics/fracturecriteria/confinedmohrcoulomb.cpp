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
#include "confinedmohrcoulomb.h"

namespace Mu {

ConfinedMohrCoulomb::ConfinedMohrCoulomb(double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down)
{
}


ConfinedMohrCoulomb::~ConfinedMohrCoulomb()
{
}

double ConfinedMohrCoulomb::grade(const ElementState &s) 
{

	if(s.getParent()->getBehaviour()->fractured())
		return 0 ;

	Vector pstress0 = s.getPrincipalStresses(Point(0, 0, 0), true) ;
	Vector pstress1 = s.getPrincipalStresses(Point(0, 1, 0), true) ;
	Vector pstress2 = s.getPrincipalStresses(Point(1, 0, 0), true) ;
	double maxStress = std::max(pstress0.max(), std::max(pstress1.max(),pstress2.max()));
	double minStress = std::min(pstress0.min(), std::min(pstress1.min(),pstress2.min()));
	
	metInCompression = false ;
	metInTension = false ;
	if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		Vector pstress3 = s.getPrincipalStresses(Point(0, 0, 1), true) ;
		maxStress = std::max(pstress3.max(),maxStress);
        minStress = std::min(pstress3.min(), minStress);
	}
	
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

FractureCriterion * ConfinedMohrCoulomb::getCopy() const
{
	return new ConfinedMohrCoulomb(*this) ;
}

Material ConfinedMohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}

}
