//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mohrcoulomb.h"

namespace Mu {

MohrCoulomb::MohrCoulomb(double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down)
{
}


MohrCoulomb::~MohrCoulomb()
{
}

double MohrCoulomb::grade(const ElementState &s) 
{

	if(s.getParent()->getBehaviour()->fractured())
		return 0 ;

	Vector pstress0 = s.getPrincipalStresses(Point(0, 0, 0), true) ;
	Vector pstress1 = s.getPrincipalStresses(Point(0, 1, 0), true) ;
	Vector pstress2 = s.getPrincipalStresses(Point(1, 0, 0), true) ;
	double maxStress = std::max(pstress0.max(), std::max(pstress1.max(),pstress2.max()));
	double minStress = std::min(pstress0.min(), std::min(pstress1.min(),pstress2.min()));
	
	if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		Vector pstress3 = s.getPrincipalStresses(Point(0, 0, 1), true) ;
		maxStress = std::max(pstress3.max(),maxStress);
                minStress = std::min(pstress3.min(), minStress);
	}
// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs(minStress/downVal) > std::abs(maxStress/upVal) ;
	metInTension = std::abs(minStress/downVal) < std::abs(maxStress/upVal) ;
	if( maxStress >= upVal )
	{
		metInTension = true;
		if( minStress <= downVal )
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

FractureCriterion * MohrCoulomb::getCopy() const
{
	return new MohrCoulomb(*this) ;
}

Material MohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}

}
