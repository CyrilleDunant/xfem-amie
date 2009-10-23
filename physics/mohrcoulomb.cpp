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
#include "../mesher/delaunay.h"

namespace Mu {

MohrCoulomb::MohrCoulomb(double up, double down)
	: upVal(up), downVal(down)
{
}


MohrCoulomb::~MohrCoulomb()
{
}

double MohrCoulomb::grade(const ElementState &s) const 
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
	if( maxStress >= upVal )
	{
		return 1. - std::abs(upVal/maxStress) ;
	}
		
	if( minStress <= downVal )
	{

		return 1. - std::abs(downVal/minStress) ;
	}

	return 0 ;
}

FractureCriterion * MohrCoulomb::getCopy() const
{
	return new MohrCoulomb(*this) ;
}

}
