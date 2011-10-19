
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "druckerprager.h"
#include "../../mesher/delaunay.h"
namespace Mu {

DruckerPrager::DruckerPrager(double thresh, double friction, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(thresh), friction(friction)
{
	setMaterialCharacteristicRadius(radius);
}


DruckerPrager::~DruckerPrager()
{
}

double DruckerPrager::grade(ElementState &s)
{
	metInCompression = true ;
	metInTension = true ;
	std::pair<Vector, Vector> str( smoothedPrincipalStressAndStrain(s) ) ;
	double maxStress = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double tr = str.first[0]+str.first[1] ;
		maxStress = tr*friction + sqrt(0.5)*sqrt((str.first[0]-tr*.5)*(str.first[0]-tr*.5)+(str.first[1]-tr*.5)*(str.first[1]-tr*.5)) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double tr = str.first[0]+str.first[1]+str.first[3] ;
		maxStress =  tr*friction + sqrt(0.5)*sqrt((str.first[0]-tr/3.)*(str.first[0]-tr/3.)+(str.first[1]-tr/3.)*(str.first[1]-tr/3.)+(str.first[2]-tr/3.)*(str.first[2]-tr/3.)) ;
	}
	if(std::abs(maxStress) > std::abs(threshold) )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}

}

FractureCriterion * DruckerPrager::getCopy() const
{
	return new DruckerPrager(threshold, friction, getMaterialCharacteristicRadius()) ;
}

Material DruckerPrager::toMaterial()
{
	Material mat ;
	return mat ;
}

}
