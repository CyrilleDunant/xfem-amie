
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
#include "nonlocalvonmises.h"
#include "../../mesher/delaunay.h"
namespace Mu {

NonLocalVonMises::NonLocalVonMises(double thresh, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(thresh)
{
	setMaterialCharacteristicRadius(radius);
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(ElementState &s)
{
	metInCompression = true ;
	metInTension = true ;
	std::pair<Vector, Vector> str( smoothedPrincipalStressAndStrain(s) ) ;
	double maxStress = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		maxStress = sqrt( ( ( str.first[0] - str.first[1] ) * ( str.first[0] - str.first[1] ) + str.first[0] * str.first[0] + str.first[1] * str.first[1] ) / 2. ) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		maxStress = sqrt( ( str.first[0] - str.first[1] ) * ( str.first[0] - str.first[1] ) + ( str.first[0] - str.first[2] ) * ( str.first[0] - str.first[2] ) + ( str.first[1] - str.first[2] ) * ( str.first[1] - str.first[2] ) ) / 6 ;
	}
	if(maxStress > threshold )
	{
		return 1. - threshold/maxStress ;
	}
	else 
	{
		return -1.+ maxStress/threshold;
	}

}

FractureCriterion * NonLocalVonMises::getCopy() const
{
	return new NonLocalVonMises(threshold, getMaterialCharacteristicRadius()) ;
}

Material NonLocalVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}

}
