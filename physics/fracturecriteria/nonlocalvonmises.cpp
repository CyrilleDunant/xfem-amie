
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
#include "../damagemodels/damagemodel.h"
namespace Amie {

NonLocalVonMises::NonLocalVonMises(double thresh, double E, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(std::abs(thresh)), E(E)
{
	setMaterialCharacteristicRadius(radius);
	met = false ;
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(ElementState &s)
{
	met = false ;
	Vector str( getSmoothedField(s, PRINCIPAL_REAL_STRESS_FIELD ) ) ;
	
	double maxStress = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{           
		maxStress = sqrt( ( ( str[0] - str[1] ) * ( str[0] - str[1] ) + str[0] * str[0] + str[1] * str[1] ) * .5 ) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		maxStress = sqrt( ( str[0] - str[1] ) * ( str[0] - str[1] ) + ( str[0] - str[2] ) * ( str[0] - str[2] ) + ( str[1] - str[2] ) * ( str[1] - str[2] ) ) / 6 ;
	}
	
	if( maxStress >= threshold )
	{
		met = true ;
		 return 1. - std::abs( threshold / maxStress );
	}
	
	return -1. + std::abs( maxStress / threshold );

}

FractureCriterion * NonLocalVonMises::getCopy() const
{
	return new NonLocalVonMises(threshold,E,  getMaterialCharacteristicRadius()) ;
}

Material NonLocalVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}

}
