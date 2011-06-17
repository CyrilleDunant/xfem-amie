
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
	Vector str( smoothedPrincipalStress(s) ) ;
// 	Vector stra( smoothedPrincipalStrain(s) ) ;
	double maxStress = 0 ;
// 	double maxStrain = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
// 		maxStrain =  sqrt( 2. / 3 * ( stra[0] * stra[0] + stra[1] * stra[1] ) ) ;
		maxStress = sqrt( ( ( str[0] - str[1] ) * ( str[0] - str[1] ) + str[0] * str[0] + str[1] * str[1] ) / 2. ) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
// 		maxStrain = sqrt( 2. / 3 * ( stra[0] * stra[0] + stra[1] * stra[1] + stra[2] * stra[2] ) ) ;
		maxStress = sqrt( ( str[0] - str[1] ) * ( str[0] - str[1] ) + ( str[0] - str[2] ) * ( str[0] - str[2] ) + ( str[1] - str[2] ) * ( str[1] - str[2] ) ) / 6 ;
	}
	
// 	if(maxStrain < POINT_TOLERANCE_2D)
// 		return -1 ;
// 	double modulus =std::abs( maxStress/maxStrain) ;
// 	double straincrit = threshold/modulus ;
// 	double c = sqrt(straincrit*straincrit*modulus*modulus+straincrit*straincrit) ;
	
// 	if(c < POINT_TOLERANCE_2D)
// 		return -1 ;
	
// 	double ss = sqrt(maxStress*maxStress+maxStrain*maxStrain) ;
	
// 	return std::abs(maxStrain/threshold)-1.;
	

	if(maxStress > threshold )
	{
		return 1. - threshold/maxStress ;
	}
	else 
	{
		return -1.+ maxStress/threshold;
	}

// 	if(ss > c )
// 	{
// 		return 1. - std::abs( c / ss )  ;
// 	}
// 	else 
// 	{
// 		return -1.+ std::abs( ss / c );
// 	}
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
