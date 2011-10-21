
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
#include "../damagemodels/plasticstrain.h"
namespace Mu {

DruckerPrager::DruckerPrager(double downthres,double upthres, double friction, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upthreshold(upthres), downthreshold(downthres), friction(friction)
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
	Vector str( smoothedStress(s, EFFECTIVE_STRESS) ) ;
	double maxStress = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double tr = str[0]+str[1] ;
		maxStress = tr*friction + sqrt(0.5)*sqrt((str[0]-tr*.5)*(str[0]-tr*.5)+(str[1]-tr*.5)*(str[1]-tr*.5)+2.*(str[2])*(str[2])) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double tr = str[0]+str[1]+str[3] ;
		maxStress =  tr*friction - sqrt((str[0]-tr*.333333333333)*(str[0]-tr*.333333333333)+(str[1]-tr*.333333333333)*(str[1]-tr*.333333333333)+(str[2]-tr*.333333333333)*(str[2]-tr*.333333333333)+2.*(str[3])*(str[3])+2.*(str[4])*(str[4])+2.*(str[5])*(str[5])) ;
	}
	if(maxStress > upthreshold && maxStress > 0)
	{
		return 1. - upthreshold/maxStress ;
	}
	else if(maxStress > 0)
	{
		return -1.+ maxStress/downthreshold;
	}
	else if(maxStress < downthreshold && maxStress <= 0)
	{
		return 1. - downthreshold/maxStress ;
	}
	else
		return -1.+ maxStress/(downthreshold);
	
	return -1 ;

}

FractureCriterion * DruckerPrager::getCopy() const
{
	return new DruckerPrager(downthreshold, upthreshold, friction, getMaterialCharacteristicRadius()) ;
}

Material DruckerPrager::toMaterial()
{
	Material mat ;
	return mat ;
}

}
