
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

DruckerPrager::DruckerPrager(double downthres,double upthres, double modulus,  double friction, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), upthreshold(upthres), downthreshold(downthres), friction(friction), modulus(modulus)
{
	setMaterialCharacteristicRadius(radius);
	metInCompression = false ;
	metInTension = false ;
	inTension = true ;
}


DruckerPrager::~DruckerPrager()
{
}

double DruckerPrager::grade(ElementState &s)
{
	double factor = 1 ;
	std::pair<Vector, Vector> stressstrain( smoothedPrincipalStressAndStrain(s, REAL_STRESS) ) ;
	Vector str = stressstrain.first ;
	Vector stra = stressstrain.second ;
	double maxStress = 0 ;
	double maxStrain = 0 ;
	double pseudomodulus = modulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	//hardening function from Jirasek et al.
	if(dynamic_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()))
	{
		pseudomodulus = modulus ;
	
		PlasticStrain* ps = static_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()) ;
		
		double kappa_0 = ps->kappa_0 ;
		double kappa_p = ps->getPlasticity() ;
		if(kappa_0 > 0)
		{
			(kappa_p < kappa_0 )?factor = ((kappa_p)*(kappa_p)-3.*(kappa_p)*kappa_0+3.*kappa_0*kappa_0)*(kappa_p)/(kappa_0*kappa_0*kappa_0):factor = 1. ;
		}
		factor *= 1.-ps->getDamage() ;
		pseudomodulus = modulus*(1.-ps->getDamage()) ;
	}

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double tr = str[0]+str[1] ;
		maxStress = tr*friction - sqrt(0.5)*sqrt((str[0]-tr*.5)*(str[0]-tr*.5)+(str[1]-tr*.5)*(str[1]-tr*.5)+2.*(str[0]-str[1])*(str[0]-str[1])) ;
		maxStrain = maxStress/(modulus*factor) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double tr = str[0]+str[1]+str[3] ;
		maxStress =  tr*friction - sqrt((str[0]-tr*.333333333333)*(str[0]-tr*.333333333333)+(str[1]-tr*.333333333333)*(str[1]-tr*.333333333333)+(str[2]-tr*.333333333333)*(str[2]-tr*.333333333333)+2.*(str[3])*(str[3])+2.*(str[4])*(str[4])+2.*(str[5])*(str[5])) ;
		maxStrain = maxStress/(modulus*factor) ;
	}

// 	if(maxStress > upthreshold && maxStress > POINT_TOLERANCE_2D)
// 	{
// // 		std::cout << factor << ", "<< maxStress<< std::endl ;
// 		metInTension = true ;
// 		metInCompression = false ;
// 		inTension = true ;
// 		return 1. - std::abs(upthreshold/maxStress) ;
// 	}
// 	else if(maxStress >= 0 && std::abs(upthreshold) > POINT_TOLERANCE_2D)
// 	{
// 		metInTension = false ;
// 		metInCompression = false ;
// 		inTension = true ;
// 		return -1.+ std::abs(maxStress/upthreshold) ;
// 	}
// 	else if(maxStress < downthreshold && maxStress < -POINT_TOLERANCE_2D)
// 	{
// 		metInTension = false ;
// 		metInCompression = true ;
// 		inTension = false ;
// 		return 1. - std::abs(downthreshold/maxStress) ;
// 	}
// 	else if(std::abs(downthreshold) > POINT_TOLERANCE_2D)
// 	{
// 		metInTension = false ;
// 		metInCompression = false ;
// 		inTension = false ;
// 		return -1.+std::abs( maxStress/downthreshold) ;
// 	}

	
	double effectiveUp = upthreshold*factor/pseudomodulus ;
	double effectiveDown = downthreshold*factor/pseudomodulus ;
	if(maxStrain > effectiveUp && maxStrain > POINT_TOLERANCE_2D)
	{
// 		std::cout << factor << ", "<< maxStress<< std::endl ;
		metInTension = true ;
		metInCompression = false ;
		inTension = true ;
		return 1. - std::abs(effectiveUp/maxStrain) ;
	}
	else if(maxStrain >= 0 && std::abs(effectiveUp) > POINT_TOLERANCE_2D)
	{
		metInTension = false ;
		metInCompression = false ;
		inTension = true ;
		return -1.+ std::abs(maxStrain/effectiveUp) ;
	}
	else if(maxStrain < effectiveDown && maxStrain < -POINT_TOLERANCE_2D)
	{
		metInTension = false ;
		metInCompression = true ;
		inTension = false ;
		return 1. - std::abs(effectiveDown/maxStrain) ;
	}
	else if(std::abs(effectiveDown) > POINT_TOLERANCE_2D)
	{
		metInTension = false ;
		metInCompression = false ;
		inTension = false ;
		return -1.+std::abs( maxStrain/effectiveDown) ;
	}
	metInTension = false ;
	metInCompression = false ;
	inTension = false ;
	
	return -1 ;

}

FractureCriterion * DruckerPrager::getCopy() const
{
	return new DruckerPrager(downthreshold, upthreshold, modulus, friction, getMaterialCharacteristicRadius()) ;
}

Material DruckerPrager::toMaterial()
{
	Material mat ;
	return mat ;
}

}
