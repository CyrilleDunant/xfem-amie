
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
namespace Amie {

DruckerPrager::DruckerPrager(double downthres,double upthres, double modulus,  double friction, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z), upthreshold(upthres), downthreshold(downthres), friction(friction), modulus(modulus)
{
	setMaterialCharacteristicRadius(radius);
	metInCompression = false ;
	metInTension = false ;
	inTension = true ;
	cap = 1 ;
	smoothingType = QUARTIC_COMPACT ;
}


DruckerPrager::~DruckerPrager()
{
}

double DruckerPrager::grade(ElementState &s)
{
	double factor = 1 ;
	
	Vector stra = getSmoothedField(PRINCIPAL_STRAIN_FIELD, s) ; //+((double)random()/RAND_MAX*2.-1.)*.0001*stressstrain.first ;
	Vector str =  (stra-s.getParent()->getBehaviour()->getImposedStrain(Point(1./3., 1./3.)))*s.getParent()->getBehaviour()->getTensor(Point(1./3., 1./3.)) ; //+((double)random()/RAND_MAX*2.-1.)*.0001*stressstrain.second ;
	double maxStress = 0 ;
	double maxStrain = 0 ;
	double pseudomodulus = modulus*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	double dfactor = 1 ;
	//hardening function from Jirasek et al.
	if(dynamic_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()))
	{
		pseudomodulus = modulus ;
	
		PlasticStrain* ps = static_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()) ;
		
// 		double kappa_0 = ps->kappa_0 ;
// 		double kappa_p = ps->getPlasticity() ;
// 		if(kappa_0 > 0 )
// 		{
// 			(kappa_p <= kappa_0 )?factor = ((kappa_p)*(kappa_p)-3.*(kappa_p)*kappa_0+3.*kappa_0*kappa_0)*(kappa_p)/(kappa_0*kappa_0*kappa_0):factor = 1. ;
// 		}
// // 		std::cout << kappa_p << " , " << kappa_0 << " :: " << std::flush ;
// 
// 		factor = std::min(factor, factor*cap) ;
		dfactor = 1.-ps->getDamage() ;
		pseudomodulus = modulus*dfactor ;
	}
	
	if(pseudomodulus < POINT_TOLERANCE)
	  return -1 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double tr = str[0]+str[1] ;
		maxStress = tr*friction - sqrt(0.5)*sqrt((str[0]-tr/3.)*(str[0]-tr/3.)+(str[1]-tr/3.)*(str[1]-tr/3.)+2.*(str[0]-str[1])*(str[0]-str[1])) ;
		maxStrain = maxStress/pseudomodulus ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double tr = str[0]+str[1]+str[3] ;
		maxStress =  tr*friction - sqrt((str[0]-tr*.333333333333)*(str[0]-tr*.333333333333)+(str[1]-tr*.333333333333)*(str[1]-tr*.333333333333)+(str[2]-tr*.333333333333)*(str[2]-tr*.333333333333)+2.*(str[3])*(str[3])+2.*(str[4])*(str[4])+2.*(str[5])*(str[5])) ;
		maxStrain = maxStress/pseudomodulus ;
	}

// 	if(maxStress > upthreshold && maxStress > POINT_TOLERANCE)
// 	{
// // 		std::cout << factor << ", "<< maxStress<< std::endl ;
// 		metInTension = true ;
// 		metInCompression = false ;
// 		inTension = true ;
// 		return 1. - std::abs(upthreshold/maxStress) ;
// 	}
// 	else if(maxStress >= 0 && std::abs(upthreshold) > POINT_TOLERANCE)
// 	{
// 		metInTension = false ;
// 		metInCompression = false ;
// 		inTension = true ;
// 		return -1.+ std::abs(maxStress/upthreshold) ;
// 	}
// 	else if(maxStress < downthreshold && maxStress < -POINT_TOLERANCE)
// 	{
// 		metInTension = false ;
// 		metInCompression = true ;
// 		inTension = false ;
// 		return 1. - std::abs(downthreshold/maxStress) ;
// 	}
// 	else if(std::abs(downthreshold) > POINT_TOLERANCE)
// 	{
// 		metInTension = false ;
// 		metInCompression = false ;
// 		inTension = false ;
// 		return -1.+std::abs( maxStress/downthreshold) ;
// 	}

// 	std::cout << factor << ", "<< dfactor << ", "<< maxStress<< std::endl ;

	
	double effectiveUp = upthreshold*factor*dfactor ;
	double effectiveDown = downthreshold*factor*dfactor ;
	if(maxStress > effectiveUp && maxStrain > POINT_TOLERANCE)
	{
// 		std::cout << factor << ", "<< maxStress<< std::endl ;
		metInTension = true ;
		metInCompression = false ;
		inTension = true ;
		return 1. - std::abs(effectiveUp/maxStress) ;
	}
	else if(maxStress >= 0 && std::abs(effectiveUp) > POINT_TOLERANCE)
	{
		metInTension = false ;
		metInCompression = false ;
		inTension = true ;
		return -1.+ std::abs(maxStress/effectiveUp) ;
	}
	else if(maxStress < effectiveDown && maxStress < -POINT_TOLERANCE)
	{
		metInTension = false ;
		metInCompression = true ;
		inTension = false ;
		return 1. - std::abs(effectiveDown/maxStress) ;
	}
	else if(std::abs(effectiveDown) > POINT_TOLERANCE)
	{
		metInTension = false ;
		metInCompression = false ;
		inTension = false ;
		return -1.+std::abs( maxStress/effectiveDown) ;
	}
	metInTension = false ;
	metInCompression = false ;
	inTension = false ;
	
	return -1 ;

}

FractureCriterion * DruckerPrager::getCopy() const
{
	DruckerPrager * ret = new DruckerPrager(downthreshold, upthreshold, modulus, friction, getMaterialCharacteristicRadius()) ;
	ret->cap = cap ;
	ret->copyEssentialParameters( this ) ;
	return ret  ;
}

}
