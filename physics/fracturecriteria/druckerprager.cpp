
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
	double factor = 1 ;
	metInCompression = true ;
	metInTension = true ;
	std::pair<Vector, Vector> stressstrain( smoothedStressAndStrain(s, EFFECTIVE_STRESS) ) ;
	Vector str = stressstrain.first ;
	Vector stra = stressstrain.second ;
	double maxStress = 0 ;
	
	//hardening function from Jirasek et al.
	if(dynamic_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()))
	{
		PlasticStrain* ps = static_cast<PlasticStrain*>(s.getParent()->getBehaviour()->getDamageModel()) ;
		double kappa_0 = ps->kappa_0 ;
		double sigma_y = upthreshold ;
		Vector istrain = ps->imposedStrain*ps->getState()[0] ;
		double kappa_p = ps->plasticVariable + sqrt(2./3.)*sqrt(istrain[0]*istrain[0]+istrain[1]*istrain[1]+istrain[2]*istrain[2]) ;
		
		if(str.min() < .05*downthreshold || str.max() > .05*upthreshold || ps->plasticVariable > POINT_TOLERANCE_2D)
		{
			if(kappa_p < kappa_0 )
				factor = std::max((kappa_p*kappa_p-3.*kappa_p*kappa_0+3.*kappa_0*kappa_0)*kappa_p/(kappa_0*kappa_0*kappa_0),0.001) ;
			else
				factor = 1. ;
		}
// 		std::cout << kappa_0 << " vs " << kappa_p << " :  " <<  (kappa_p < kappa_0) << " : "<< factor << std::endl ;
		
	}
	
	

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double tr = str[0]+str[1] ;
		maxStress = (tr*friction + sqrt(0.5)*sqrt((str[0]-tr*.5)*(str[0]-tr*.5)+(str[1]-tr*.5)*(str[1]-tr*.5)+2.*(str[2])*(str[2])))*1.5 ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double tr = str[0]+str[1]+str[3] ;
		maxStress =  tr*friction - sqrt((str[0]-tr*.333333333333)*(str[0]-tr*.333333333333)+(str[1]-tr*.333333333333)*(str[1]-tr*.333333333333)+(str[2]-tr*.333333333333)*(str[2]-tr*.333333333333)+2.*(str[3])*(str[3])+2.*(str[4])*(str[4])+2.*(str[5])*(str[5])) ;
	}
	if(maxStress > upthreshold*factor && maxStress > 0)
	{
		return 1. - std::abs(factor*upthreshold/maxStress) ;
	}
	else if(maxStress >= 0)
	{
		return -1.+ std::abs(maxStress/(factor*upthreshold));
	}
	else if(maxStress < factor*downthreshold && maxStress < 0)
	{
		return 1. - std::abs(factor*downthreshold/maxStress) ;
	}
	else
		return -1.+std::abs( maxStress/(factor*downthreshold));
	
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
