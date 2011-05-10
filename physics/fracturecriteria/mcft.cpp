//
// C++ Implementation: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mcft.h"
#include "../damagemodels/damagemodel.h"
#include <set>

namespace Mu {

MCFT::MCFT(double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down)
{
}


MCFT::~MCFT()
{
}

double MCFT::grade(ElementState &s)
{
	Vector pstrain = s.getPrincipalStrains(s.getParent()->getCenter()) ;
	Vector pstress = s.getPrincipalStresses(s.getParent()->getCenter()) ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double tensionCritStrain = 2e6/37e9 ;
	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain/critStrain ;
	
	double mcftFactor = (2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)/(0.8-0.34*tstrain/critStrain) ;
	double maxCompression = -std::abs(downVal)*mcftFactor ;
	
	if(mcftFactor > 1 || mcftFactor < 0)
		maxCompression = -std::abs(downVal) ;
	
	double maxTension = upVal ;
	
//	if(tstrain > tensionCritStrain)
//	{
		//Yamamoto model 
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;
		
		//MCFT model 
 		maxTension = upVal/(1.+sqrt(500.*tstrain)) ;
		
		//perfectly brittle
// 		maxTension = 0 ;
//	}
	
	

	metInCompression = cstrain <= 0 && std::abs(cstress/maxCompression) > std::abs(tstress/maxTension) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs(cstress/maxCompression) < std::abs(tstress/maxTension) || cstrain >= 0;
	
	
	std::vector<double> crits ;
	crits.push_back(-1) ;
	
	if( cstress <= 0 && std::abs(cstress) >= std::abs(maxCompression) )
	{
		metInCompression = true ;
		crits.push_back(1. - std::abs(maxCompression/cstress)) ;
	}
	
	if(tstress >= 0 && std::abs(tstress) >= std::abs(maxTension))
	{
		metInTension = true ;
		crits.push_back(1. - std::abs(maxTension/tstress)) ;
	}
	
	
	if(tstress >= 0 && std::abs(tstress)  < std::abs(maxTension))
	{
		crits.push_back(-1. + std::abs(tstress/maxTension)) ;
	}

	if(cstress <= 0 && std::abs(cstress) < std::abs(downVal))
	{
		crits.push_back(-1. + std::abs(cstress/downVal)) ;
	}
	
	std::sort(crits.begin(), crits.end());
	return crits.back() ;
}

FractureCriterion * MCFT::getCopy() const
{
	return new MCFT(*this) ;
}

Material MCFT::toMaterial()
{
	Material mat ;
	return mat ;
}

}
