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

double MCFT::grade(const ElementState &s) 
{
	Vector pstrain = s.getPrincipalStrains(s.getParent()->getCenter()) ;
	Vector pstress = s.getPrincipalStresses(s.getParent()->getCenter()) ;

	double tstrain = pstrain[0];
	double cstrain = pstrain[pstrain.size()-1];
	double tstress = pstress[0];
	double cstress = pstress[pstrain.size()-1];
// 	if(tstrain < 0)
// 		tstrain = 0 ;
// 	if(cstrain > 0)
// 		cstrain = 0 ;
// 	if(tstress < 0)
// 		tstress = 0 ;
// 	if(cstress > 0)
// 		cstress = 0 ;
	

	double critStrain = -0.002 ;//-0.002
	double renormCompressionStrain = -cstrain/critStrain ;
	
	double maxCompression = -std::abs(downVal)*(2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)/(0.8-0.34*tstrain/critStrain) ;
	
// 	if(std::abs((2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)/(0.8-0.34*tstrain/critStrain)) > 1)
// 		maxCompression = -std::abs(downVal) ;
	
	double maxTension = upVal ;
	if(tstrain > -critStrain )
	{
		//Yamamoto model 
// 		maxTension = upVal/(1.+sqrt(200000000.*(tstrain+critStrain))) ;
		
		//MCFT model 
		maxTension = upVal/(1.+sqrt(500.*tstrain)) ;
		
		//perfectly brittle
// 		maxTension = 0 ;
	}

	metInCompression = std::abs(cstress/maxCompression) > std::abs(tstress/maxTension) ;
	metInTension = std::abs(cstress/maxCompression) < std::abs(tstress/maxTension) ;
	
	
	std::vector<double> crits ;
	crits.push_back(-1) ;
	if( cstress < 0 && std::abs(cstress) >= std::abs(maxCompression) )
	{
		metInCompression = true ;
		crits.push_back(1. - std::abs(maxCompression/cstress)) ;
	}
	
	if(tstress > 0 && std::abs(tstress) >= std::abs(maxTension))
	{
		metInTension = true ;
		crits.push_back(1. - std::abs(maxTension/tstress)) ;
	}
	
	
	if(tstress > 0 && std::abs(tstress)  < std::abs(maxTension))
	{
		crits.push_back(-1. + std::abs(tstress/maxTension)) ;
	}

	if(cstress < 0 && std::abs(cstress) < std::abs(maxCompression))
	{
		crits.push_back(-1. + std::abs(cstress/maxCompression)) ;
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
