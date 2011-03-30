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
#include "fractionmcft.h"
#include "../damagemodels/damagemodel.h"

namespace Mu {

FractionMCFT::FractionMCFT(double up, double down, Matrix concreteCGTensor, double phi, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down), concreteCGTensor(concreteCGTensor), phi(phi)
{
}


FractionMCFT::~FractionMCFT()
{
}

double FractionMCFT::grade(const ElementState &s) 
{
	
	
	if(s.getParent()->getBehaviour()->getDamageModel()->fractured() == true)
		return -1 ;
	
	Vector pstrain = s.getPrincipalStrains(s.getParent()->getCenter()) ;
	Vector strains = s.getStrain(Point(1./3., 1./3.),true) ;
	Vector stresses = strains*concreteCGTensor*(1.-phi); //assume stresses are scaled by the fraction of steel
	
	if(s.getParent()->getBehaviour()->hasInducedForces())
			stresses -= s.getParent()->getBehaviour()->getImposedStress(s.getParent()->getCenter()) ;
	Vector pstress(2) ;
	pstress[0] = 0.5*(stresses[0]+stresses[1]) - 
		0.5*sqrt(
			(stresses[0]-stresses[1])*(stresses[0]-stresses[1]) + 
			(stresses[2]*stresses[2])
			) ;
	pstress[1] = 0.5*(stresses[0]+stresses[1]) + 
		0.5*sqrt(
			(stresses[0]-stresses[1])*(stresses[0]-stresses[1]) + 
			(stresses[2]*stresses[2])
			) ;
		

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	metInCompression = false ;
	metInTension = false ;
	
	double critStrain = -0.002 ;//-0.002
	double renormCompressionStrain = cstrain/critStrain ;
	
	double mcftFactor = (2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)/(0.8-0.34*tstrain/critStrain) ;
	double maxCompression = -std::abs(downVal)*mcftFactor ;
	
	if(mcftFactor > 1 || mcftFactor < 0)
		maxCompression = -std::abs(downVal) ;
	
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

	if(cstress < 0 && std::abs(cstress) < std::abs(downVal))
	{
		crits.push_back(-1. + std::abs(cstress/downVal)) ;
	}
	
	std::sort(crits.begin(), crits.end());
	return crits.back() ;
}

FractureCriterion * FractionMCFT::getCopy() const
{
	return new FractionMCFT(*this) ;
}

Material FractionMCFT::toMaterial()
{
	Material mat ;
	return mat ;
}

}
