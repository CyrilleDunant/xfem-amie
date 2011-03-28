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

FractionMCFT::FractionMCFT(double up, double down, Matrix concreteCGTensor, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down), concreteCGTensor(concreteCGTensor)
{
}


FractionMCFT::~FractionMCFT()
{
}

double FractionMCFT::grade(const ElementState &s) 
{
	
	
	if(s.getParent()->getBehaviour()->getDamageModel()->fractured() == true)
		return -1 ;
	Vector pstrain = -s.getPrincipalStrains(s.getParent()->getCenter()) ;
	
	Vector strains = s.getStrain(Point(1./3., 1./3.),true) ;
	Vector stresses = strains*concreteCGTensor;
	
	if(s.getParent()->getBehaviour()->hasInducedForces())
			stresses -= s.getParent()->getBehaviour()->getImposedStress(s.getParent()->getCenter()) ;
	Vector pstress(2) ;
	pstress[0] = -(stresses[0]+stresses[1])/2. - 
		sqrt(
			(stresses[0]-stresses[1])*(stresses[0]-stresses[1])/4. + 
			(stresses[2]*stresses[2])
			) ;
	pstress[1] = -(stresses[0]+stresses[1])/2. + 
		sqrt(
			(stresses[0]-stresses[1])*(stresses[0]-stresses[1])/4. + 
			(stresses[2]*stresses[2])
			) ;
		

	

	metInCompression = false ;
	metInTension = false ;
	
	double critStrain = -0.002 ;
	double renormCompressionStrain = pstrain.max()/critStrain ;
	
	double maxCompression = -(2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)*downVal/(0.8-0.34*pstrain.min()/critStrain) ;
	
	double maxTension = upVal/(1.+sqrt(200.*std::abs(pstrain.min()))) ;
	
	if( -pstress.min() >= maxTension )
	{
		metInTension = true ;
		if(-pstress.max() <= maxCompression)
		{
			metInCompression = true ;

			return std::max(1. - std::abs(maxTension/-pstress.min() ), 1. - std::abs(maxCompression/-pstress.max())) ;
		}

		return 1. - std::abs(maxTension/-pstress.min() ) ;
	}
		
	if( -pstress.max() <= maxCompression )
	{
		metInCompression = true ;

		return 1. - std::abs(maxCompression/-pstress.max()) ;
	}
	
	double s0 = -1. + std::abs(-pstress.min()/maxTension);
	double s1 = -1. + std::abs(-pstress.max()/maxCompression) ;
	
	if(-pstress.max() > 0)
	{
		return s0 ;
	}
	
	if(-pstress.min() < 0)
	{
		return s1 ;
	}
	
	
	if(std::abs(s0) > std::abs(s1))
	{

		return s0 ;
	}

	return s1;
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
