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
	Vector pstrain = -s.getPrincipalStrains(s.getParent()->getCenter()) ;
	Vector pstress = -s.getPrincipalStresses(s.getParent()->getCenter()) ;

	metInCompression = false ;
	metInTension = false ;
	
	double critStrain = -0.002 ;
	double renormCompressionStrain = pstrain.max()/critStrain ;
	
	double maxCompression = -(2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)*downVal/(0.8-0.34*pstrain.min()/critStrain) ;
	
	double maxTension = upVal/(1.+sqrt(200.*std::abs(pstrain.min()))) ;
	
	if( -pstress.min() >= maxTension && pstrain.min() < critStrain)
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
		return s0 ;

	return s1;
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
