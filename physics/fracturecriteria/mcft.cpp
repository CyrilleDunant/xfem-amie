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
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include <set>

namespace Mu
{

MCFT::MCFT( double up, double down, double youngModulus, double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	physicalCharacteristicRadius = charRad ;
	tensionCritStrain = up / youngModulus ;
}


MCFT::~MCFT()
{
}

double MCFT::grade( ElementState &s )
{
	Vector pstrain = s.getPrincipalStrains( s.getParent()->getCenter() ) ;
	Vector pstress = s.getPrincipalStresses( s.getParent()->getCenter() ) ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain / critStrain ;
	double mcftFactor = 1 ;
	double maxCompression = -std::abs( downVal ) ;
	if(cstrain < critStrain)
	{
		mcftFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
		double maxCompression = -std::abs( downVal ) * mcftFactor ;
	}

	if( mcftFactor > 1 || mcftFactor < 0 )
		maxCompression = -std::abs( downVal ) ;

	double maxTension = upVal ;

	if(tstrain > tensionCritStrain)
	{
	//Yamamoto model
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;

	//MCFT model
	maxTension = upVal / ( 1. + sqrt( 200.*tstrain ) ) ;

	//perfectly brittle
// 		maxTension = 0 ;
	}



	metInCompression = cstrain <= 0 && std::abs( cstress / maxCompression ) > std::abs( tstress / maxTension ) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs( cstress / maxCompression ) < std::abs( tstress / maxTension ) || cstrain >= 0;


	std::vector<double> crits ;
	crits.push_back( -1 ) ;

	if( cstress <= 0 && std::abs( cstress ) >= std::abs( maxCompression ) )
	{
		metInCompression = true ;
		crits.push_back( 1. - std::abs( maxCompression / cstress ) ) ;
	}

	if( tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) )
	{
		metInTension = true ;
		crits.push_back( 1. - std::abs( maxTension / tstress ) ) ;
	}


	if( tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) )
	{
		crits.push_back( -1. + std::abs( tstress / maxTension ) ) ;
	}

	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) )
	{
		crits.push_back( -1. + std::abs( cstress / downVal ) ) ;
	}

	std::sort( crits.begin(), crits.end() );
	return crits.back() ;
}

FractureCriterion *MCFT::getCopy() const
{
	return new MCFT( *this ) ;
}

Material MCFT::toMaterial()
{
	Material mat ;
	return mat ;
}


NonLocalMCFT::NonLocalMCFT( double up, double down, double youngModulus,  double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	physicalCharacteristicRadius = charRad ;
	tensionCritStrain = up/ youngModulus ;
	strainBroken = false ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::grade( ElementState &s )
{
	
	Vector pstrain(smoothedPrincipalStrain(s)) ;
// 	Vector localpstrain = s.getPrincipalStrains(s.getParent()->getCenter()) ;
	Vector pstress(smoothedPrincipalStress(s)) ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
// 	double localtstrain = localpstrain.max();
// 	double localcstrain = localpstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();
;
	double critStrain = -0.002 ;
	double mcftFactor = 1 ;
// 	double renormCompressionStrain = cstrain / critStrain ;
	double maxCompression = -std::abs( downVal ) ;

	double tensionModulus  = std::abs(tstress/tstrain) ;
	double compressionModulus = std::abs(cstress/cstrain) ;
	
	double maxTension = upVal ;
	
	
	double tstraincrit = 0 ; 
	double prevtstraincrit = tstraincrit ;
	do
	{
		prevtstraincrit = tstraincrit ;
		if(tstraincrit > tensionCritStrain|| (strainBroken && tstrain > 0))
		{
			tstraincrit = upVal / (( 1. + sqrt( 200.*tstraincrit ) )*tensionModulus) ;
		}
		else
		{
			tstraincrit = tensionCritStrain ;
		}
		
	} while (std::abs(prevtstraincrit - tstraincrit) > POINT_TOLERANCE_2D) ;
	
	maxTension = tstraincrit*tensionModulus ;
	
// 	if(tstrain > tensionCritStrain || (strainBroken && tstrain > 0))
// 	{
// 		maxTension = upVal / ( 1. + sqrt( 200.*tstrain ) ) ;
		
// 	}
	if(tstrain > tstraincrit)
		strainBroken = true ;
	
	
	
	double cstraincrit = 0 ; 
	double prevtcstraincrit = cstraincrit ;
	do
	{
		prevtcstraincrit = cstraincrit ;
		if(cstraincrit < critStrain)
		{
			double renormCompressionStrain = cstraincrit / critStrain ;
			cstraincrit = -std::abs( downVal ) * ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / (( 0.8 - 0.34 * tstraincrit / critStrain )*compressionModulus) ;
		}
		else
		{
			cstraincrit = critStrain ;
		}
		
	} while (std::abs(prevtcstraincrit - cstraincrit) > POINT_TOLERANCE_2D) ;
	maxCompression = cstraincrit*compressionModulus ;
	
// 	if(cstrain < critStrain)
// 	{
// 		mcftFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
// 		maxCompression = -std::abs( downVal ) * mcftFactor ;
// 	}


// 	if( mcftFactor > 1 || mcftFactor < 0 )
// 		maxCompression = -std::abs( downVal ) ;


	metInCompression = cstrain <= 0 && std::abs( cstress / maxCompression ) > std::abs( tstress / maxTension ) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs( cstress / maxCompression ) < std::abs( tstress / maxTension ) || cstrain >= 0;


	std::vector<double> crits ;
	crits.push_back( -1 ) ;

	if( cstress <= 0 && std::abs( cstress ) >= std::abs( maxCompression ) )
	{
		metInCompression = true ;
// 		crits.push_back( 1. - std::abs( maxCompression / cstress ) ) ;
		double c = sqrt(cstraincrit*cstraincrit*compressionModulus*compressionModulus+cstraincrit*cstraincrit) ;
		double s = sqrt(cstress*cstress+cstrain*cstrain) ;
		crits.push_back( 1. - std::abs( c / s ) ) ;
	}
	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) )
	{
		double c = sqrt(cstraincrit*cstraincrit*compressionModulus*compressionModulus+cstraincrit*cstraincrit) ;
		double s = sqrt(cstress*cstress+cstrain*cstrain) ;
		crits.push_back( -1. + std::abs( s / c ) ) ;
	}

	if( strainBroken && tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) )
	{
		metInTension = true ;
		double c = sqrt(tstraincrit*tstraincrit*tensionModulus*tensionModulus+tstraincrit*tstraincrit) ;
		double s = sqrt(tstress*tstress+tstrain*tstrain) ;
		crits.push_back( 1. - std::abs( c / s ) ) ;
	}

	if( tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) )
	{
		double c = sqrt(tstraincrit*tstraincrit*tensionModulus*tensionModulus+tstraincrit*tstraincrit) ;
		double s = sqrt(tstress*tstress+tstrain*tstrain) ;
		crits.push_back( -1. + std::abs( s / c ) ) ;
	}



	std::sort( crits.begin(), crits.end() );

	return crits.back() ;
}

FractureCriterion *NonLocalMCFT::getCopy() const
{
	return new NonLocalMCFT( *this ) ;
}

Material NonLocalMCFT::toMaterial()
{
	Material mat ;
	return mat ;
}

}
