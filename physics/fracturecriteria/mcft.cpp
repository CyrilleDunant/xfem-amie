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
	double compressionFactor = 1 ;
	double maxCompression = -std::abs( downVal ) ;
	if(cstrain < critStrain)
	{
		compressionFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
		double maxCompression = -std::abs( downVal ) * compressionFactor ;
	}

	if( compressionFactor > 1 || compressionFactor < 0 )
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
	tensionCritStrain = up / youngModulus ;
// 	upVal *= 1.87 ;
	strainBroken = false ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::grade( ElementState &s )
{
	
	Vector pstrain(smoothedPrincipalStrain(s)) ;
	Vector pstress(smoothedPrincipalStress(s)) ;
	
// 	if(s.getParent()->getBehaviour()->getDamageModel()->getState().max() < POINT_TOLERANCE_2D && pstrain.max() > tensionCritStrain && pstress.max() > POINT_TOLERANCE_2D)
// 	{
// 		std::cout << "before : " << tensionCritStrain << std::flush ;
// 		tensionCritStrain = upVal / (pstress.max()/pstrain.max()) ;
// 		std::cout << " ; after : " << tensionCritStrain << std::endl ;
// 		exit(0) ;
// 	}
// 	

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain / critStrain ;
	double compressionFactor = 1 ;
	double maxCompression = -std::abs( downVal ) ;
	if(cstrain < critStrain)
	{
// 		compressionFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
// 		double maxCompression = -std::abs( downVal ) * compressionFactor ;
		double n = 0.8 + downVal/17e6 ;
		double k = std::max(0.67 + downVal/62e6, 1.) ;
		compressionFactor = n*renormCompressionStrain / (n - 1. + pow(renormCompressionStrain, n*k)) ;
	}
	maxCompression *= compressionFactor ;
	if( compressionFactor > 1 || compressionFactor < 0 )
		maxCompression = -std::abs( downVal ) ;

	
	double maxTension = upVal;
	
	if(tstrain >= tensionCritStrain)
	{
	//Yamamoto model
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;

	//MCFT model
		maxTension = upVal / ( 1. + sqrt( 200.*tstrain ) ) ;
// 		if(s.getParent()->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE_2D && !strainBroken)
// 		{
// 			strainBroken = true ;
// 			tensionCritStrain/= 1.87 ;
// 			upVal /= 1.87 ;
// 		}
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

	if( tstrain >= tensionCritStrain && tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) )
	{
		metInTension = true ;
		crits.push_back( 1. - std::abs( maxTension / tstress ) ) ;
	}


	if( tstrain >= tensionCritStrain && tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) )
	{
		crits.push_back( -1. + std::abs( tstress / maxTension ) ) ;
	}
	
	if(tstrain < tensionCritStrain && tstrain > 0)
	{
		crits.push_back( -1. + std::abs( tstrain / tensionCritStrain ) ) ;
	}

	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) )
	{
		crits.push_back( -1. + std::abs( cstress / downVal ) ) ;
	}

	std::sort( crits.begin(), crits.end() );
// 	if(crits.back() < -.999 &&  tstrain > 0)
// 		std::cout << " + "<< tstrain << "  " << std::abs( tstrain / tensionCritStrain ) << std::endl ;
// 	if(crits.back() < -.999 &&  cstrain < 0)
// 		std::cout << " + "<< cstrain << "  " << std::abs( cstrain / critStrain ) << std::endl ;
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
