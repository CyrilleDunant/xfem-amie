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

MCFT::MCFT( double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
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

	double tensionCritStrain = 2e6 / 37e9 ;
	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain / critStrain ;

	double mcftFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
	double maxCompression = -std::abs( downVal ) * mcftFactor ;

	if( mcftFactor > 1 || mcftFactor < 0 )
		maxCompression = -std::abs( downVal ) ;

	double maxTension = upVal ;

//	if(tstrain > tensionCritStrain)
//	{
	//Yamamoto model
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;

	//MCFT model
	maxTension = upVal / ( 1. + sqrt( 500.*tstrain ) ) ;

	//perfectly brittle
// 		maxTension = 0 ;
//	}



	metInCompression = cstrain <= 0 && std::abs( cstress / maxCompression ) > std::abs( tstress / maxTension ) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs( cstress / maxCompression ) < std::abs( tstress / maxTension ) || cstrain >= 0;


	std::vector<double> crits ;
	crits.push_back( -1 ) ;

	if( cstress <= 0 && std::abs( cstress ) >= std::abs( maxCompression ) ) {
		metInCompression = true ;
		crits.push_back( 1. - std::abs( maxCompression / cstress ) ) ;
	}

	if( tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) ) {
		metInTension = true ;
		crits.push_back( 1. - std::abs( maxTension / tstress ) ) ;
	}


	if( tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) ) {
		crits.push_back( -1. + std::abs( tstress / maxTension ) ) ;
	}

	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) ) {
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


NonLocalMCFT::NonLocalMCFT( double up, double down, double charRad, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
	physicalCharacteristicRadius = charRad ;
}


NonLocalMCFT::~NonLocalMCFT()
{
}

double NonLocalMCFT::grade( ElementState &s )
{

	Vector str( s.getPrincipalStressAtNodes() ) ;
	Vector stra( s.getPrincipalStrainAtNodes() ) ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL ) {
		str *= s.getParent()->area() ;
		stra *= s.getParent()->area() ;
		double fact = s.getParent()->area() ;

		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ ) {
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;

			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() ) {
				continue ;
			}

			double dc = squareDist2D( s.getParent()->getCenter(), ci->getCenter() ) ;

			if( ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured() ) {
				double d =  exp( -dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) );
				double a = ci->area() ;
				Vector pstress( ci->getState().getPrincipalStressAtNodes() ) ;
				Vector pstrain( ci->getState().getPrincipalStrainAtNodes() ) ;

				if( !ci->getBehaviour()->fractured() ) {
					str += pstress * a * d ;
					stra += pstrain * a * d ;
					fact += a * d ;

					if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius ) { // MIRROR_X
						str += pstress * a * d ;
						stra += pstrain * a * d ;
						fact += a * d ;
					}

					if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius ) { // MIRROR_Y
						str += pstress * a * d ;
						stra += pstrain * a * d ;
						fact += a * d ;
					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius ) { // MIRROR_XY
						str += pstress * a * d ;
						stra += pstrain * a * d ;
						fact += a * d ;
					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius ) { // MIRROR_XY
						str += pstress * a * d ;
						stra += pstrain * a * d ;
						fact += a * d ;
					}
				}
			}
		}

		str /= fact ;
		stra /=fact ;

	} else
		if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) {
			str *= s.getParent()->volume() ;
			stra *= s.getParent()->volume() ;
			double fact = s.getParent()->volume() ;

			// gaussian smooth
			for( size_t i = 0 ; i < cache.size() ; i++ ) {
				DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;

				if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() ) {
					continue ;
				}

				double dc = squareDist3D( s.getParent()->getCenter(), ci->getCenter() ) ;

				if( ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured() ) {
					double d =  exp( -dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) );
					double a = ci->volume() ;
					Vector pstress = ci->getState().getPrincipalStressAtNodes() ;
					Vector pstrain = ci->getState().getPrincipalStrainAtNodes() ;

					if( !ci->getBehaviour()->fractured() ) {
						str += pstress * a * d ;
						stra += pstrain * a * d ;
						fact += a * d ;

						if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius ) { // MIRROR_X
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius ) { // MIRROR_Y
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius ) { // MIRROR_Y
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;;
						}

						if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}

						if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius ) { // MIRROR_XY
							str += pstress * a * d ;
							stra += pstrain * a * d ;
							fact += a * d ;
						}
					}
				}
			}

			str /= fact ;
			stra /= fact ;
		}



	Vector pstrain = stra ;
	Vector pstress = str ;

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	double tensionCritStrain = 2e6 / 37e9 ;
	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain / critStrain ;

	double mcftFactor = ( 2.*renormCompressionStrain - renormCompressionStrain * renormCompressionStrain ) / ( 0.8 - 0.34 * tstrain / critStrain ) ;
	double maxCompression = -std::abs( downVal ) * mcftFactor ;

	if( mcftFactor > 1 || mcftFactor < 0 )
		maxCompression = -std::abs( downVal ) ;

	double maxTension = upVal ;

//	if(tstrain > tensionCritStrain)
//	{
	//Yamamoto model
//		maxTension = upVal/(1.+sqrt(2e6*(tstrain+tensionCritStrain))) ;

	//MCFT model
	maxTension = upVal / ( 1. + sqrt( 500.*tstrain ) ) ;

	//perfectly brittle
// 		maxTension = 0 ;
//	}



	metInCompression = cstrain <= 0 && std::abs( cstress / maxCompression ) > std::abs( tstress / maxTension ) || tstrain <= 0;
	metInTension = tstrain >= 0 && std::abs( cstress / maxCompression ) < std::abs( tstress / maxTension ) || cstrain >= 0;


	std::vector<double> crits ;
	crits.push_back( -1 ) ;

	if( cstress <= 0 && std::abs( cstress ) >= std::abs( maxCompression ) ) {
		metInCompression = true ;
		crits.push_back( 1. - std::abs( maxCompression / cstress ) ) ;
	}

	if( tstress >= 0 && std::abs( tstress ) >= std::abs( maxTension ) ) {
		metInTension = true ;
		crits.push_back( 1. - std::abs( maxTension / tstress ) ) ;
	}


	if( tstress >= 0 && std::abs( tstress )  < std::abs( maxTension ) ) {
		crits.push_back( -1. + std::abs( tstress / maxTension ) ) ;
	}

	if( cstress <= 0 && std::abs( cstress ) < std::abs( downVal ) ) {
		crits.push_back( -1. + std::abs( cstress / downVal ) ) ;
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
