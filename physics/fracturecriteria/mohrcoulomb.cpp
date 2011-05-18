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
#include "mohrcoulomb.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"

namespace Mu
{

MohrCoulomb::MohrCoulomb( double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
}


MohrCoulomb::~MohrCoulomb()
{
	for( size_t i = 0 ; i < testPoints.size() ; i++ )
		delete testPoints[i] ;
}

double MohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return 0 ;

	if( testPoints.size() == 0 ) {
		if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL ) {
			testPoints.resize( 4 );
			testPoints[0] = new Point( 0, 0, 0 ) ;
			testPoints[1] = new Point( 1, 0, 0 ) ;
			testPoints[2] = new Point( 0, 1, 0 ) ;
			testPoints[3] = new Point( 0, 0, 1 ) ;
		} else {
			testPoints.resize( 3 );
			testPoints[0] = new Point( 0, 0, 0 ) ;
			testPoints[1] = new Point( 1, 0, 0 ) ;
			testPoints[2] = new Point( 0, 1, 0 ) ;
		}
	}

	Vector pstress = s.getPrincipalStresses( testPoints, true ) ;

	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	if( maxStress >= upVal ) {
		metInTension = true;

		if( minStress <= downVal )
			metInCompression = true ;

		return 1. - std::abs( upVal / maxStress ) ;
	}

	if( minStress <= downVal ) {
		metInCompression = true ;
		return 1. - std::abs( downVal / minStress ) ;
	}

	double s0 = -1. + std::abs( maxStress / upVal );
	double s1 = -1. + std::abs( minStress / downVal ) ;

	if( minStress > 0 ) {
		return s0 ;
	}

	if( maxStress < 0 ) {
		return s1 ;
	}



	if( std::abs( s0 ) > std::abs( s1 ) )
		return s0 ;

	return s1;
}

FractureCriterion *MohrCoulomb::getCopy() const
{
	return new MohrCoulomb( *this ) ;
}

Material MohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}

NonLocalMohrCoulomb::NonLocalMohrCoulomb( double up, double down, MirrorState mirroring, double delta_x, double delta_y, double delta_z ) : FractureCriterion( mirroring, delta_x, delta_y, delta_z )
	, upVal( up ), downVal( down )
{
}


NonLocalMohrCoulomb::~NonLocalMohrCoulomb()
{
	for( size_t i = 0 ; i < testPoints.size() ; i++ )
		delete testPoints[i] ;
}

double NonLocalMohrCoulomb::grade( ElementState &s )
{

	if( s.getParent()->getBehaviour()->fractured() )
		return 0 ;

	Vector str( s.getPrincipalStressAtNodes() ) ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		std::vector<double> fact ;
		double area = s.getParent()->area() ;
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
			fact.push_back(area);
			
		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;

			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() )
			{
				continue ;
			}

			if( ci->getBehaviour()->getFractureCriterion() &&  !ci->getBehaviour()->fractured() )
			{
				std::vector<double> dc ;
				for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					dc.push_back( squareDist2D( ci->getCenter(), s.getParent()->getCenter() )) ;
				
				std::vector<double> d ;
				for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					d.push_back( exp( -dc[j] / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) ));
				
				Vector pstress( ci->getState().getPrincipalStressAtNodes() ) ;
				Vector pstrain( ci->getState().getPrincipalStrainAtNodes() ) ;

				if( !ci->getBehaviour()->fractured() )
				{
					for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < 2 ; k++)
						{
							str[j*2+k] += pstress[j*2+k] * d[j] * area;
						}
						fact[j] += area ;
					}

					if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 2 ; k++)
							{
								str[j*2+k] += pstress[j*2+k] * d[j] * area;
							}
							fact[j] += area ;
						}
					}

					if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 2 ; k++)
							{
								str[j*2+k] += pstress[j*2+k] * d[j] * area;
							}
							fact[j] += area ;
						}
					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 2 ; k++)
							{
								str[j*2+k] += pstress[j*2+k] * d[j] * area;
							}
							fact[j] += area ;
						}
					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 2 ; k++)
							{
								str[j*2+k] += pstress[j*2+k] * d[j] * area;
							}
							fact[j] += area ;
						}
					}
					
					for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
							fact[j] += area ;
				}
			}
		}

		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < 2 ; k++)
			{
				str[j*2+k] /= fact[j] ;
			}
		}

	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		std::vector<double> fact ;
		double volume = s.getParent()->volume() ;
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
			fact.push_back(volume);

		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;

			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() )
			{
				continue ;
			}



			if( ci->getBehaviour()->getFractureCriterion() &&  !ci->getBehaviour()->fractured() )
			{
				std::vector<double> dc ;
				for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					dc.push_back( squareDist3D( ci->getCenter(), s.getParent()->getCenter() )) ;
				
				std::vector<double> d ;
				volume = ci->volume() ;
				
				for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					d.push_back( exp( -dc[j] / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) ));
				Vector pstress = ci->getState().getPrincipalStressAtNodes() ;
				Vector pstrain = ci->getState().getPrincipalStrainAtNodes() ;

				if( !ci->getBehaviour()->fractured() )
				{
					for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < 3 ; k++)
						{
							str[j*3+k] += pstress[j*3+k] * d[j] * volume;
						}
						fact[j] += volume ;
					}

					if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}

					}

					if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}

					if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
					{
						for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						{
							for(size_t k = 0 ; k < 3 ; k++)
							{
								str[j*3+k] += pstress[j*3+k] * d[j] * volume;
							}
							fact[j] += volume ;
						}
					}
					
					for(size_t j = 0 ; j < ci->getBoundingPoints().size() ; j++)
						fact[j] += volume ;
				}
			}
		}

		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < 3 ; k++)
			{
				str[j*3+k] /= fact[j] ;
			}
		}
	}

	Vector pstress(0., s.getParent()->spaceDimensions()) ;
	for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
	{
		for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
		{
			pstress[k] += str[j*s.getParent()->spaceDimensions()+k] / s.getParent()->getBoundingPoints().size() ;
		}
	}
	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	if( maxStress >= upVal ) {
		metInTension = true;

		if( minStress <= downVal )
			metInCompression = true ;

		return 1. - std::abs( upVal / maxStress ) ;
	}

	if( minStress <= downVal ) {
		metInCompression = true ;
		return 1. - std::abs( downVal / minStress ) ;
	}

	double s0 = -1. + std::abs( maxStress / upVal );
	double s1 = -1. + std::abs( minStress / downVal ) ;

	if( minStress > 0 ) {
		return s0 ;
	}

	if( maxStress < 0 ) {
		return s1 ;
	}



	if( std::abs( s0 ) > std::abs( s1 ) )
		return s0 ;

	return s1;
}

FractureCriterion *NonLocalMohrCoulomb::getCopy() const
{
	return new NonLocalMohrCoulomb( *this ) ;
}

Material NonLocalMohrCoulomb::toMaterial()
{
	Material mat ;
	return mat ;
}

}
