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

	if( testPoints.size() == 0 )
	{
		if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		{
			testPoints.resize( 4 );
			testPoints[0] = new Point( 0, 0, 0 ) ;
			testPoints[1] = new Point( 1, 0, 0 ) ;
			testPoints[2] = new Point( 0, 1, 0 ) ;
			testPoints[3] = new Point( 0, 0, 1 ) ;
		}
		else
		{
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

	if( maxStress >= upVal )
	{
		metInTension = true;

		if( minStress <= downVal )
			metInCompression = true ;

		return 1. - std::abs( upVal / maxStress ) ;
	}

	if( minStress <= downVal )
	{
		metInCompression = true ;
		return 1. - std::abs( downVal / minStress ) ;
	}

	double s0 = -1. + std::abs( maxStress / upVal );
	double s1 = -1. + std::abs( minStress / downVal ) ;

	if( minStress > 0 )
	{
		return s0 ;
	}

	if( maxStress < 0 )
	{
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
		double area = s.getParent()->area() ;
		str *= area ;
		double fact = area;
			
		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| !ci->getBehaviour()->getFractureCriterion() 
				|| (!(ci->getBehaviour()->getTensor(ci->getCenter()).isNull()) && ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->fractured()
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3. * physicalCharacteristicRadius * physicalCharacteristicRadius)
			{
				continue ;
			}

			double d = exp( -dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) );

			Vector pstress( ci->getState().getPrincipalStressAtNodes() ) ;
			
			area = ci->area() ;

			str += pstress * d * area;
			fact += area ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				str += pstress * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				str += pstress * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				str += pstress * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				str += pstress * d * area;
				fact += area ;
			}
		}
		str /= fact ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double fact ;
		double volume = s.getParent()->volume() ;
		fact = volume;

		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()  
				|| ci->getBehaviour()->getFractureCriterion() 
				|| (!(ci->getBehaviour()->getTensor(ci->getCenter()).isNull()) && ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3.* physicalCharacteristicRadius * physicalCharacteristicRadius
			)
			{
				continue ;
			}

			

			volume = ci->volume() ;
			double d =  exp(-dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius )) ;
			Vector pstress = ci->getState().getPrincipalStressAtNodes() ;

			if( !ci->getBehaviour()->fractured() )
			{
				str += pstress * d * volume;
				fact += volume ;

				if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					fact += volume ;
				}
			}
		}
		str /= fact ;
	}

	Vector pstress( 0., s.getParent()->spaceDimensions() ) ;

	for( size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++ )
	{
		for( size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++ )
		{
			pstress[k] += str[j * s.getParent()->spaceDimensions() + k] / s.getParent()->getBoundingPoints().size() ;
		}
	}

	double maxStress = pstress.max() ;
	double minStress = pstress.min() ;

// 	std::cout << pstress0[0] << ", " << pstress0[1] << ", "<< pstress0[2] << std::endl ;
	metInTension = false ;
	metInCompression = false ;
	metInCompression = std::abs( minStress / downVal ) > std::abs( maxStress / upVal ) ;
	metInTension = std::abs( minStress / downVal ) < std::abs( maxStress / upVal ) ;

	if( maxStress >= upVal )
	{
		metInTension = true;

		if( minStress <= downVal )
			metInCompression = true ;

		return 1. - std::abs( upVal / maxStress ) ;
	}

	if( minStress <= downVal )
	{
		metInCompression = true ;
		return 1. - std::abs( downVal / minStress ) ;
	}

	double s0 = -1. + std::abs( maxStress / upVal );
	double s1 = -1. + std::abs( minStress / downVal ) ;

	if( minStress > 0 )
	{
		return s0 ;
	}

	if( maxStress < 0 )
	{
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
