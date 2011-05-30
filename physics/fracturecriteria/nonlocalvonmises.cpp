
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "nonlocalvonmises.h"
#include "../../mesher/delaunay.h"
namespace Mu {

NonLocalVonMises::NonLocalVonMises(double thresh, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(thresh)
{
	setMaterialCharacteristicRadius(radius);
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(ElementState &s)
{
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
	
	double maxStress = 0 ;
		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
		{
			maxStress = sqrt( ( ( str[0] - str[1] ) * ( str[0] - str[1] ) + str[0] * str[0] + str[1] * str[1] ) / 2. ) ;
		}
		else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		{
			maxStress = sqrt( ( str[0] - str[1] ) * ( str[0] - str[1] ) + ( str[0] - str[2] ) * ( str[0] - str[2] ) + ( str[1] - str[2] ) * ( str[1] - str[2] ) ) / 6 ;
		}
	
	if(maxStress > threshold )
	{
		return 1. - std::abs(threshold/maxStress) ;
	}
	else 
	{
		return -1.+ std::abs(maxStress/threshold);
	}
}

FractureCriterion * NonLocalVonMises::getCopy() const
{
	return new NonLocalVonMises(threshold, radius) ;
}

Material NonLocalVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}

}
