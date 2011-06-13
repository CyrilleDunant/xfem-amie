//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
using namespace Mu ;

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) : neighbourhoodradius(.0005), 
neighbourhoodvolume(-1), 
physicalCharacteristicRadius(.008), 
scoreAtState(-1), 
nonLocalScoreAtState(-1),
metInTension(false), metInCompression(false), metAtStep(false), 
mirroring(mirroring), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z), 
deltaScoreAtState(0), 
deltaEnergyAtState(0), 
energyDamageDifferential(0), 
criterionDamageDifferential(0), 
energyIndexed(false), 
noEnergyUpdate(true), 
mesh2d(NULL), mesh3d(NULL), 
stable(true), checkpoint(false), inset(false),
scoreTolerance(1e-4)
{
}

Vector FractureCriterion::smoothedPrincipalStrain(ElementState &s) const
{
	Vector stra(s.getParent()->spaceDimensions()) ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double area = s.getParent()->area() ;

		Vector strainAtNodes(s.getPrincipalStrainAtNodes()) ;
		
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
			{
				stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area/s.getParent()->getBoundingPoints().size() ;
			}
		}

		double fact = area ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| !ci->getBehaviour()->getFractureCriterion() 
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() && ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->fractured()
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 4. * physicalCharacteristicRadius * physicalCharacteristicRadius)
			{
				continue ;
			}
			
			//this is to eliminate scaling effects ;
			double factor = 1 ;
			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
				factor = std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]) ;
			
			double d = exp( -dc / (2.* physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;

			Vector strainAtNodes(ci->getState().getPrincipalStrainAtNodes()) ;
			
			for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
			{
				for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
				{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
				}
			}
			fact += area*d ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}
		}
		stra /= fact ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double volume = s.getParent()->volume() ;
		Vector strainAtNodes(s.getPrincipalStrainAtNodes()) ;
		
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
			{
				stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume/s.getParent()->getBoundingPoints().size() ;
			}
		}

		double fact = volume ;
		
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;
			
			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()  
				|| ci->getBehaviour()->getFractureCriterion() 
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() &&ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3.* physicalCharacteristicRadius * physicalCharacteristicRadius
			)
			{
				continue ;
			}

			double volume = ci->volume() ;
			double factor = 1 ;
			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
				factor = std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]) ;
			
			double d = exp( -dc / ( 2.*physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;

			Vector strainAtNodes(ci->getState().getPrincipalStrainAtNodes()) ;
			
			if( !ci->getBehaviour()->fractured() )
			{

				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
				
				if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
					{
						for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
						{
							stra[k] += strainAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
						}
					}
					fact += volume*d ;
				}
			}
		}
		stra /= fact ;
	}
	
	return stra ;
}

Vector FractureCriterion::smoothedPrincipalStress( ElementState &s) const
{
	Vector str(s.getParent()->spaceDimensions() ) ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double area = s.getParent()->area() ;

		Vector stressAtNodes(s.getPrincipalStressAtNodes()) ;
		
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
			{
				str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area/s.getParent()->getBoundingPoints().size() ;
			}
		}

		double fact = area ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| !ci->getBehaviour()->getFractureCriterion() 
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() && ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->fractured()
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 4. * physicalCharacteristicRadius * physicalCharacteristicRadius)
			{
				continue ;
			}
			
			area = ci->area() ;
			//this is to eliminate scaling effects ;
			double factor = 1 ;
			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
				factor = std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]) ;
			
			double d = exp( -dc / (2.* physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;

			Vector stressAtNodes(ci->getState().getPrincipalStressAtNodes()) ;
			
			for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
			{
				for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
				{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
				}
			}
			fact += area*d ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*area*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += area*d ;
			}
		}
		str /= fact ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double volume = s.getParent()->volume() ;
		Vector stressAtNodes(s.getPrincipalStressAtNodes()) ;
		
		for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
		{
			for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
			{
				str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume/s.getParent()->getBoundingPoints().size() ;
			}
		}

		double fact = volume ;
		
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;
			
			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()  
				|| ci->getBehaviour()->getFractureCriterion() 
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() &&ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3.* physicalCharacteristicRadius * physicalCharacteristicRadius
				|| ci->getBehaviour()->fractured()
			)
			{
				continue ;
			}

			volume = ci->volume() ;
			double factor = 1 ;
			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
				factor = std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]) ;
			
			double d = exp( -dc / ( 2.*physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;
			
			Vector stressAtNodes(ci->getState().getPrincipalStressAtNodes()) ;
			

			for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
			{
				for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
				{
					str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
				}
			}
			fact += volume*d ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}

			if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
				{
					for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
					{
						str[k] += stressAtNodes[j*s.getParent()->spaceDimensions()+k]*volume*d/s.getParent()->getBoundingPoints().size() ;
					}
				}
				fact += volume*d ;
			}
		}
		str /= fact ;
	}
	
	return str ;
}

double FractureCriterion::getDeltaEnergy(const ElementState & s, double delta_d)
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = static_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron *>((*mesh3d)[ cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[ cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * >((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return (originalenergy-energy)/(delta_d) ;
}

void FractureCriterion::initialiseCache(const ElementState & s)
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	if(testedTri)
	{
		if(!cache.empty())
		{
			cache.clear();
			area.clear();
		}
		Circle epsilon(neighbourhoodradius,testedTri->getCenter()) ;
		if(!testedTri->tree)
			return ;
		mesh2d = &testedTri->tree->getTree() ;
		std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getNeighbouringElementsInGeometry(testedTri, &epsilon);
		std::vector<DelaunayTriangle *> neighbourhood ;
		for(size_t i = 0 ; i < testedTri->neighbourhood.size() ; i++)
		{
			if(testedTri->getNeighbourhood(i)->getBehaviour() && 
				testedTri->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR && 
				testedTri->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTri->getNeighbourhood(i));
				cache.push_back(testedTri->getNeighbourhood(i)->index);
				area.push_back(testedTri->getNeighbourhood(i)->area());
			}
		}
		
		
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour() && 
				tempcache[i]->getBehaviour()->type != VOID_BEHAVIOUR && 
				tempcache[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache[i]->index);
					area.push_back(tempcache[i]->area());
				}
			}
		}
		
		if(cache.empty())
			cache.push_back(testedTri->index);
	}
	else if(testedTet)
	{
		if(!cache.empty())
			cache.clear();
		
		Sphere epsilon(neighbourhoodradius,testedTet->getCenter()) ;
		if(!testedTet->tree)
			return ;
		mesh3d = &testedTet->tree->getTree() ;
		std::vector<DelaunayTetrahedron *> tempcache3d = testedTet->tree->getConflictingElements(&epsilon);
		std::vector<DelaunayTetrahedron *> neighbourhood ;
		for(size_t i = 0 ; i < testedTet->neighbourhood.size() ; i++)
		{
			if(testedTet->getNeighbourhood(i)->getBehaviour()
				&& testedTet->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR  
				&& testedTet->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTet->getNeighbourhood(i));
				cache.push_back(testedTet->getNeighbourhood(i)->index);
				area.push_back(testedTet->getNeighbourhood(i)->volume());
			}
		}
		
		for(size_t i = 0 ; i < tempcache3d.size() ; i++)
		{
			if(tempcache3d[i]->getBehaviour()
				&& tempcache3d[i]->getBehaviour()->type != VOID_BEHAVIOUR 
				&& tempcache3d[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache3d[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache3d[i]->index);
					area.push_back(tempcache3d[i]->volume());
				}
			}
		}
		if(cache.empty())
			cache.push_back(testedTet->index);
	}
	
	std::sort(cache.begin(), cache.end()) ;
}

std::pair<double, double> FractureCriterion::getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,   tric->second,    tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j <  tric->neighbour.size() ; j++)
			{
				if( tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>( tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					score += 1./(1.-sc);
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * > ((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					score += 1./(1.-sc) ;
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return std::make_pair( (energy-originalenergy)/(delta_d),(score-originalscore)/(delta_d)) ;
}

std::pair<double, double> FractureCriterion::setChange(const ElementState &s) 
{
	stable = true ;

	if( !s.getParent()->getBehaviour()->getDamageModel())
		return std::make_pair(0.,0.) ;
	if(cache.size() == 0)
	{
		initialiseCache(s);
	}
	if(!inset && !checkpoint)
		return std::make_pair(0.,0.) ;
	
	if(mesh2d)
	{
		
		// outside of the checkpoints, we only care about the order of the elements in 
		// term of their score. At the checkpoint, we consider the elements which
		// have met their criterion
		if(checkpoint) //new iteration
		{
			if(!metAtStep)
			{
				inset = false ;
				damagingSet.clear();
				return std::make_pair(0.,0.) ;
			}
			
			std::vector<unsigned int> newSet ;
			std::multimap<double, DelaunayTriangle *> sortedElements ;
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
				if(ci->getBehaviour()->getFractureCriterion())
				{
					double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
					sortedElements.insert( std::make_pair(renormScore, ci)) ;
				}
			}
			double thresholdScore = 0 ;
			if(!sortedElements.empty())
				thresholdScore = sortedElements.rbegin()->first ;
			double minscore = thresholdScore ;
			Point maxLocus = sortedElements.rbegin()->second->getCenter() ;
			
			if(!sortedElements.empty() && thresholdScore > 0 )
			{
				for(auto i = sortedElements.rbegin() ; i != sortedElements.rend() ; i++ )
				{
					if(i->second->getBehaviour()->getFractureCriterion()->met() && std::abs(i->first-thresholdScore) < scoreTolerance/*squareDist2D(maxLocus, i->second->getCenter()) <= physicalCharacteristicRadius*physicalCharacteristicRadius*/)
					{
						newSet.push_back(i->second->index);
						minscore = i->second->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
					}
					else
						break ;
				}
			}

			if(std::abs(nonLocalScoreAtState - thresholdScore) >= scoreTolerance)
			{
				proximitySet.clear() ;
				return std::make_pair(0.,0.) ;
			}
			inset = true ;
			if(!newSet.empty())
				std::stable_sort(newSet.begin(), newSet.end());
			damagingSet = newSet ;
			std::set<unsigned int> newProximity ;
			
			if(!sortedElements.empty()&& thresholdScore > 0)
			{
				for(auto i = sortedElements.rbegin() ; i != sortedElements.rend() ; i++ )
				{
					if(i->second->getBehaviour()->getFractureCriterion()->met() && std::abs(i->first-thresholdScore) < scoreTolerance/*squareDist2D(maxLocus, i->second->getCenter()) <= physicalCharacteristicRadius*physicalCharacteristicRadius*/)
					{
						continue ;
					}
					else
						newProximity.insert(i->second->index) ;
				}
			}
			
			if(newProximity.empty())
			{
				std::cerr << "element too small!" << std::endl ;
				DelaunayTriangle * ci = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
				for(size_t i = 0 ; i < ci->neighbourhood.size() ; i++)
				{
					if(static_cast<DelaunayTriangle *>(ci->getNeighbourhood(i))->getBehaviour()->getFractureCriterion())
						newProximity.insert(ci->neighbourhood[i]);
				}
			}
			proximitySet.clear() ;
			proximitySet.insert(proximitySet.end(), newProximity.begin(), newProximity.end()) ;
			
			
			double maxscore = 0 ;

			if(!proximitySet.empty())
			{
				maxscore = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[0]])->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
			}
			for(size_t i = 1 ; i < proximitySet.size() ; i++)
			{
				double nls = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[i]])->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
				if(nls > maxscore)
					maxscore = nls ;
			}

			return std::make_pair(maxscore - minscore, thresholdScore-minscore) ;
		}
		else
		{

			if(damagingSet.empty())
				return std::make_pair(0.,0.) ;
			
			DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[0]]) ;
			double maxscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

			double thresholdScore = maxscore ;
			for(size_t i = 1 ; i < damagingSet.size() ; i++)
			{
				double nls = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[i]])->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

				if(nls < maxscore)
				{
					maxscore = nls ;
				}
				if(nls > thresholdScore)
				{
					thresholdScore = nls ;
				}
			}
			double minscore = 0 ;
			if(!proximitySet.empty())
			{
				minscore = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[0]])->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

				for(size_t i = 1 ; i < proximitySet.size() ; i++)
				{
					double nls = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[i]])->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

					if(nls > minscore)
						minscore = nls ;
				}
			}

			return std::make_pair(maxscore - minscore, thresholdScore-maxscore) ;
		}
	}
	else
	{
		DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
		// outside of the checkpoints, we only care about the order of the elements in 
		// term of their score. At the checkpoint, we consider the elements which
		// have met their criterion
		if(checkpoint) //new iteration
		{
			std::vector<unsigned int> newSet ;
			std::multimap<double, DelaunayTriangle *> sortedElements ;
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
				if(ci->getBehaviour()->getFractureCriterion())
				{
					sortedElements.insert( std::make_pair(ci->getBehaviour()->getFractureCriterion()->getScoreAtState(), ci)) ;
				}
			}
			
			double maxscore =  sortedElements.rbegin()->first ;
			for(auto i = sortedElements.rbegin() ; i != sortedElements.rend() ; i++ )
			{
				if(std::abs(i->second->getBehaviour()->getFractureCriterion()->getScoreAtState()-maxscore) < scoreTolerance && i->second->getBehaviour()->getFractureCriterion()->met())
				{
					newSet.push_back(i->second->index);
				}
				else
					break ;
			}
			std::stable_sort(newSet.begin(), newSet.end());
		
			inset = std::binary_search(newSet.begin(), newSet.end(), testedTet->index) ;
			if(!inset)
			{
				damagingSet.clear();
				return std::make_pair(0., 0.) ;
			}

			damagingSet = newSet ;
			inset = true ;
			
			return std::make_pair(0., 0.) ;
		}
		else
		{
			DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
			checkpoint = false ;
			
			if(damagingSet.empty())
				return std::make_pair(0., 0.) ;
			if(!metAtStep)
				return std::make_pair(0., 0.) ;
			
			DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[damagingSet[0]]) ;
			double minscore = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
			if(!ci->getBehaviour()->getFractureCriterion()->met())
				return std::make_pair(0., 0.) ;
			for(size_t i = 1 ; i < damagingSet.size() ; i++)
			{
				ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[damagingSet[i]]) ;
				if(!ci->getBehaviour()->getFractureCriterion()->met())
					return std::make_pair(0., 0.) ;
				if(ci->getBehaviour()->getFractureCriterion()->getScoreAtState() < minscore)
					minscore = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
			}
			
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
				if(ci->getBehaviour()->getFractureCriterion() )
				{
					if(ci->getBehaviour()->getFractureCriterion()->getScoreAtState() > minscore && !std::binary_search(damagingSet.begin(), damagingSet.end(), ci->index))
						return std::make_pair(0., 0.) ;
				}
			}
			
			return std::make_pair(0., 0.) ;
		}
	}
	
	//shut up the compiler
	return std::make_pair(0., 0.) ;
}

void FractureCriterion::step(ElementState &s)
{
		if(cache.empty() )
				initialiseCache(s) ;

	if(energyIndexed && s.getDeltaTime() > POINT_TOLERANCE_2D )
		noEnergyUpdate = true ;
	
	if(energyIndexed && noEnergyUpdate /*&& scoreAtState > -.5*/)
	{
		noEnergyUpdate = false ;
		currentEnergy = 0 ;
		bool buildvolume = false ;
		if(neighbourhoodvolume < 0)
		{
			neighbourhoodvolume = 0 ;
			buildvolume = true ;
		}
	
		if(mesh2d)
		{
			Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTriangle *tri = static_cast<DelaunayTriangle *>( (*mesh2d)[getCache()[i]] ) ;
				if(c.in(*tri->first) || c.in(*tri->second) || c.in(*tri->third))
				{
					currentEnergy += tri->getState().elasticEnergy()*tri->area() ;
					if(buildvolume)
						neighbourhoodvolume += tri->area() ;
				}
			}
		}
		else if(mesh3d)
		{
			Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTetrahedron *tet = static_cast<DelaunayTetrahedron *>( (*mesh3d)[getCache()[i]] ) ;
				if(c.in(*tet->first) || c.in(*tet->second) || c.in(*tet->third) || c.in(*tet->fourth))
				{
					currentEnergy += tet->getState().elasticEnergy()*tet->volume() ;
					if(buildvolume)
						neighbourhoodvolume += tet->volume() ;
				}
			}
		}
		deltaEnergyAtState = (currentEnergy-previousEnergy) ;
		
		previousEnergy = currentEnergy ;
		
	}
	
	if(energyIndexed /*&& scoreAtState > -.5*/)
	{
		std::pair<double, double> dedc = getDeltaEnergyDeltaCriterion(s, 0.0001) ;
		energyDamageDifferential = dedc.first ;
		criterionDamageDifferential = dedc.second ;
	}
	
	if(s.getParent()->getBehaviour()->fractured())
	{
		scoreAtState = -1 ;
		return ;
	}
	
	scoreAtState = grade(s) ;
}

void FractureCriterion::computeNonLocalState(ElementState &s, NonLocalSmoothingType st)
{
	metAtStep = false ;

	if( s.getParent()->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel()->fractured())
	{
		return  ;
	}
	
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	
	switch (st)
	{
		case NULL_SMOOTH :
		{
			metAtStep = scoreAtState > 0 ;
			nonLocalScoreAtState = scoreAtState ;
			return ;
		}
		case MAX_PROXIMITY_SMOOTH :
		{
			if(testedTri)
			{
				nonLocalScoreAtState = scoreAtState ;
				if (scoreAtState < 0)
				{
					metAtStep = false ;
					return  ;
				}
				double maxNeighbourhoodScore = 0 ;
				double matchedArea = 0 ;
// 				std::map<double, DelaunayTriangle *> scores ;
				std::vector<double> scores ;
				std::vector<DelaunayTriangle *> trisToTest ;
				DelaunayTriangle * maxLocus = NULL;
// 				double areamax = 0 ;
				if(!cache.empty())
				{
					for(size_t i = 0 ; i< cache.size() ; i++)
					{
						DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;

						double s = 0. ;
						if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
						{
							s = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
							if(s > 0)
								scores.push_back(s);
						}
						
						if(s > 0)
						{
							trisToTest.push_back(ci) ;
						}

						if(s > maxNeighbourhoodScore)
						{
							maxNeighbourhoodScore = s ;
							maxLocus = ci ;
						}
					}
				}
				
				if(maxNeighbourhoodScore < -scoreTolerance)
				{
					metAtStep = false ;
					return  ;
				}
				
				if(maxLocus)
				{
					std::vector<DelaunayTriangle *> toTest ;
					std::set<DelaunayTriangle *> matchingElements ;
					matchingElements.insert(maxLocus) ;
					toTest.push_back(maxLocus);
					while(!toTest.empty())
					{
						std::vector<DelaunayTriangle *> newToTest ;
						for(size_t i = 0 ; i < toTest.size() ; i++)
						{
							for(size_t j = 0 ; j < toTest[i]->neighbourhood.size() ; j++)
							{
								if(matchingElements.find(toTest[i]->getNeighbourhood(j)) == matchingElements.end()
									&& toTest[i]->getNeighbourhood(j)->getBehaviour()->getFractureCriterion()
									&& toTest[i]->getNeighbourhood(j)->getBehaviour()->getFractureCriterion()->getScoreAtState() > 0
								)
								{
									newToTest.push_back(toTest[i]->getNeighbourhood(j));
									matchingElements.insert(toTest[i]->getNeighbourhood(j)) ;
									double a = toTest[i]->getNeighbourhood(j)->area() ;
									matchedArea += a ;
									if(mirroring == MIRROR_X && std::abs(toTest[i]->getNeighbourhood(j)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
									{
										matchedArea += a ;
									}
									if(mirroring == MIRROR_Y &&  std::abs(toTest[i]->getNeighbourhood(j)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
									{
										matchedArea += a ;
									}
									if(mirroring == MIRROR_XY &&  std::abs(toTest[i]->getNeighbourhood(j)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
									{
										matchedArea += a ;
									}
									if(mirroring == MIRROR_XY &&  std::abs(toTest[i]->getNeighbourhood(j)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
									{
										matchedArea += a ;
									}
								}
								if(matchedArea >= 2.*M_PI*physicalCharacteristicRadius*physicalCharacteristicRadius)
									goto endloop ;
							}
						}

						toTest = newToTest ;
					}
				}
endloop:

// 				if(matchedArea < 2.*M_PI*physicalCharacteristicRadius*physicalCharacteristicRadius)
// 				{
// 					metAtStep = false ;
// 					return ;
// 				}

				std::sort(scores.begin(), scores.end()) ;
				double threshold = scores[round((scores.size()-1)*.9)] ;
				bool nearmaxlocus = false;
				
				for(size_t i = 0 ; i< trisToTest.size() ; i++)
				{
					if(trisToTest[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() > threshold)
					{
						if(squareDist2D(trisToTest[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
						{
							nearmaxlocus = true ;
							break ;
						}
					}
				}

				if (nearmaxlocus)
				{
					metAtStep = true ;
					return  ;
				}
				metAtStep = false ;
				return  ;

			}
			
			if(testedTet)
			{
				nonLocalScoreAtState = scoreAtState ;
				if(testedTet->visited())
				{
					metAtStep = false ;
					return  ;
				}
						
				if (scoreAtState <= 0)
				{
					metAtStep = false ;
					return  ;
				}

				double maxNeighbourhoodScore = 0 ;
// 				double matchedArea = 0 ;
// 				std::map<double, DelaunayTetrahedron *> scores ;
// 				std::map<DelaunayTetrahedron *, double> areatemp ;
// 				DelaunayTetrahedron * maxLocus = NULL;
				
				if(!cache.empty())
				{
					for(size_t i = 0 ; i< cache.size() ; i++)
					{
						DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;

						if( !ci->getBehaviour()->fractured())
						{
							double s = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
// 							scores[-s] =  ci;
							if(s > maxNeighbourhoodScore)
							{
								maxNeighbourhoodScore = s ;
// 								maxLocus = ci ;
							}
						}
// 						else if(ci->getBehaviour()->fractured())
// 						{
// 							double s = POINT_TOLERANCE_2D ;
// // 							scores[-s] =  ci;
// 						}
// 						areatemp[ci] = area[i] ;

					}
				}
				
				if(maxNeighbourhoodScore < 0)
				{
					metAtStep = false ;
					return  ;
				}
				
				std::vector<DelaunayTetrahedron *> maxloci ;
				
				for(size_t i = 0 ; i< cache.size() ; i++)
				{
					DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
					if(ci->getBehaviour()->getFractureCriterion())
						if(std::abs(ci->getBehaviour()->getFractureCriterion()->getScoreAtState()-maxNeighbourhoodScore) < scoreTolerance)
							maxloci.push_back(ci) ;
				}
				
// 				bool foundcutoff = false ;
// 				double thresholdscore = maxNeighbourhoodScore ;
// 				
// 				for(auto i = scores.begin() ; i != scores.end() ; ++i)
// 				{
// 					
// 					if(!foundcutoff)
// 					{
// 						if(-i->first > 0 )
// 						{
// 							matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_X && std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_Y &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_Z &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_XZ &&  std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_XZ &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_YZ &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 							if(mirroring == MIRROR_YZ &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
// 								matchedArea += areatemp[i->second] ;
// 						}
// 						if (matchedArea > 1.333333333333333333*physicalCharacteristicRadius*physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI)
// 						{
// 							thresholdscore = -i->first ;
// 							foundcutoff  = true ;
// 							break ;
// 						}
// 					}
// 		// 			else
// 		// 				i->second->visited = true ;
// 				}
// 				if (!foundcutoff )
// 				{
// 					metAtStep = false ;
// 					return  ;
// 				}

				for(size_t i = 0 ; i < maxloci.size() ; i++)
					if(squareDist3D(maxloci[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
				{
					metAtStep = true ;
					return  ;
				}
				metAtStep = false ;
				return  ;
			}
			else if(testedHex)
			{
				std::set<HexahedralElement *> neighbourhood ;
				std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
				for(size_t i = 0 ; i < neighbours.size() ; i++)
				{
					for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
					{
						if(neighbours[i]->neighbourhood[j] != testedHex 
							&& !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
							neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
					}
				}
				double score = grade(s) ;
				double maxNeighbourhoodScore = 0 ;
				if(!neighbourhood.empty())
				{
					for(auto i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
					{
						if((*i)->getBehaviour()->getFractureCriterion() 
							&& !(*i)->getBehaviour()->fractured())
							maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
																							(*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
						if((*i)->getBehaviour()->changed())
						{
							maxNeighbourhoodScore = 10.*score ;
							break ;
						}
						
						if (maxNeighbourhoodScore > score)
							break ;
						
					}
				}
				
				if( score > 0 )
				{
					if(score > maxNeighbourhoodScore)
					{
						metAtStep = true ;
						return  ;
					}
				}
				metAtStep = false ;
				return  ;
			}
			else
			{
				std::cout << " criterion not implemented for this kind of element" << std::endl ;
				metAtStep = false ;
				return  ;
			}
		}
		case GAUSSIAN_SMOOTH :
		{
			if(testedTri)
			{
				if(testedTri->visited)
				{
					metAtStep = false ;
					return  ;
				}
				
				double str = 0 ;
				double fact = 0 ;
				
				for(size_t i = 0 ; i< cache.size() ; i++)
				{
					DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
					double dc = squareDist2D(s.getParent()->getCenter(), ci->getCenter()) ;						
					if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
					{
						double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
						double a = ci->area() ;
						double s = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
						if(!ci->getBehaviour()->fractured())
						{
							str += s*a*d ;
							fact+=a*d ;
							if(mirroring == MIRROR_X && std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
							{
								str +=s*a*d ;
								fact+=a*d ;
							}
							if(mirroring == MIRROR_Y &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
							{
								str +=s*a*d ;
								fact+=a*d ;
							}
							if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
							{
								str +=s*a*d ;
								fact+=a*d ;
							}
							if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
							{
								str +=s*a*d ;
								fact+=a*d ;
							}
						}
					}
				}
				
				double smoothscore = str/fact ;
				
				metAtStep =  (smoothscore > scoreTolerance) ;
				if(scoreAtState >= 0)
				  std::cout << scoreAtState << "  " << smoothscore << std::endl ;
				nonLocalScoreAtState = smoothscore ;
				return ;

			}
			if(testedTet)
			{

				if(testedTet->visited())
				{
					metAtStep = false ;
					return  ;
				}
				
				double str = 0 ;
				double fact = 0 ;
				
				for(size_t i = 0 ; i< cache.size() ; i++)
				{
					DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
					double dc = squareDist3D(s.getParent()->getCenter(), ci->getCenter()) ;
					if(ci->getBehaviour()->getFractureCriterion())
					{
						double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
						double a = ci->volume() ;
						str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
						fact+=d*a ;
						if(mirroring == MIRROR_X && std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Y &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Z &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
					}
				}
				
				double smoothscore = str/fact ;
				metAtStep =  smoothscore > 0 ;
				nonLocalScoreAtState = smoothscore ;
				return ;
			}
			else if(testedHex)
			{
				std::set<HexahedralElement *> neighbourhood ;
				std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
				for(size_t i = 0 ; i < neighbours.size() ; i++)
				{
					for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
					{
							neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
					}
				}
				double str = 0 ;
				double fact = 0 ;
				
				for(auto ci = neighbourhood.begin() ; ci != neighbourhood.end() ; ++ci)
				{
					double dc = squareDist3D(s.getParent()->getCenter(), (*ci)->getCenter()) ;
					if((*ci)->getBehaviour()->getFractureCriterion())
					{
						double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
						double a = (*ci)->volume() ;
						str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
						fact+=d*a ;
						if(mirroring == MIRROR_X && std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Y &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Z &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
					}
				}
				
				double smoothscore = str/fact ;
				metAtStep =  smoothscore > 0 ;
				nonLocalScoreAtState = smoothscore ;
				return ;
			}
			else
			{
				std::cout << " criterion not implemented for this kind of element" << std::endl ;
				metAtStep = false ;
				return  ;
			}
		}
	}
	//shut up the compiler
		metAtStep = false ;
		return  ;
}

FractureCriterion::~FractureCriterion()
{
}

void FractureCriterion::setNeighbourhoodRadius(double r)
{
	neighbourhoodradius = r ;
	cache.clear() ;
}

void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
	physicalCharacteristicRadius = r ;
}

bool FractureCriterion::met() const
{
	return metAtStep ;
}

Material FractureCriterion::toMaterial()
{
	Material mat ;
	return mat ;
}


