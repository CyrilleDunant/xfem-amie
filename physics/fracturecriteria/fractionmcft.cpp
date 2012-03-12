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
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"

namespace Mu {

FractionMCFT::FractionMCFT(double up, double down , Matrix steelCGTensor, double youngModulus, double phi, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, upVal(up), downVal(down), steelCGTensor(steelCGTensor), phi(phi)
{
	tensionCritStrain = up / youngModulus ;
	metInCompression = false ;
	metInTension = false ;
}


FractionMCFT::~FractionMCFT()
{
}

double FractionMCFT::grade(ElementState &s)
{
	
	Vector str( s.getPrincipalStressAtNodes() ) ;
	Vector stra( s.getPrincipalStrainAtNodes() ) ;
	
	Vector totalStrain = s.getStrain(s.getParent()->getCenter()) ;
	Vector totalStress = s.getStress(s.getParent()->getCenter()) ;
	Vector concreteStress = (s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter()) - steelCGTensor*phi)*totalStrain ;
	
	double factor = sqrt(std::inner_product(&concreteStress[0], &concreteStress[concreteStress.size()], &concreteStress[0], double(0)))/sqrt(std::inner_product(&totalStress[0], &totalStress[totalStress.size()], &totalStress[0], double(0))) ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double area = s.getParent()->area() ;
		str *= area ;
		stra *= area ;
		double fact = area;
			
		// gaussian smooth
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| !ci->getBehaviour()->getFractureCriterion() 
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() &&ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->fractured()
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3. * physicalCharacteristicRadius * physicalCharacteristicRadius)
			{
				continue ;
			}

			double d = exp( -dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius ) );

			Vector pstress( ci->getState().getPrincipalStressAtNodes() ) ;
			pstress *= factor ;
			Vector pstrain( ci->getState().getPrincipalStrainAtNodes() ) ;
			
			area = ci->area() ;

			str += pstress * d * area;
			stra += pstrain * d * area;
			fact += area ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				str += pstress * d * area;
				stra += pstrain * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				str += pstress * d * area;
				stra += pstrain * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				str += pstress * d * area;
				stra += pstrain * d * area;
				fact += area ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				str += pstress * d * area;
				stra += pstrain * d * area;
				fact += area ;
			}
		}
		str /= fact ;
		stra /= fact ;
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
				|| (!ci->getBehaviour()->getTensor(ci->getCenter()).isNull() &&ci->getBehaviour()->getTensor(ci->getCenter())[0][0] < POINT_TOLERANCE_3D)
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 3.* physicalCharacteristicRadius * physicalCharacteristicRadius
			)
			{
				continue ;
			}

			

			volume = ci->volume() ;
			double d =  exp(-dc / ( physicalCharacteristicRadius * physicalCharacteristicRadius )) ;
			Vector pstress = ci->getState().getPrincipalStressAtNodes() ;
			pstress *= factor ;
			Vector pstrain = ci->getState().getPrincipalStrainAtNodes() ;

			if( !ci->getBehaviour()->fractured() )
			{
				str += pstress * d * volume;
				stra += pstrain * d * volume;
				fact += volume ;

				if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					str += pstress * d * volume;
					stra += pstrain * d * volume;
					fact += volume ;
				}
			}
		}
		str /= fact ;
		stra /= fact ;
	}

	Vector pstrain(0., s.getParent()->spaceDimensions()) ;
	Vector pstress(0., s.getParent()->spaceDimensions()) ;
	for(size_t j = 0 ; j < s.getParent()->getBoundingPoints().size() ; j++)
	{
		for(size_t k = 0 ; k < s.getParent()->spaceDimensions() ; k++)
		{
			pstrain[k] += stra[j*s.getParent()->spaceDimensions()+k]  ;
			pstress[k] += str[j*s.getParent()->spaceDimensions()+k]  ;
		}
	}
	pstrain /= s.getParent()->getBoundingPoints().size() ;
	pstress /= s.getParent()->getBoundingPoints().size() ;
	

	double tstrain = pstrain.max();
	double cstrain = pstrain.min();
	double tstress = pstress.max();
	double cstress = pstress.min();

	metInCompression = false ;
	metInTension = false ;
	
	double critStrain = -0.002 ;
	double renormCompressionStrain = cstrain/critStrain ;
	
	double mcftFactor = (2.*renormCompressionStrain-renormCompressionStrain*renormCompressionStrain)/(0.8-0.34*tstrain/critStrain) ;
	double maxCompression = -std::abs(downVal)*mcftFactor ;
	
	if(mcftFactor > 1 || mcftFactor < 0)
		maxCompression = -std::abs(downVal) ;
	
	double maxTension = upVal ;
	
//	if(tstrain > tensionCritStrain )
//	{
		//Yamamoto model 
//		maxTension = upVal/(1.+sqrt(200000000.*(tstrain+tensionCritStrain))) ;
		
		//MCFT model 
 		maxTension = upVal/(1.+sqrt(500.*tstrain)) ;
		
		//perfectly brittle
// 		maxTension = 0 ;
//	}
	
	std::vector<double> crits ;
	crits.push_back(-1) ;
	
	if( cstress <= 0 && std::abs(cstress) >= std::abs(maxCompression) )
	{
		metInCompression = true ;
		crits.push_back(1. - std::abs(maxCompression/cstress)) ;
	}
	
	if(tstress >= 0 && std::abs(tstress) >= std::abs(maxTension))
	{
		metInTension = true ;
		crits.push_back(1. - std::abs(maxTension/tstress)) ;
	}
	
	
	if(tstress >= 0 && std::abs(tstress)  < std::abs(maxTension))
	{
		crits.push_back(-1. + std::abs(tstress/maxTension)) ;
	}

	if(cstress <= 0 && std::abs(cstress) < std::abs(downVal))
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
