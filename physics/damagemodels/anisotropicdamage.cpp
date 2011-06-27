//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "anisotropicdamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

AnisotropicLinearDamage::AnisotropicLinearDamage() 
{
	getState(true).resize(4, 0.) ;
	getPreviousState().resize(4, 0.) ;
	isNull = false ;

	inCompression = false ;
	inTension = false ;
}

std::pair<Vector, Vector> AnisotropicLinearDamage::computeDamageIncrement(ElementState & s)
{
	inCompression = false ;
	inTension = false ;
	Vector ret(0., 4) ;
	
// 	double E_2 = s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())[0][0] ; E_2*=E_2 ;
// 	double l_2 = s.getParent()->area() ; 
// 	double maxincrement = std::abs((l_2*E_2-1.)/(l_2+l_2*E_2)) ;
	double compressionDamage = 0;
	double tensionDamagex = 0;
	double tensionDamagey = 0;
	double tensionDamagez = 0;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		compressionDamage = 1 ; 
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
	{
		inTension = true ;
		
		Vector angle = s.getStrain(s.getParent()->getCenter()) ;
		double n = sqrt(angle[0]*angle[0]+angle[1]*angle[1]) ;
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
			n = sqrt(angle[0]*angle[0]+angle[1]*angle[1]+angle[2]*angle[2]) ;
		
		if(angle[0] > 0)
		{
			tensionDamagex = angle[0] ; 
		}
		if(angle[1] > 0)
		{
			tensionDamagey = angle[1] ; 
		}
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && angle[2] > 0)
		{
			tensionDamagez = angle[2] ; 
		}

	}
	ret[0] = compressionDamage ;
	ret[1] = tensionDamagex ;
	ret[2] = tensionDamagey ;
	ret[3] = tensionDamagez ;
	if( std::abs(ret).max() > POINT_TOLERANCE_2D)
		ret /= std::abs(ret).max() ;
	return std::make_pair(state,ret) ;
// 	std::cout << getState().sum() << std::flush ;
}

void AnisotropicLinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < getState().size() -1 ; i++)
		getState(true)[i] = std::min(getState()[i]+d,0.9999) ;
}

Matrix AnisotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.;

	if(inTension)
	{
// 		return m*(1.-tensionDamage) ;
		if(getState()[1] < secondaryThresholdDamageDensity)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 1.-getState()[1] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 0. ;
		}
		
		if(getState()[2] < secondaryThresholdDamageDensity)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 1.-getState()[2] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 0. ;
		}
		
		if(getState()[3] < secondaryThresholdDamageDensity && ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 1.-getState()[3] ;
		}
		else if( ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 0. ;
		}
		
		
		if(ret.numRows() <= 3)
		{
			if(getState()[0] < secondaryThresholdDamageDensity )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 1.-getState()[0] ;
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 0. ;
			}
		}
		else
		{
			if(getState()[0] < secondaryThresholdDamageDensity )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 1.-getState()[0] ;
					ret[4][j]*= 1.-getState()[0] ;
					ret[5][j]*= 1.-getState()[0] ;
				}
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 0. ;
					ret[4][j]*= 0. ;
					ret[5][j]*= 0. ;
				}
			}
		}
		
		return ret ;
	}
	if(inCompression)
	{
		if(getState()[0] >= thresholdDamageDensity)
			return m*0. ;
		
		return m*(1.-getState()[0]) ;
	}
	return ret*(1.-std::max(std::max(getState()[1], std::max(getState()[2],getState()[3])), getState()[0])) ;
}

Matrix AnisotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.;

	if(inTension)
	{
// 		return m*(1.-tensionDamage) ;
		if(getPreviousState()[1] < secondaryThresholdDamageDensity)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 1.-getPreviousState()[1] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 0. ;
		}
		
		if(getPreviousState()[2] < secondaryThresholdDamageDensity)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 1.-getPreviousState()[2] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 0. ;
		}
		
		if(getPreviousState()[3] < secondaryThresholdDamageDensity && ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 1.-getPreviousState()[3] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 0. ;
		}
		
		if(ret.numRows() <= 3)
		{
			if(getPreviousState()[0] < secondaryThresholdDamageDensity )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 1.-getPreviousState()[0] ;
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 0. ;
			}
		}
		else
		{
			if(getPreviousState()[0] < secondaryThresholdDamageDensity )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 1.-getPreviousState()[0] ;
					ret[4][j]*= 1.-getPreviousState()[0] ;
					ret[5][j]*= 1.-getPreviousState()[0] ;
				}
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 0. ;
					ret[4][j]*= 0. ;
					ret[5][j]*= 0. ;
				}
			}
		}
		
		return ret ;
	}
	if(inCompression)
	{
		if(getPreviousState()[0] >= thresholdDamageDensity)
			return m*0. ;
		
		return m*(1.-getPreviousState()[0]) ;
	}
	
	return ret*(1.-std::max(std::max(getPreviousState()[1], std::max(getPreviousState()[2],getPreviousState()[3])), getPreviousState()[0])) ;
}

bool AnisotropicLinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
	if(inTension)
		if(std::min(getState()[2], std::min(getState()[1],getState()[3])) >= secondaryThresholdDamageDensity)
			return true ;
		
	if(inCompression)
		if(getState()[0] >= thresholdDamageDensity)
			return true ;
		
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
	return getState()[2]  >= secondaryThresholdDamageDensity || getState()[1]  >= secondaryThresholdDamageDensity  || getState()[3]  >= secondaryThresholdDamageDensity || getState()[0] >= thresholdDamageDensity;
}

AnisotropicLinearDamage::~AnisotropicLinearDamage()
{
}


}
