//
// C++ Implementation: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "anisotropicdamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

AnisotropicLinearDamage::AnisotropicLinearDamage(double characteristicRadius) : DamageModel(characteristicRadius)
{
	state.resize(4, 0.) ;
	previousstate.resize(4, 0.) ;
	isNull = false ;
	state = 0 ;

	inCompression = false ;
	inTension = false ;
}

Vector AnisotropicLinearDamage::computeDamageIncrement(ElementState & s)
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
		compressionDamage = (1.-state[0]) ; 
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
			tensionDamagex = angle[0]/n*(1.-state[1]) ; 
		}
		if(angle[1] > 0)
		{
			tensionDamagey = angle[1]/n*(1.-state[2]) ; 
		}
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL && angle[2] > 0)
		{
			tensionDamagez = angle[2]/n*(1.-state[3]) ; 
		}

	}
	ret[0] = compressionDamage ;
	ret[1] = tensionDamagex ;
	ret[2] = tensionDamagey ;
	ret[3] = tensionDamagez ;
	return ret ;
// 	std::cout << state.sum() << std::flush ;
}

void AnisotropicLinearDamage::artificialDamageStep(double d)
{
	for(size_t i = 0 ; i < state.size() -1 ; i++)
		state[i] = std::min(state[i]+d,0.9999) ;
}

Matrix AnisotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.;

	if(inTension)
	{
// 		return m*(1.-tensionDamage) ;
		if(state[1] < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 1.-state[1] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 0. ;
		}
		
		if(state[2] < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 1.-state[2] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 0. ;
		}
		
		if(state[3] < secondaryThresholdDamageDensity/fraction && ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 1.-state[3] ;
		}
		else if( ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 0. ;
		}
		
		
		if(ret.numRows() <= 3)
		{
			if(state[0] < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 1.-state[0] ;
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 0. ;
			}
		}
		else
		{
			if(state[0] < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 1.-state[0] ;
					ret[4][j]*= 1.-state[0] ;
					ret[5][j]*= 1.-state[0] ;
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
		if(state[0] >= thresholdDamageDensity/fraction)
			return m*0. ;
		
		return m*(1.-state[0]) ;
	}
	return ret*(1.-std::max(std::max(state[1], std::max(state[2],state[3])), state[0])) ;
}

Matrix AnisotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return m*0.;

	if(inTension)
	{
// 		return m*(1.-tensionDamage) ;
		if(previousstate[1] < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 1.-previousstate[1] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 0. ;
		}
		
		if(previousstate[2] < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 1.-previousstate[2] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 0. ;
		}
		
		if(previousstate[3] < secondaryThresholdDamageDensity/fraction && ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 1.-previousstate[3] ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 0. ;
		}
		
		if(ret.numRows() <= 3)
		{
			if(previousstate[0] < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 1.-previousstate[0] ;
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 0. ;
			}
		}
		else
		{
			if(previousstate[0] < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 1.-previousstate[0] ;
					ret[4][j]*= 1.-previousstate[0] ;
					ret[5][j]*= 1.-previousstate[0] ;
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
		if(previousstate[0] >= thresholdDamageDensity/fraction)
			return m*0. ;
		
		return m*(1.-previousstate[0]) ;
	}
	
	return ret*(1.-std::max(std::max(previousstate[1], std::max(previousstate[2],previousstate[3])), previousstate[0])) ;
}

bool AnisotropicLinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
// 	if(inTension)
// 		if(std::min(tensionDamagey, std::min(tensionDamagex,tensionDamagez)) >= secondaryThresholdDamageDensity/fraction)
// 			return true ;
		
	if(inCompression)
		if(state[0] >= thresholdDamageDensity/fraction)
			return true ;
		
		return false ;
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
// 	return tensionDamagey >= secondaryThresholdDamageDensity/fraction || tensionDamagex >= secondaryThresholdDamageDensity/fraction || compressionDamage >= thresholdDamageDensity/fraction  ;
}

AnisotropicLinearDamage::~AnisotropicLinearDamage()
{
}


}
