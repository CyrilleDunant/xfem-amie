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

AnisotropicLinearDamage::AnisotropicLinearDamage(int numDof, double characteristicRadius) : DamageModel(characteristicRadius)
{
	state.resize(4, 0.) ;
	previousstate.resize(4, 0.) ;
	isNull = false ;
	state = 0 ;
	tensionDamagex = 0 ;
	tensionDamagey = 0 ;
	tensionDamagez = 0. ;
	shearDamage = 0 ;
	compressionDamage = 0 ;
	inCompression = false ;
	inTension = false ;
}

const Vector & AnisotropicLinearDamage::damageState() const
{
	return state ;
}


Vector & AnisotropicLinearDamage::damageState() 
{
	return state ;
}
void AnisotropicLinearDamage::step(ElementState & s)
{
	previousstate = state ;
	inCompression = false ;
	inTension = false ;
	if(fraction < 0)
	{
		double volume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			volume = sqrt(s.getParent()->area()) ;
		else
			volume = pow(s.getParent()->volume(), 2./3.) ;
		
		double charVolume ;
		if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			charVolume = sqrt(M_PI*characteristicRadius*characteristicRadius) ;
		else
			charVolume = pow(4./3.*M_PI*characteristicRadius*characteristicRadius*characteristicRadius, 2./3.) ;
		fraction = volume/charVolume ;
		if(fraction > 1)
			std::cout << "elements too large for damage characteristic radius!" << std::endl ;
		fraction = std::min(fraction, 1.) ;
	}
	
// 	double E_2 = s.getParent()->getBehaviour()->getTensor(s.getParent()->getCenter())[0][0] ; E_2*=E_2 ;
// 	double l_2 = s.getParent()->area() ; 
// 	double maxincrement = std::abs((l_2*E_2-1.)/(l_2+l_2*E_2)) ;
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInCompression)
	{
		inCompression = true ;
		compressionDamage +=/* std::min(*/damageDensityIncrement*fraction/*, maxincrement )*/ ; 
		compressionDamage = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE, compressionDamage) ;
		compressionDamage = std::min(.99999, compressionDamage) ;
		compressionDamage = std::max(0., compressionDamage) ;
	}
	
	if(s.getParent()->getBehaviour()->getFractureCriterion()->metInTension)
	{
		inTension = true ;
		
		Vector stress = s.getStrain(s.getParent()->getCenter()) ;
		double norm = sqrt(stress[0]*stress[0]+stress[1]*stress[1]) ;
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
			norm = sqrt(stress[0]*stress[0]+stress[1]*stress[1]+stress[2]*stress[2]) ;
		double factorx= std::abs(stress[0])/norm ;
		double factory= std::abs(stress[1])/norm ;
		double factorz= std::abs(stress[2])/norm ;
		
		tensionDamagex += factorx*/*std::min(*/damageDensityIncrement*fraction/*, maxincrement )*/ ; 
		tensionDamagex = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamagex) ;
		tensionDamagex = std::min(.99999, tensionDamagex) ;
		tensionDamagex = std::max(0., tensionDamagex) ;
		
		tensionDamagey += factory*/*std::min(*/damageDensityIncrement*fraction/*, maxincrement )*/ ; 
		tensionDamagey = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamagey) ;
		tensionDamagey = std::min(.99999, tensionDamagey) ;
		tensionDamagey = std::max(0., tensionDamagey) ;
		if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
		{
			tensionDamagez += factorz*/*std::min(*/damageDensityIncrement*fraction/*, maxincrement )*/ ; 
			tensionDamagez = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, tensionDamagez) ;
			tensionDamagez = std::min(.99999, tensionDamagez) ;
			tensionDamagez = std::max(0., tensionDamagez) ;
		}
		
		shearDamage += /*std::min(*/damageDensityIncrement*fraction/*, maxincrement )*/ ; 
		shearDamage = std::min(secondaryThresholdDamageDensity/fraction+POINT_TOLERANCE, shearDamage) ;
		shearDamage = std::min(.99999, shearDamage) ;
		shearDamage = std::max(0., shearDamage) ;
	}
	state[0] = compressionDamage ;
	state[1] = tensionDamagex ;
	state[2] = tensionDamagey ;
	state[3] = tensionDamagez ;
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
		if(tensionDamagex < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 1.-tensionDamagex ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[0][j]*= 0. ;
		}
		
		if(tensionDamagey < secondaryThresholdDamageDensity/fraction)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 1.-tensionDamagey ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[1][j]*= 0. ;
		}
		
		if(tensionDamagez < secondaryThresholdDamageDensity/fraction && ret.numRows() > 3)
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 1.-tensionDamagez ;
		}
		else
		{
			for(size_t j = 0 ; j < m.numCols() ;j++)
				ret[2][j]*= 0. ;
		}
		
		if(ret.numRows() <= 3)
		{
			if(shearDamage < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 1.-shearDamage ;
			}
			else
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
					ret[2][j]*= 0. ;
			}
		}
		else
		{
			if(shearDamage < secondaryThresholdDamageDensity/fraction )
			{
				for(size_t j = 0 ; j < m.numCols() ;j++)
				{
					ret[3][j]*= 1.-shearDamage ;
					ret[4][j]*= 1.-shearDamage ;
					ret[5][j]*= 1.-shearDamage ;
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
// 		for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
// 		for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
// 		{
// 			for(size_t j = 0 ; j < m.numCols() ;j++)
// 			{
// 				ret[i][j] *= 1.-tensionDamage ;
// 			}
// 		}
	}
	if(inCompression)
	{
		if(compressionDamage >= thresholdDamageDensity/fraction)
			return m*0. ;
		
		return m*(1.-compressionDamage) ;
// 		for(size_t i = 0 ; i < (m.numRows()+1)/2 ;i++)
// 		{
// 			for(size_t j = 0 ; j < m.numCols() ;j++)
// 			{
// 				ret[i][j] *= 1.-compressionDamage ;
// 			}
// 		}
	}
	return ret*(1.-std::max(std::max(tensionDamagex, std::max(tensionDamagey,tensionDamagez)), compressionDamage)) ;
}

Matrix AnisotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;
	
	if(fractured())
		return ret*0.0001	;
	//this is a silly way of distinguishing between 2D and 3D
	for(size_t i = 0 ; i < (m.numRows()+1)/2 ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(previousstate[0]) ;
		}
	}

	for(size_t i = (m.numRows()+1)/2 ; i < m.numRows() ;i++)
	{
		for(size_t j = 0 ; j < m.numCols() ;j++)
		{
			ret[i][j] *= 1.-(previousstate[i]) ;
		}
	}

	return ret ;
}

bool AnisotropicLinearDamage::fractured() const
{
	if (fraction < 0)
		return false ;
	
	if(inTension)
		if(std::min(tensionDamagey, std::min(tensionDamagex,tensionDamagez)) >= secondaryThresholdDamageDensity/fraction)
			return true ;
		
	if(inCompression)
		if(compressionDamage >= thresholdDamageDensity/fraction)
			return true ;
		
		return false ;
// 	std::cout << std::max(tensionDamage, compressionDamage) <<  " " << thresholdDamageDensity/**fraction*/ << std::endl ;
// 	return tensionDamagey >= secondaryThresholdDamageDensity/fraction || tensionDamagex >= secondaryThresholdDamageDensity/fraction || compressionDamage >= thresholdDamageDensity/fraction  ;
}

AnisotropicLinearDamage::~AnisotropicLinearDamage()
{
}


}
