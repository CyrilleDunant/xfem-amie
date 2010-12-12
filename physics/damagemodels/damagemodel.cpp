//
// C++ Interface: damagemodel
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "damagemodel.h"

namespace Mu
{

/** \brief Damage model interface */

	
	DamageModel::DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		isNull = true ; 
		thresholdDamageDensity = .999999999999 ;
		secondaryThresholdDamageDensity = .99999999999999 ;
		damageDensityIncrement = .0001 ;
		fraction = -1 ;
	} ;
	
	double DamageModel::getThresholdDamageDensity() const
	{
		return thresholdDamageDensity ;
	}
	
	double DamageModel::getSecondaryThresholdDamageDensity() const
	{
		return secondaryThresholdDamageDensity ;
	}
	
	double DamageModel::getCharacteristicRadius() const
	{
		return characteristicRadius ;
	}
	
	double DamageModel::getDamageDensityIncrement() const
	{
		return damageDensityIncrement ;
	}
	void DamageModel::setThresholdDamageDensity(double d)
	{
		thresholdDamageDensity = d ;
	}
	
	void DamageModel::setSecondaryThresholdDamageDensity(double d)
	{
		secondaryThresholdDamageDensity = d ;
	}
	
	void DamageModel::setMaterialCharacteristicRadius(double d)
	{
		characteristicRadius = d ;
	}
	
	void DamageModel::setDamageDensityIncrement(double d)
	{
		damageDensityIncrement = d ;
	}
} ;	
