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

#ifndef __DAMAGE_MODEL_H__
#define __DAMAGE_MODEL_H__

#include "../../elements/integrable_entity.h"
#include "../../utilities/matrixops.h"

namespace Mu
{

/** \brief Damage model interface */
class DamageModel
{
protected:
	double characteristicRadius ;
	double thresholdDamageDensity ;
	double damageDensityIncrement ;
	double fraction ;
public:
	DamageModel(double characteristicRadius) : characteristicRadius(characteristicRadius)
	{ 
		isNull = true ; 
		thresholdDamageDensity = .2 ;
		damageDensityIncrement = .1 ;
		fraction = -1 ;
	} ;
	
	void setThresholdDamageDensity(double d)
	{
		thresholdDamageDensity = d ;
	}
	
	void setDamageDensityIncrement(double d)
	{
		damageDensityIncrement = d ;
	}
	
	/** \brief Return a vector of values describing the damage stage of the material
	 * 
	 * @return a Vector
	 */
	virtual const Vector & damageState() const = 0 ;
	virtual Vector & damageState() = 0 ;

	/** \brief Increment the damage from the current state of the element considered
	 * 
	 * @param s ElementState
	 */
	virtual void step(ElementState & s) = 0 ;

	/** \brief Increment the damage from an external value
	 * 
	 * @param d damage
	 */
	virtual void artificialDamageStep(double d) = 0 ;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const = 0 ;

	/** \brief Modify the rigidity matrix according to the damage model
	 * 
	 * @param m Matrix to modify
	 * @return a new Matrix
	 */
	virtual Matrix apply(const Matrix & m) const = 0 ;

	/** \brief Modify the rigidity matrix according to the damage model function (for space-time behaviour)
	 * 
	 * @param m Matrix to modify
	 * @return a new FunctionMatrix
	 */
// 	virtual FunctionMatrix apply(const Matrix & m) const = 0 ;
	bool isNull ;
} ;

/** \brief Null damage model. A material with this damage model can never be damaged.*/
class NullDamage : public DamageModel
{
protected:
	Vector state ;
public:

	NullDamage() : DamageModel(0), state(0) { } ;

	/** \brief Return an empty vector
	 * 
	 * @return an empty vector
	 */
	virtual const Vector & damageState() const { return state ;} ;

	/** \brief Do nothing
	 * 
	 * @param s ElementState
	 */
	virtual void step(ElementState & s) { } ;

	/** \brief Do nothing
	 * 
	 * @param d damage
	 */
	virtual void artificialStep(double d) { } ;

	/** \brief return a copy of the Matrix given as argument
	 * 
	 * @param m Matrix
	 * @return m
	 */
	virtual Matrix apply(const Matrix & m) const {return m ;} ;

// 	virtual FunctionMatrix apply(const Matrix & m) const
// 	{
// 		FunctionMatrix ret(m.numRows(), m.numCols()) ;
// 		ret += m ;
// 		return ret ;
// 	}
} ;

}

#endif
