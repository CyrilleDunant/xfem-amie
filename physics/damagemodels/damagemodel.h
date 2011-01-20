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
	double secondaryThresholdDamageDensity ;
	double damageDensityIncrement ;
	double fraction ;
public:
	
	Vector state ;
	Vector previousstate ;
	bool isNull ;
    
	
	DamageModel(double characteristicRadius);
	
	double getThresholdDamageDensity() const;
	
	double getSecondaryThresholdDamageDensity() const;
	
	double getCharacteristicRadius() const ;
	
	double getDamageDensityIncrement() const ;

	void setThresholdDamageDensity(double d);
	
	void setSecondaryThresholdDamageDensity(double d) ;
	
	void setMaterialCharacteristicRadius(double d) ;
	
	void setDamageDensityIncrement(double d) ;
	
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

	/** \brief Get previous damage value (if stored in the model) */
	virtual Vector getPreviousDamage() = 0 ;

	/** \brief Get previous previous damage value (if stored in the model) */
	virtual Vector getPreviousPreviousDamage() = 0 ;

	/** \brief Impose previous and previous previous damage (if stored in the model) */
	virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) = 0 ;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const = 0 ;

	/** \brief Modify the rigidity matrix according to the damage model
	 * 
	 * @param m Matrix to modify
	 * @return a new Matrix
	 */
	virtual Matrix apply(const Matrix & m) const = 0 ;
	
	virtual Matrix applyPrevious(const Matrix & m) const = 0 ;

	virtual DamageModel * getCopy() const = 0 ;
	/** \brief Modify the rigidity matrix according to the damage model function (for space-time behaviour)
	 * 
	 * @param m Matrix to modify
	 * @return a new FunctionMatrix
	 */
	
	virtual bool changed() const ;

} ;

/** \brief Null damage model. A material with this damage model can never be damaged.*/
class NullDamage : public DamageModel
{
public:

	NullDamage() : DamageModel(0) { state.resize(0); previousstate.resize(0);} ;

	/** \brief Return the vector of variables describing the damage state
	 * 
	 * @return damage state vector
	 */
	virtual const Vector & damageState() const { return state ;} ;
	virtual Vector & damageState() { return state ;} ;

	/** \brief Do nothing
	 * 
	 * @param s ElementState
	 */
	virtual void step(ElementState & s) { } ;

	/** \brief Do nothing
	 * 
	 * @param d damage
	 */
	virtual void artificialDamageStep(double d) { } ;

	/** \brief returns 0 */
	virtual Vector getPreviousDamage() {return Vector(0) ; } ;

	/** \brief returns 0 */
	virtual Vector getPreviousPreviousDamage() {return Vector(0) ; } ;

	/** \brief Do nothing */
	virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) { } ;

	/** \brief return a copy of the Matrix given as argument
	 * 
	 * @param m Matrix
	 * @return m
	 */
	virtual Matrix apply(const Matrix & m) const {return m ;} ;
	virtual Matrix applyPrevious(const Matrix & m) const {return m ;} ;
	
	virtual bool fractured() const {return false ; } ;
	
	virtual DamageModel * getCopy() const { return new NullDamage() ;}

// 	virtual FunctionMatrix apply(const Matrix & m) const
// 	{
// 		FunctionMatrix ret(m.numRows(), m.numCols()) ;
// 		ret += m ;
// 		return ret ;
// 	}
} ;

}

#endif
