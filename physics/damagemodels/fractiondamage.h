//
// C++ Interface: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef FRACTIONDAMAGE_H
#define FRACTIONDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/** \brief Linear Damage model, with different compressive and tensile damage effects
 * This damage allows to simulate a mix of two materials, one which fails and one which doesn't
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class FractionLinearDamage : public DamageModel
{
public:

	bool inCompression ;
	bool inTension ;
	Matrix remnant ;
	double phi ;
	/** \brief Constructor, set the number of degrees of freedom and a strain limit for failure
	 * 
	 * @param numDof 
	 * @param threshold 
	 */
	FractionLinearDamage(int numDof, double characteristicRadius, Matrix remnant, double phi) ;

	virtual ~FractionLinearDamage();

	/** \brief return the damage state.
	 * 
	 * @return two-element damage state
	 */
	virtual const Vector & damageState() const ;
	virtual Vector & damageState() ;


	/** \brief Increment the damage from an external value.
	 * 
	 * @param d damage
	 */
	virtual void artificialDamageStep(double d) ;

	/** \brief returns 0 */
	virtual Vector getPreviousDamage() {return Vector(0) ; } ;

	/** \brief returns 0 */
	virtual Vector getPreviousPreviousDamage() {return Vector(0) ; } ;

	/** \brief Do nothing */
	virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) { } ;

	/** \brief Increment the damage. 
	 * The formula used for the increment varies depending on whether the element is 
	 * in tension or compression. Once it has been determined whether we are in tension or 
	 * compression, the corresponding damage value is modified.
	 * 
	 * @param s ElementState to consider
	 */
	virtual Vector computeDamageIncrement(ElementState & s) ;

	/** \brief compute the stiffness matrix from the damage state
	 * 
	 * \f$ K_{ij} = K_{ij}*d_i \f$
	 * @param m original Matrix
	 * @return modified Matrix
	 */
	virtual Matrix apply(const Matrix & m) const;
	
	virtual Matrix applyPrevious(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const ;
	
	virtual DamageModel * getCopy() const { return new FractionLinearDamage(state.size()-1, getCharacteristicRadius(), remnant, phi) ;}
	
};

}

#endif
