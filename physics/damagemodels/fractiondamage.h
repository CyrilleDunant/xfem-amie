//
// C++ Interface: lineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
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
	FractionLinearDamage(Matrix remnant, double phi) ;

	virtual ~FractionLinearDamage();

	/** \brief Increment the damage from an external value.
	 * 
	 * @param d damage
	 */
	virtual void artificialDamageStep(double d) ;

	/** \brief Do nothing */
	virtual void artificialPreviousDamage(Vector previous, Vector previousprevious) { } ;

	/** \brief Increment the damage. 
	 * The formula used for the increment varies depending on whether the element is 
	 * in tension or compression. Once it has been determined whether we are in tension or 
	 * compression, the corresponding damage value is modified.
	 * 
	 * @param s ElementState to consider
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

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
	
	virtual DamageModel * getCopy() const { return new FractionLinearDamage(remnant, phi) ;}
	
};

}

#endif
