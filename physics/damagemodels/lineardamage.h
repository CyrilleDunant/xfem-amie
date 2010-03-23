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
#ifndef MULINEARDAMAGE_H
#define MULINEARDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/** \brief Linear Damage model, with different compressive and tensile damage effects
 * This damage allows a better simulation of materials exhibitng different failure behaviour in traction and compression.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class LinearDamage : public DamageModel
{
protected:
	Vector state ;
public:

	/** \brief Constructor, set the number of degrees of freedom and a strain limit for failure
	 * 
	 * @param numDof 
	 * @param threshold 
	 */
	LinearDamage(int numDof) ;

	virtual ~LinearDamage();

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

	/** \brief Increment the damage. 
	 * The formula used for the increment varies depending on whether the element is 
	 * in tension or compression. Once it has been determined whether we are in tension or 
	 * compression, the corresponding damage value is modified.
	 * 
	 * @param s ElementState to consider
	 */
	virtual void step(ElementState & s) ;

	/** \brief compute the stiffness matrix from the damage state
	 * 
	 * \f$ K_{ij} = K_{ij}*d_i \f$
	 * @param m original Matrix
	 * @return modified Matrix
	 */
	virtual Matrix apply(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const ;
	
};

}

#endif
