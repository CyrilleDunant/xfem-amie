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
#ifndef MUANISOLINEARDAMAGE_H
#define MUANISOLINEARDAMAGE_H

#include "damagemodel.h"

namespace Amie {

/** \brief Linear Damage model, with different compressive and tensile damage effects
 * This damage allows a better simulation of materials exhibitng different failure behaviour in traction and compression.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class AnisotropicLinearDamage : public DamageModel
{
public:

	/** \brief Constructor, set the number of degrees of freedom and a strain limit for failure
	 * 
	 * @param characteristicRadius material characteristic radius
	 */
	AnisotropicLinearDamage() ;

	virtual ~AnisotropicLinearDamage();

	/** \brief Increment the damage. 
	 * The formula used for the increment varies depending on whether the element is 
	 * in tension or compression. Once it has been determined whether we are in tension or 
	 * compression, the corresponding damage value is modified.
	 * 
	 * @param s ElementState to consider
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;
	
	virtual void computeDelta(ElementState & s) ;

	/** \brief compute the stiffness matrix from the damage state
	 * 
	 * \f$ K_{ij} = K_{ij}*d_i \f$
	 * @param m original Matrix
	 * @return modified Matrix
	 */
	virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const ;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured(int direction = -1) const ;
	
	virtual DamageModel * getCopy() const ; 

};

}

#endif
