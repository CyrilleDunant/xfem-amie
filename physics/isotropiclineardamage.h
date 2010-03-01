//
// C++ Interface: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_ISOTROPIC_LINEARDAMAGE_H
#define MU_ISOTROPIC_LINEARDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class IsotropicLinearDamage : public DamageModel
{
protected:
	Vector state ;
	double eps ;
public:
	/** \brief Constructor. Set the number of degrees of freedom
	 * 
	 * @param numDof number of degrees of freedom
	 */
	IsotropicLinearDamage(int numDof) ;

	virtual ~IsotropicLinearDamage();

	/** \brief Return the single-element Vector in which is stored the scaling factor
	 * 
	 * @return single-element Vector in which is stored the scaling factor
	 */
	virtual const Vector & damageState() const ;

	/** \brief Increment the damage
	 * 
	 * The formula used is \f$ d += .1e^{\frac{A_e}{\pi\epsilon^2}} ; \f$
	 * @param s ElementState passed as a parameter
	 */
	virtual void step(ElementState & s) ;

	/** \brief compute the new stifness matrix after damage
	 * 
	 * \f$ K' = K(1-d) \f$
	 * @param m Matrix to modify
	 * @return the new Matrix
	 */
	virtual Matrix apply(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const  ;
};

}

#endif
