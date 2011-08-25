//
// C++ Interface: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_ISOTROPIC_NONLOCALLINEARDAMAGE_H
#define MU_ISOTROPIC_NONLOCALLINEARDAMAGE_H

#include "damagemodel.h"

namespace Mu {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class NonLocalIsotropicLinearDamage : public DamageModel
{
protected:
	double eps ;
public:
	/** \brief Constructor. Set the number of degrees of freedom
	 * 
	 * @param numDof number of degrees of freedom
	 */
	NonLocalIsotropicLinearDamage() ;

	virtual ~NonLocalIsotropicLinearDamage() { };
	
	/** \brief Increment the damage
	 * 
	 * @param s ElementState passed as a parameter
	 */
	virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

	/** \brief compute the new stifness matrix after damage
	 * 
	 * \f$ K' = K(1-d) \f$
	 * @param m Matrix to modify
	 * @return the new Matrix
	 */
	virtual Matrix apply(const Matrix & m) const;
	
	virtual Matrix applyPrevious(const Matrix & m) const;

	/** \brief return true is the element concerned is fractured 
		*/
	virtual bool fractured() const  ;
	
	virtual DamageModel * getCopy() const { return new NonLocalIsotropicLinearDamage() ;}
	
	virtual void postProcess() ;
};

}

#endif
