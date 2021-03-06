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
#ifndef MU_STRAIN_BROKEN_ISOTROPIC_LINEARDAMAGE_H
#define MU_STRAIN_BROKEN_ISOTROPIC_LINEARDAMAGE_H

#include "damagemodel.h"

namespace Amie {

/** \brief Isotropic linear damage model. The stifness of an affected element is scaled by a factor between 1 and 1 - .9999999
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class StrainBrokenIsotropicLinearDamage final: public DamageModel
{
protected:
    double eps ;
    double limitStrain ;
public:
    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    StrainBrokenIsotropicLinearDamage(int numDof, double limitStrain) ;

    virtual ~StrainBrokenIsotropicLinearDamage();

    /** \brief Increment the damage
     *
     * The formula used is \f$ d += .1e^{\frac{A_e}{\pi\epsilon^2}} ; \f$
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;
    virtual void computeDelta(ElementState & s) ;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;
    virtual DamageModel * getCopy() const ;

};

}

#endif
