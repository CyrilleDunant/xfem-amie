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
#ifndef MUFIBERLINEARDAMAGE_H
#define MUFIBERLINEARDAMAGE_H

#include "damagemodel.h"

namespace Amie {

/** \brief Linear Damage model, with different compressive and tensile damage effects
 * This damage allows a better simulation of materials exhibitng different failure behaviour in traction and compression.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
/*PARSE FiberBasedIsotropic : DamageModel
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[maximum_damage] 0.6 // damage above which an element is considered broken
*/
class FiberBasedIsotropicLinearDamage : public DamageModel
{
protected:
    double fibreFraction ;
    bool alt ; 
public:
    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    FiberBasedIsotropicLinearDamage(double f = 0.1, double cutoff = 0.6) ;

    virtual ~FiberBasedIsotropicLinearDamage();

    /** \brief Increment the damage
     *
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;
    virtual void computeDelta(ElementState & s) ;

    /** \brief Increment the damage from the current state of the element considered
     *
     * @param s ElementState
     */
    virtual void step(ElementState & s, double maxscore)  ;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured() const  ;

    virtual DamageModel * getCopy() const ; 

};

}

#endif
