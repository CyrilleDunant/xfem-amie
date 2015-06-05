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
#ifndef MULINEARDAMAGE_H
#define MULINEARDAMAGE_H

#include "damagemodel.h"

namespace Amie {

/** \brief Linear Damage model, with different compressive and tensile damage effects
 * This damage allows a better simulation of materials exhibitng different failure behaviour in traction and compression.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class LinearDamage final: public DamageModel
{
public:

    bool inCompression ;
    bool inTension ;
    /** \brief Constructor, set the number of degrees of freedom and a strain limit for failure
     *
     * @param numDof
     * @param threshold
     */
    LinearDamage() ;

    virtual ~LinearDamage();

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
     * \f$ K_{ij} = K_{ij}*(1-d_i) \f$
     * @param m original Matrix
     * @return modified Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured() const ;

    virtual DamageModel * getCopy() const ;

    virtual void computeDelta(ElementState & s) ;

};

}

#endif
