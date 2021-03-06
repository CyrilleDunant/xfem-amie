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
#ifndef MU_SPACETIME_FIBERLINEARDAMAGE_H
#define MU_SPACETIME_FIBERLINEARDAMAGE_H

#include "damagemodel.h"

namespace Amie {

/*PARSE SpaceTimeFiberBasedIsotropic DamageModel
    @value[damage_increment] 0.1 // damage increment which is applied at each step of the damage algorithm
    @value[time_tolerance] 0.001 // minimum time between two successive damage events
    @value[maximum_damage] 0.6 // damage above which an element is considered broken
*/
class SpaceTimeFiberBasedIsotropicLinearDamage : public DamageModel
{
protected:
    double fibreFraction ;
    double timeTolerance ;
    Function visc ;

public:
    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    SpaceTimeFiberBasedIsotropicLinearDamage(double f = 0.1, double tol = 0.001, double cutoff = 0.6) ;

    virtual ~SpaceTimeFiberBasedIsotropicLinearDamage();

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

    void setLogitViscousDamageLaw(double a, double b, double c) ;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    virtual Matrix applyViscous(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    virtual void postProcess() ;

    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;

    virtual DamageModel * getCopy() const ;
};

}

#endif
