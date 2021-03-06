//
// C++ Interface: vonmises
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUBoundVONMISES_H
#define MUBoundVONMISES_H

#include "fracturecriterion.h"
#include "../damagemodels/damagemodel.h"

namespace Amie {

/** \brief The von Mises fracture criterion is met when the vonMises stress reaches a threshold level

	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
/*PARSE BoundedVonMises FractureCriterion
    @value[tensile_strength] // value of the tensile strength of the material
    @value[maximum_damage] // value of damage above which an element is considered broken
*/
class BoundedVonMises : public FractureCriterion
{
public:
    double threshold ;
    double damageThreshold ;
    DamageModel * dmodel ;
public:
    /** \brief Constructor
     * @param thres Set the maximum stress.
     */
    BoundedVonMises(double thres, double damageThreshold);

    virtual ~BoundedVonMises();

    /** \brief Return a copy of this criterion
     */
    virtual FractureCriterion * getCopy() const;

    /** \brief Return normalised distance to the fracture surface
     *
     * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
     * @param s ElementState to consider
    */
    virtual double grade(ElementState &s)  ;

};


}

#endif
