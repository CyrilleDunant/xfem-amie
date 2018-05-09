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
#ifndef MUNLDRUCKERPAGER_H
#define MUNLDRUCKERPAGER_H

#include "fracturecriterion.h"

namespace Amie {

/** \brief The von Mises fracture criterion is met when the vonMises stress reaches a threshold level

	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class DruckerPrager : public FractureCriterion
{
    bool metInTension ;
    bool metInCompression ;
    bool inTension ;
public:
    double upthreshold ;
    double downthreshold ;
    double friction ;
    double modulus ;
    double cap ;
    virtual bool directionInTension(size_t direction, double t = 0)
    {
        return inTension ;
    }
    virtual bool directionInCompression(size_t direction, double t = 0)
    {
        return !inTension ;

    }
    virtual bool directionMet(size_t direction, double t = 0)
    {
        if(direction == 0)
            return metInTension ;
        return metInCompression ;
    }
public:
    /** \brief Constructor
     * @param thres Set the maximum stress.
     */
    DruckerPrager(double downthres,double upthres,double modulus,  double friction, double radius);

    virtual ~DruckerPrager();

    /** \brief Return a copy of this criterion
     */
    virtual FractureCriterion * getCopy() const;

    /** \brief Return normalised distance to the fracture surface
     *
     * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
     * @param s ElementState to consider
    */
    virtual double grade(ElementState &s)  ;

    virtual bool fractured() {
        return false ;
    }

};

}

#endif
