//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUCONFINEDMOHRCOULOMBSTRAIN_H
#define MUCONFINEDMOHRCOULOMBSTRAIN_H

#include "fracturecriterion.h"

namespace Amie {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when the maximum principal stresses is below or above the prescribed limits 
	
*/
/*PARSE ConfinedMohrCoulombWithStrainLimit FractureCriterion
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
    @value[tensile_strain] // maximum strain in tension (positive)
*/
class ConfinedMohrCoulombWithStrainLimit : public FractureCriterion
{

public:
	double upVal ;
	double downVal ;
	double strainLimit ;
	bool metInCompression  ;
	bool metInTension  ;
	virtual bool directionInTension(size_t direction, double t = 0) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return metInTension ;}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	ConfinedMohrCoulombWithStrainLimit(double up, double down, double strainLimit);

	virtual ~ConfinedMohrCoulombWithStrainLimit();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual void scale(double d) {upVal *= d ; downVal *= d ; };
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ; } ;
};

}

#endif
