//
// C++ Interface: ruptureenergy
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ENR_CRIT_H
#define ENR_CRIT_H

#include "fracturecriterion.h"

namespace Mu {

/** \brief This FractureCriterion is met when a certain (elastic) Energy density is reached
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	
*/
class RuptureEnergy : public FractureCriterion
{
	

public:

	double energy ;

/** \brief Constructor. Set the threshold energy density
*
* @param energy Threshold energy density
*/
	RuptureEnergy(double energy);

	virtual ~RuptureEnergy();

/** \brief return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief Return the normalised distance to the fracture surface.
*
* The distance is computed thus: \f$ 1 - \frac{E_c A_e}{\int_e \sigma \epsilon de} \f$
* @param s ElementState to consider
*/
	virtual double grade(const ElementState &s) const  ;


	virtual Material toMaterial() ;
};

}

#endif
