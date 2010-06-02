//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MAX_STRAIN_H__
#define MAX_STRAIN_H__

#include "fracturecriterion.h"
#include "../../mesher/delaunay_3d.h"

namespace Mu {

/** \brief Maximum strain fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The maximum (tensile) strain criterion is met when a strain limit is reached.
	
*/
class MaximumStrain : public FractureCriterion
{
	double upVal ;
public:
	/** \brief Constructor 
	 * @param up Set the maximum strain. 
	 */
	MaximumStrain(double up);

	virtual ~MaximumStrain();

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; strain}{max\; Mises; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
	virtual double grade(const ElementState &s) const  ;
 	
	/** \brief Return a copy of this criterion
	 */
	virtual FractureCriterion * getCopy() const;

	virtual Material toMaterial() ;
};

}

#endif
