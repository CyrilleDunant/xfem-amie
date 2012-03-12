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
	bool metInCompression  ;
	bool metInTension  ;
public:
	
	virtual bool directionInTension(size_t direction) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction) {return metInTension ;}

	/** \brief Constructor 
	 * @param up Set the maximum strain. 
	 */
	MaximumStrain(double up, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

	virtual ~MaximumStrain();

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; strain}{max\; Mises; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
	virtual double grade(ElementState &s)  ;
 	
	/** \brief Return a copy of this criterion
	 */
	virtual FractureCriterion * getCopy() const;

	virtual Material toMaterial() ;
	
	virtual void scale(double d) {upVal *= d ;};
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal*20e9 ; } ;
};

}

#endif
