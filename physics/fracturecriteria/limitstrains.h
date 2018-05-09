//
// C++ Interface: limit strains
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LIMIT_STRAIN_H__
#define LIMIT_STRAIN_H__

#include "fracturecriterion.h"
#include "../../mesher/delaunay_3d.h"

namespace Amie {

/** \brief Maximum strain fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The maximum (tensile) strain criterion is met when a strain limit is reached.
	
*/
/*PARSE LimitStrains FractureCriterion
    @value[tensile_strain] // maximum strain in tension (positive)
    @value[compressive_strain] // maximum strain in compression (negative)
*/
class LimitStrains : public FractureCriterion
{
	double maxUpVal ;
	double maxDownVal ;
public:
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction, double t = 0) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return metInTension ;}
	/** \brief Constructor 
	 * @param up Set the maximum strain. 
	 */
	LimitStrains(double maxdown, double maxup);

	virtual ~LimitStrains();

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; strain}{max\; Mises; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
	virtual double grade(ElementState &s)  ;
 	
	/** \brief Return a copy of this criterion
	 */
	virtual FractureCriterion * getCopy() const;


};

}

#endif
