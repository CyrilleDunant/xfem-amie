//
// C++ Interface: limit strains
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DELTA_ENERGY_H__
#define DELTA_ENERGY_H__

#include "fracturecriterion.h"
#include "../../mesher/delaunay_3d.h"
#include "../../mesher/delaunay.h"

namespace Mu {

/** \brief Maximum strain fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	This criterion computes the energy difference from damaging the element.
	
*/
class DeltaEnergy : public FractureCriterion
{
	double radius ;
	FractureCriterion * criterion ;
	std::vector<DelaunayTriangle *> cache2d ;
	std::vector<DelaunayTetrahedron *> cache3d ;
public:
	/** \brief Constructor 
	 * @param radius Radius in which the change of energy is considered
	 * @param criterion underlying fracture criterion considered for this element
	 */
	DeltaEnergy(double radius, FractureCriterion * criterion);

	virtual ~DeltaEnergy();

	/** \brief Return delta energy. 
	 *
	 * Because the fracture criterion is used as a proxy for energy
	 * this value can be normalised.
	 * @param s ElementState to consider
	*/
	virtual double grade(const ElementState &s)  ;
 	
	/** \brief Return a copy of this criterion
	 */
	virtual FractureCriterion * getCopy() const;

	virtual Material toMaterial() ;
	
	virtual void step(const ElementState &s) ;

};

}

#endif
