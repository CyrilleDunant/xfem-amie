//
// C++ Interface: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ARTIFICIAL_FRACTURE_H
#define ARTIFICIAL_FRACTURE_H

#include "fracturecriterion.h"

namespace Amie {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

	/**
	Artificial only fracture criterion fracture criterion
	
	*/
	class ArtificialFracture : public FractureCriterion
	{
	public:
		ArtificialFracture(MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);
		virtual void initialiseCache(const ElementState & s) { } ;
	
		virtual ~ArtificialFracture();
		
		/** \brief Return false: artificial damage only
		 * 
		 * @param s ElementState ton consider
		 */
		virtual bool met(const ElementState & s) {return false ; } ;

		/** \brief Return a normalised distance to the fracture surface, 
		 * 
		 * The returned value lies between 0 and 1
		 * @param  ElementState to consider
		 * @return a value between 0 and 1
		 */
		virtual double grade(const ElementState & s) = 0 ;
		
		/** \brief Produce a copy of the fracture criterion
		 * 
		 * @return a new FractureCriterion
		 */
		virtual FractureCriterion * getCopy() const = 0;

		/** \brief set the neighbourhood in which to consider other elements to check for failure.
		 * 
		 * @param r new radius
		 */
		virtual void setNeighbourhoodRadius(double r) { } ;
		
		virtual double getTensileLimit(const ElementState & s) const {return 2e6 ; } ;
	
	};

} ;

#endif
