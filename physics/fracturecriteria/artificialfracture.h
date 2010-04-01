//
// C++ Interface: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ARTIFICIAL_FRACTURE_H
#define ARTIFICIAL_FRACTURE_H

#include "fracturecriterion.h"

namespace Mu {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

	/**
	Artificial only fracture criterion fracture criterion
	
	*/
	class ArtificialFracture : public FractureCriterion
	{
	public:
		ArtificialFracture() ;
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
		virtual double grade(const ElementState & s) const = 0 ;
		
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
	
	};

} ;

#endif