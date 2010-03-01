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
#ifndef MUFRACTURECRITERION_H
#define MUFRACTURECRITERION_H

#include "../elements/integrable_entity.h"

namespace Mu {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

	/**
	Abstract definition of a fracture criterion
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class FractureCriterion
	{
		std::vector<DelaunayTriangle *> cache ;
		std::vector<DelaunayTetrahedron *> cache3d ;
		double eps ;
	public:
		FractureCriterion() ;
		virtual void initialiseCache(const ElementState & s) ;
	
		virtual ~FractureCriterion();
		
		/** \brief Return true if the fracture criterion is met
		 * 
		 * @param s ElementState ton consider
		 * @return true if the fracture criterion is met
		 */
		virtual bool met(const ElementState & s) ;

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
		virtual void setNeighbourhoodRadius(double r) ;
	
	};

} ;

#endif