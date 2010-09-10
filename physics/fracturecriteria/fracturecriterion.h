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

#include "../../elements/integrable_entity.h"

namespace Mu {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

	/**
	Abstract definition of a fracture criterion
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class FractureCriterion
	{
	protected:
		std::vector<DelaunayTriangle *> cache ;
		std::vector<DelaunayTetrahedron *> cache3d ;
		std::vector<double> area ;
		double neighbourhoodradius ;
		double physicalCharacteristicRadius ;
		double scoreAtState ;
		
	public:
		
		
		
		FractureCriterion() ;
		virtual void initialiseCache(const ElementState & s) ;
	
		virtual ~FractureCriterion();
		
		void step(const Mu::ElementState& s) ;
		
		/** \brief Return true if the fracture criterion is met
		 * 
		 * @param s ElementState ton consider
		 * @return true if the fracture criterion is met
		 */
		virtual bool met(const ElementState & s) ;

		/** \brief Return a normalised distance to the fracture surface, 
		 * 
		 * The returned value lies between -1 and 1
		 * @param  ElementState to consider
		 * @return a value between -1 and 1
		 */
		virtual double grade(const ElementState & s) = 0 ;
		virtual double getSteppedScore() const {return scoreAtState ;} ;
		
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
		double getNeighbourhoodRadius() const { return neighbourhoodradius ;} ;

		virtual void setMaterialCharacteristicRadius(double r) ;
		double getMaterialCharacteristicRadius() const { return physicalCharacteristicRadius ;} ;
		
		virtual Material toMaterial() ;
		
		const std::vector<DelaunayTriangle *> & getCache() const { return cache ; } ;
		const std::vector<DelaunayTetrahedron *> & getCache3d() const { return cache3d ; } ;
	
	};

} ;

#endif
