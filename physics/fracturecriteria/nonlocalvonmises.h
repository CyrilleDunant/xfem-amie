//
// C++ Interface: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUNLVONMISES_H
#define MUNLVONMISES_H

#include "fracturecriterion.h"

namespace Amie {

	/** \brief The von Mises fracture criterion is met when the vonMises stress reaches a threshold level
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/

	/*PARSE NonLocalVonMises FractureCriterion
		@value[tensile_strength] // maximum stress of the material
		@value[young_modulus] // Young modulus of the material
		@value[material_characteristic_radius] // characteristic length for the non-local damage
	*/
	class NonLocalVonMises : public FractureCriterion
	{
		bool met ;
	public:
		double threshold ;
		double E ;
		virtual bool directionInTension(size_t direction, double t = 0) {return met ;}
		virtual bool directionInCompression(size_t direction, double t = 0) {return met ;}
		virtual bool directionMet(size_t direction, double t = 0) {return met;}
	public:
	/** \brief Constructor 
	 * @param thres Set the maximum stress. 
	 */
		NonLocalVonMises(double thres, double E, double radius);
	
		virtual ~NonLocalVonMises();

	/** \brief Return a copy of this criterion
	 */
		virtual FractureCriterion * getCopy() const;

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
		virtual double grade(ElementState &s)  ;
		
		virtual void scale(double d ) { threshold *= d ;}
		
		virtual double getTensileLimit(const ElementState & s) const {return threshold ;};
	};

} 

#endif
