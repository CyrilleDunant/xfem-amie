//
// C++ Interface: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUCONFINEDVONMISES_H
#define MUCONFINEDVONMISES_H

#include "fracturecriterion.h"

namespace Mu {

	/** \brief The von Mises fracture criterion is met when the vonMises stress reaches a threshold level
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class ConfinedVonMises : public FractureCriterion
	{
	public:
		double thresholdup ;
		double thresholddown ;
	public:
	/** \brief Constructor 
	 * @param thres Set the maximum stress. 
	 */
		ConfinedVonMises(double threshup, double theshdown);
	
		virtual ~ConfinedVonMises();

	/** \brief Return a copy of this criterion
	 */
		virtual FractureCriterion * getCopy() const;

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
		virtual double grade(const ElementState &s)  ;

		virtual Material toMaterial() ;
	};
	

} ;

#endif