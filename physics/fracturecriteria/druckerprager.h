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
#ifndef MUNLDRUCKERPAGER_H
#define MUNLDRUCKERPAGER_H

#include "fracturecriterion.h"

namespace Mu {

	/** \brief The von Mises fracture criterion is met when the vonMises stress reaches a threshold level
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class DruckerPrager : public FractureCriterion
	{
	public:
		double threshold ;
		double friction ;
	public:
	/** \brief Constructor 
	 * @param thres Set the maximum stress. 
	 */
		DruckerPrager(double thres, double friction, double radius, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);
	
		virtual ~DruckerPrager();

	/** \brief Return a copy of this criterion
	 */
		virtual FractureCriterion * getCopy() const;

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
		virtual double grade(ElementState &s)  ;
		
		virtual bool fractured() {return false ; }

		virtual Material toMaterial() ;
		
		virtual void scale(double d ) { threshold *= d ;}
	};

} ;

#endif
