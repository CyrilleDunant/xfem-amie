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
#ifndef GEO_BEASED_CONTACT_H
#define GEO_BEASED_CONTACT_H

#include "collisiondetector.h"

namespace Amie {

	/** \brief Detect collisions from penetration to a geometry
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/


	class GeometryBasedContact : public CollisionDetector
	{
    public:    
        Geometry *geo ;
        
	public:
	/** \brief Constructor 
	 * @param thres Set the maximum stress. 
	 */
		GeometryBasedContact(Geometry *geo);
	
		virtual ~GeometryBasedContact();

		virtual double grade(ElementState &s)  ;
		
		virtual void scale(double d ) { }
		
		virtual FractureCriterion * getCopy() const;
    };

} 

#endif
