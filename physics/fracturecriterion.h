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

	/**
	Abstract definition of a fracture criterion
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class FractureCriterion
	{
	public:
		FractureCriterion() ;
	
		virtual ~FractureCriterion();
		
		virtual bool met(const ElementState &) const = 0 ;
		
		virtual FractureCriterion * getCopy() const;
	
	};

} ;

#endif
