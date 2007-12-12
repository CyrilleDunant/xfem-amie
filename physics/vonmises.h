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
#ifndef MUVONMISES_H
#define MUVONMISES_H

#include "fracturecriterion.h"

namespace Mu {

	/**
	The von Mises fracture criterion is met when the vonMises stress reaches a threshold level
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class VonMises : public FractureCriterion
	{
		double threshold ;
	public:
		VonMises(double thres);
	
		virtual ~VonMises();
		
		virtual bool met(const ElementState & s) const ;
	
	};

} ;

#endif
