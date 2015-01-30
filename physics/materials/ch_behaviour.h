// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CH_BEHAVIOUR_H
#define CH_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
	struct CHBehaviour : public Stiffness
	{
        // per http://www.sciencedirect.com/science/article/pii/S1359645408008720
		CHBehaviour(double E=35.4e9, double nu=0.28, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // CH_BEHAVIOUR_H
