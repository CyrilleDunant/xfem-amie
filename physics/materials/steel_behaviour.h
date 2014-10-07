// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STEEL_BEHAVIOUR_H
#define STEEL_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
	struct SteelBehaviour : public WeibullDistributedStiffness
	{
		SteelBehaviour(double E=210e9, double nu=0.27, double tensile = 276e6, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // PASTE_BEHAVIOUR
