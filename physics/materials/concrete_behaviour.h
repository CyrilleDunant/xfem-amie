// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONCRETE_BEHAVIOUR_H
#define CONCRETE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct ConcreteBehaviour : public WeibullDistributedStiffness
	{
		bool reinforced ;
		ConcreteBehaviour(double E=37e9, double nu=0.3, double tensile=2e6, double compressive = -37e6, planeType pt = PLANE_STRESS, bool reinforced = false,  SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, MirrorState mirroring = NO_MIRROR, double dx = 0, double dy = 0, double dz = 0) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // PASTE_BEHAVIOUR
