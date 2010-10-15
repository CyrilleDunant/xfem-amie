#ifndef PASTE_BEHAVIOUR_H
#define PASTE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct PasteBehaviour : public WeibullDistributedStiffness
	{
		PasteBehaviour(double E=12e9, double nu=0.3, double tensile=135000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // PASTE_BEHAVIOUR
