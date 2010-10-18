#ifndef REBAR_BEHAVIOUR_H
#define REBAR_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct RebarBehaviour : public WeibullDistributedStiffness
	{
		RebarBehaviour(double E=59e9, double nu=0.3, double tensile=570000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // PASTE_BEHAVIOUR