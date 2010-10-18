#ifndef STEEL_BEHAVIOUR_H
#define STEEL_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct SteelBehaviour : public WeibullDistributedStiffness
	{
		SteelBehaviour(double E=59e9, double nu=0.3, double tensile=570000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;

} ;

#endif // PASTE_BEHAVIOUR