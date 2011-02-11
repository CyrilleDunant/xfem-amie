#ifndef AGGREGATE_BEHAVIOUR_H
#define AGGREGATE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct AggregateBehaviour : public WeibullDistributedStiffness
	{
		AggregateBehaviour(double E=59e9, double nu=0.3, double tensile=5.7e6, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;
	
	struct ElasticOnlyAggregateBehaviour : public AggregateBehaviour
	{
		ElasticOnlyAggregateBehaviour(double E=59e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	
	} ;

} ;

#endif // PASTE_BEHAVIOUR
