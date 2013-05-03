// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef AGGREGATE_BEHAVIOUR_H
#define AGGREGATE_BEHAVIOUR_H

#include "../weibull_distributed_stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct AggregateBehaviour : public WeibullDistributedStiffness
	{
		double up ;
		double yield ;
		double c ;
		AggregateBehaviour(double E=59e9, double nu=0.3, double up = 0.00075, double yield = 0.001, double c = 12000., SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
		
	} ;
	
	struct ElasticOnlyAggregateBehaviour : public AggregateBehaviour
	{
		ElasticOnlyAggregateBehaviour(double E=59e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	
	} ;

	struct ViscoElasticOnlyAggregateBehaviour : public AggregateBehaviour
	{
		ViscoElasticOnlyAggregateBehaviour(double E=59e9, double nu=0.3, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		virtual Form * getCopy() const ;
	} ;

} ;

#endif // PASTE_BEHAVIOUR