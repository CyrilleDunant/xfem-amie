#include "aggregate_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

AggregateBehaviour::AggregateBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim), -8*tensile,tensile)
{
	materialRadius = 0.0002 ;
	neighbourhoodRadius = 0.0003 ;
}

Form * AggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new MohrCoulomb(up*factor,down*factor)) ;
	return ret ;
}



