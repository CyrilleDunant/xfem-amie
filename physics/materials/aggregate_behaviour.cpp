#include "aggregate_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

AggregateBehaviour::AggregateBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim), -8*tensile,tensile)
{
	materialRadius = 0.0005 ;
	neighbourhoodRadius = 0.002 ;
}

Form * AggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	return new Stiffness(param*factor) ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new MohrCoulomb(up*factor,down*factor), materialRadius) ;
	ret->setNeighbourhoodRadius(neighbourhoodRadius) ;
	ret->criterion->setNeighbourhoodRadius(neighbourhoodRadius);
	ret->criterion->setMaterialCharacteristicRadius(materialRadius);
	ret->dfunc->setDamageDensityIncrement(.01);
	ret->dfunc->setMaterialCharacteristicRadius(materialRadius);
	ret->dfunc->setThresholdDamageDensity(.5);
	return ret ;
}



