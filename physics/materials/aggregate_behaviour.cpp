// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "aggregate_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

AggregateBehaviour::AggregateBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, -8*tensile,tensile)
{
	materialRadius = 0.0005 ;
}

Form * AggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalExponentiallyDecreasingMohrCoulomb(up*factor,down*factor, 3.*up*factor/E, 2.*down*factor/E, E)) ;
	ret->criterion->setMaterialCharacteristicRadius(materialRadius);
// 	ret->dfunc->setThresholdDamageDensity(1.);
	return ret ;
}

ElasticOnlyAggregateBehaviour::ElasticOnlyAggregateBehaviour(double E, double nu, SpaceDimensionality dim) : AggregateBehaviour(E,nu,0.,dim)
{

}

Form * ElasticOnlyAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	return new Stiffness(param*factor) ;
}



