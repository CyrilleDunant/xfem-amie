// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "aggregate_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../viscoelasticity.h"
#include "../viscoelasticity_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../damagemodels/rotatingcrack.h"
#include "../damagemodels/fiberbasedisotropiclineardamage.h"
#include "../fracturecriteria/maxstrain.h"
#include "../fracturecriteria/creeprupture.h"
#include "../damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../utilities/random.h"

using namespace Amie ;

AggregateBehaviour::AggregateBehaviour(double E, double nu, double up_, double yield, double c, SpaceDimensionality dim, double var) : WeibullDistributedStiffness(E,nu, dim, 0.,0.), up(up_), yield(yield), c(c)
{
	materialRadius = 0.0004 ;
    variability = var ;
}

Form * AggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
    StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalLinearlyDecreasingMohrCoulomb(up*factor,-8000.*up*factor, 3.*factor*up/E, -factor*24000.*up/E,E), new FiberBasedIsotropicLinearDamage(0.2,0.99)) ;
//     StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalLinearlyDecreasingMohrCoulomb(up*factor,-8000.*up*factor, 3.*factor*up/E, -factor*24000.*up/E,E), new FiberBasedIsotropicLinearDamage(0.05,0.95)) ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);

	return copy ;
}

ElasticOnlyAggregateBehaviour::ElasticOnlyAggregateBehaviour(double E, double nu, SpaceDimensionality dim, double var) : AggregateBehaviour(E,nu,0.,0.,0.,dim, var)
{
//	variability = 0.5 ;
}

Form * ElasticOnlyAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	Stiffness * copy =  new Stiffness(param*factor) ;

	return copy ;
}


ViscoElasticOnlyAggregateBehaviour::ViscoElasticOnlyAggregateBehaviour(double E, double nu, SpaceDimensionality dim, double var) : AggregateBehaviour(E,nu,0.,0.,0.,dim, var), freeblocks(0)
{

}

Form * ViscoElasticOnlyAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	return new Viscoelasticity( PURE_ELASTICITY, param*factor, 2+freeblocks )  ;

}

ViscoDamageAggregateBehaviour::ViscoDamageAggregateBehaviour(double E, double nu, double up, double r, SpaceDimensionality dim, double var) : AggregateBehaviour(E,nu,up,0.,0.,dim, var), rad(r), freeblocks(0)
{
	materialRadius = r ;
}

Form * ViscoDamageAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	ViscoelasticityAndFracture * copy = new ViscoelasticityAndFracture( PURE_ELASTICITY, param*factor, new SpaceTimeNonLocalMaximumStrain(up), new IsotropicLinearDamage(), 2+freeblocks )  ;
	copy->getFractureCriterion()->setScoreTolerance(1e-4) ;
	copy->getDamageModel()->setDamageDensityTolerance(0.05) ;
	copy->getDamageModel()->setThresholdDamageDensity(0.8) ;
	copy->criterion->setMaterialCharacteristicRadius(rad) ;
	copy->getDamageModel()->setNeedGlobalMaximumScore(true) ;
	return copy ;
}
