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

using namespace Amie ;

AggregateBehaviour::AggregateBehaviour(bool e, bool st, double E, double nu, double up, double r, SpaceDimensionality dim, planeType pt, double var, int blocks)  : WeibullDistributedStiffness(E,nu, dim, -8000.*up,up, pt, var, r), freeblocks(blocks), spaceTime(st), elastic(e)
{

}

Form * AggregateBehaviour::getCopy() const 
{
        std::default_random_engine generator(std::rand());
        std::weibull_distribution< double > distribution(5, 1);
        double weib = distribution(generator) ;
	double factor = 1. - variability + variability*weib ;

	Form * copy ;

	if(spaceTime)
	{
		if(!elastic)
		{
			copy = new ViscoelasticityAndFracture( PURE_ELASTICITY, param*factor, new SpaceTimeNonLocalMaximumStrain(up/E), new IsotropicLinearDamage(), 2+freeblocks )  ;
			copy->getFractureCriterion()->setScoreTolerance(1e-4) ;
			copy->getDamageModel()->setDamageDensityTolerance(0.05) ;
			copy->getDamageModel()->setThresholdDamageDensity(0.8) ;
			copy->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius) ;
			copy->getDamageModel()->setNeedGlobalMaximumScore(true) ;
		}
		else
		{
			copy = new Viscoelasticity( PURE_ELASTICITY, param*factor, 2+freeblocks )  ;
		}
	}
	else
	{
		if(!elastic)
		{
			copy = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,-8000.*up*factor, E), new FiberBasedIsotropicLinearDamage(0.1,0.699)) ;
			copy->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
		}
		else
		{
			copy = new Stiffness(param*factor) ;
		}
	}

//	return new Stiffness(param*factor) ;
//     StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalLinearlyDecreasingMohrCoulomb(up*factor,-8000.*up*factor, 3.*factor*up/E, -factor*24000.*up/E,E), new FiberBasedIsotropicLinearDamage(0.1,0.95)) ;

	return copy ;
}


