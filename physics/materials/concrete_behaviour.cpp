// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "concrete_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../fracturecriteria/mcft.h"
#include "../damagemodels/rotatingcrack.h"

using namespace Amie ;

ConcreteBehaviour::ConcreteBehaviour(double E, double nu, double compressive, planeType pt, RedistributionType rtype, SpaceDimensionality dim, MirrorState mirroring , double dx ,double  dy, double dz) : WeibullDistributedStiffness(E,nu, dim, compressive,0, pt, 0.05, 0.064, mirroring, dx , dy , dz ), rtype(rtype)
{

}

Form * ConcreteBehaviour::getCopy() const 
{
// 	return new Stiffness(param) ;
        std::default_random_engine generator(std::rand());
        std::weibull_distribution< double > distribution(5, 1);
        double weib = distribution(generator) ;
	double factor = 1. - variability + variability*weib ;
//	weib = distribution(generator) ;
// 	double upFactor = factor ; //1 -.7+.7*weib ; 
	NonLocalMCFT * fcrit = new NonLocalMCFT(down*factor,E*factor, materialRadius,rtype, mirroring , dx, dy, dz) ;
	fcrit->rebarLocationsAndDiameters = rebarLocationsAndDiameters ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, fcrit, /*new IsotropicLinearDamage()*/new RotatingCrack(E*factor, nu)) ;
// 	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMCFT(up, down,E, materialRadius, mirroring , dx, dy, dz), new NonLocalIsotropicLinearDamage()) ;
// 	copy->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
// 	ret->getDamageModel()->setThresholdDamageDensity(1);
	
	return copy ;
}



