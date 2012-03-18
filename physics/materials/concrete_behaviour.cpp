// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "concrete_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../fracturecriteria/mcft.h"
#include "../../utilities/random.h"
#include "../damagemodels/rotatingcrack.h"

using namespace Mu ;

ConcreteBehaviour::ConcreteBehaviour(double E, double nu, double compressive, planeType pt, RedistributionType rtype, SpaceDimensionality dim, MirrorState mirroring , double dx ,double  dy, double dz) : WeibullDistributedStiffness(E,nu, dim, compressive,0, pt, mirroring, dx , dy , dz ), rtype(rtype)
{
	materialRadius = 0.015 ;
	neighbourhoodRadius = materialRadius*1.5 ;
	variability = 0. ;
}

Form * ConcreteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	weib = RandomNumber().weibull(1,5) ;
	double upFactor = factor ; //1 -.7+.7*weib ; 
	NonLocalMCFT * fcrit = new NonLocalMCFT(down*factor,E*factor, materialRadius,rtype, mirroring , dx, dy, dz) ;
	fcrit->rebarLocationsAndDiameters = rebarLocationsAndDiameters ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, fcrit, /*new IsotropicLinearDamage()*/new RotatingCrack(E, nu)) ;
// 	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMCFT(up, down,E, materialRadius, mirroring , dx, dy, dz), new NonLocalIsotropicLinearDamage()) ;
	ret->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
	return ret ;
}



