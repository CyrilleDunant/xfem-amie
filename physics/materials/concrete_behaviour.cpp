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

using namespace Mu ;

ConcreteBehaviour::ConcreteBehaviour(double E, double nu, double tensile, double compressive, planeType pt, bool reinforced, SpaceDimensionality dim, MirrorState mirroring , double dx ,double  dy, double dz) : WeibullDistributedStiffness(E,nu, dim, compressive,tensile, pt, mirroring, dx , dy , dz ), reinforced(reinforced)
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
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMCFT(up*upFactor, down*factor,E*factor, materialRadius,reinforced, mirroring , dx, dy, dz)/*, new AnisotropicLinearDamage()*/) ;
// 	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMCFT(up, down,E, materialRadius, mirroring , dx, dy, dz), new NonLocalIsotropicLinearDamage()) ;
	ret->getFractureCriterion()->setMaterialCharacteristicRadius(materialRadius);
	return ret ;
}



