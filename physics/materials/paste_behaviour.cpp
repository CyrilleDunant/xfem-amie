// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "paste_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

PasteBehaviour::PasteBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, -8*tensile,tensile)
{
	materialRadius = 0.001 ;
	neighbourhoodRadius = 0.002 ;
}

Form * PasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,down*factor), materialRadius) ;
	ret->setNeighbourhoodRadius(neighbourhoodRadius) ;
	ret->criterion->setNeighbourhoodRadius(neighbourhoodRadius);
	ret->criterion->setMaterialCharacteristicRadius(materialRadius);
	ret->dfunc->setThresholdDamageDensity(.9);
	return ret ;
}

ElasticOnlyPasteBehaviour::ElasticOnlyPasteBehaviour(double E, double nu, SpaceDimensionality dim) : PasteBehaviour(E,nu,0.,dim)
{

}

Form * ElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	return new Stiffness(param*factor) ;
}



