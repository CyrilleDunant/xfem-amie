// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "paste_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../generalized_fd_maxwell.h"
#include "../maxwell.h"
#include "../parallel_behaviour.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

PasteBehaviour::PasteBehaviour(double E, double nu, double up, double yield, double c, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, 0.,0.), up(up), yield(yield), c(c)
{
	materialRadius = 0.00025 ;
}

Form * PasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalInverseRootMohrCoulomb(up*factor,yield*factor, E, c)) ;
	ret->criterion->setMaterialCharacteristicRadius(materialRadius);
// 	ret->dfunc->setThresholdDamageDensity(.999);
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

ViscoElasticOnlyPasteBehaviour::ViscoElasticOnlyPasteBehaviour(double E, double nu, double Evisc, double nuvisc, double tau, SpaceDimensionality dim) : PasteBehaviour(E, nu, 0., dim), Evisc(Evisc), nuvisc(nuvisc), tau(tau) 
{

}

Form * ViscoElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. ;//- variability + variability*weib ;
	Matrix Kvisc = Material::cauchyGreen(std::make_pair(Evisc, nuvisc), true, param.numRows() == 3 ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) ;
	Vector decay(1./(tau*24*60*60), param.numRows()) ;
	
	Stiffness * stiff = new Stiffness(param*factor) ;
	NewmarkNumeroffMaxwell * max = new NewmarkNumeroffMaxwell(Kvisc, decay) ;
	
	return new ParallelBehaviour(stiff, max) ;
}

