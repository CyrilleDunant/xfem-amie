// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "aggregate_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../viscoelasticity.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../damagemodels/rotatingcrack.h"
#include "../damagemodels/fiberbasedisotropiclineardamage.h"
#include "../../utilities/random.h"

using namespace Mu ;

AggregateBehaviour::AggregateBehaviour(double E, double nu, double up, double yield, double c, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, 0.,0.), up(up), yield(yield), c(c)
{
	materialRadius = 0.0001 ;
}

Form * AggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,-8.*up*factor, E*factor), new FiberBasedIsotropicLinearDamage(0.1,0.6)) ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);
// 	ret->dfunc->setThresholdDamageDensity(1.);
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ;
}

ElasticOnlyAggregateBehaviour::ElasticOnlyAggregateBehaviour(double E, double nu, SpaceDimensionality dim) : AggregateBehaviour(E,nu,0.,dim)
{

}

Form * ElasticOnlyAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	Stiffness * copy =  new Stiffness(param*factor) ;
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ;
}


ViscoElasticOnlyAggregateBehaviour::ViscoElasticOnlyAggregateBehaviour(double E, double nu, SpaceDimensionality dim) : AggregateBehaviour(E,nu,0.,dim)
{

}

Form * ViscoElasticOnlyAggregateBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	return new Viscoelasticity( PURE_ELASTICITY, param*factor, 2 )  ;
}
