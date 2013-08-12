// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "paste_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
#include "../viscoelasticity.h"
#include "../viscoelasticity_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../damagemodels/fiberbasedisotropiclineardamage.h"
#include "../../utilities/random.h"
#include "../fracturecriteria/maxstrain.h"
#include "../fracturecriteria/creeprupture.h"
#include "../damagemodels/spacetimefiberbasedisotropiclineardamage.h"

using namespace Mu ;

PasteBehaviour::PasteBehaviour(double E, double nu, double up, double yield, double c, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, 0.,0.), up(up), yield(yield), c(c)
{
	materialRadius = 0.00018 ;
}

Form * PasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,-8.*up*factor, E*factor), new FiberBasedIsotropicLinearDamage(0.1,0.6)) ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);
 	copy->dfunc->setThresholdDamageDensity(.6);
	
/*	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}*/
	return copy ;
}

ElasticOnlyPasteBehaviour::ElasticOnlyPasteBehaviour(double E, double nu, SpaceDimensionality dim) : PasteBehaviour(E,nu,0.,0.,0.,dim)
{

}

Form * ElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	return new Stiffness(param*factor) ;
}

ViscoElasticOnlyPasteBehaviour::ViscoElasticOnlyPasteBehaviour(double E, double nu, double e1, double e2, SpaceDimensionality dim) : PasteBehaviour(E, nu, 0.,0.,0., dim), e_1(e1), e_2(e2)
{

}

Form * ViscoElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;

	Matrix C0 = param*factor ;
	Matrix C1 = C0*e_1 ;
	Matrix C2 = C0*e_2 ;
	Matrix E1 = C1*10. ;
	Matrix E2 = C2*300. ;
	
	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(C1,E1));
	branches.push_back(std::make_pair(C2,E2));
	
	return new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C0, branches) ;
}

ViscoDamagePasteBehaviour::ViscoDamagePasteBehaviour(double E, double nu, double e1, double e2 , double up_, double r, SpaceDimensionality dim) : PasteBehaviour(E, nu, up_,0.,0., dim), e_1(e1), e_2(e2)
{
	materialRadius = r ;
}

Form * ViscoDamagePasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;

	Matrix C0 = param*factor ;
	Matrix C1 = C0*e_1 ;
	Matrix C2 = C0*e_2 ;
	Matrix E1 = C1*10. ;
	Matrix E2 = C2*300. ;
	
	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(C1,E1));
	branches.push_back(std::make_pair(C2,E2));
	
	ViscoelasticityAndFracture * copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStrain(up, up*param[0][0]*factor), new SpaceTimeFiberBasedIsotropicLinearDamage(0.1,1e-9)) ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius) ;
	return copy ;

}
