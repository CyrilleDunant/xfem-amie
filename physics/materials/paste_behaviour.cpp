// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
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
#include <cstdio> 

using namespace Mu ;

PasteBehaviour::PasteBehaviour(double E, double nu, double up, double yield, double c, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, 0.,0.), up(up), yield(yield), c(c)
{
	materialRadius = 0.0005 ;
}

Form * PasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	return new Stiffness(param*factor) ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new NonLocalMohrCoulomb(up*factor,-8.*up*factor, E*factor), new FiberBasedIsotropicLinearDamage(0.05,0.9)) ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);
//  	copy->dfunc->setThresholdDamageDensity(.8);
	
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

ViscoElasticOnlyPasteBehaviour::ViscoElasticOnlyPasteBehaviour(double E, double nu, double e1, double e2, SpaceDimensionality dim) : PasteBehaviour(E, nu, 0.,0.,0., dim), e_1(e1), e_2(e2), freeblocks(0)
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
	
	return new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C0, branches,0,freeblocks) ;
}

ViscoDamagePasteBehaviour::ViscoDamagePasteBehaviour(double E, double nu, double e1, double e2 , double up_, double r, SpaceDimensionality dim) : PasteBehaviour(E, nu, up_,0.,0., dim), e_1(e1), e_2(e2), freeblocks(0), ctype(STRAIN_CRITERION)
{
	materialRadius = r ;
	stressFraction = 1. ;
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
	
	ViscoelasticityAndFracture * copy ;
	IsotropicLinearDamage * dampaste = new IsotropicLinearDamage() ;
//	dampaste->setLogitViscousDamageLaw(0.025, 0.3, 2.5) ;

	switch(ctype)
	{
		case STRAIN_CRITERION:
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStrain(up, up*param[0][0]*factor), dampaste, 0, freeblocks) ;
			break ;
		case STRESS_CRITERION:
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStress(up, up*param[0][0]*factor), dampaste, 0, freeblocks) ;
			break ;
		case MIXED_CRITERION:
			double k = 1. + C0[0][0]/C1[0][0] + C0[0][0]/C2[0][0] ;
//			std::cout << k << std::endl ;
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalEllipsoidalMixedCriterion(up, up*param[0][0]*stressFraction, param[0][0], param[0][0]/k), dampaste, 0, freeblocks) ;
			break ;
	}
	copy->criterion->setMaterialCharacteristicRadius(materialRadius) ;
	copy->getFractureCriterion()->setScoreTolerance(1e-4) ;
	copy->getDamageModel()->setDamageDensityTolerance(0.05) ;
	copy->getDamageModel()->setThresholdDamageDensity(0.8) ;
	copy->getDamageModel()->setNeedGlobalMaximumScore(true) ;
	return copy ;

}

PseudoBurgerViscoElasticOnlyPasteBehaviour::PseudoBurgerViscoElasticOnlyPasteBehaviour(double E, double nu, double e1, double t2, SpaceDimensionality dim) : PasteBehaviour(E, nu, 0.,0.,0., dim), e_1(e1), t_2(t2)
{
	variability = 0. ;//0.001 ;
}

Form * PseudoBurgerViscoElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;

	Matrix C0 = param*factor ;
	Matrix C1 = C0*e_1 ;
	Matrix C2 = C0*0.1 ;
	Matrix E1 = C1*2. ;
	Matrix E2 = C2*t_2 ;
	
	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(C1,E1));
	branches.push_back(std::make_pair(C2,E2));
	
	return new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C0, branches) ;
}

PseudoBurgerViscoDamagePasteBehaviour::PseudoBurgerViscoDamagePasteBehaviour(double E, double nu, double e1, double t2, double up, double r, SpaceDimensionality dim) : PasteBehaviour(E, nu, up,0.,0., dim), e_1(e1), t_2(t2), freeblocks(0), ctype(STRAIN_CRITERION)
{
	stressFraction = 0.85 ;
	materialRadius = r ;
	variability = 0. ;//0.001 ;
}

Form * PseudoBurgerViscoDamagePasteBehaviour::getCopy() const 
{
	
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
//	std::cout << factor << std::endl ;

	Matrix C0 = param*factor ;
	Matrix C1 = C0*e_1 ;
	Matrix C2 = C0*0.1 ;
	Matrix E1 = C1*2. ;
	Matrix E2 = C2*t_2 ;
// 	C0[2][2] *= .9 ;
// 	C1[2][2] *= .9 ;
// 	C2[2][2] *= .9 ;
// 	E1[2][2] *= .9 ;
// 	E2[2][2] *= .9 ;
	
	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(C1,E1));
	branches.push_back(std::make_pair(C2,E2));

	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.05, 1e-9, 0.8 ) ;
	dampaste->setLogitViscousDamageLaw(0.025, 0.3, 2.5) ;

	ViscoelasticityAndFracture * copy ;
	switch(ctype)
	{
		case STRAIN_CRITERION:
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStrain(up, up*param[0][0]*factor), dampaste, 0, freeblocks) ;
			break ;
		case STRESS_CRITERION:
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalMaximumStress(up, up*param[0][0]*factor), dampaste, 0, freeblocks) ;
			break ;
		case MIXED_CRITERION:
			double k = 1. + C0[0][0]/C1[0][0] + C0[0][0]/C2[0][0] ;
//			std::cout << k << std::endl ;
			copy = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, C0, branches, new SpaceTimeNonLocalEllipsoidalMixedCriterion(up, up*param[0][0]*stressFraction, param[0][0], param[0][0]/k), dampaste, 0, freeblocks) ;
			break ;
	}
	copy->criterion->setMaterialCharacteristicRadius(materialRadius) ;
	return copy ;

}


