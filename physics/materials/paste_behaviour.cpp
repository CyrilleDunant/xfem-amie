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
	time_d = true ;
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
// 	IsotropicLinearDamageRate * dampaste = new IsotropicLinearDamageRate() ;
	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.025, 1e-9, 0.8 ) ;
	dampaste->setLogitViscousDamageLaw(0.025, 0.3, 2.5) ;

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
	copy->getDamageModel()->setDamageDensityTolerance(0.01) ;
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

	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.05, 1e-9, 0.7 ) ;
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



HydratingMechanicalCementPaste::HydratingMechanicalCementPaste(FeatureTree * diffusionTree) :LinearForm(Material::cauchyGreen(std::make_pair(1,0.2), true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS), true, true, 2),  diffusionTree(diffusionTree)
{
	v.push_back(XI);
	v.push_back(ETA);
}

Form * HydratingMechanicalCementPaste::getCopy() const
{
	return new HydratingMechanicalCementPaste(diffusionTree) ;
}

void HydratingMechanicalCementPaste::step(double timestep, ElementState & currentState, double maxscore)
{
	DelaunayTriangle * diffusionElement = diffusionTree->getElements2D()[dynamic_cast<DelaunayTriangle *>(currentState.getParent())->index] ;
	Vector saturation = diffusionElement->getState().getDisplacements() ;
	double effectiveSaturation = (saturation[0]+saturation[1]+saturation[2])/3 ;
	double doh = dynamic_cast<HydratingDiffusionCementPaste *>(diffusionElement->getBehaviour())->getDegreeOfHydration() ;
	
	param = getMechanicalProperties(effectiveSaturation,doh) ;
}

Matrix HydratingMechanicalCementPaste::getMechanicalProperties(double effectiveSaturation, double doh) 
{
	double C [10]= { 6328.4269, - 5426.419, 20602.939, - 13932.055, 19700.509, 105229.09, 25708.722, - 52055.068, 47843.694, - 136323.96 };

	//in MPa
	double E = C[0] + C[1]*effectiveSaturation + C[2]*doh + C[3]*effectiveSaturation*effectiveSaturation + C[4]*doh*effectiveSaturation + C[5]*doh*doh + C[6]*effectiveSaturation*effectiveSaturation*effectiveSaturation + C[7]*doh*effectiveSaturation*effectiveSaturation + C[8]*doh*doh*effectiveSaturation + C[9]*doh*doh*doh;

	if (E<0){ E=0; }
	E *= 1e6 ;
		
	double nu = 0.2 ;
	
	return Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS) ;
}

Vector HydratingMechanicalCementPaste::getAutogeneousForce(double saturation, double doh) 
{

}


void HydratingMechanicalCementPaste::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

bool HydratingMechanicalCementPaste::fractured() const {return false ; }

HydratingMechanicalCementPaste::~HydratingMechanicalCementPaste() {} ;

bool HydratingMechanicalCementPaste::changed() const {return false ; }

Matrix HydratingMechanicalCementPaste::getTensor(const Point & p, IntegrableEntity * e, int g) const 
{
	return param ;
}


HydratingDiffusionCementPaste::HydratingDiffusionCementPaste(): LinearForm(Matrix(2,2), true, false, 1)  
{
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(TIME_VARIABLE);
	doh = 0.3 ;
	change = false ;
	double D1_bz= (-5.412*doh*doh - 0.987*doh + 3.713)*1e-10;     //[m2/s], Duffison coefficient in function of DoH

	double a = 0.04 ;
	int n = 6;
	double Abz=0.9002*doh-0.1503;    //Bolzmann function parameter in function of DoH
	double Bbz=0.2815*doh+0.0297;   
	double h_b=1.-1./(1.+exp((1.-Abz)/Bbz));   //Bolzmann function for inversed desorption isotherm
	
	param[0][0]= D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
	param[1][1]= D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
}

Form * HydratingDiffusionCementPaste::getCopy() const
{
	HydratingDiffusionCementPaste * ret = new HydratingDiffusionCementPaste() ;
	ret->doh = doh ;
	return ret ;
}

void HydratingDiffusionCementPaste::step(double timestep, ElementState & currentState, double maxscore) 
{
	if(saturation.size() == 0)
		saturation.resize(6, 1.);
	else
		saturation = currentState.getDisplacements() ;
	double effectiveSaturation = (saturation[0]+saturation[1]+saturation[2])/3. ;
	effectiveSaturation = std::min(.999, effectiveSaturation) ;
	double deltaH = std::max(0.,getDeltaDoH(effectiveSaturation, currentState)) ;
	doh = std::min(doh+deltaH, 1.) ;
	double D = std::abs(getDiffusionCoefficient(effectiveSaturation, currentState) ) ;
	param[0][0]=D ;
	param[1][1]=D ;
	currentState.getParent()->behaviourUpdated = true ;
}

double HydratingDiffusionCementPaste::getDeltaDoH(double saturation, ElementState & currentState) 
{
	double thresholdRelativeHumidity = 1.199*doh-0.1503 ;
	double maxRelativeHumidity = 0.999 ;
	double factor = (saturation-thresholdRelativeHumidity)/(maxRelativeHumidity-thresholdRelativeHumidity) ;
	
	if(factor < POINT_TOLERANCE_2D)
		return 0. ;

	DelaunayTriangle * diffusionElement = dynamic_cast<DelaunayTriangle *>(currentState.getParent()) ;
	
	double deltaTime = currentState.getNodalDeltaTime() ;
// 	double currentTime = diffusionElement->getBoundingPoint(0).t + deltaTime*.5;
	
	double A1 =-0.347;
	double A2 = 278.8;
	double A3 =-0.677;
	double A4 = 1.37;

	return ((-(A1/A2)*exp(-doh/A2) - (A3/A4)*exp(-doh/A4) )* deltaTime)*factor;
	
}

double HydratingDiffusionCementPaste::getDiffusionCoefficient(double saturation, ElementState & currentState) 
{
	double D1_bz= (-5.412*doh*doh - 0.987*doh + 3.713)*1e-10;     //[m2/s], Duffison coefficient in function of DoH

	double a = 0.04 ;
	int n = 6;
	double Abz=0.9002*doh-0.1503;    //Bolzmann function parameter in function of DoH
	double Bbz=0.2815*doh+0.0297;   
	double h_b=1.-1./(1.+exp((saturation-Abz)/Bbz));   //Bolzmann function for inversed desorption isotherm
	return D1_bz*(a + (1.-a)*1./(1.+pow((1.-h_b)/0.25,n)))*24*3600;
}


void HydratingDiffusionCementPaste::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{

	ret[0][0] = (vm->ieval(VectorGradientDot(p_i) * param * VectorGradient(p_j, true),  gp, Jinv, v)
		  + vm->ieval(VectorGradient(p_i) * param * VectorGradientDot(p_j, true),  gp, Jinv, v)) ;
}

void HydratingDiffusionCementPaste::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine *vm) const
{
	ret[0][0] = (vm->ieval(Differential(p_i, TIME_VARIABLE) * Differential(p_j, TIME_VARIABLE),  gp, Jinv, v)
		  + vm->ieval(DoubleDifferential(p_j, TIME_VARIABLE, TIME_VARIABLE) * p_i,  gp, Jinv, v)) ;
}

HydratingDiffusionCementPaste::~HydratingDiffusionCementPaste() { } ;

/** \brief return true if the damage state has been modfied*/
bool HydratingDiffusionCementPaste::changed() const { return change ; }

/** \brief Return the (damaged) Stifness tensor*/
Matrix HydratingDiffusionCementPaste::getTensor(const Point & p, IntegrableEntity * e , int g ) const { return param ; }



