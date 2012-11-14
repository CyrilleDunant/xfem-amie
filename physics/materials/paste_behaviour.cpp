// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "paste_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../stiffness.h"
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

ViscoElasticOnlyPasteBehaviour::ViscoElasticOnlyPasteBehaviour(SpaceDimensionality dim) : PasteBehaviour(12e9, 0.3, 0, dim)
{
	Kinf = 1.6752e6 ;
	Ginf = 1.1173e6 ;
	
	addBranch(3.9268e9, 2.6365e9, 0.003) ;
	addBranch(3.8072e9, 2.5542e9, 0.03) ;
	addBranch(206.88e6, 119.67e6, 0.33) ;
	addBranch(249.50e6, 191.00e6, 3) ;
	addBranch(289.68e6, 199.45e6, 27) ;
	addBranch(336.62e6, 231.77e6, 285) ;
	addBranch(238.57e6, 135.30e6, 3000) ;
}

void ViscoElasticOnlyPasteBehaviour::addBranch(double k, double g, double eta) 
{
	branches.push_back(std::make_pair(std::make_pair(k,g), eta)) ;
}

void ViscoElasticOnlyPasteBehaviour::clearBranches() 
{
	branches.clear() ;
}

Form * ViscoElasticOnlyPasteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1. - variability + variability*weib ;
	
	Matrix Kvisc = Material::cauchyGreen(std::make_pair(Kinf, Ginf), false, param.numRows() == 3 ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) * factor;
	
	std::vector<Matrix> rig ;
	std::vector<double> eta ;
	
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		  rig.push_back( Material::cauchyGreen(std::make_pair(branches[i].first.first, branches[i].first.second), false, param.numRows() == 3 ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL) * factor )  ;
		  eta.push_back( branches[i].second * factor ) ;
	}
	
	return new GeneralizedIterativeMaxwell(Kvisc, rig, eta) ;
}

