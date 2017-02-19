// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "steel_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"

using namespace Amie ;

SteelBehaviour::SteelBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu,dim, -8*tensile,tensile, PLANE_STRESS, 0.2, 0.0002, YOUNG_POISSON)
{
	materialRadius = 0.0002 ;
}

Form * SteelBehaviour::getCopy() const 
{
        std::default_random_engine generator(std::rand());
        std::weibull_distribution< double > distribution(5, 1);
        double weib = distribution(generator) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new MohrCoulomb(up*factor,down*factor, E)) ;
	
	return copy ;
}



