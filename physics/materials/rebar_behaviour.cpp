// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "rebar_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Mu ;

RebarBehaviour::RebarBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu, dim, -8*tensile,tensile)
{
	materialRadius = 0.0002 ;
}

Form * RebarBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new MohrCoulomb(up*factor,down*factor)) ;
	return ret ;
}



