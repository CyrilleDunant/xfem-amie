// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "caco3_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"



using namespace Amie ;

CaCO3Behaviour::CaCO3Behaviour(SpaceDimensionality dim, double E, double nu) : Stiffness(E,nu, dim, PLANE_STRESS, YOUNG_POISSON)
{
}



Form * CaCO3Behaviour::getCopy() const 
{
	return new DerivedStiffness(param) ;
}



