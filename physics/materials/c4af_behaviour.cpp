// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 017
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "c4af_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"

using namespace Amie ;

C4AFBehaviour::C4AFBehaviour(double E, double nu, SpaceDimensionality dim) : Stiffness(E, nu ,dim, PLANE_STRESS, YOUNG_POISSON)
{
}

Form * C4AFBehaviour::getCopy() const 
{
	return new DerivedStiffness(param) ;
}