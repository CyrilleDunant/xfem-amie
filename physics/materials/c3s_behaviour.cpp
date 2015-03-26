// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "c3s_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Amie ;

C3SBehaviour::C3SBehaviour(double E, double nu, SpaceDimensionality dim) : Stiffness(Tensor::cauchyGreen(std::make_pair(E,nu), true,dim))
{
}

Form * C3SBehaviour::getCopy() const 
{
	return new DerivedStiffness(param) ;
}



