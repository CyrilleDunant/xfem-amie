// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "ch_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Amie ;

CHBehaviour::CHBehaviour(double E, double nu, SpaceDimensionality dim) : Stiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim))
{
}

Form * CHBehaviour::getCopy() const 
{
	return new Stiffness(param) ;
}



