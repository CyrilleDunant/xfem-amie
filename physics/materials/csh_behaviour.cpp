// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "csh_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Amie ;

CSHBehaviour::CSHBehaviour(CSHType type, double E, double nu, SpaceDimensionality dim) : Stiffness(Material::cauchyGreen(std::make_pair(E,nu), true,dim))
{
    if(type == INNER_CSH)
        param *= 1.8 ;
}

Form * CSHBehaviour::getCopy() const 
{
	return new Stiffness(Matrix(param.numRows(), param.numCols(), const_cast<Vector *>(&param.array()))) ;
}


