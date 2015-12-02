// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "gel_behaviour.h"
#include "../viscoelasticity_and_imposed_deformation.h"

using namespace Amie ;

GelBehaviour::GelBehaviour(bool st, double E, double nu, double alpha, SpaceDimensionality dim, planeType pt, int b) : StiffnessWithImposedDeformation(E,nu, alpha, dim, pt), spaceTime(st), freeblocks(b)
{

}

Form * GelBehaviour::getCopy() const
{
	if(!spaceTime)
		return new  StiffnessWithImposedDeformation(param, imposed);

	Matrix c = param ;
	Vector a = imposed ;
	return new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, c, a, 2 + freeblocks ) ;
}
