// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "gel_behaviour.h"
#include "../viscoelasticity_and_imposed_deformation.h"

using namespace Mu ;

GelBehaviour::GelBehaviour(double E, double nu, double alpha, SpaceDimensionality dim) : StiffnessWithImposedDeformation(E,nu, alpha, dim)
{

}

ViscoElasticOnlyGelBehaviour::ViscoElasticOnlyGelBehaviour(double E, double nu, double alpha, SpaceDimensionality dim): GelBehaviour(E,nu,alpha,dim)
{

}

Form * ViscoElasticOnlyGelBehaviour::getCopy() const
{
	Matrix c = param ;
	Vector a = imposed ;
	return new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, c, a ) ;
}