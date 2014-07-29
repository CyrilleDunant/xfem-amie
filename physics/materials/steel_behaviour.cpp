// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "steel_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../../utilities/random.h"

using namespace Amie ;

SteelBehaviour::SteelBehaviour(double E, double nu, double tensile, SpaceDimensionality dim) : WeibullDistributedStiffness(E,nu,dim, -8*tensile,tensile)
{
	materialRadius = 0.0002 ;
}

Form * SteelBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * copy = new StiffnessAndFracture(param*factor, new MohrCoulomb(up*factor,down*factor)) ;
	
	if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ;
}



