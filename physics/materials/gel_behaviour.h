// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GEL_BEHAVIOUR_H
#define GEL_BEHAVIOUR_H

#include "../stiffness_with_imposed_deformation.h"
#include "../../geometry/geometry_base.h"

namespace Mu
{
	struct GelBehaviour : public StiffnessWithImposedDeformation
	{
		GelBehaviour(double E=22e9, double nu=0.28, double alpha=0.22, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
	} ;
	
	struct ViscoElasticOnlyGelBehaviour : public GelBehaviour
	{
		ViscoElasticOnlyGelBehaviour(double E=22e9, double nu=0.28, double alpha=0.22, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;
		
		Form * getCopy() const ;
	} ;
	
} ;

#endif // GEL_BEHAVIOUR
