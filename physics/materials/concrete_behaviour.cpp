// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "concrete_behaviour.h"
#include "../stiffness_and_fracture.h"
#include "../homogenization/homogenization_base.h"
#include "../fracturecriteria/mohrcoulomb.h"
#include "../fracturecriteria/mcft.h"
#include "../../utilities/random.h"

using namespace Mu ;

ConcreteBehaviour::ConcreteBehaviour(double E, double nu, double tensile, double compressive, SpaceDimensionality dim, MirrorState mirroring , double dx ,double  dy, double dz) : WeibullDistributedStiffness(E,nu, dim, compressive,tensile, mirroring, dx , dy , dz )
{
	materialRadius = 0.06 ;
	neighbourhoodRadius = 0.18 ;
}

Form * ConcreteBehaviour::getCopy() const 
{
	double weib = RandomNumber().weibull(1,5) ;
	double factor = 1 - variability + variability*weib ;
	StiffnessAndFracture * ret = new StiffnessAndFracture(param*factor, new NonLocalMCFT(up, down,E, materialRadius, mirroring , dx, dy, dz)) ;
	return ret ;
}



