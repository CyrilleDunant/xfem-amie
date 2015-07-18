// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GEL_BEHAVIOUR_H
#define GEL_BEHAVIOUR_H

#include "../stiffness_with_imposed_deformation.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct GelBehaviour : public StiffnessWithImposedDeformation
{
    //volumic expansion is 22%, so linear is 7%,            cube root 0.028
    // variation in density 2.06 to 2.22 of amorphous silica suggests 0.025257
    // bulteel suggests 0.5 expansion,                             so 0.0367
    GelBehaviour(double E=22e9, double nu=0.18, double alpha=0.025, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

} ;

struct ViscoElasticOnlyGelBehaviour : public GelBehaviour
{
    int freeblocks ;

    ViscoElasticOnlyGelBehaviour(double E=22e9, double nu=0.18, double alpha=0.028, SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL) ;

    Form * getCopy() const ;
} ;

}

#endif // GEL_BEHAVIOUR
