// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 2017
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GYP_BEHAVIOUR_H
#define GYP_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct GypBehaviour : public Stiffness
{
    // from https://www.researchgate.net/profile/EJ_Garboczi/publication/303751822_Hydrate_dissolution_influence_on_the_Young's_modulus_of_cement_pastes/links/575acdde08ae9a9c95518e1a/Hydrate-dissolution-influence-on-the-Youngs-modulus-of-cement-pastes.pdf
    GypBehaviour(double E=30e9, double nu=0.3, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // C4AF_BEHAVIOUR_H