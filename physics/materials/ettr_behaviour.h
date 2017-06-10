// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 2017
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ETTR_BEHAVIOUR_H
#define ETTR_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct EttrBehaviour : public Stiffness
{
    // from https://www.researchgate.net/profile/EJ_Garboczi/publication/303751822_Hydrate_dissolution_influence_on_the_Young's_modulus_of_cement_pastes/links/575acdde08ae9a9c95518e1a/Hydrate-dissolution-influence-on-the-Youngs-modulus-of-cement-pastes.pdf
    EttrBehaviour(double E=22.4e9, double nu=0.25, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // Ettr_BEHAVIOUR_H