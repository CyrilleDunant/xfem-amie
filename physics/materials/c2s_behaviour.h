// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 2017
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef C2S_BEHAVIOUR_H
#define C2S_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct C2SBehaviour : public Stiffness
{
    // http://www.sciencedirect.com/science/article/pii/S0008884600005056
    C2SBehaviour(double E=130e9, double nu=0.3, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // C3S_BEHAVIOUR_H
