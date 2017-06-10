// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 2017
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef C4AF_BEHAVIOUR_H
#define C4AF_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct C4AFBehaviour : public Stiffness
{
    // per http://www.sciencedirect.com/science/article/pii/S1359645408008720
    C4AFBehaviour(double E=125e9, double nu=0.3, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // C4AF_BEHAVIOUR_H