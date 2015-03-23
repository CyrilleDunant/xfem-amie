// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef C3S_BEHAVIOUR_H
#define C3S_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct C3SBehaviour : public Stiffness
{
    // per http://www.sciencedirect.com/science/article/pii/S1359645408008720
    C3SBehaviour(double E=138e9, double nu=0.28, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // C3S_BEHAVIOUR_H
