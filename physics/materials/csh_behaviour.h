// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CSH_BEHAVIOUR_H
#define CSH_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
typedef enum
{
    INNER_CSH,
    OUTER_CSH
} CSHType ;

struct CSHBehaviour : public Stiffness
{
    // per http://www.sciencedirect.com/science/article/pii/S1359645408008720
    CSHBehaviour(CSHType type, double E=25e9, double nu=0.25,  SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

} ;

#endif // PASTE_BEHAVIOUR
