// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
// Author: Zhangli Hu <zhangli.hu@epfl.ch>, (C) 2017
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FILLER_BEHAVIOUR_H
#define FILLER_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct FillerBehaviour : public Stiffness
{
    
    FillerBehaviour(double E=138e9, double nu=0.28, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // Filler_BEHAVIOUR_H