// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CaCO3_BEHAVIOUR_H
#define CaCO3_BEHAVIOUR_H

#include "../stiffness.h"
#include "../../geometry/geometry_base.h"

namespace Amie
{
struct CaCO3Behaviour : public Stiffness
{
    // per http://ciks.cbt.nist.gov/~garbocz/paper148/node10.htm#
    CaCO3Behaviour(SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL, double E=79.6e9, double nu=0.31) ;

    virtual Form * getCopy() const ;

} ;

}

#endif // CaCO3_BEHAVIOUR_H
