//
// C++ Interface: damagemodel
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CONTACT_MODEL_H
#define CONTACT_MODEL_H

#include "../../elements/integrable_entity.h"
#include "../damagemodels/damagemodel.h"
#include "../../utilities/matrixops.h"

namespace Amie
{


class ContactModel : public DamageModel
{

public:

    ContactModel();
    virtual ~ContactModel() {};

    virtual void step(ElementState & s, double maxscore)  ;


} ;


}

#endif
