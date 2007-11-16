//
// C++ Interface: crackinitiation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUCRACKINITIATION_H
#define MUCRACKINITIATION_H

#include "enrichmentbehaviour.h"

namespace Mu {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class CrackInitiation : public EnrichmentBehaviour
{
public:
    CrackInitiation();

    ~CrackInitiation();

	virtual std::vector<EnrichmentFeature *> step(double dt, DelaunayTree * dt) const;
};

}

#endif
