//
// C++ Interface: crackinitiation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUCRACKINITIATION_H
#define MUCRACKINITIATION_H

#include "enrichmentbehaviour.h"

namespace Mu {

/** \brief Initiate Enrichment cracks from a fatigue Law
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class CrackInitiation : public EnrichmentBehaviour
{
public:
    CrackInitiation();

    ~CrackInitiation();

/** \brief return the new cracks to add*/
	virtual std::vector<EnrichmentFeature *> step(double t, DelaunayTree * dt) const;
};

}

#endif
