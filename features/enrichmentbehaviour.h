//
// C++ Interface: enrichmentbehaviour
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUENRICHMENTBEHAVIOUR_H
#define MUENRICHMENTBEHAVIOUR_H

#include "features.h"

namespace Mu {

/** Enrichment-generating behaviour interface.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class EnrichmentBehaviour{
public:
    EnrichmentBehaviour();
    ~EnrichmentBehaviour();

	virtual std::vector<EnrichmentFeature *> step(double dt, DelaunayTree * dt) const = 0;

};

}

#endif
