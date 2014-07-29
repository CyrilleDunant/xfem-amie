//
// C++ Interface: enrichmentbehaviour
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUENRICHMENTBEHAVIOUR_H
#define MUENRICHMENTBEHAVIOUR_H

#include "features.h"

namespace Amie {

/** \brief Enrichment-generating behaviour interface.
* 
* Such a behaviour adds enrichments according to some law.
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class EnrichmentBehaviour{
public:
    EnrichmentBehaviour();
    ~EnrichmentBehaviour();

	virtual std::vector<EnrichmentFeature *> step(double t, DelaunayTree * dt) const = 0;

};

}

#endif
