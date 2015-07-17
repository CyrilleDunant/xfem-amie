
// Author: Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __PLACEMENT_H__
#define __PLACEMENT_H__

#include <iostream>
#include <vector>

#include "../features/features.h"
#include "../features/inclusion.h"
namespace Amie
{

std::vector<Feature *> placement2D(const Geometry* box, std::vector<Feature *> inclusions, double minDist = 0., int placedAggregates = 0, int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;

std::vector<Feature *> placement3D(const Geometry* box, std::vector<Feature *> inclusions, double minDist = 0., int placedAggregates = 0, int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;

std::vector<Feature *> placement2DInInclusions(const Geometry* box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, double minDist = 0., int placedAggregates = 0, int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;

std::vector<Feature *> placement2DOnEdge( const Geometry * box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, bool onVertex = false, double minDist = 0., int placedAggregates = 0., int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;

}

#endif
