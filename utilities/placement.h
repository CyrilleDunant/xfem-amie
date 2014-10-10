
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
	std::vector<Feature *> placement(const Amie::Geometry* box, std::vector< Amie::Feature* > inclusions, int* nombreGranulatsPlaces, int nombreGranulatsDejaPlaces = 0, int triesMax = 6400, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>(), bool verbose = true) ;
	std::vector<Amie::EllipsoidalInclusion *> placement_with_rotation(const Amie::Geometry* box, std::vector< Amie::EllipsoidalInclusion* > inclusions, int* nombreGranulatsPlaces, int triesMax, bool verbose = true) ;
	
	std::vector<Feature *> placement2D(const Geometry* box, std::vector<Feature *> inclusions, double minDist = 0., int placedAggregates = 0, int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;
	std::vector<Feature *> placement2DInInclusions(const Geometry* box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, double minDist = 0., int placedAggregates = 0, int tries = 6400, double orientation = M_PI, std::vector<Geometry *> exclusionZones = std::vector<Geometry *>()) ;
	
// 	double masseInitiale;
// 	double densite;

	double chiffreAleatoire(double );
} ;

#endif
