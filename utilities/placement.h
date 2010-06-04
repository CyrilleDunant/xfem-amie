
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PLACEMENT_H__
#define __PLACEMENT_H__

#include <iostream>
#include <vector>

#include "../features/features.h"
#include "../features/inclusion.h"
namespace Mu
{
	std::vector<Feature *> placement(const Mu::Geometry* box, std::vector< Mu::Feature* > inclusions, int* nombreGranulatsPlaces, int triesMax, bool verbose = true) ;
	std::vector<Mu::EllipsoidalInclusion *> placement_with_rotation(const Mu::Geometry* box, std::vector< Mu::EllipsoidalInclusion* > inclusions, int* nombreGranulatsPlaces, int triesMax, bool verbose = true) ;
// 	double masseInitiale;
// 	double densite;


	double chiffreAleatoire(double );
} ;

#endif
