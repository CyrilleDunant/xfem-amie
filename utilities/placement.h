
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
namespace Mu
{
	std::vector<Feature *> placement(const Mu::Geometry* box, std::vector< Mu::Feature* > inclusions, int* nombreGranulatsPlaces, int nombreGranulatsDejaPlaces = 0, int triesMax = 6400, bool verbose = true) ;
	std::vector<Mu::EllipsoidalInclusion *> placement_with_rotation(const Mu::Geometry* box, std::vector< Mu::EllipsoidalInclusion* > inclusions, int* nombreGranulatsPlaces, int triesMax, bool verbose = true) ;
// 	double masseInitiale;
// 	double densite;

	double chiffreAleatoire(double );
} ;

#endif
