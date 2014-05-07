//
// C++ Implementation: delaunay
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parallel_delaunay_3d.h"
#include <omp.h>
#include <limits>

// #define DEBUG
// #undef DEBUG

using namespace Mu ;


ParallelDelaunayTree3D::ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2, Point *p3, const std::vector<Geometry *> & domains) : domains(domains)
{
	for(size_t i = 0 ; i < domains.size() ; i++)
	{
		meshes.push_back(new DelaunayTree3D(p0, p1, p2, p3));
	}
}

void ParallelDelaunayTree3D::insert(Point * p)
{
	
}
