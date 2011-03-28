// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_base.h"
#include "level_set.h"
#include <limits>

#include "../mesher/delaunay.h"
#include "../polynomial/vm_function_base.h"

using namespace Mu ;


LevelSet::LevelSet()
{
	gType = LEVEL_SET ;
	this->center = Point(0,0) ;
}

LevelSet::LevelSet(Point p)
{
	gType = LEVEL_SET ;
	this->center = p ;
}

LevelSet::LevelSet(double x, double y)
{
	gType = LEVEL_SET ;
	this->center = Point(x,y) ;
}

LevelSet::LevelSet(Geometry * g)
{
	gType = LEVEL_SET ;
	this->source = g ;
	this->center = g->getCenter() ;
//	this->signedDistance = new Function(g->getCenter()) ;
}

LevelSet::LevelSet(Function f)
{
	gType = LEVEL_SET ;
	distanceFunction.push_back(f) ;
}

LevelSet::LevelSet(std::vector<Function> f)
{
	gType = LEVEL_SET ;
	for(size_t i = 0 ; i < f.size() ; i++)
		distanceFunction.push_back(f[i]) ;
}

std::vector<Point> LevelSet::getBoundingBox() const 
{
	std::vector<Point> ret ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	return ret ;
}
