// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009 (added: ellipses)
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_base.h"
#include "level_set.h"
#include <limits>
#include <iomanip>

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

std::vector<Point> LevelSet::getBoundingBox() const 
{
	std::vector<Point> ret ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	ret.push_back(getCenter()) ;
	return ret ;
}