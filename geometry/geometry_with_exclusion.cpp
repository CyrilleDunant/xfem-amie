
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_with_exclusion.h"

using namespace Mu ;

GeometryWithExclusion::GeometryWithExclusion(Geometry * f, Geometry * e)
{
	father = f ;
	exclusions.push_back(e) ;
}

GeometryWithExclusion::GeometryWithExclusion::GeometryWithExclusion(Geometry * f, std::vector<Geometry *> e)
{
	father = f ;
	for(size_t i = 0 ; i < e.size() ; i++)
		exclusions.push_back(e[i]) ;
}
	
void GeometryWithExclusion::project(Point * p) const 
{
	Point base(p->getX(), p->getY(), p->getZ()) ;
	Point projOnFather(base) ;
	father->project(&projOnFather) ;
	Point proj(projOnFather) ;
	double distance = dist(base, proj) ;
	for(size_t i = 0 ; i < exclusions.size() ; i++)
	{
		Point projOnExclusion(base) ;
		exclusions[i]->project(&projOnExclusion) ;
		double distanceToExclusion = dist(base,projOnExclusion) ;
		if(distanceToExclusion < distance)
		{
			proj = projOnExclusion ;
			distance = distanceToExclusion ;
		}
	}
}

bool GeometryWithExclusion::in(const Point & p) const  
{
	if(!father->in(p))
		return false ;
		
	for(size_t i = 0 ; i < exclusions.size() ; i++)
	{
		if(exclusions[i]->in(p))
			return false ;
	}
	
	return true ;
}

double GeometryWithExclusion::area() const
{
	double a = father->area() ;
	for(size_t i = 0 ; i < exclusions.size() ; i++)
		a -= exclusions[i]->area() ;
	return a ;
}

double GeometryWithExclusion::volume() const 
{
	double v = father->volume() ;
	for(size_t i = 0 ; i < exclusions.size() ; i++)
		v -= exclusions[i]->volume() ;
	return v ;

}

bool GeometryWithExclusion::intersects(const Geometry * g) const 
{
	bool intersectsFather = father->intersects(g) ;
	if(!intersectsFather)
		return false ;
	for(size_t i = 0 ; i < exclusions.size() ; i++)
	{
		if(!exclusions[i]->intersects(g) && exclusions[i]->in(g->getCenter()))
			return false ;
	}
	
	return true ;
}

std::vector<Point> GeometryWithExclusion::intersection(const Geometry * g) const
{
	std::vector<Point> ret ;
	std::vector<Point> tmp = father->intersection(g) ;
	for(size_t i = 0 ; i < tmp.size() ; i++)
	{
		bool in = true ;
		size_t j = 0 ;
		while(in && j < exclusions.size())
		{
			if(exclusions[j]->in(tmp[i]))
				in = false ;
		}
		if(in)
			ret.push_back(tmp[i]) ;
	}
	for(size_t i = 0 ; i < exclusions.size() ; i++)
	{
		tmp = exclusions[i]->intersection(g) ;
		for(size_t i = 0 ; i < tmp.size() ; i++)
			ret.push_back(tmp[i]) ;		
	}
	
	return ret ;
}

