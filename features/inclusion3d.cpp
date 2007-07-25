//
// C++ Implementation: inclusion3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "inclusion3d.h"

using namespace Mu ;

Inclusion3D::Inclusion3D(Feature *father, double r, double x, double y, double z) : Sphere(r, x, y,z ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Sphere(r*1.07, x, y,z) ;
	this->boundary2 = new Sphere(r*0.93, x, y,z) ;
}
	
Inclusion3D::Inclusion3D(Feature *father, double r,  Point center) : Sphere(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Sphere(r*1.07, center) ;
	this->boundary2 = new Sphere(r*0.93, center) ;
}
	
Inclusion3D::Inclusion3D(double r, double x, double y, double z) : Sphere(r, x,y,z ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Sphere(r*1.02, x,y,z) ;
	this->boundary2 = new Sphere(r*0.93, x, y,z) ;
}
	
Inclusion3D::Inclusion3D(double r, Point center) : Sphere(r,center ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Sphere(r*1.02, center) ;
	this->boundary2 = new Sphere(r*0.93, center) ;
}
	
bool Inclusion3D::interacts(Feature * f) const 	{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
	if(f->inBoundary((*i)))
		return true ;
return false ;
}
	
std::vector<Geometry *> Inclusion3D::getRefinementZones(size_t level) const 
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Sphere(sqrt(radius)*(1.2), Sphere::getCenter())) ;
	if(level > 1)
		ret.push_back(new Sphere(sqrt(radius)*(1.15), Sphere::getCenter())) ;
	if(level > 2)
		ret.push_back(new Sphere(sqrt(radius)*(1.08), Sphere::getCenter())) ;
	return ret ;
}

std::vector<DelaunayTriangle *> Inclusion3D::getTriangles( DelaunayTree * dt)  { 
	return std::vector<DelaunayTriangle *> (0) ;
}
	
std::vector<DelaunayTetrahedron *> Inclusion3D::getTetrahedrons(const DelaunayTree_3D * dt)  {
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->conflicts(this->boundary) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren()->size() ;  j++)
		{
			if(this->getChild(j)->in(temp[i]->getCenter()))
			{
				inChild = true ; 
				break ;
			}
		}
		if(this->boundary->in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}
	
	
Point * Inclusion3D::pointAfter(size_t i) { 
	return NULL ; 
}
	

	
void Inclusion3D::sample(size_t n)
{
	this->Sphere::sampleSurface(n) ;
	
// 	size_t numberOfRadii = static_cast<size_t>(sqrt(std::max(4*n, (size_t)37))) ;
// 	double r = sqrt(radius) ;
// 	double radiusExt = (r + (r*.5/(double)numberOfRadii) ) ;
// 	
// 	double radiusInt = (r - (r*.5/(double)numberOfRadii) ) ;
// 	
// 	delete this->boundary ;
// 	delete this->boundary2 ;
// 	
// 	this->boundary = new Sphere(radiusExt, *getCenter()) ;
// 	this->boundary2 = new Sphere(radiusInt, *getCenter()) ;
}


VirtualInclusion3D::VirtualInclusion3D(Feature *father, double r, double x, double y, double z) : Sphere(r, x, y,z ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = NULL ; //new Sphere(r*1.07, x, y,z) ;
	this->boundary2 =  NULL ; //new Sphere(r*0.93, x, y,z) ;
}

VirtualInclusion3D::VirtualInclusion3D(Feature *father, double r,  Point center) : Sphere(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary =  NULL ; //new Sphere(r*1.07, center) ;
	this->boundary2 =  NULL ; //new Sphere(r*0.93, center) ;
}

VirtualInclusion3D::VirtualInclusion3D(double r, double x, double y, double z) : Sphere(r, x,y,z ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary =  NULL ; //new Sphere(r*1.02, x,y,z) ;
	this->boundary2 =  NULL ; //new Sphere(r*0.93, x, y,z) ;
}

VirtualInclusion3D::VirtualInclusion3D(double r, Point center) : Sphere(r,center ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Sphere(r*1.02, center) ;
	this->boundary2 = new Sphere(r*0.93, center) ;
}

bool VirtualInclusion3D::interacts(Feature * f) const 	{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;
}

bool VirtualInclusion3D::inBoundary(const Point *) const
{
	return false ;
}

bool VirtualInclusion3D::inBoundary(const Point) const
{
	return false ;
}

std::vector<Geometry *> VirtualInclusion3D::getRefinementZones(size_t level) const 
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Sphere(sqrt(radius)*(1.2), Sphere::getCenter())) ;
	if(level > 1)
		ret.push_back(new Sphere(sqrt(radius)*(1.15), Sphere::getCenter())) ;
	if(level > 2)
		ret.push_back(new Sphere(sqrt(radius)*(1.08), Sphere::getCenter())) ;
	return ret ;
}

std::vector<DelaunayTriangle *> VirtualInclusion3D::getTriangles( DelaunayTree * dt)  { 
	return std::vector<DelaunayTriangle *> (0) ;
}

std::vector<DelaunayTetrahedron *> VirtualInclusion3D::getTetrahedrons(const DelaunayTree_3D * dt)  {
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->conflicts(dynamic_cast<Sphere *>(this)) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren()->size() ;  j++)
		{
			if(this->getChild(j)->in(temp[i]->getCenter()))
			{
				inChild = true ; 
				break ;
			}
		}
		if(this->in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}


Point * VirtualInclusion3D::pointAfter(size_t i) { 
	return NULL ; 
}



void VirtualInclusion3D::sample(size_t n)
{
}
