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
#include "../physics/physics_base.h"
#include "../physics/void_form.h"
#include "../utilities/xml.h"


using namespace Mu ;

Inclusion3D::Inclusion3D(Feature *father, double r, double x, double y, double z) : Sphere(r, x, y,z ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}
	
Inclusion3D::Inclusion3D(Feature *father, double r,  Point center) : Sphere(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}
	
Inclusion3D::Inclusion3D(double r, double x, double y, double z) : Sphere(r, x,y,z ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}
	
Inclusion3D::Inclusion3D(double r, Point center) : Sphere(r,center ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

XMLTree * Inclusion3D::toXML()
{
	XMLTree * inc = new XMLTree("inclusion3d") ;
	XMLTree * geom = new XMLTree("geometry") ;
	XMLTree * c = static_cast<Sphere *>(this)->toXML() ;
	geom->addChild(c) ;
	inc->addChild(geom) ;
	return inc ;
}


bool Inclusion3D::interacts(Feature * f, double d) const 	
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
	if(f->inBoundary(*(*i), d))
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

std::vector<DelaunayTriangle *> Inclusion3D::getElements( Mesh<DelaunayTriangle, DelaunayTreeItem> * dt)  { 
	return std::vector<DelaunayTriangle *> (0) ;
}
	
std::vector<DelaunayTetrahedron *> Inclusion3D::getElements( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dt)  {
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->getConflictingElements(this->getPrimitive()) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren().size() ;  j++)
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
	
Point * Inclusion3D::pointAfter(size_t i) 
{ 
	return NULL ; 
}
	
void Inclusion3D::sample(size_t n)
{
	this->Sphere::sampleSurface(n) ;
}

OctahedralInclusion::OctahedralInclusion(Feature *father, double r, double x, double y, double z) :RegularOctahedron(r, x, y,z ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

OctahedralInclusion::OctahedralInclusion(Feature *father, double r,  Point center):RegularOctahedron(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

OctahedralInclusion::OctahedralInclusion(double r, double x, double y, double z) : RegularOctahedron(r, x,y,z ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

OctahedralInclusion::OctahedralInclusion(double r, Point center) : RegularOctahedron(r,center ), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}
	
bool OctahedralInclusion::interacts(Feature * f, double d) const 
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
	if(f->inBoundary(*(*i), d))
		return true ;
return false ;
}
	
std::vector<Geometry *> OctahedralInclusion::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new RegularOctahedron(getRadius()*(1.2), RegularOctahedron::getCenter())) ;
	if(level > 1)
		ret.push_back(new RegularOctahedron(getRadius()*(1.15), RegularOctahedron::getCenter())) ;
	if(level > 2)
		ret.push_back(new RegularOctahedron(getRadius()*(1.08), RegularOctahedron::getCenter())) ;
	return ret ;
}
	
std::vector<DelaunayTriangle *> OctahedralInclusion::getTriangles( DelaunayTree * dt) 
{
	return std::vector<DelaunayTriangle *> (0) ;
}
	
std::vector<DelaunayTetrahedron *> OctahedralInclusion::getTetrahedrons( DelaunayTree3D * dt)
{
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->conflicts(this->getPrimitive()) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren().size() ;  j++)
		{
			if(this->getChild(j)->in(temp[i]->getCenter()))
			{
				inChild = true ; 
				break ;
			}
		}
		if(in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

Point * OctahedralInclusion::pointAfter(size_t i)
{
	return NULL ; 
}

void OctahedralInclusion::sample(size_t n)
{
	this->RegularOctahedron::sampleSurface(n) ;
}




VirtualInclusion3D::VirtualInclusion3D(Feature *father, double r, double x, double y, double z) : Sphere(r, x, y,z ), VirtualFeature(father)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

VirtualInclusion3D::VirtualInclusion3D(Feature *father, double r,  Point center) : Sphere(r, center ), VirtualFeature(father)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

VirtualInclusion3D::VirtualInclusion3D(double r, double x, double y, double z) : Sphere(r, x,y,z ), VirtualFeature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

VirtualInclusion3D::VirtualInclusion3D(double r, Point center) : Sphere(r,center ), VirtualFeature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

bool VirtualInclusion3D::interacts(Feature * f, double d) const 	{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
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

Feature * VirtualInclusion3D::getSource()
{
	return this ;
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

std::vector<DelaunayTetrahedron *> VirtualInclusion3D::getTetrahedrons( DelaunayTree3D * dt)  {
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->conflicts(dynamic_cast<Sphere *>(this)) ;
	
	for(size_t i = 0 ; i < temp.size() ; i++)
	{
		bool inChild = false ;
		for(size_t j = 0 ;  j< this->getChildren().size() ;  j++)
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
