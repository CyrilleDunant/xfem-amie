//
// C++ Implementation: pore3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "pore3d.h"

using namespace Mu ;


Pore3D::Pore3D(Feature *father, double r, double x, double y, double z) : Feature(father), Sphere(r, x, y, z)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Sphere(r*1.07, x, y,z) ;
}

Pore3D::Pore3D(Feature *father, double r, Point center) : Feature(father), Sphere(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Sphere(r*1.07, center) ;
}

Pore3D::Pore3D(double r, double x, double y, double z): Feature(NULL), Sphere(r, x,y,z)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Sphere(r*1.07, x,y,z) ;
}

Pore3D::Pore3D(double r, Point center) : Feature(NULL), Sphere(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Sphere(r*1.07, center) ;
}

bool Pore3D::interacts(Feature * f) const
{
	return this->boundary->intersects(f->getBoundary()) ;
}
	
std::vector<Geometry *> Pore3D::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Sphere(sqrt(radius)*(1.15), Sphere::getCenter())) ;
	if(level > 1)
		ret.push_back(new Sphere(sqrt(radius)*(1.07), Sphere::getCenter())) ;
	if(level > 2)
		ret.push_back(new Sphere(sqrt(radius)*(1.027), Sphere::getCenter())) ;
	return ret ;
}


std::vector<DelaunayTriangle *> Pore3D::getTriangles( DelaunayTree * dt){
	return std::vector<DelaunayTriangle *>(0) ;
}

std::vector<DelaunayTetrahedron *> Pore3D::getTetrahedrons( DelaunayTree_3D * dt){
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->conflicts(this->boundary) ;
	
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
		if(this->boundary->in(temp[i]->getCenter()) || inChild )
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

Point * Pore3D::pointAfter(size_t i) {
	return NULL ;
}

void Pore3D::sample(size_t n) {
	this->Sphere::sampleSurface(2*n) ;
	
	size_t numberOfRadii = static_cast<size_t>(sqrt(std::max(4*n, (size_t)27))) ;
	double r = sqrt(radius) ;
	double radiusExt = (r + (r*(2)/(double)numberOfRadii) ) ;
	
	delete this->boundary ;
	
	this->boundary = new Sphere(radiusExt, getCenter()) ;
	
// 	this->Sphere::sampleBoundingSurface(4*n) ;/*BoundingSurface ;*/
// 	this->Sphere::sampleSurface(n/4) ;
// // // 	this->Sphere::sampleSurface((size_t)round(sqrt(n))) ;
// 	this->inPoints->resize(1) ;
// 	(*this->inPoints)[0] = new Point(this->center.x, this->center.y, this->center.z) ;
}
