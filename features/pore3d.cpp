//
// C++ Implementation: pore3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "pore3d.h"
#include "../physics/physics_base.h"
#include "../physics/void_form.h"
#include "pore.h"

 namespace Amie {


Pore3D::Pore3D(Feature *father, double r, double x, double y, double z) : Feature(father), Sphere(r, x, y, z)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

Pore3D::Pore3D(Feature *father, double r, Point center) : Feature(father), Sphere(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

Pore3D::Pore3D(double r, double x, double y, double z): Feature(nullptr), Sphere(r, x,y,z)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

Pore3D::Pore3D(double r, Point center) : Feature(nullptr), Sphere(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
}

bool Pore3D::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
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

void Pore3D::setRadius(double newR)
{
	Sphere::setRadius(newR);
	isUpdated = true ;
}

std::vector<DelaunayTriangle *> Pore3D::getElements2D( FeatureTree * dt){
	return std::vector<DelaunayTriangle *>(0) ;
}

std::vector<DelaunayTetrahedron *> Pore3D::getElements3D( FeatureTree * dt){
	std::vector<DelaunayTetrahedron *> ret  ;
	
	std::vector<DelaunayTetrahedron *> temp = dt->get3DMesh()->getConflictingElements(this->getPrimitive()) ;
	
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


void Pore3D::sample(double linearDensity, double surfaceDensityFactor, Sampler * sampler) {
	this->Sphere::sampleSurface(linearDensity, surfaceDensityFactor, sampler) ;
	
// 	size_t numberOfRadii = static_cast<size_t>(sqrt(std::max(4*n, (size_t)27))) ;
// 	double r = sqrt(radius) ;
// 	double radiusExt = (r + (r*(2)/(double)numberOfRadii) ) ;
// 	
// 	delete this->boundary ;
// 	
// 	this->boundary = new Sphere(radiusExt, getCenter()) ;
	
// 	this->Sphere::sampleBoundingSurface(4*n) ;/*BoundingSurface ;*/
// 	this->Sphere::sampleSurface(n/4) ;
// // // 	this->Sphere::sampleSurface((size_t)round(sqrt(n))) ;
// 	this->inPoints->resize(1) ;
// 	(*this->inPoints)[0] = new Point(this->center.getX(), this->center.getY(), this->center.getZ()) ;
}

 }
