//
// C++ Implementation: pore
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pore.h"
#include "../physics/physics_base.h"
#include "../physics/void_form.h"

using namespace Mu ;

std::vector<DelaunayTriangle *> Pore::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *> ret;
	
	std::vector<DelaunayTriangle *> temp = dt->getElements2D(this->getPrimitive()) ;
	
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
		if(this->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

Pore::Pore(Feature *father, double r, double x, double y): Feature(father), Circle(r, x, y)
{
	this->isEnrichmentFeature = false ;
	delete this->behaviour ; 
	this->behaviour = new VoidForm() ;
}

Pore::Pore(Feature *father, double r, Point center): Feature(father), Circle(r, center)
{
	this->isEnrichmentFeature = false ;
	delete this->behaviour ; 
	this->behaviour = new VoidForm() ;
}

Pore::Pore(double r, double x, double y) : Feature(NULL), Circle(r, x, y)
{
	this->isEnrichmentFeature = false ;
	delete this->behaviour ; 
	this->behaviour = new VoidForm() ;
}

Pore::Pore(double r, Point center) :  Feature(NULL), Circle(r, center)
{
	this->isEnrichmentFeature = false ;
	delete this->behaviour ; 
	this->behaviour = new VoidForm() ;
}

void Pore::setRadius(double newR)
{
	Circle::setRadius(newR);
	isUpdated = true ;
}

bool Pore::interacts(Feature * f, double d) const 
{
	for(PointSet::const_iterator i = this->begin() ; i != this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}	

void Pore::sample(size_t n)
{
	this->Circle::sampleSurface(n) ;
// 	this->Circle::sampleBoundingSurface(n) ;
// 	this->inPoints->resize(1) ;
// 	(*this->inPoints)[0] = new Point(this->center.x, this->center.y) ;
}

std::vector<Geometry *> Pore::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
// 	ret.push_back(new Circle(radius + radius/(2.*boundingPoints.size()), this->Circle::getCenter())) ;
// 	ret.push_back(new Circle(radius + radius/(1.5*boundingPoints.size()), this->Circle::getCenter())) ;
	if(level > 0)
		ret.push_back(new Circle(radius + radius*(0.2), this->Circle::getCenter())) ;
	if(level > 1)
		ret.push_back(new Circle(radius + radius*(0.1), this->Circle::getCenter())) ;
	if(level > 2)
		ret.push_back(new Circle(radius + radius*(0.05), this->Circle::getCenter())) ;
	//ret.push_back(new Circle(radius + radius*(0.1), this->Circle::getCenter())) ;
	//	ret.push_back(new Circle(radius + radius/(0.1*boundingPoints.size()), this->Circle::getCenter())) ;
	
	return ret ;
}


std::vector<DelaunayTriangle *> TriangularPore::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *> ret ;
	std::vector<DelaunayTriangle *> temp = dt->getElements2D(this->getPrimitive()) ;
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
		if(this->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

TriangularPore::TriangularPore(Feature * father,const Point & a, const Point & b, const Point & c) : Feature(father), Triangle(a, b, c)
{
	this->behaviour = new VoidForm() ;
	this->isEnrichmentFeature = false ;
}

TriangularPore::TriangularPore(const Point & a, const Point & b, const Point & c) : Feature(NULL) ,Triangle(a, b, c)
{
	this->behaviour = new VoidForm() ;
	this->isEnrichmentFeature = false ;
}


void TriangularPore::sample(size_t n)
{
// 	n = std::max(3*n, (size_t)32) ;

	this->sampleSurface(n*1.5) ;
}


std::vector<Geometry *> TriangularPore::getRefinementZones(size_t level) const
{
	Point a = getBoundingPoint(0) ;
	Point b = getBoundingPoint(boundingPoints.size()/3) ;
	Point c = getBoundingPoint(2*boundingPoints.size()/3) ;
	Point va = a-getCenter() ;
	Point vb = b-getCenter() ;
	Point vc = c-getCenter() ;
	std::vector<Geometry *> ret ;
	
	if(level > 0)
		ret.push_back(new Triangle(a+va*0.2,b+vb*0.2,c+vc*0.2)) ;
	if(level > 1)
		ret.push_back(new Triangle(a+va*0.15,b+vb*0.15,c+vc*0.15)) ;
	if(level > 2)
		ret.push_back(new Triangle(a+va*0.1,b+vb*0.1,c+vc*0.1)) ;
	return ret ;
}

bool TriangularPore::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;

}

