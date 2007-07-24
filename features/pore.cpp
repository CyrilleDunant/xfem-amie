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

using namespace Mu ;

std::vector<DelaunayTriangle *> Pore::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *> ret;
	
	std::vector<DelaunayTriangle *> temp = dt->conflicts(this->boundary) ;
	
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
		if(this->boundary->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

Pore::Pore(Feature *father, double r, double x, double y): Feature(father), Circle(r, x, y)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Circle(r+r*0.07, x, y) ;
}

Pore::Pore(Feature *father, double r, Point center): Feature(father), Circle(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Circle(r+r*0.07, center) ;
}

Pore::Pore(double r, double x, double y) : Feature(NULL), Circle(r, x, y)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Circle(r+r*0.07, x, y) ;
}

Pore::Pore(double r, Point center) :  Feature(NULL), Circle(r, center)
{
	this->isEnrichmentFeature = false ;
	this->behaviour = new VoidForm() ;
	this->boundary = new Circle(r+r*0.07, center) ;
}

bool Pore::interacts(Feature * f) const 
{
	return this->boundary->intersects(f->getBoundary()) ;
// 	for(Point ** i =this->begin() ; i < this->end() ; i++)
// 		if(f->inBoundary((*i)))
// 			return true ;
// 	return false ;
}	

void Pore::sample(size_t n)
{
	delete this->boundary ;
	double numberOfRings =round((double)n/(2. * M_PI )) ;
	double r = getRadius()*(1. + .6/(numberOfRings+1)) ;
	this->boundary = new Circle(r, this->getCenter()) ;
	this->Circle::sampleSurface(n) ;
// 	this->Circle::sampleBoundingSurface(n) ;
// 	this->inPoints->resize(1) ;
// 	(*this->inPoints)[0] = new Point(this->center.x, this->center.y) ;
}

Point * Pore::pointAfter(size_t i)
{
	double theta_i = atan2((*this->boundingPoints)[i]->y, (*this->boundingPoints)[i]->x) ;
	double theta_ip = atan2((*this->boundingPoints)[(i+1)%this->boundingPoints->size()]->y, (*this->boundingPoints)[(i+1)%this->boundingPoints->size()]->x) ;
	double theta = 0.5*theta_i + 0.5*theta_ip ;
	
	Point * to_insert = new Point(cos(theta)*this->getRadius()+ this->Circle::getCenter().x, sin(theta)*this->getRadius()+ this->Circle::getCenter().y) ;
	std::valarray<Point *> temp(this->boundingPoints->size()+1) ;
	std::copy(&(*this->boundingPoints)[0], &(*this->boundingPoints)[i], &temp[0]) ;
	temp[i+1] = to_insert ;
	std::copy(&(*this->boundingPoints)[i+1], &(*this->boundingPoints)[this->boundingPoints->size()], &temp[i+2]) ;
	this->boundingPoints->resize(temp.size()) ;
	std::copy(&temp[0],&temp[temp.size()] , &(*this->boundingPoints)[0]) ;
	return to_insert ;
}

std::vector<Geometry *> Pore::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
// 	ret.push_back(new Circle(radius + radius/(2.*this->boundingPoints->size()), this->Circle::getCenter())) ;
// 	ret.push_back(new Circle(radius + radius/(1.5*this->boundingPoints->size()), this->Circle::getCenter())) ;
	if(level > 0)
		ret.push_back(new Circle(radius + radius*(0.2), this->Circle::getCenter())) ;
	if(level > 1)
		ret.push_back(new Circle(radius + radius*(0.1), this->Circle::getCenter())) ;
	if(level > 2)
		ret.push_back(new Circle(radius + radius*(0.05), this->Circle::getCenter())) ;
	//ret.push_back(new Circle(radius + radius*(0.1), this->Circle::getCenter())) ;
	//	ret.push_back(new Circle(radius + radius/(0.1*this->boundingPoints->size()), this->Circle::getCenter())) ;
	
	return ret ;
}


std::vector<DelaunayTriangle *> TriangularPore::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *> ret ;
	
	std::vector<DelaunayTriangle *> temp = dt->conflicts(this->boundary) ;
	
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
		if(this->boundary->in(temp[i]->getCenter()) && inChild == false)
			ret.push_back(temp[i]) ;
	}
	return ret ;
}

TriangularPore::TriangularPore(Feature * father,const Point & a, const Point & b, const Point & c) : Feature(father), Triangle(a, b, c)
{
	this->behaviour = new VoidForm() ;
	this->isEnrichmentFeature = false ;
	Point va = a-getCenter() ;
	Point vb = a-getCenter() ;
	Point vc = a-getCenter() ;
	this->boundary = new Triangle(a+va*0.03,b+vb*0.03,c+vc*0.03) ;
	this->boundary2 = new Triangle(a-va*0.001,b-vb*0.001,c-vc*0.001) ;
}

TriangularPore::TriangularPore(const Point & a, const Point & b, const Point & c) : Feature(NULL) ,Triangle(a, b, c)
{
	this->behaviour = new VoidForm() ;
	this->isEnrichmentFeature = false ;
	Point va = a-getCenter() ;
	Point vb = b-getCenter() ;
	Point vc = c-getCenter() ;
	this->boundary = new Triangle(a+va*0.03,b+vb*0.03,c+vc*0.03) ;
	this->boundary2 = new Triangle(a-va*0.001,b-vb*0.001,c-vc*0.001) ;
}


Point * TriangularPore::pointAfter(size_t i)
{
	return NULL ;
}


void TriangularPore::sample(size_t n)
{
// 	delete this->boundary ;
// 	this->boundary = new Circle(radius + radius/(0.25*n), this->Circle::getCenter()) ;
	this->sampleBoundingSurface(n*3) ;
}


std::vector<Geometry *> TriangularPore::getRefinementZones(size_t level) const
{
	Point a = getBoundingPoint(0) ;
	Point b = getBoundingPoint(boundingPoints->size()/3) ;
	Point c = getBoundingPoint(2*boundingPoints->size()/3) ;
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

bool TriangularPore::interacts(Feature * f) const
{
	for(Point ** i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;

}

