//
// C++ Implementation: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "inclusion.h"

using namespace Mu ;

std::vector<DelaunayTriangle *> Inclusion::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->conflicts(this->boundary) ;
	
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

Inclusion::Inclusion(Feature * father,double r, double x, double y) : Circle(r, x, y ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(r+r*0.1, x, y) ;
	this->boundary2 = new Circle(r-r*0.1, x, y) ;
}

Inclusion::Inclusion(Feature * father,double r, Point center) : Circle(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(r+r*0.1, center) ;
	this->boundary2 = new Circle(r-r*0.1, center) ;
}

Inclusion::Inclusion(double r, double x, double y) : Circle(r, x, y), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(r+r*0.1, x, y) ;
	this->boundary2 = new Circle(r-r*0.1, x, y) ;
}

Inclusion::Inclusion(double r, Point center) : Circle(r, center),  Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	this->boundary = new Circle(r+r*0.1, center) ;
	this->boundary2 = new Circle(r-r*0.1, center) ;
}

Point * Inclusion::pointAfter(size_t i)
{
	double theta_i = atan2(boundingPoints[i]->y, boundingPoints[i]->x) ;
	double theta_ip = atan2(boundingPoints[(i+1)%this->boundingPoints.size()]->y, boundingPoints[(i+1)%this->boundingPoints.size()]->x) ;
	double theta = 0.5*theta_i + 0.5*theta_ip ;
	
	Point * to_insert = new Point(cos(theta)*this->getRadius()+ this->Circle::getCenter().x, sin(theta)*this->getRadius()+ this->Circle::getCenter().y) ;
	std::valarray<Point *> temp(this->boundingPoints.size()+1) ;
	std::copy(&boundingPoints[0], &boundingPoints[i], &temp[0]) ;
	temp[i+1] = to_insert ;
	std::copy(&boundingPoints[i+1], &boundingPoints[this->boundingPoints.size()], &temp[i+2]) ;
	this->boundingPoints.resize(temp.size()) ;
	std::copy(&temp[0],&temp[temp.size()] , &boundingPoints[0]) ;
	return to_insert ;
}


void Inclusion::sample(size_t n)
{
	delete this->boundary ;
	double numberOfRings =round((double)n/(2. * M_PI )) ;
	double r = getRadius()*(1.+ .6/(numberOfRings+1)) ;
	this->boundary = new Circle(r, getCenter()) ;
	std::cout << r << std::endl; 
	this->sampleSurface(n) ;
}

std::vector<Geometry *> Inclusion::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Circle(getRadius()*2., this->Circle::getCenter())) ;
	if(level > 1)
		ret.push_back(new Circle(getRadius() * 1.5, this->Circle::getCenter())) ;
	if(level > 2)
		ret.push_back(new Circle(getRadius() * 1.1, this->Circle::getCenter())) ;
	return ret ;
}

bool Inclusion::interacts(Feature * f) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;
}


std::vector<DelaunayTriangle *> TriangularInclusion::getTriangles( DelaunayTree * dt) 
{
	std::vector<DelaunayTriangle *> ret ;
	
	std::vector<DelaunayTriangle *>  temp = dt->conflicts(this->boundary) ;
	
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

TriangularInclusion::TriangularInclusion(Feature * father,const Point & a, const Point & b, const Point & c) : Triangle(a, b, c), Feature(father)
{
	this->isEnrichmentFeature = false ;
	Point va = a-getCenter() ;
	Point vb = a-getCenter() ;
	Point vc = a-getCenter() ;
	this->boundary = new Triangle(a+va*0.02,b+vb*0.02,c+vc*0.02) ;
	this->boundary2 = new Triangle(a-va*0.02,b-vb*0.02,c-vc*0.02) ;
}

TriangularInclusion::TriangularInclusion(const Point & a, const Point & b, const Point & c) : Triangle(a, b, c), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
	Point va = a-getCenter() ;
	Point vb = b-getCenter() ;
	Point vc = c-getCenter() ;
	this->boundary = new Triangle(a+va*0.02,b+vb*0.02,c+vc*0.02) ;
	this->boundary2 = new Triangle(a-va*0.02,b-vb*0.02,c-vc*0.02) ;
}


Point * TriangularInclusion::pointAfter(size_t i)
{
	return NULL ;
}


void TriangularInclusion::sample(size_t n)
{
// 	delete this->boundary ;
// 	this->boundary = new Circle(radius + radius/(0.25*n), this->Circle::getCenter()) ;
	this->sampleSurface(2*n) ;
}

std::vector<Geometry *> TriangularInclusion::getRefinementZones(size_t level) const
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

bool TriangularInclusion::interacts(Feature * f) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary((*i)))
			return true ;
	return false ;

}
