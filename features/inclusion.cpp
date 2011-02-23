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
#include "../physics/spatially_distributed_stiffness.h"
#include "../polynomial/vm_base.h"
#include "../utilities/xml.h"
#include "../geometry/geometry_base.h"


using namespace Mu ;

std::vector<DelaunayTriangle *> Inclusion::getElements2D( FeatureTree* dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->getElements2D(getPrimitive()) ;
	
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

void Inclusion::setRadius(double newR)
{
	Circle::setRadius(newR);
	isUpdated = true ;
}

Inclusion::Inclusion(Feature * father,double r, double x, double y) : Circle(r, x, y ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(Feature * father,double r, Point center) : Circle(r, center ), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(double r, double x, double y) : Circle(r, x, y), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(double r, Point center) : Circle(r, center),  Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(Circle c) : Circle(c.getRadius(), c.getCenter()), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(XMLTree *  xml) : Circle(xml->getChild(0)), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}


XMLTree * Inclusion::toXML()
{
	XMLTree * inc = new XMLTree("inclusion") ;
	XMLTree * c = this->Circle::toXML() ;
	XMLTree * b = behaviour->toXML() ;
	inc->addChild(c) ;
	inc->addChild(b) ;
	return inc ;
}


void Inclusion::sample(size_t n)
{
	this->sampleSurface(n*5/4) ;
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

bool Inclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}


std::vector<DelaunayTriangle *> TriangularInclusion::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *> ret ;
	
	std::vector<DelaunayTriangle *>  temp = dt->getElements2D(this->getPrimitive()) ;
	
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

TriangularInclusion::TriangularInclusion(Feature * father,const Point & a, const Point & b, const Point & c) : Triangle(a, b, c), Feature(father)
{
	this->isEnrichmentFeature = false ;

}

TriangularInclusion::TriangularInclusion(const Point & a, const Point & b, const Point & c) : Triangle(a, b, c), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}


void TriangularInclusion::sample(size_t n)
{
// 	delete this->boundary ;
// 	this->boundary = new Circle(radius + radius/(0.25*n), this->Circle::getCenter()) ;
	this->sampleSurface(n/2) ;
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

bool TriangularInclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}


EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, Point center, Point a, Point b) : Ellipse(center, a, b), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, Point center, Point a, double b) : Ellipse(center, a, b), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Point center, Point a, Point b) : Ellipse(center, a, b), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Point center, Point a, double b) : Ellipse(center, a, b), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}



bool EllipsoidalInclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}
	
std::vector<Geometry *> EllipsoidalInclusion::getRefinementZones(size_t level) const
{
	std::vector<Geometry *> ret ;
	if(level > 0)
		ret.push_back(new Ellipse(this->getCenter(), this->getMajorAxis() * 2., this->getMinorAxis() * 2.)) ;
	if(level > 1)
		ret.push_back(new Ellipse(this->getCenter(), this->getMajorAxis() * 1.5, this->getMinorAxis() * 1.5)) ;
	if(level > 2)
		ret.push_back(new Ellipse(this->getCenter(), this->getMajorAxis() * 1.1, this->getMinorAxis() * 1.1)) ;
	return ret ;
}
	
std::vector<DelaunayTriangle *> EllipsoidalInclusion::getElements2D( FeatureTree * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->getElements2D(this->getPrimitive()) ;
	
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

void EllipsoidalInclusion::sample(size_t n)
{
	this->sampleSurface(n) ;

}

ITZFeature::ITZFeature(Feature *father, Feature * g, const Matrix & m, const Matrix & p, double l, double ca,double cb) : LevelSet(g), VirtualFeature(father)
{
	behaviour = new SpatiallyDistributedStiffness(m,p,l,ca,cb) ;
	source = g ;
	length = l ;
	this->isEnrichmentFeature = false ;
}

bool ITZFeature::in(const Point & p) const
{
        Point proj(p.x,p.y) ;
        source->project(&proj) ;
        return dist(p,proj) < getLength() ;

}

Form * ITZFeature::getBehaviour( const Point & p)
{
        Point proj(p.x,p.y) ;
        source->project(&proj) ;
	static_cast<SpatiallyDistributedStiffness *>(behaviour)->setDistance(dist(p,proj)) ;
	
	return this->behaviour->getCopy() ;
}



