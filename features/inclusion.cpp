//
// C++ Implementation: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
			if(this->getChild(j)->in(temp[i]->getCenter()) /*&& !(this->getChild(j)->isVirtualFeature)*/)
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

Inclusion::Inclusion(double r, double x, double y) : Circle(r, x, y), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(double r, Point center) : Circle(r, center),  Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(Circle c) : Circle(c.getRadius(), c.getCenter()), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}

Inclusion::Inclusion(XMLTree *  xml) : Circle(xml->getChild(0)), Feature(nullptr)
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

TriangularInclusion::TriangularInclusion(const Point & a, const Point & b, const Point & c) : Triangle(a, b, c), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}


void TriangularInclusion::sample(size_t n)
{
// 	delete this->boundary ;
// 	this->boundary = new Circle(radius + radius/(0.25*n), this->Circle::getCenter()) ;
	this->sampleSurface(n*5/4) ;
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

std::vector<DelaunayTriangle *> RectangularInclusion::getElements2D( FeatureTree * dt) 
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

RectangularInclusion::RectangularInclusion(Feature * father,const Point & a, const Point & b, const Point & c, const Point & d) : OrientedRectangle(a, b, c, d), Feature(father)
{
	this->isEnrichmentFeature = false ;
	
}

RectangularInclusion::RectangularInclusion(const Point & a, const Point & b, const Point & c, const Point & d) : OrientedRectangle(a, b, c, d), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}


void RectangularInclusion::sample(size_t n)
{
	this->sampleSurface(n*3/2) ;
}

std::vector<Geometry *> RectangularInclusion::getRefinementZones(size_t level) const
{
	Point a = getBoundingPoint(0) ;
	Point b = getBoundingPoint(boundingPoints.size()/4) ;
	Point c = getBoundingPoint(2*boundingPoints.size()/4) ;
	Point d = getBoundingPoint(3*boundingPoints.size()/4) ;
	Point va = a-getCenter() ;
	Point vb = b-getCenter() ;
	Point vc = c-getCenter() ;
	Point vd = d-getCenter() ;
	std::vector<Geometry *> ret ;
	
	if(level > 0)
		ret.push_back(new OrientedRectangle(a+va*0.2,b+vb*0.2,c+vc*0.2,d+vd*0.2)) ;
	if(level > 1)
		ret.push_back(new OrientedRectangle(a+va*0.15,b+vb*0.15,c+vc*0.15,d+vd*0.15)) ;
	if(level > 2)
		ret.push_back(new OrientedRectangle(a+va*0.1,b+vb*0.1,c+vc*0.1,d+vd*0.1)) ;
	return ret ;
}

bool RectangularInclusion::interacts(Feature * f, double d) const
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

EllipsoidalInclusion::EllipsoidalInclusion(Point center, Point a, Point b) : Ellipse(center, a, b), Feature(nullptr)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Point center, Point a, double b) : Ellipse(center, a, b), Feature(nullptr)
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
	this->sampleSurface(n*7/4) ;

}

ITZFeature::ITZFeature(Feature *father, Feature * g, double l, double a, double b) : LevelSet(g), VirtualFeature(father)
{
	min = a ;
	max = b ;
//	behaviour = new SpatiallyDistributedStiffness(m,p,l,ca,cb) ;
	sources.push_back(g) ;
	length = l ;
	this->isEnrichmentFeature = false ;
}

ITZFeature::ITZFeature(Feature *father, std::vector<Inclusion *> & g, double l, double a, double b) : LevelSet(dynamic_cast<Circle *>(g[0])), VirtualFeature(father)
{
	min = a ;
	max = b ;
//	behaviour = new SpatiallyDistributedStiffness(m,p,l,ca,cb) ;
	for(size_t i = 0 ; i < g.size() ; i++)
		sources.push_back(g[i]) ;
	length = l ;
	this->isEnrichmentFeature = false ;
}

bool ITZFeature::in(const Point & p) const
{
	double d = -1 ;
	for(size_t i = 0 ; i < sources.size() ; i++)
	{
		Point proj(p.x,p.y) ;
		Circle r(sources[i]->getRadius(), sources[i]->getCenter().x, sources[i]->getCenter().y-getLength()*0.5) ;
		r.project(&proj) ;
		if(sources[i]->in(p))
			return false ;
		if(d < 0 || dist(p,proj) < d)
			d = dist(p,proj) ;
	}
        return d < getLength() ;

}

Form * ITZFeature::getBehaviour( const Point & p)
{
	double d = -1 ;
	for(size_t i = 0 ; i < sources.size() ; i++)
	{
		Point proj(p.x,p.y) ;
		Circle r(sources[i]->getRadius(), sources[i]->getCenter().x, sources[i]->getCenter().y-getLength()*0.5) ;
		r.project(&proj) ;
		if(sources[i]->in(p))
			continue ;
		if(d < 0 || dist(p,proj) < d)
			d = dist(p,proj) ;
	}
	
	Form * ret = this->getFather()->getBehaviour(p)->getCopy() ;
	
//	static_cast<SpatiallyDistributedStiffness *>(behaviour)->setDistance(dist(p,proj)) ;
	
 	return ret ;
}



