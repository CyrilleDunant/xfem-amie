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

using namespace Mu ;

std::vector<DelaunayTriangle *> Inclusion::getElements( Mesh<DelaunayTriangle> * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->getConflictingElements(this->getPrimitive()) ;
	
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

bool Inclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}


std::vector<DelaunayTriangle *> TriangularInclusion::getElements( Mesh<DelaunayTriangle> * dt) 
{
	std::vector<DelaunayTriangle *> ret ;
	
	std::vector<DelaunayTriangle *>  temp = dt->getConflictingElements(this->getPrimitive()) ;
	
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

bool TriangularInclusion::interacts(Feature * f, double d) const
{
	for(PointSet::const_iterator i =this->begin() ; i < this->end() ; i++)
		if(f->inBoundary(*(*i), d))
			return true ;
	return false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, double a, double b, double originX, double originY, double axisX, double axisY) : Ellipse(a,b,originX,originY,axisX,axisY), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, double a, double b, double originX, double originY) : Ellipse(a,b,originX,originY), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, double a, double b, const Point center, const Point axis) : Ellipse(a,b,center,axis), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(Feature *father, double a, double b, const Point center) : Ellipse(a,b,center), Feature(father)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(double a, double b, double originX, double originY, double axisX, double axisY) : Ellipse(a,b,originX,originY,axisX,axisY), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(double a, double b, double originX, double originY) : Ellipse(a,b,originX,originY), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}

EllipsoidalInclusion::EllipsoidalInclusion(double a, double b, const Point center, const Point axis)  : Ellipse(a,b,center,axis), Feature(NULL)
{
	this->isEnrichmentFeature = false ;
}
	
EllipsoidalInclusion::EllipsoidalInclusion(double a, double b, const Point center)  : Ellipse(a,b,center), Feature(NULL)
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
		ret.push_back(new Ellipse(getRadius() * 2., this->Ellipse::getMinorRadius(), this->Ellipse::getCenter(),this->Ellipse::getMajorAxis())) ;
	if(level > 1)
		ret.push_back(new Ellipse(getRadius() * 1.5, this->Ellipse::getMinorRadius(), this->Ellipse::getCenter(),this->Ellipse::getMajorAxis())) ;
	if(level > 2)
		ret.push_back(new Ellipse(getRadius() * 1.1, this->Ellipse::getMinorRadius(), this->Ellipse::getCenter(),this->Ellipse::getMajorAxis())) ;
	return ret ;
}
	
std::vector<DelaunayTriangle *> EllipsoidalInclusion::getElements( Mesh<DelaunayTriangle> * dt) 
{
	std::vector<DelaunayTriangle *>ret;
	
	std::vector<DelaunayTriangle *>temp = dt->getConflictingElements(this->getPrimitive()) ;
	
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

Point * EllipsoidalInclusion::pointAfter(size_t i)
{
	Point bi = getBoundingPoint(i) ;
	bi = bi - this->Ellipse::getCenter() ;
	Point bip = getBoundingPoint((i+1)%this->boundingPoints.size()) ;
	bip = bip - this->Ellipse::getCenter() ;
	double theta_i = atan2(bi * this->Ellipse::getMinorAxis(), bi * this->Ellipse::getMajorAxis()) ;
	double theta_ip = atan2(bip * this->Ellipse::getMinorAxis(), bip * this->Ellipse::getMajorAxis()) ;
//	double theta_i = atan2((boundingPoints[i] - this->Ellipse::getCenter()) * this->Ellipse::getMinorAxis(), (boundingPoints[i] - this->Ellipse::getCenter()) * this->Ellipse::getMajorAxis()) ;
//	double theta_ip = atan2((boundingPoints[(i+1)%this->boundingPoints.size()]-this->Ellipse::getCenter()) * this->Ellipse::getMinorAxis(), (boundingPoints[(i+1)%this->boundingPoints.size()]-this->Ellipse::getCenter()) * this->Ellipse::getMajorAxis()) ;
	double theta = 0.5*theta_i + 0.5*theta_ip ;
	
	Point * to_insert = new Point(this->Ellipse::getCenter() + this->Ellipse::getMajorAxis() * (cos(theta)*this->Ellipse::getMajorRadius()) + this->Ellipse::getMinorAxis() * (cos(theta)*this->Ellipse::getMinorRadius())) ;
	std::valarray<Point *> temp(this->boundingPoints.size()+1) ;
	std::copy(&boundingPoints[0], &boundingPoints[i], &temp[0]) ;
	temp[i+1] = to_insert ;
	std::copy(&boundingPoints[i+1], &boundingPoints[this->boundingPoints.size()], &temp[i+2]) ;
	this->boundingPoints.resize(temp.size()) ;
	std::copy(&temp[0],&temp[temp.size()] , &boundingPoints[0]) ;
	return to_insert ;
}

void EllipsoidalInclusion::sample(size_t n)
{
	this->sampleSurface(n) ;

//	std::cout << "number of points => " << n << std::endl ;
//	std::cout << "sample ellipsoidal inclusion => done" << std::endl ;
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
  
	if(std::sqrt((squareDist2D(p,source->getCenter()))) > (source->getRadius()+getLength()))
	  return false ;
  	double alpha = dynamic_cast<Ellipse *>(source)->getMajorAxis().angle() ;
	Point prot((p-source->getCenter()).x*cos(-alpha)-(p-source->getCenter()).y*sin(-alpha), 
		   +(p-source->getCenter()).x*sin(-alpha)+(p-source->getCenter()).y*cos(-alpha)) ;
	
	Function x("x") ;
	Function y("y") ;
	x = x/(dynamic_cast<Ellipse *>(source)->getMajorRadius()+getLength()) ;
	y = y/(dynamic_cast<Ellipse *>(source)->getMinorRadius()+getLength()) ;
	Function ellipse = (x^2)+(y^2) - 1 ;
	
	double dx = VirtualMachine().eval(ellipse, prot) ;
	return dx < 0 ;
	
// 	Point temp(p) ;
// 	getSource()->project(&temp) ;
// /*	if((p-temp).norm() < getLength())
// 	  std::cout << getSource()->in(p) << std::endl ;*/
// 	return (((p - temp).norm() < getLength()) && !getSource()->in(p)) ;
}

Form * ITZFeature::getBehaviour( const Point & p)
{
	double alpha = dynamic_cast<Ellipse *>(source)->getMajorAxis().angle() ;
	Point prot((p-source->getCenter()).x*cos(-alpha)-(p-source->getCenter()).y*sin(-alpha), 
		   +(p-source->getCenter()).x*sin(-alpha)+(p-source->getCenter()).y*cos(-alpha)) ;
	
	Function x("x") ;
	Function y("y") ;
	x = x/dynamic_cast<Ellipse *>(source)->getMajorRadius() ;
	y = y/dynamic_cast<Ellipse *>(source)->getMinorRadius() ;
	Function ellipse = (x^2)+(y^2) - 1 ;
	
	
	double dx = VirtualMachine().eval(ellipse, prot) ;
	
	Ellipse ell(dynamic_cast<Ellipse *>(source)->getMajorRadius(),
		    dynamic_cast<Ellipse *>(source)->getMinorRadius(),
		    0,0,1,0) ;
	double theta = prot.angle() ;
	Point pell(std::cos(theta),std::sin(theta)) ;
	double r = ell.getRadiusOnEllipse(theta) ;
	Point ptry = pell * r ;
	double dist = squareDist2D(ptry,prot) ;

	int signn = + 1 ;
	double thetan = theta + 0.01 ;
	Point pelln(std::cos(thetan),std::sin(thetan)) ;
	double rn = ell.getRadiusOnEllipse(thetan) ;
	Point ptryn = pelln * rn ;
	double distn = squareDist2D(ptryn,prot) ;
	if(distn > dist)
	{
	  signn = - 1 ;
	  thetan = theta - 0.001 ;
	  pelln.x = std::cos(thetan) ;
	  pelln.y = std::sin(thetan) ;
	  rn = ell.getRadiusOnEllipse(thetan) ;
	  ptryn = pelln * rn ;
	  distn = squareDist2D(ptryn,prot) ;
	}
	int i = 0 ;
	
	while(distn < dist)
	{
	  i++ ;
	  dist = distn ;
	  thetan = thetan + signn * 0.01 ;
	  pelln.x = std::cos(thetan) ;
	  pelln.y = std::sin(thetan) ;
	  rn = ell.getRadiusOnEllipse(thetan) ;
	  ptryn = pelln * rn ;
	  distn = squareDist2D(ptryn,prot) ;
	}
		
//	if(dx>0)
//	  std::cout << sqrt((distn+dist)/2) << " ; " << getLength() << std::endl ;
	
	static_cast<SpatiallyDistributedStiffness *>(behaviour)->setDistance(sqrt((distn+dist)/2)) ;
	
	return this->behaviour ;
}



