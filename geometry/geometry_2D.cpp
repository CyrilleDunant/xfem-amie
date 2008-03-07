// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_2D.h"

using namespace Mu ;

Triangle::Triangle() : ConvexGeometry(3)
{
	gType = TRIANGLE ;
	assert(this->size() == 3) ;
	
	boundingPoints[0] = new Point(0, 1) ;
	boundingPoints[1] = new Point(0, 0) ;
	boundingPoints[2] = new Point(1, 0) ;
	
	computeCenter() ;
	computeCircumCenter() ;
	radius = 0.471405;
	sqradius = radius*radius ;
}

Triangle::Triangle( const Point & p0,  const Point & p1,  const Point & p2) : ConvexGeometry(3)
{
	gType = TRIANGLE ;
	
	assert(this->size() == 3) ;
	
	boundingPoints[0] = new Point(p0) ;
	boundingPoints[1] = new Point(p1) ;
	boundingPoints[2] = new Point(p2) ;
	
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
			if(!isTrigoOriented())
				std::cout << "arrrgh !" << std::endl ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = sqrt((squareDist2D(p1, circumCenter)+ squareDist2D(p0, circumCenter))/2.);
	sqradius = radius*radius ;
	
}

Parallelogramme::Parallelogramme() : ConvexGeometry(4)
{
	gType =PARALLELOGRAMME;
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(-1,1);
	boundingPoints[1] = new Point(-1,-1);
	boundingPoints[2] = new Point(1, -1);
	boundingPoints[3] = new Point(1, 1);
	computeCenter() ;
}

Parallelogramme::Parallelogramme( const Point & p0,  const Point & p1,  const Point & p2,  const Point & p3): ConvexGeometry(4)
{
	gType = PARALLELOGRAMME ;
	
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(p0) ;
	boundingPoints[1] = new Point(p1) ;
	boundingPoints[2] = new Point(p2) ;
	boundingPoints[3] = new Point(p3) ;
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0.z == p1.z == p2.z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = (squareDist2D(p1, circumCenter) + squareDist2D(p0, circumCenter) + squareDist2D(p2, circumCenter)+squareDist2D(p3, circumCenter))/4.;
	
}

Parallelogramme::Parallelogramme( const Point *p0,  const Point *p1,  const Point *p2,  const Point *p3): ConvexGeometry(4)
{
	gType = PARALLELOGRAMME ;
	
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(*p0) ;
	boundingPoints[1] = new Point(*p1) ;
	boundingPoints[2] = new Point(*p2) ;
	boundingPoints[3] = new Point(*p3) ;
	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = (squareDist2D(*p1, circumCenter) + squareDist2D(*p0, circumCenter) + squareDist2D(*p2, circumCenter)+squareDist2D(*p3, circumCenter))/4.;
	
}


void Parallelogramme::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		this->center += this->getPoint(i) ;
	}
	
	this->center = this->center/this->size() ;
}

double Parallelogramme::getRadius() const
{
	return sqrt(radius) ;
}

Point Parallelogramme::getCircumCenter() const
{
	return this->circumCenter ;
}

void Parallelogramme::computeCircumCenter()
{	
	if (fabs(boundingPoints[1]->y-boundingPoints[0]->y) < 20*std::numeric_limits<double>::epsilon()) 
	{
		double m2 = - (boundingPoints[2]->x-boundingPoints[1]->x) / (boundingPoints[2]->y-boundingPoints[1]->y);
		double mx2 = (boundingPoints[1]->x + boundingPoints[2]->x) / 2.0;
		double my2 = (boundingPoints[1]->y + boundingPoints[2]->y) / 2.0;
		double xc = (boundingPoints[1]->x + boundingPoints[0]->x) / 2.0;
		double yc = fma(m2, (xc - mx2), my2);
		
		circumCenter = Point(xc, yc) ;
	} 
	else if (fabs(boundingPoints[2]->y-boundingPoints[1]->y) < 20*std::numeric_limits<double>::epsilon()) 
	{
		double m1 = - (boundingPoints[1]->x-boundingPoints[0]->x) / (boundingPoints[1]->y-boundingPoints[0]->y);
		double mx1 = (boundingPoints[0]->x + boundingPoints[1]->x) / 2.0;
		double my1 = (boundingPoints[0]->y + boundingPoints[1]->y) / 2.0;
		double xc = (boundingPoints[2]->x + boundingPoints[1]->x) / 2.0;
		double yc = fma(m1, (xc - mx1), my1);
		
		circumCenter = Point(xc, yc) ;
	} 
	else 
	{
		double m1 = - (boundingPoints[1]->x-boundingPoints[0]->x) / (boundingPoints[1]->y-boundingPoints[0]->y);
		double m2 = - (boundingPoints[2]->x-boundingPoints[1]->x) / (boundingPoints[2]->y-boundingPoints[1]->y);
		double mx1 = (boundingPoints[0]->x + boundingPoints[1]->x) / 2.0;
		double mx2 = (boundingPoints[1]->x + boundingPoints[2]->x) / 2.0;
		double my1 = (boundingPoints[0]->y + boundingPoints[1]->y) / 2.0;
		double my2 = (boundingPoints[1]->y + boundingPoints[2]->y) / 2.0;
		double xc = fma(m1, mx1, fma(- m2 , mx2, my2 - my1)) / (m1 - m2);
		double yc = fma(m1, (xc - mx1), my1);
		
		circumCenter = Point(xc, yc) ;
	}
}

bool Parallelogramme::inCircumCircle(const Point p) const
{
	return  squareDist2D(circumCenter, p) < radius  ;
}

bool Parallelogramme::inCircumCircle(const Point *p) const
{
	return  squareDist2D(circumCenter, (*p)) < radius  ;
}

bool Parallelogramme::is1D() const
{ 
	return false ; 
} 

double Parallelogramme::area() const
{
	assert(this->boundingPoints.size() == 4) ;
	
	Segment s0((*boundingPoints[0]), (*boundingPoints[1])) ;
	Segment s1((*boundingPoints[0]), (*boundingPoints[2])) ;

	return ((s0.vector())^(s1.vector())).norm() ;
}


void Parallelogramme::project(Point * p) const
{
	Segment s(getCenter(), *p) ;
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		if(s.on(getPoint(i)) && !Segment(getPoint(i), *p).on(getCenter()))
		{
			p->x = getPoint(i).x ;
			p->y = getPoint(i).y ;
			return ;
		}
	}
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		Segment seg(getPoint(i), getPoint((i+1)%boundingPoints.size())) ;
		if(s.intersects(seg))
		{
			Point t = s.intersection(seg);
			p->x = t.x ;
			p->y = t.y ;
			return ;
		}
	}
	
}


bool Parallelogramme::in(const Point &p) const
{
	return ConvexPolygon::in(p) ;
}

void Parallelogramme::sampleBoundingSurface(size_t num_points)
{
	assert(num_points%4== 0) ;
	
	Point v0(*boundingPoints[0]) ;
	Point v1(*boundingPoints[1]) ;
	Point v2(*boundingPoints[2]) ;
	Point v3(*boundingPoints[3]) ;
	this->boundingPoints.resize(num_points) ;
	
	boundingPoints[0] = new Point(v0) ;
	
	for(size_t i = 1 ; i < num_points/4; i++)
	{
		boundingPoints[i] = new Point(v0*3.*i/num_points + v1*(1.-3.*i/num_points)) ;
	}
	
	boundingPoints[num_points/3] = new Point(v1) ;
	
	for(size_t i = num_points/3+1 ; i < 2*num_points/3 ; i++)
	{
		boundingPoints[i] = new Point(v1*3.*(i-num_points/3)/num_points + v2*(1.-3.*(i-num_points/3)/num_points)) ;
	}
	
	boundingPoints[2*num_points/3] = new Point(v2) ;
	
	for(size_t i = 2*num_points/3+1 ; i < num_points ; i++)
	{
		boundingPoints[i] = new Point(v2*3.*(i-2*num_points/3)/num_points + v0*(1.-3.*(i-2*num_points/3)/num_points)) ;
	}
}

void Parallelogramme::sampleSurface(size_t num_points)
{
	if(this->size() == 4)
		this->sampleBoundingSurface(num_points) ;
}

Triangle::Triangle( Point *p0,  Point *p1,  Point *p2): ConvexGeometry(3)
{
	gType = TRIANGLE ;
	
	assert(this->size() == 3) ;
	

	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;

	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
			if(!isTrigoOriented())
				std::cout << "arrrgh !" << std::endl ;
		}
	}
	
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if(p0->z == p1->z == p2->z == 0)
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
		
	}
	
	radius = (dist(p1, &circumCenter)+ dist(p0, &circumCenter))/2.;
	sqradius = radius*radius ;
}


void Triangle::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		this->center += this->getPoint(i) ;
	}
	
	this->center = this->center/this->size() ;
}

double Triangle::getRadius() const
{
	return radius ;
}

Point Triangle::getCircumCenter() const
{
	return this->circumCenter ;
}

std::vector<Point> Triangle::getBoundingBox() const
{
	std::vector<Point> box ;
	box.push_back(getCircumCenter()+Point(getRadius(), -getRadius())) ;
	box.push_back(getCircumCenter()+Point(getRadius(), getRadius())) ;
	box.push_back(getCircumCenter()+Point(-getRadius(), getRadius())) ;
	box.push_back(getCircumCenter()+Point(-getRadius(), -getRadius())) ;
	
	return box ;
}

void Triangle::computeCircumCenter()
{	
	
	if (std::abs(getBoundingPoint(1).y-getBoundingPoint(0).y) < POINT_TOLERANCE) 
	{
		double m2  =  (getBoundingPoint(1).x-getBoundingPoint(2).x ) / (getBoundingPoint(2).y-getBoundingPoint(1).y);
		double mx2 = (getBoundingPoint(1).x + getBoundingPoint(2).x) ;
		double my2 = (getBoundingPoint(1).y + getBoundingPoint(2).y) ;
		double xc  = (getBoundingPoint(1).x + getBoundingPoint(0).x) ;
		double yc  = m2 * (xc - mx2) + my2;
		
		circumCenter.set(xc/2., yc/2.) ;
	} 
	else if (std::abs(getBoundingPoint(2).y-getBoundingPoint(1).y) < POINT_TOLERANCE) 
	{
		double m1  =  (getBoundingPoint(0).x - getBoundingPoint(1).x ) / (getBoundingPoint(1).y-getBoundingPoint(0).y);
		double mx1 = (getBoundingPoint(0).x + getBoundingPoint(1).x) ;
		double my1 = (getBoundingPoint(0).y + getBoundingPoint(1).y) ;
		double xc  = (getBoundingPoint(2).x + getBoundingPoint(1).x) ;
		double yc  = m1 * (xc - mx1) + my1;
		
		circumCenter.set(xc/2., yc/2.) ;
	} 
	else 
	{
		double m1  = (getBoundingPoint(0).x-getBoundingPoint(1).x) / (getBoundingPoint(1).y-getBoundingPoint(0).y);
		double m2  = (getBoundingPoint(1).x-getBoundingPoint(2).x) / (getBoundingPoint(2).y-getBoundingPoint(1).y);
		double mx1 = (getBoundingPoint(0).x + getBoundingPoint(1).x) ;
		double mx2 = (getBoundingPoint(1).x + getBoundingPoint(2).x) ;
		double my1 = (getBoundingPoint(0).y + getBoundingPoint(1).y) ;
		double my2 = (getBoundingPoint(1).y + getBoundingPoint(2).y) ;
		double xc  = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		double yc  = m1 * (xc - mx1) + my1;
		
		circumCenter.set(xc/2., yc/2.) ;
	}
}

bool Triangle::inCircumCircle(const Point & p) const
{
	double x = circumCenter.x -p.x ;
	double y = circumCenter.y -p.y ;
	return  fma(x, x, y*y)< sqradius - 2.*POINT_TOLERANCE  ;
}

bool Triangle::inCircumCircle(const Point *p) const
{
	double x = circumCenter.x -p->x ;
	double y = circumCenter.y -p->y ;
	return  fma(x, x, y*y) < sqradius - 2.*POINT_TOLERANCE  ;
}


double Triangle::area() const
{
// 	assert(this->boundingPoints.size() == 3) ;
	
	Segment s0(getBoundingPoint(0), getBoundingPoint(this->boundingPoints.size()/3)) ;
	Segment s1(getBoundingPoint(0), getBoundingPoint(this->boundingPoints.size()/3 * 2)) ;

	return 0.5*(s0.vector()^s1.vector()).norm() ;
}


void Triangle::project(Point * p) const
{
	Segment s(getCenter(), *p) ;
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		if(s.on(getBoundingPoint(i)) && !Segment(getBoundingPoint(i), *p).on(getCenter()))
		{
			p->x = getBoundingPoint(i).x ;
			p->y = getBoundingPoint(i).y ;
			return ;
		}
	}
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		Segment seg(getBoundingPoint(i), getBoundingPoint((i+1)%boundingPoints.size())) ;
		if(s.intersects(seg))
		{
			Point t = s.intersection(seg);
			p->x = t.x ;
			p->y = t.y ;
			return ;
		}
	}
	
}

bool Triangle::in(const Point &p) const
{
	return ConvexPolygon::in(p) ;
	Segment s0(getBoundingPoint(0), getBoundingPoint(getBoundingPoints().size()/3)) ;
	Segment s1(getBoundingPoint(getBoundingPoints().size()/3), getBoundingPoint(2*getBoundingPoints().size()/3)) ;
	Segment s2(getBoundingPoint(0), getBoundingPoint(2*getBoundingPoints().size()/3)) ;
	
	Segment test(p, getCenter()) ; 
	
	return !(s0.intersects(test) || s1.intersects(test)|| s2.intersects(test)) ;
}


void Triangle::sampleBoundingSurface(size_t num_points)
{
	assert(num_points%3 == 0) ;
	
	Point * v0 =boundingPoints[0] ;
	Point * v1 =boundingPoints[1] ;
	Point * v2 =boundingPoints[2] ;
	
	this->boundingPoints.resize(num_points) ;
	
	boundingPoints[0] = v0 ;
	
	for(size_t i = 1 ; i < num_points/3 ; i++)
	{
		boundingPoints[i] = new Point((*v0)*3.*i/num_points + (*v1)*(1.-3.*i/num_points)) ;
	}
	
	boundingPoints[num_points/3] = v1 ;
	
	for(size_t i = num_points/3+1 ; i < 2*num_points/3 ; i++)
	{
		boundingPoints[i] = new Point((*v1)*3.*(i-num_points/3)/num_points + (*v2)*(1.-3.*(i-num_points/3)/num_points)) ;
	}
	
	boundingPoints[2*num_points/3] = v2 ;
	
	for(size_t i = 2*num_points/3+1 ; i < num_points ; i++)
	{
		boundingPoints[i] = new Point((*v2)*3.*(i-2*num_points/3)/num_points + (*v0)*(1.-3.*(i-2*num_points/3)/num_points)) ;
	}
}

void Triangle::sampleSurface(size_t num_points)
{
	num_points -= num_points%3 ;
	
	if(this->size() == 3)
		this->sampleBoundingSurface(num_points) ;
	
	std::vector<Point> newPoints ;
	
	size_t end_i = 2*boundingPoints.size()/3 ;
	
	for(int i = 0 ; i < (int)num_points/3-1 ; i++)
	{
		for(int j = 0 ; j < i+1 ; j++)
		{
			double fact = (double)(j+1)/(double)(i+1) ;
			double d = dist(getBoundingPoint(i+1), getBoundingPoint(end_i-1-i))*0.17/(i+1.) ;
			double xrand = ((double)random()/(double)(RAND_MAX+1)*2.-1.)*d ;
			double yrand = ((double)random()/(double)(RAND_MAX+1)*2.-1.)*d ;
			newPoints.push_back(getBoundingPoint(i+1)*(1.-fact) + getBoundingPoint(end_i-1-i)*fact+Point(xrand, yrand)) ;
		}
	}
	
	this->inPoints.resize(newPoints.size()) ;
	
	for(size_t i = 0 ; i < inPoints.size() ; i++)
	{
		inPoints[i] = new Point(newPoints[i]) ;
	}
	
}

void Triangle::print() const 
{
	std::cout << "inPoints" << std::endl ;
	
	for(size_t i = 0 ; i < inPoints.size() ; i++ )
	{
		getInPoint(i).print() ;
	}
	std::cout << std::endl ;
	
		std::cout << "boundingPoints" << std::endl ;
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++ )
	{
		getBoundingPoint(i).print() ;
	}
	std::cout << std::endl ;
	
}

Rectangle::Rectangle(double x, double y, double originX, double originY) : ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE ;
	this->center = Point(originX,originY) ;
	topLeft = Point(originX-0.5*x, originY+0.5*y) ; boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(originX+0.5*x, originY+0.5*y) ; boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(originX+0.5*x, originY-0.5*y) ; boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft =  Point(originX-0.5*x, originY-0.5*y) ; boundingPoints[3] = new Point(bottomLeft) ;
	
	
} 

Rectangle::Rectangle(double x, double y, Point &center) :  ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE ;
	this->center = center ;
	topLeft = Point(center.x-0.5*x, center.y+0.5*y) ; boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(center.x+0.5*x, center.y+0.5*y) ; boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(center.x+0.5*x, center.y-0.5*y) ; boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft =  Point(center.x-0.5*x, center.y-0.5*y) ; boundingPoints[3] = new Point(bottomLeft) ;
	
	
}

Rectangle::Rectangle() :  ConvexGeometry(4), size_y(2), size_x(2)
{
	gType = RECTANGLE ;
	this->center = Point(0,0) ; 
	topLeft = Point(-1, 1) ;boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(1,1) ;boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(1, -1) ;boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft = Point(-1, -1) ;boundingPoints[3] = new Point(bottomLeft) ;
	
	
}

std::vector<Point> Rectangle::getBoundingBox() const
{
	std::vector<Point> box ;
	box.push_back(topLeft) ;
	box.push_back(topRight) ;
	box.push_back(bottomRight) ;
	box.push_back(bottomLeft) ;
	
	return box ;
}


void Rectangle::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
		this->center += this->getPoint(i) ;
	
	if(this->size() != 0)
		this->center = this->center/this->size() ;
}

double  Rectangle::getRadius() const
{
	return sqrt(width()*width()*0.25 + height()*height()*0.25) ;
}

bool Rectangle::in(const Point & p) const 
{
	if(p.x < getCenter().x - 0.5*size_x)
		return false ;
	if(p.x  > getCenter().x + 0.5*size_x)
		return false ;
	if(p.y > getCenter().y + 0.5*size_y)
		return false ;
	if(p.y  < getCenter().y - 0.5*size_y)
		return false ;
	
	return true ;
}

double Rectangle::area() const
{
	return this->size_x*this->size_y ;
}

double Rectangle::width() const
{
	return size_x ;
}

double Rectangle::height() const
{
	return size_y ;
}


void Rectangle::project(Point * p) const
{
	
	Line l(getCenter(), *p) ;
	
	Segment A(topLeft, bottomLeft) ;
	if(l.intersects(A))
	{
		Point inter = l.intersection(A) ;
		p->x = inter.x ;
		p->y = inter.y ;
		return ;
	}
	
	Segment B(topLeft, topRight) ;
	
	if(l.intersects(B))
	{
		Point inter = l.intersection(B) ;
		p->x = inter.x ;
		p->y = inter.y ;
		return ;
	}
	
	Segment C(topRight, bottomRight) ;
	
	if(l.intersects(C))
	{
		Point inter = l.intersection(C) ;
		p->x = inter.x ;
		p->y = inter.y ;
		return ;
	}
	
	Segment D(bottomRight, bottomLeft) ;
	
	if(l.intersects(D))
	{
		Point inter = l.intersection(D) ;
		p->x = inter.x ;
		p->y = inter.y ;
		return ;
	}
}

void Rectangle::sampleBoundingSurface(size_t num_points)
{
	assert(num_points%4 == 0) ;
	double perimeter = 2*(size_x+size_y) ;
	
	double distanceBetweenPoints = perimeter/num_points ;
	
	this->numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongX = size_x/(this->numberOfPointsAlongX-1) ;
	
	this->numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongY = size_y/(this->numberOfPointsAlongY-1) ;
	
	num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;
	
	boundingPoints.resize(num_points) ;
	
	for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
	{
// 		double randx=((2.*random()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
		double randy=((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == 0 || i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[i] = new Point(center.x-0.5*size_x, center.y + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
			double randx=((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongX-1)
			randx = 0 ;
		boundingPoints[numberOfPointsAlongY+i-1] = new Point( center.x-0.5*size_x+i*distanceBetweenPointsAlongX+ randx, 
			getCenter().y-0.5*size_y);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		double randy=((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[numberOfPointsAlongX+numberOfPointsAlongY+i-2] = new Point(center.x+0.5*size_x,
			center.y-0.5*size_y+i*distanceBetweenPointsAlongY+ randy);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		double randx= ((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
		boundingPoints[2*numberOfPointsAlongY+numberOfPointsAlongX+i-3] = new Point(
			center.x + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
			center.y + 0.5*size_y) ;
		
	}
}

void Rectangle::sampleSurface(size_t num_points)
{
	if(this->size() == 4)
	{
		this->Rectangle::sampleBoundingSurface((size_t)round((double)num_points*std::max(size_x/size_y, size_y/size_x)/M_PI*2.)) ;
	}
	
	size_t nip = static_cast<size_t>((numberOfPointsAlongX-2)*(numberOfPointsAlongY-2)) ;
	
	inPoints.resize(nip) ;
		
	double distanceBetweenPointsAlongX = size_x/(this->numberOfPointsAlongX-1) ;
	double distanceBetweenPointsAlongY = size_y/(this->numberOfPointsAlongY-1) ;
	
	if(nip > 0)
	{
		for(size_t i = 0 ; i < this->numberOfPointsAlongX-2 ; i++)
		{
			for(size_t j = 0 ; j < this->numberOfPointsAlongY-2 ; j++)
			{	
				double randx=((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_x/numberOfPointsAlongX) ;
				double randy=((2.*random()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
				
				inPoints[i*(numberOfPointsAlongY-2)+j] = new Point(center.x - 0.5*size_x + (double)(i+1)*distanceBetweenPointsAlongX+randx,
					center.y - 0.5*size_y + (double)(j+1)*distanceBetweenPointsAlongY+ randy) ;
			}
		}
	}
}

Circle::Circle(double r, double originX, double originY)
{
	gType = CIRCLE ;
	this->center = Point(originX, originY) ;
	this->radius = r ;
	this->sqradius = r*r ;
}

Circle::Circle(double r, const Point *center)
{
	gType = CIRCLE ;
	this->center = Point(*center) ;
	this->radius = r ; 
	this->sqradius = r*r ;
}

Circle::Circle(double r, const Point center)
{
	gType = CIRCLE ;
	this->center = center ;
	this->radius = r ; 
	this->sqradius = r*r ;
}

void Circle::setRadius(double newr)
{
	this->radius = newr ;
	
	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		project(&getBoundingPoint(i)) ;
	}
	this->sqradius = newr*newr ;
}

void Circle::computeCenter()
{
}

std::vector<Point> Circle::getBoundingBox() const
{
	std::vector<Point> box ;
	box.push_back(getCenter()+Point(getRadius(), -getRadius())) ;
	box.push_back(getCenter()+Point(getRadius(), getRadius())) ;
	box.push_back(getCenter()+Point(-getRadius(), getRadius())) ;
	box.push_back(getCenter()+Point(-getRadius(), -getRadius())) ;
	
	return box ;
}

void Circle::project(Point * p) const
{	
	if(*p == center)
	{
		p->x = center.x + radius ;
		p->y = center.y ;
	}
	
	Segment seg(center, *p) ;
	
	Point proj = center + seg.vector()*(radius/seg.vector().norm()) ;
	
	p->x = proj.x ;
	p->y = proj.y ;
	
	return ;
	
}


void Circle::sampleBoundingSurface(size_t num_points)
{
	this->boundingPoints.resize(num_points) ;
	
	double angle = 2*M_PI/ (num_points) ;
	
	for (size_t i = 0 ; i< num_points ; i++)
	{
		double randa= 0;//((2.*(double)random()/(RAND_MAX+1.0))-1.)*0.15*(M_PI/num_points) ;
		boundingPoints[i] = new Point(getRadius()*cos((double)i*angle+randa) + getCenter().x, getRadius()*sin((double)i*angle+randa) + getCenter().y);
// 		std::cout << "x = " << boundingPoints[i]->x << ", y = " << boundingPoints[i]->y << std::endl ;
	}
}

void Circle::sampleSurface(size_t num_points)
{
	assert(!sampled) ;
	if(boundingPoints.size() == 0)
		this->sampleBoundingSurface(num_points) ;
	sampled = true ;
	size_t numberOfRings = static_cast<size_t>((double)num_points/(2. * M_PI )) ;
// 	if(numberOfRings > 0)
// 		numberOfRings-- ;
	assert(numberOfRings >= 0) ;
	double angle = 2.*M_PI/ (num_points) ;
	double offset = 0 ;
	
	//std::cout << "we have " << numberOfRings<< " rings" << std::endl ;
	size_t num_points_start = num_points ;
	
	std::vector<Point*> temp ;
	
	for (size_t i = 0 ; i< numberOfRings ; ++i)
	{
		double r = getRadius()*(1. - (double)(i + 1)/(numberOfRings+1)) ;
		//std::cout << "radius is " << r << std::endl ;
		
		for (size_t j = 0 ; j< num_points ; ++j)
		{
			double randa= 0 ; //((2.*(double)random()/(RAND_MAX+1.0))-1.)*0.2*(M_PI/num_points) ;
			double randr= 0 ; //(.2*r/(numberOfRings+1))*((double)random()/RAND_MAX*2.-1.0) ;
			temp.push_back(new Point((r+randr)*cos((double)(j+0.5*(i))*angle+randa+offset) + getCenter().x, (r+randr)*sin((double)(j+0.5*(i))*angle+randa) + getCenter().y));
		}
		
		num_points = (size_t)(/*std::max(*/(double)num_points_start*(r/getRadius())/*, (double)8)*/) ;
		
		angle = 2.*M_PI/ (num_points) ;
		
		offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;
	}
	
	inPoints.resize(temp.size() + 1) ;
	inPoints[0] = new Point(center) ;
	std::copy(temp.begin(), temp.end(),&inPoints[1]) ;
	//std::cout << "we have " << num_points << " sample points" << std::endl ;
	
}

bool Circle::in(const Point & v) const 
{
	double val[2] = {getCenter().x-v.x, getCenter().y-v.y} ;
	return val[0]*val[0] + val[1]*val[1] < sqradius ;
}

double Circle::getRadius() const
{
	return radius ;
}

double Circle::area() const
{
	return M_PI*sqradius ;
}


LayeredCircle::LayeredCircle(std::vector<double> radii,double originX,double originY): Circle(radii[0], originX, originY)
{
	this->gType = LAYERED_CIRCLE ;
	std::sort(radii.begin(), radii.end()) ;
	Circle::setRadius(*radii.rbegin()) ;
	radiuses = radii ;
	center = Point(originX, originY) ;
}


LayeredCircle::LayeredCircle(std::vector<double> radii, const Point c) : Circle(radii[0], c)
{
	this->gType = LAYERED_CIRCLE ;
	std::sort(radii.begin(), radii.end()) ;
	Circle::setRadius(*radii.rbegin()) ;
	radiuses = radii ;
	center = c ;
}

LayeredCircle::LayeredCircle(double r, const Point center) : Circle(r, center)
{
	radiuses.push_back(r) ;
}

void LayeredCircle::sampleSurface(size_t num_points)
{
	assert(!sampled) ;
	if(boundingPoints.size() == 0)
		this->sampleBoundingSurface(num_points) ;
	sampled = true ;
	size_t numberOfRings = static_cast<size_t>((double)num_points/(2. * M_PI )) ;
	assert(numberOfRings >= 0) ;
	double angle = 2.*M_PI/ (num_points) ;
	double offset = 0 ;
	
	size_t num_points_start = num_points ;
	std::vector<Point*> temp ;
	std::vector<double> originalSamplingRadiuses ;
	double meanDelta = getRadius()/(numberOfRings+1);
	
	for (size_t i = 0 ; i< numberOfRings ; ++i)
	{
		originalSamplingRadiuses.push_back(getRadius()*(1. - (double)(i + 1)/(numberOfRings+1))) ;
	}
	
	originalSamplingRadiuses.insert(originalSamplingRadiuses.end(), radiuses.begin(),  radiuses.end()-1) ;
	
	std::sort(originalSamplingRadiuses.begin(),originalSamplingRadiuses.end()) ;
	
	for(size_t i = 0 ;  i < originalSamplingRadiuses.size() ; i++)
	{
		std::cout << originalSamplingRadiuses[i] << std::endl ;
	}
	bool deleted = true ;
// 	while(deleted)
// 	{
// 		deleted = false ;
// 		for(size_t i = 0 ; i < originalSamplingRadiuses.size() ;++i)
// 		{
// 			if(originalSamplingRadiuses[i+1] - originalSamplingRadiuses[i] < .2*meanDelta)
// 			{
// 				bool foundI = false ;
// 				bool foundIp1 = false ;
// 				for(size_t j = 0 ; j < radiuses.size() ; j++)
// 				{
// 					if(std::abs(radiuses[j]-originalSamplingRadiuses[i]) < .2*meanDelta)
// 						foundI = true ;
// 					if(std::abs(radiuses[j]-originalSamplingRadiuses[i+1]) < .2*meanDelta)
// 						foundIp1 = true ;
// 				}
// 	
// 				if(!foundI && foundIp1 )
// 				{
// 					originalSamplingRadiuses.erase(i+1) ;
// 					deleted = true ;
// 					break ;
// 				}
// 				else if (foundI && !foundIp1)
// 				{
// 					originalSamplingRadiuses.erase(i) ;
// 					deleted = true ;
// 					break ;
// 				}
// 			}
// 		}
// 	}
	
	for(size_t i = 0 ;  i < originalSamplingRadiuses.size() ; i++)
	{
		std::cout << originalSamplingRadiuses[i] << std::endl ;
	}
	
	for (size_t i = originalSamplingRadiuses.size()-1 ; i != 0 ; --i)
	{
		double r = originalSamplingRadiuses[i] ;
		for (size_t j = 0 ; j< num_points ; ++j)
		{
			double randa= 0 ; //((2.*(double)random()/(RAND_MAX+1.0))-1.)*0.2*(M_PI/num_points) ;
			double randr= 0 ; //(.2*r/(numberOfRings+1))*((double)random()/RAND_MAX*2.-1.0) ;
			temp.push_back(new Point((r+randr)*cos((double)(j+0.5*(i))*angle+randa+offset) + getCenter().x, (r+randr)*sin((double)(j+0.5*(i))*angle+randa) + getCenter().y));
		}
		
		num_points = (size_t)((double)num_points_start*(r/getRadius())) ;
		
		angle = 2.*M_PI/ (num_points) ;
		
		offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;
	}
	
	inPoints.resize(temp.size() + 1) ;
	inPoints[0] = new Point(center) ;
	std::copy(temp.begin(), temp.end(),&inPoints[1]) ;
}

void LayeredCircle::setRadius(double newr)
{
	double ratio = newr/(*radiuses.rbegin()) ;
	radiuses.pop_back() ;
	radiuses.push_back(newr) ;
	std::sort(radiuses.begin(), radiuses.end()) ;
	
	
	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		getBoundingPoint(i).x = (getBoundingPoint(i).x - center.x)*ratio + center.x ;
		getBoundingPoint(i).y = (getBoundingPoint(i).y - center.y)*ratio + center.y ;
	}
	
	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		getInPoint(i).x = (getInPoint(i).x - center.x)*ratio + center.x ;
		getInPoint(i).y = (getInPoint(i).y - center.y)*ratio + center.y ;
	}
	
	this->sqradius = newr*newr ;
}

void LayeredCircle::addRadius(double newr)
{
	radiuses.push_back(newr) ;
	std::sort(radiuses.begin(), radiuses.end()) ;
	this->sqradius = (*radiuses.rbegin())*(*radiuses.rbegin()) ;
}


SegmentedLine::SegmentedLine(const std::valarray<Point *> & points) : NonConvexGeometry(0)
{
	gType = SEGMENTED_LINE ;
	
	this->center = Point() ;
	this->boundingPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
		boundingPoints[i] = points[i] ;
}

void SegmentedLine::computeCenter()
{
}

double SegmentedLine::getRadius() const
{
	//! \todo make it do something
	return 0 ;
}

Point * SegmentedLine::getHead() const 
{
	return  boundingPoints[0];
}

Point * SegmentedLine::getTail() const 
{
	return boundingPoints[this->boundingPoints.size()-1] ;
}

void SegmentedLine::sampleBoundingSurface(size_t num_points) 
{ 
}

void SegmentedLine::project(Point *p) const
{

	std::map<double, Point> projections ;
	
	for(size_t i = 1 ; i < boundingPoints.size()-1 ; i++)
	{
		Point proj = Segment(getBoundingPoint(i), getBoundingPoint(i+1)).project(*p) ;
		projections[squareDist2D(*p, proj)] = proj ;
	}
	
	p->x = projections.begin()->second.x ;
	p->y = projections.begin()->second.y ;
}

void SegmentedLine::sampleSurface(size_t num_points) 
{ 
	
} 

bool SegmentedLine::in(const Point & v) const
{
	for (size_t i = 0 ; i < boundingPoints.size() ;i++)
		if(*boundingPoints[i] == v)
			return true ;
	
	return false ;
}

