// Authors: 	Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//		Alain Giorla <alain.giorla@epfl.ch>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_2D.h"
#include "../polynomial/vm_base.h"

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
	
	
	if((p0.z == p1.z)  && (p0.z == p2.z) &&  (p0.z == 0))
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
			if(!isTrigoOriented())
			{
				std::cout << "arrrgh !" << std::endl ;
				boundingPoints[0]->print() ;
				boundingPoints[1]->print() ;
				boundingPoints[2]->print() ;
			}
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if((p0.z == p1.z)  && (p0.z == p2.z) &&  (p0.z == 0))
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = sqrt((squareDist2D(p1, circumCenter)+ squareDist2D(p0, circumCenter))/2.);
	sqradius = radius*radius ;
	
}

Triangle::Triangle(XMLTree * xml)
{
	gType = TRIANGLE ;
	assert(this->size() == 3) ;
	
	Point p0(0,1) ;
	Point p1(0,0) ;
	Point p2(1,0) ;

	boundingPoints[0] = new Point(0, 1) ;
	boundingPoints[1] = new Point(0, 0) ;
	boundingPoints[2] = new Point(1, 0) ;

	if(xml->match("triangle"))
	{
		p0 = Point(xml->getChild(0)) ;
		p1 = Point(xml->getChild(1)) ;
		p2 = Point(xml->getChild(2)) ;

		boundingPoints[0] = new Point(p0) ;
		boundingPoints[1] = new Point(p1) ;
		boundingPoints[2] = new Point(p2) ;
	}
	
	computeCenter() ;
	computeCircumCenter() ;

	if((p0.z == p1.z)  && (p0.z == p2.z) &&  (p0.z == 0))
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = sqrt((squareDist2D(p1, circumCenter)+ squareDist2D(p0, circumCenter))/2.);
	sqradius = radius*radius ;
	
}


XMLTree * Triangle::toXML()
{
	XMLTree * tri = new XMLTree("triangle") ;
	tri->addChild(boundingPoints[0]->toXML()) ;
	tri->addChild(boundingPoints[1]->toXML()) ;
	tri->addChild(boundingPoints[2]->toXML()) ;
	return tri ;
}


OrientedRectangle::OrientedRectangle() : ConvexGeometry(4)
{
	gType =PARALLELOGRAMME;
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(-1,1);
	boundingPoints[1] = new Point(-1,-1);
	boundingPoints[2] = new Point(1, -1);
	boundingPoints[3] = new Point(1, 1);
	computeCenter() ;
}

OrientedRectangle::OrientedRectangle( const Point & p0,  const Point & p1,  const Point & p2,  const Point & p3): ConvexGeometry(4)
{
	gType = PARALLELOGRAMME ;
	
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(p0) ;
	boundingPoints[1] = new Point(p1) ;
	boundingPoints[2] = new Point(p2) ;
	boundingPoints[3] = new Point(p3) ;
	
	if((p0.z == p1.z)  && (p0.z == p2.z) &&  (p0.z == 0))
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if((p0.z == p1.z)  && (p0.z == p2.z) &&  (p0.z == 0))
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = (squareDist2D(p1, circumCenter) + squareDist2D(p0, circumCenter) + squareDist2D(p2, circumCenter)+squareDist2D(p3, circumCenter))/4.;
	
}

OrientedRectangle::OrientedRectangle( const Point *p0,  const Point *p1,  const Point *p2,  const Point *p3): ConvexGeometry(4)
{
	gType = PARALLELOGRAMME ;
	
	assert(this->size() == 4);
	
	boundingPoints[0] = new Point(*p0) ;
	boundingPoints[1] = new Point(*p1) ;
	boundingPoints[2] = new Point(*p2) ;
	boundingPoints[3] = new Point(*p3) ;
	
	if((p0->z == p1->z)  && (p0->z == p2->z) &&  (p0->z == 0))
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}
	
	computeCircumCenter() ;
	computeCenter() ;
	
	if((p0->z == p1->z)  && (p0->z == p2->z) &&  (p0->z == 0))
	{
		if(!this->in(this->getCenter()))
		{
			assert(false) ;
		}
	}
	
	radius = (squareDist2D(*p1, circumCenter) + squareDist2D(*p0, circumCenter) + squareDist2D(*p2, circumCenter)+squareDist2D(*p3, circumCenter))/4.;
	
}

XMLTree * OrientedRectangle::toXML()
{
	XMLTree * rect = new XMLTree("oriented rectangle") ;
	rect->addChild(boundingPoints[0]->toXML()) ;
	rect->addChild(boundingPoints[1]->toXML()) ;
	rect->addChild(boundingPoints[2]->toXML()) ;
	rect->addChild(boundingPoints[3]->toXML()) ;
	return rect ;
}


void OrientedRectangle::computeCenter()
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		this->center += this->getPoint(i) ;
	}
	
	this->center = this->center/this->size() ;
}

double OrientedRectangle::getRadius() const
{
	return sqrt(radius) ;
}

Point OrientedRectangle::getCircumCenter() const
{
	return this->circumCenter ;
}

void OrientedRectangle::computeCircumCenter()
{	
	if (fabs(boundingPoints[1]->y-boundingPoints[0]->y) < 20*POINT_TOLERANCE) 
	{
		double m2 = - (boundingPoints[2]->x-boundingPoints[1]->x) / (boundingPoints[2]->y-boundingPoints[1]->y);
		double mx2 = (boundingPoints[1]->x + boundingPoints[2]->x) / 2.0;
		double my2 = (boundingPoints[1]->y + boundingPoints[2]->y) / 2.0;
		double xc = (boundingPoints[1]->x + boundingPoints[0]->x) / 2.0;
		double yc = fma(m2, (xc - mx2), my2);
		
		circumCenter = Point(xc, yc) ;
	} 
	else if (fabs(boundingPoints[2]->y-boundingPoints[1]->y) < 20*POINT_TOLERANCE) 
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

bool OrientedRectangle::inCircumCircle(const Point p) const
{
	return  squareDist2D(circumCenter, p) < radius  ;
}

bool OrientedRectangle::inCircumCircle(const Point *p) const
{
	return  squareDist2D(circumCenter, (*p)) < radius  ;
}

bool OrientedRectangle::is1D() const
{ 
	return false ; 
} 

double OrientedRectangle::area() const
{
	assert(this->boundingPoints.size() == 4) ;
	
	Segment s0((*boundingPoints[0]), (*boundingPoints[1])) ;
	Segment s1((*boundingPoints[0]), (*boundingPoints[2])) ;

	return ((s0.vector())^(s1.vector())).norm() ;
}


void OrientedRectangle::project(Point * p) const
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


bool OrientedRectangle::in(const Point &p) const
{
	
	bool in = false ;
	
	for (size_t i = 0, j  =  boundingPoints.size()-1; i <  boundingPoints.size(); j = i++)
	{
		if (
			(((boundingPoints[i]->y <= p.y ) 
				&& (p.y<boundingPoints[j]->y)) 
				|| ((boundingPoints[j]->y <= p.y) 
				&& (p.y<boundingPoints[i]->y))) 
				&& (p.x < (boundingPoints[j]->x - boundingPoints[i]->x) * (p.y - boundingPoints[i]->y) / (boundingPoints[j]->y - boundingPoints[i]->y) + boundingPoints[i]->x))
			in = !in;
	}
	
	return in ;
	
}
std::vector<Point> OrientedRectangle::getSamplingBoundingPoints(size_t num_points) const 
{
	std::vector<Point> ret ;
	num_points = num_points + num_points%4 ;
	
	Point  v0(*boundingPoints[0]) ;
	Point  v1(*boundingPoints[1]) ;
	Point  v2(*boundingPoints[2]) ;
	Point  v3(*boundingPoints[3]) ;
	
	ret.push_back(v0) ;
	
	for(size_t i = 1 ; i < num_points/4; i++)
	{
		ret.push_back(v0*4.*i/num_points + v1*(1.-4.*i/num_points)) ;
	}
	
	ret.push_back(v1) ;
	
	for(size_t i = num_points/4+1 ; i < 2*num_points/4 ; i++)
	{
		ret.push_back(v1*4.*(i-num_points/4)/num_points + v2*(1.-4.*(i-num_points/4)/num_points)) ;
	}
	
	ret.push_back(v2) ;
	
	for(size_t i = 2*num_points/4+1 ; i < num_points ; i++)
	{
		ret.push_back(v2*4.*(i-2*num_points/4)/num_points + v3*(1.-4.*(i-2*num_points/4)/num_points)) ;
	}
	
	ret.push_back(v3) ;
	
	for(size_t i = 3*num_points/4+1 ; i < num_points ; i++)
	{
		ret.push_back(v3*4.*(i-3*num_points/4)/num_points + v0*(1.-4.*(i-3*num_points/4)/num_points)) ;
	}
	
	return ret ;
}


void OrientedRectangle::sampleBoundingSurface(size_t num_points)
{
	num_points = num_points + num_points%4 ;
	
	Point * v0(boundingPoints[0]) ;
	Point * v1(boundingPoints[1]) ;
	Point * v2(boundingPoints[2]) ;
	Point * v3(boundingPoints[3]) ;
	this->boundingPoints.resize(num_points) ;
	
	boundingPoints[0] = v0 ;
	
	for(size_t i = 1 ; i < num_points/4; i++)
	{
		boundingPoints[i] = new Point(*v0*4.*i/num_points + *v1*(1.-4.*i/num_points)) ;
	}
	
	boundingPoints[num_points/4] = v1 ;
	
	for(size_t i = num_points/4+1 ; i < 2*num_points/4 ; i++)
	{
		boundingPoints[i] = new Point(*v1*4.*(i-num_points/4)/num_points + *v2*(1.-4.*(i-num_points/4)/num_points)) ;
	}
	
	boundingPoints[2*num_points/4] = v2 ;
	
	for(size_t i = 2*num_points/4+1 ; i < num_points ; i++)
	{
		boundingPoints[i] = new Point(*v2*4.*(i-2*num_points/4)/num_points + *v3*(1.-4.*(i-2*num_points/4)/num_points)) ;
	}
	
	boundingPoints[3*num_points/4] = v3 ;
	
	for(size_t i = 3*num_points/4+1 ; i < num_points ; i++)
	{
		boundingPoints[i] = new Point(*v3*4.*(i-3*num_points/4)/num_points + *v0*(1.-4.*(i-3*num_points/4)/num_points)) ;
	}
	
}

void OrientedRectangle::sampleSurface(size_t num_points)
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

	
	if((p0->z == p1->z)  && (p0->z == p2->z) &&  (p0->z == 0))
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
	
	if((p0->z == p1->z)  && (p0->z == p2->z) &&  (p0->z == 0))
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

const Point & Triangle::getCircumCenter() const
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
	
	if (std::abs(getBoundingPoint(1).y-getBoundingPoint(0).y) < 100.*POINT_TOLERANCE) 
	{
		double m2  =  (getBoundingPoint(1).x-getBoundingPoint(2).x ) / (getBoundingPoint(2).y-getBoundingPoint(1).y);
		double mx2 = (getBoundingPoint(1).x + getBoundingPoint(2).x) ;
		double my2 = (getBoundingPoint(1).y + getBoundingPoint(2).y) ;
		double xc  = (getBoundingPoint(1).x + getBoundingPoint(0).x) ;
		double yc  = m2 * (xc - mx2) + my2;
		
		circumCenter.set(xc/2., yc/2.) ;
	} 
	else if (std::abs(getBoundingPoint(2).y-getBoundingPoint(1).y) < 100.*POINT_TOLERANCE) 
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
	if(p.x > circumCenter.x+1.01*radius)
		return false ;
	if(p.x < circumCenter.x-1.01*radius)
		return false ;
	if(p.y > circumCenter.y+1.01*radius)
		return false ;
	if(p.y < circumCenter.y-1.01*radius)
		return false ;

	if(squareDist2D(circumCenter, p) < .99*sqradius)
		return true ;
	
	double delta = POINT_TOLERANCE ;
	Point a(p) ; a.x += delta ; a.y += delta ; 
	Point c(p) ; c.x += delta ; c.y -= delta ; 
	Point e(p) ; e.x -= delta ; e.y += delta ; 
	Point g(p) ; g.x -= delta ; g.y -= delta ; 

	return  squareDist2D(circumCenter, a) < sqradius 
		&&  squareDist2D(circumCenter, c) < sqradius
		&&  squareDist2D(circumCenter, e) < sqradius
		&&  squareDist2D(circumCenter, g) < sqradius;
	double x = circumCenter.x -p.x ;
	double y = circumCenter.y -p.y ;
	return  fma(x, x, y*y)< sqradius*(1. - 100.*POINT_TOLERANCE)  ;
}

bool Triangle::inCircumCircle(const Point *p) const
{
	if(p->x > circumCenter.x+1.01*radius)
		return false ;
	if(p->x < circumCenter.x-1.01*radius)
		return false ;
	if(p->y > circumCenter.y+1.01*radius)
		return false ;
	if(p->y < circumCenter.y-1.01*radius)
		return false ;

	if(squareDist2D(circumCenter, *p) < .99*sqradius)
		return true ;
	
	double delta = POINT_TOLERANCE ;
	Point a(*p) ; a.x += delta ; a.y += delta ; 
	Point c(*p) ; c.x += delta ; c.y -= delta ; 
	Point e(*p) ; e.x -= delta ; e.y += delta ; 
	Point g(*p) ; g.x -= delta ; g.y -= delta ; 

	return  squareDist2D(circumCenter, a) < sqradius 
		&&  squareDist2D(circumCenter, c) < sqradius
		&&  squareDist2D(circumCenter, e) < sqradius
		&&  squareDist2D(circumCenter, g) < sqradius;
	double x = circumCenter.x -p->x ;
	double y = circumCenter.y -p->y ;
	return  fma(x, x, y*y) < sqradius*(1. - 100.*POINT_TOLERANCE)  ;
}


double Triangle::area() const
{
// 	assert(this->boundingPoints.size() == 3) ;
	int pointsInTimePlane = this->boundingPoints.size() ;
	if(getBoundingPoint(0).t != 0)
	{
		pointsInTimePlane = 0 ;
		double init = getBoundingPoint(0).t ;
		int counter = 0 ;
		while(std::abs(getBoundingPoint(counter++).t-init) < 1e-8)
			pointsInTimePlane++ ;
	}
	Segment s0(getBoundingPoint(0), getBoundingPoint(pointsInTimePlane/3)) ;
	Segment s1(getBoundingPoint(0), getBoundingPoint(pointsInTimePlane/3 * 2)) ;

	return 0.5*std::abs((s0.vector()^s1.vector()).z) ;
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
	bool isAPoint = false ;
	for (int i = 0; i <  boundingPoints.size(); i++)
	{
		if(p == *boundingPoints[i])
		{
			isAPoint = true ;
			break ;
		}
	}
	
	Point proj(p) ; project(&proj) ;
	bool isOnSurface = dist(p, proj) < POINT_TOLERANCE ;
	
	Segment s(p, getCenter()) ;
	return !s.intersects(this) || isAPoint || isOnSurface;
		
	bool in = false ;
	
	for (int i = 0, j  =  boundingPoints.size()-1; i <  boundingPoints.size(); j = i++)
	{
		if( std::abs(boundingPoints[j]->y - boundingPoints[i]->y) > 2.*POINT_TOLERANCE)
		{
			if (
				(((boundingPoints[i]->y < p.y + 2.*POINT_TOLERANCE) 
					&& (p.y-2.*POINT_TOLERANCE < boundingPoints[j]->y)) 
					|| ((boundingPoints[j]->y < p.y+2.*POINT_TOLERANCE) 
					&& (p.y < boundingPoints[i]->y+2.*POINT_TOLERANCE))) 
					&& (p.x < (boundingPoints[j]->x - boundingPoints[i]->x) 
					 * (p.y - boundingPoints[i]->y) 
					 / (boundingPoints[j]->y - boundingPoints[i]->y) 
					 + boundingPoints[i]->x + 2.*POINT_TOLERANCE))
				in = !in;
		}
	}
	
	return in ;
	
}

std::vector<Point> Triangle::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;
	assert(num_points%3 == 0) ;
	
	Point v0 = *boundingPoints[0] ;
	Point v1 = *boundingPoints[1] ;
	Point v2 = *boundingPoints[2] ;
		
	ret.push_back(v0) ;
	
	for(size_t i = 1 ; i < num_points/3 ; i++)
	{
		ret.push_back(v0*3.*i/num_points + v1*(1.-3.*i/num_points)) ;
	}
	
	ret.push_back(v1) ;
	
	for(size_t i = num_points/3+1 ; i < 2*num_points/3 ; i++)
	{
		ret.push_back(v1*3.*(i-num_points/3)/num_points + v2*(1.-3.*(i-num_points/3)/num_points)) ;
	}
	
	ret.push_back(v2) ;
	
	for(size_t i = 2*num_points/3+1 ; i < num_points ; i++)
	{
		ret.push_back(v2*3.*(i-2*num_points/3)/num_points + v0*(1.-3.*(i-2*num_points/3)/num_points)) ;
	}

	return ret ;
}

void Triangle::sampleBoundingSurface(size_t num_points)
{
	num_points -= num_points%3 ;
	if(num_points == 0)
		num_points = 3 ;
	Point * v0 = &getBoundingPoint(0) ;
	Point * v1 = &getBoundingPoint(1) ;
	Point * v2 = &getBoundingPoint(2) ;
	
	getBoundingPoints().resize(num_points) ;
	
	boundingPoints[0] = v0 ;
	
	for(size_t i = 1 ; i < num_points/3 ; i++)
	{
		getBoundingPoints()[i] = new Point((*v1)*3.*i/num_points + (*v0)*(1.-3.*i/num_points)) ;
	}
	
	getBoundingPoints()[num_points/3] = v1 ;
	
	for(size_t i = num_points/3+1 ; i < 2*num_points/3 ; i++)
	{
		getBoundingPoints()[i] = new Point((*v2)*3.*(i-num_points/3)/num_points + (*v1)*(1.-3.*(i-num_points/3)/num_points)) ;
	}
	
	getBoundingPoints()[2*num_points/3] = v2 ;
	
	for(size_t i = 2*num_points/3+1 ; i < num_points ; i++)
	{
		getBoundingPoints()[i] = new Point((*v0)*3.*(i-2*num_points/3)/num_points + (*v2)*(1.-3.*(i-2*num_points/3)/num_points)) ;
	}
	
}

void Triangle::sampleSurface(size_t num_points)
{
	num_points += 3-num_points%3 ;
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
			double xrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
			double yrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
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

Rectangle::Rectangle(double x, double y, const Point &center) :  ConvexGeometry(4), size_y(y), size_x(x)
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

Rectangle::Rectangle(XMLTree * xml) : ConvexGeometry(4)
{
	gType = RECTANGLE ;
	this->center = Point(0,0) ; 
	topLeft = Point(-1, 1) ;
	topRight = Point(1, 1) ;
	bottomRight = Point(1, -1) ;
	bottomLeft = Point(-1, -1) ;

	if(xml->match("rectangle")) ;
	{
		this->center = Point(xml->getChild(0)->getChild(0)) ;
		size_x = xml->getChild(1)->buildDouble().second ;
		size_y = xml->getChild(2)->buildDouble().second ;
		topLeft = Point(center.x-0.5*size_x, center.y+0.5*size_y) ;
		topRight = Point(center.x+0.5*size_x, center.y+0.5*size_y) ;
		bottomRight = Point(center.x+0.5*size_x, center.y-0.5*size_y) ;
		bottomLeft =  Point(center.x-0.5*size_x, center.y-0.5*size_y) ;
	}

	boundingPoints[0] = new Point(topLeft) ;
	boundingPoints[1] = new Point(topRight) ;
	boundingPoints[2] = new Point(bottomRight) ;
	boundingPoints[3] = new Point(bottomLeft) ;

}


XMLTree * Rectangle::toXML()
{
	XMLTree * rect = new XMLTree("rectangle") ;
	XMLTree * c = new XMLTree("center") ;
	c->addChild(this->getCenter().toXML()) ;
	rect->addChild(c) ;
	rect->addChild(new XMLTree("x",size_x)) ;
	rect->addChild(new XMLTree("y",size_y)) ;
	return rect ;
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
	if(p.x < getCenter().x - 0.5*width())
		return false ;
	if(p.x  > getCenter().x + 0.5*width())
		return false ;
	if(p.y > getCenter().y + 0.5*height())
		return false ;
	if(p.y  < getCenter().y - 0.5*height())
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
	Segment A(topLeft, bottomLeft) ;
	
	Segment B(topLeft, topRight) ;
	
	Segment C(topRight, bottomRight) ;

	Segment D(bottomRight, bottomLeft) ;

	std::map<double, Point> tries ;
	Point p0 = A.project(*p) ;
	Point p1 = B.project(*p) ;
	Point p2 = C.project(*p) ;
	Point p3 = D.project(*p) ;
	tries[dist(p, &p0)] = p0 ;
	tries[dist(p, &p1)] = p1 ;
	tries[dist(p, &p2)] = p2 ;
	tries[dist(p, &p3)] = p3 ;
	*p = tries.begin()->second ;
}

std::vector<Point> Rectangle::getSamplingBoundingPoints(size_t num_points) const
{
	double perimeter = 2*(size_x+size_y) ;
	
	double distanceBetweenPoints = perimeter/num_points ;
	std::vector<Point> ret ;
	
	double numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;
	
	double numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y/distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;
	
	num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;
	
	for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
	{
// 		double randx=((2.*rand()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == 0 || i == numberOfPointsAlongY-1)
			randy = 0 ;
		ret.push_back(Point(center.x-0.5*size_x, center.y + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy)) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
		double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongX-1)
			randx = 0 ;
		ret.push_back(Point( center.x-0.5*size_x+i*distanceBetweenPointsAlongX+ randx, 
		                     getCenter().y-0.5*size_y));
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongY-1)
			randy = 0 ;
		ret.push_back(Point(center.x+0.5*size_x,
		                    center.y-0.5*size_y+i*distanceBetweenPointsAlongY+ randy));
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
		ret.push_back(Point(
			center.x + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
		                     center.y + 0.5*size_y)) ;
	}
	
	return ret ;
}

void Rectangle::sampleBoundingSurface(size_t num_points)
{
  //	assert(num_points%4 == 0) ;
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
// 		double randx=((2.*rand()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == 0 || i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[i] = new Point(center.x-0.5*size_x, center.y + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
			double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongX-1)
			randx = 0 ;
		boundingPoints[numberOfPointsAlongY+i-1] = new Point( center.x-0.5*size_x+i*distanceBetweenPointsAlongX+ randx, 
			getCenter().y-0.5*size_y);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[numberOfPointsAlongX+numberOfPointsAlongY+i-2] = new Point(center.x+0.5*size_x,
			center.y-0.5*size_y+i*distanceBetweenPointsAlongY+ randy);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
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
				double randx= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.1*(size_x/numberOfPointsAlongX) ;
				double randy= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.1*(size_y/numberOfPointsAlongY) ;
				
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

Circle::Circle(XMLTree * xml)
{
	gType = CIRCLE ;
	this->center = Point(0,0) ;
	this->radius = 1 ; 

	if(xml->match("circle"))
	{
		this->center = Point(xml->getChild(0)->getChild(0)) ;
		this->radius = xml->getChild(1)->buildDouble().second ;
	}

	this->sqradius = this->radius * this->radius ;
}


XMLTree * Circle::toXML()
{
	XMLTree * cir = new XMLTree("circle") ;
	XMLTree * c = new XMLTree("center") ;
	c->addChild(this->getCenter().toXML()) ;
	cir->addChild(c) ;
	cir->addChild(new XMLTree("radius",radius)) ;
	return cir ;
}


void Circle::setRadius(double newr)
{
	double ratio = newr/(radius) ;
	
	
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
	
	this->radius = newr ;
	this->sqradius = newr*newr ;
}

void Circle::computeCenter()
{
}

std::vector<Point> Circle::getBoundingBox() const
{
	std::vector<Point> box ;
	box.push_back(getCenter()+Point(-getRadius(), getRadius())) ;
	box.push_back(getCenter()+Point(getRadius(), getRadius())) ;
	box.push_back(getCenter()+Point(getRadius(), -getRadius())) ;
	box.push_back(getCenter()+Point(-getRadius(), -getRadius())) ;
	
	return box ;
}

void Circle::project(Point * p) const
{
	if(squareDist2D(p, &center ) < POINT_TOLERANCE*POINT_TOLERANCE)
	{
		p->x +=getRadius() ;
		return ;
	}
	
	Line l(*p, *p-getCenter()) ;
	
	std::vector<Point> inter = l.intersection(this) ;
	if(inter.empty() || inter.size() == 1)
	{
		p->print() ;
		getCenter().print() ;
	}
	if(squareDist3D(inter[0], *p) < squareDist3D(inter[1], *p))
	{
		*p = inter[0] ;
		return ;
	}
	*p = inter[1] ;
}

std::vector<Point> Circle::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;
	
	double angle = 2.*M_PI/ (num_points) ;
	
	for (size_t i = 0 ; i< num_points ; i++)
	{
		ret.push_back(Point(getRadius()*cos((double)i*angle) + getCenter().x, getRadius()*sin((double)i*angle) + getCenter().y));
	}
	
	return ret ;
}

std::vector<Point> Circle::getSamplingBoundingPointsOnArc(size_t num_points, const Point & start, const Point & finish) const
{
	std::vector<Point> ret ;
	Point init(start) ;
	project(&init) ;
	init -= getCenter() ;
	Point fin(finish) ;
	project(&fin) ;
	fin -= getCenter() ;
	double angle = (init.angle() -fin.angle())/num_points ;
	
	
	for (double i = 0 ; i< num_points ; i++)
	{
		Point newPoint(init.x*cos(i*angle)+init.y*sin(i*angle), -init.x*sin(i*angle)+init.y*cos(i*angle)) ;
		newPoint+= getCenter() ;
		ret.push_back(newPoint);
	}
	return ret ;
}

void Circle::sampleBoundingSurface(size_t num_points)
{
	getBoundingPoints().resize(num_points) ;
	
	double angle = 2.*M_PI/ (num_points) ;
	
	for (size_t i = 0 ; i < num_points ; i++)
	{		
		double randa= 0;//((2.*(double)rand()/(RAND_MAX+1.0))-1.)*0.15*(M_PI/num_points) ;
		getBoundingPoints()[i] = new Point(getRadius()*cos((double)i*angle) + getCenter().x, getRadius()*sin((double)i*angle+randa) + getCenter().y);
// 		std::cout << "x = " << boundingPoints[i]->x() << ", y = " << boundingPoints[i]->y << std::endl ;
	}
}

void Circle::sampleSurface(size_t num_points)
{
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
			double randa= 0 ; //((2.*(double)rand()/(RAND_MAX+1.0))-1.)*0.2*(M_PI/num_points) ;
			double randr= 0 ; //(.2*r/(numberOfRings+1))*((double)rand()/RAND_MAX*2.-1.0) ;
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
	if(v.x < center.x-getRadius())
		return false ;
	if(v.x > center.x+getRadius())
		return false ;
	if(v.y < center.y-getRadius())
		return false ;
	if(v.y > center.y+getRadius())
		return false ;
	
	return squareDist2D(v, getCenter()) < getSquareRadius() ;
}

double Circle::getRadius() const
{
	return radius ;
}

double Circle::getSquareRadius() const
{
	return sqradius ;
}

double Circle::area() const
{
	return M_PI*sqradius ;
}


const std::vector<double> & LayeredCircle::getRadii() const
{
	return radiuses ;
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

XMLTree * LayeredCircle::toXML()
{
	XMLTree * circle = new XMLTree("layered circle") ;
	XMLTree * c = new XMLTree("center") ;
	c->addChild(this->getCenter().toXML()) ;
	circle->addChild(c) ;
	circle->addChild(new XMLTree("radius",radiuses)) ;
	return circle ;
}



void LayeredCircle::sampleSurface(size_t num_points)
{

	std::vector<double> samplingRadiuses = radiuses;
	bool to_add = true ;
	std::vector<Point*> temp ;
	size_t numberOfRings = static_cast<size_t>((double)num_points/(2. * M_PI )) ;
	double meanDelta = getRadius()/(numberOfRings+1);

	while(to_add)
	{
		to_add = false ;
		std::vector<double> newRadii ;
		if(samplingRadiuses[0] > meanDelta*.25)
		{
			newRadii.push_back(samplingRadiuses[0]*.5) ;
			newRadii.push_back(samplingRadiuses[0]) ;
			to_add = true ;
// 			std::cout << "a-added " << samplingRadiuses[0]*.5  << std::endl ;
		}
		for(size_t i = 1 ; i < samplingRadiuses.size() ; i++)
		{
			newRadii.push_back(samplingRadiuses[i]) ;
			if(samplingRadiuses[i]-samplingRadiuses[i-1] > meanDelta*.5)
			{
				newRadii.push_back(samplingRadiuses[i-1]*.5+samplingRadiuses[i]*.5) ;
				to_add = true ;
// 				std::cout << "b-added " << samplingRadiuses[i-1]*.5+samplingRadiuses[i]*.5  << std::endl ;
			}
		}
		std::sort(newRadii.begin() , newRadii.end()) ;
// 		for(size_t i = 0 ; i < newRadii.size() ; i++)
// 		{
// 			std::cout << newRadii[i] << "  "<< std::flush ;
// 		}
// 		std::cout <<  std::endl ;
		samplingRadiuses = newRadii ;
	}

	
	num_points = std::max(num_points * 2.* getRadius()/meanDelta, 12.);
	size_t num_points_start = num_points ;
	if(boundingPoints.size() == 0)
		this->sampleBoundingSurface(num_points) ;
	sampled = true ;

	double angle = 2.*M_PI/ (num_points) ;
	double offset = 0 ;

	for (int i = samplingRadiuses.size()-2 ; i >= 0 ; --i)
	{
		double r = samplingRadiuses[i] ;
		for (size_t j = 0 ; j< num_points ; ++j)
		{
			double randa= 0 ; //((2.*(double)rand()/(RAND_MAX+1.0))-1.)*0.2*(M_PI/num_points) ;
			double randr= 0 ; //(.2*r/(numberOfRings+1))*((double)rand()/RAND_MAX*2.-1.0) ;
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
	this->radius = newr ;
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
	
	if(points[points.size()/2])
		this->center = *points[points.size()/2] ;
	else
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

std::vector<Point> SegmentedLine::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
		ret.push_back(*boundingPoints[i]) ;
	
	return ret ;
}

void SegmentedLine::sampleBoundingSurface(size_t num_points) 
{ 
}

void SegmentedLine::project(Point *p) const
{

	std::map<double, Point> projections ;
	
	for(size_t i = 0 ; i < boundingPoints.size()-1 ; i++)
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
	return false ;
	for (size_t i = 0 ; i < boundingPoints.size() ;i++)
		if(*boundingPoints[i] == v)
			return true ;
	
	return false ;
}

Ellipse::Ellipse(double a, double b, double originX,double originY, double axisX, double axisY)
{
	gType = ELLIPSE ;
	this->center = Point(originX, originY) ;
	double axisNorm = sqrt(axisX*axisX + axisY*axisY) ;
	if(axisNorm==0)
	{
		this->majoraxis = Point(1,0) ;
	} else {
		this->majoraxis = Point(axisX/axisNorm, axisY/axisNorm) ;
	}
	this->majorradius = std::max(a,b) ;
	this->minorradius = std::min(a,b) ;
	this->setSqRadius() ;
	this->setExcentricity() ;
}

Ellipse::Ellipse(double a, double b, double x,double y)
{
	gType = ELLIPSE ;
	this->center = Point(x, y) ;
	double angle = (((double)rand()/(double)RAND_MAX)-0.5)*M_PI/2 ;
//	std::cout << axisX << " ; " << axisY << std::endl ;
	majoraxis = Point(cos(angle), sin(angle)) ;

	this->majorradius = std::max(a,b) ;
	this->minorradius = std::min(a,b) ;
	this->setSqRadius() ;
	this->setExcentricity() ;
}

Ellipse::Ellipse(double a, double b, const Point *center, const Point *axis)
{
	gType = ELLIPSE ;
	this->center = Point(*center) ;
	if(axis->norm()==0)
	{
		this->majoraxis = Point(1,0) ;
	} else {
		this->majoraxis.x = (1./axis->norm()) * axis->x ;
		this->majoraxis.y = (1./axis->norm()) * axis->y ;
	}
	this->majorradius = std::max(a,b) ;
	this->minorradius = std::min(a,b) ;
	this->setSqRadius() ;
}

Ellipse::Ellipse(double a, double b, const Point center, const Point axis)
{
	gType = ELLIPSE ;
	this->center = Point(center) ;
	if(axis.norm()==0)
	{
		this->majoraxis = Point(1,0) ;
	} else {
		this->majoraxis.x = (1./axis.norm()) * axis.x ;
		this->majoraxis.y = (1./axis.norm()) * axis.y ;
	}
	this->majorradius = std::max(a,b) ;
	this->minorradius = std::min(a,b) ;
	this->setSqRadius() ;
}

Ellipse::Ellipse(double a, double b, const Point center)
{
	gType = ELLIPSE ;
	this->center = Point(center) ;
	double axisX = (double)rand()/(double)RAND_MAX ;
	double axisY = (double)rand()/(double)RAND_MAX ;
//	std::cout << axisX << " ; " << axisY << std::endl ;
	Point axis(axisX,axisY) ;
	if(axis.norm()==0)
	{
		this->majoraxis = Point(1,0) ;
	} else {
		this->majoraxis.x = (1./axis.norm()) * axis.x ;
		this->majoraxis.y = (1./axis.norm()) * axis.y ;
	}
	this->majorradius = std::max(a,b) ;
	this->minorradius = std::min(a,b) ;
	this->setSqRadius() ;
}

Ellipse::Ellipse(const Ellipse &e)
{
	gType = ELLIPSE ;
	this->center = e.getCenter() ;
	this->majoraxis = e.getMajorAxis() ;
	this->majorradius = e.getMajorRadius() ; 
	this->minorradius = e.getMinorRadius() ;
	this->setSqRadius() ;
}

XMLTree * Ellipse::toXML()
{
	XMLTree * ell = new XMLTree("ellipse") ;
	XMLTree * c = new XMLTree("center") ;
	c->addChild(this->getCenter().toXML()) ;
	XMLTree * ax = new XMLTree("axis") ;
	ax->addChild(this->getMajorAxis().toXML()) ;
	ell->addChild(c) ;
	ell->addChild(ax) ;
	ell->addChild(new XMLTree("ra",majorradius)) ;
	ell->addChild(new XMLTree("rb",minorradius)) ;
	return ell ;
}


double Ellipse::getAxisAngle() const
{
	// normally, axis is a normalized vector, so that sin(theta) = axis.y and cos(theta) = axis.x
/*	double theta = 0 ;
	if (majoraxis.x>0)
	{
		theta = asin(majoraxis.y) ;
	} else {
		theta = asin(-majoraxis.y) ;
	}*/
	return majoraxis.angle() ;	
}

Point Ellipse::getFocus(bool dir) const
{
	int dirsign = 1 ;
	if(dir)
		{ dirsign = - 1 ; }
	Point focus(center + majoraxis * (excentricity * majorradius * dirsign)) ;
	return focus ;
}

const std::pair<Point, Point> Ellipse::getBothFocus()
{
	std::pair<Point, Point> focus(getFocus(true),getFocus(false)) ;
	return focus ;
}

Point Ellipse::getMinorAxis() const
{
	return Point (-majoraxis.y, majoraxis.x) ;
}

void Ellipse::setRadius(double newa, double newb)
{
	double ratioa = newa/(majorradius) ;
	double ratiob = newb/(minorradius) ;

	Point minoraxis = this->getMinorAxis() ;
	Point inlocalreferencial(center) ;
	double inlocalreferencialmajor = 1. ;
	double inlocalreferencialminor = 0. ;

	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		inlocalreferencial = getBoundingPoint(i) - center ;
		inlocalreferencialmajor = this->majoraxis * inlocalreferencial ;
		inlocalreferencialminor = minoraxis * inlocalreferencial ;
		inlocalreferencial = (this->majoraxis) * (inlocalreferencialmajor * ratioa) ;
		inlocalreferencial = inlocalreferencial + minoraxis * (inlocalreferencialminor * ratiob) ;
		getBoundingPoint(i) = center + inlocalreferencial ;
	}
	
	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		inlocalreferencial = getInPoint(i) - center ;
		inlocalreferencialmajor = this->majoraxis * inlocalreferencial ;
		inlocalreferencialminor = minoraxis * inlocalreferencial ;
		inlocalreferencial = (this->majoraxis) * (inlocalreferencialmajor * ratioa) ;
		inlocalreferencial = inlocalreferencial + minoraxis * (inlocalreferencialminor * ratiob) ;
		getInPoint(i) = center + inlocalreferencial ;
	}


	if (newa<newb) 
	{
		this->majorradius = newb ;
		this->minorradius = newa ;
		this->majoraxis = getMinorAxis() ;
	} else {
		this->majorradius = newa ;
		this->minorradius = newb ;
	}
	this->setSqRadius() ;
}

void Ellipse::setAxis(Point newaxis)
{
	Point minoraxis = this->getMinorAxis() ;
	Point normnewaxis(newaxis/newaxis.norm()) ;
	Point orthoaxis(normnewaxis.y, normnewaxis.x) ;
	Point inlocalreferencial(center) ;
	Point rotinlocalreferencial(center) ;
	double inlocalreferencialmajor = 1. ;
	double inlocalreferencialminor = 0. ;

	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		inlocalreferencial = (getBoundingPoint(i) - center) ;
		inlocalreferencialmajor = (inlocalreferencial * majoraxis) ;
		inlocalreferencialminor = (inlocalreferencial * minoraxis) ;
		rotinlocalreferencial = (normnewaxis * inlocalreferencialmajor) + (orthoaxis * inlocalreferencialminor) ;
		getBoundingPoint(i) = center + rotinlocalreferencial ;
	}
	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		inlocalreferencial = (getInPoint(i) - center) ;
		inlocalreferencialmajor = (inlocalreferencial * majoraxis) ;
		inlocalreferencialminor = (inlocalreferencial * minoraxis) ;
		rotinlocalreferencial = (normnewaxis * inlocalreferencialmajor) + (orthoaxis * inlocalreferencialminor) ;
		getInPoint(i) = center + rotinlocalreferencial ;
	}
	
	this->majoraxis = normnewaxis ;
}

void Ellipse::computeCenter()
{
}

Point Ellipse::project(Point p) const
{	

	Point minoraxis = this->getMinorAxis() ;
	Point proj(p) ;
	if(p == center)
	{
		proj.x = center.x + majorradius*majoraxis.x ;
		proj.y = center.y + majorradius*majoraxis.y ;
	} else {
		Segment seg(center, p) ;
		proj = center - majoraxis * (seg.vector()*majoraxis) * (majorradius/seg.vector().norm()) - minoraxis * (seg.vector()*minoraxis) * (minorradius/seg.vector().norm()) ;
	}	
//	p->x = proj.x ;
//	p->y = proj.y ;
	
	return proj ;
	
}

void Ellipse::project(Point * p) const
{ 
	double alpha = majoraxis.angle() ;
	Point prot((*p-center).x*cos(-alpha)-(*p-center).y*sin(-alpha), 
		   +(*p-center).x*sin(-alpha)+(*p-center).y*cos(-alpha)) ;
	

	Ellipse ell(getMajorRadius(),
		    getMinorRadius(),
		    0,0,1,0) ;
	double theta = prot.angle() ;
	Point pell(std::cos(theta),std::sin(theta)) ;
	double r = ell.getRadiusOnEllipse(theta) ;
	Point ptry = pell * r ;
	double dist = squareDist2D(ptry,prot) ;

	int signn = + 1 ;
	double thetan = theta + 0.1 ;
	Point pelln(std::cos(thetan),std::sin(thetan)) ;
	double rn = ell.getRadiusOnEllipse(thetan) ;
	Point ptryn = pelln * rn ;
	double distn = squareDist2D(ptryn,prot) ;
	if(distn > dist)
	{
	  signn = - 1 ;
	  thetan = theta - 0.1 ;
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
	  thetan = thetan + signn * 0.1 ;
	  pelln.x = std::cos(thetan) ;
	  pelln.y = std::sin(thetan) ;
	  rn = ell.getRadiusOnEllipse(thetan) ;
	  ptryn = pelln * rn ;
	  distn = squareDist2D(ptryn,prot) ;
	}
	
//	std::cout << i << std::endl ;
//	p->print() ;
	
	p->x = center.x + ptryn.x*cos(alpha) - ptryn.y*sin(alpha) ;
	p->y = center.y + ptryn.x*sin(alpha) + ptryn.y*cos(alpha) ;

		   
/*	Function x("x") ;
	Function y("y") ;
	x /= majorradius ;
	y /= minorradius ;
	Function ellipse("x 2 ^ y 2 ^ + 1 -") ;
	
	double dx = VirtualMachine().deval(ellipse, XI, prot) ;
	double dy = VirtualMachine().deval(ellipse, ETA, prot) ;
	
/*	double A = (prot.x * prot.x) / (majorradius * majorradius) + (prot.y * prot.y) / (minorradius * minorradius) - 1 ;
	double B = 2 *(prot.x * dx / (majorradius * majorradius) + prot.y * dy / (minorradius * minorradius)) ;
	double C = (dx * dx) / (majorradius * majorradius) + (dy * dy) / (minorradius * minorradius) ;
	double DELTA = B*B - 4*A*C ;
	
	if(DELTA < 0)
	{
	  std::cout << "unable to project" << std::endl ;
	  return ;
	}
	if(DELTA == 0)
	{
	  std::cout << "you may have a mistake somewhere..." << std::endl ;
	  return ;
	}
	double t1 = (- B + sqrt(DELTA)) / (2 * A) ;
	double t2 = (- B - sqrt(DELTA)) / (2 * A) ;
	double tf = t1 ;
	if(std::abs(t2) < std::abs(t1))
	  tf = t2 ;
	double xf = prot.x + tf * dx ;
	double yf = prot.y + tf * dy ;

/*	Point rmajaxis((majoraxis.x)*sin(alpha)+(majoraxis.y)*cos(alpha), -(majoraxis.x)*sin(alpha)+(majoraxis.y)*cos(alpha)) ;
	Point rminoraxis(-rmajaxis.y, rmajaxis.x) ;
	if(std::abs(rmajaxis.angle()-M_PI) < POINT_TOLERANCE || std::abs(rmajaxis.angle()) < POINT_TOLERANCE)
	{
	  prot.x *= minorradius/majorradius ;
	}
	else
	{
	  prot.y *= minorradius/majorradius ;
	}
	
	Circle(minorradius, 0,0).project(&prot) ;
	
	 if(std::abs(rmajaxis.angle()-M_PI) < POINT_TOLERANCE || std::abs(rmajaxis.angle()) < POINT_TOLERANCE)
	{
	  prot.x /= minorradius/majorradius ;
	}
	else
	{
	  prot.y /= minorradius/majorradius ;
	}
	
	Point ret ;
	ret.x = xf*cos(alpha)+yf*sin(alpha) + center.x;
	ret.y =-xf*sin(alpha)+yf*cos(alpha) + center.y;
	p->x = ret.x ;
	p->y = ret.y ;
	
	/*double theta = (ret - center).angle() ;
	std::cout << getRadiusOnEllipse(tetha) << std::endl ;
	
	return ;
	
//	std::cout << "i am your father" << std::endl ;

//	Point minoraxis = this->getMinorAxis() ;
//	Point proj(p) ;
	if(*p == center)
	{
		p->x = center.x + majorradius*majoraxis.x ;
		p->y = center.y + majorradius*majoraxis.y ;
		return ;
	}

//	std::cout << "in Ellipse::project()" << std::endl ;

	Point proj(p->x, p->y) ; 
	proj = proj - center ;
	Point origine(0,0) ;
	Point proj_(proj) ;
	proj_.x = (proj * getMajorAxis()) * std::cos(getMajorAxis().angle()) - (proj * getMinorAxis()) * std::sin(getMajorAxis().angle()) ;
	proj_.y = (proj * getMajorAxis()) * std::sin(getMajorAxis().angle()) + (proj * getMinorAxis()) * std::cos(getMajorAxis().angle()) ;
	proj_.x = proj_.x * minorradius / majorradius ;
//	proj = toSmallCircle(proj, origine, true) ;
	Circle smallcircle(minorradius, origine) ;
	smallcircle.project(&proj_) ;
//	std::cout << (std::abs(squareDist2D(proj_, origine) - minorradius*minorradius) < POINT_TOLERANCE) << std::endl ;
	proj_.x = proj_.x * majorradius / minorradius ;
//	proj = toSmallCircle(proj, origine, false) ;
	proj.x = proj_.x * std::cos(getMajorAxis().angle()) + proj_.y * std::sin(getMajorAxis().angle()) ;
	proj.y = - proj_.x * std::sin(getMajorAxis().angle()) + proj_.y * std::cos(getMajorAxis().angle()) ;
	proj = getMajorAxis() * proj.x + getMinorAxis() * proj.y ;
	proj = proj + center ;
	

//	Point test(p->x,p->y) ;
//	test.print() ;
//	test = toSmallCircle(test, true) ;
//	test = toSmallCircle(test, false) ;
//	test.print() ;
//	std::cout << (test.x - p->x < POINT_TOLERANCE) << " ; " << (test.y - p->y < POINT_TOLERANCE) << std::endl ;
	
//	double theta = (proj-center).angle() - getMajorAxis().angle() ;
//	double r = getRadiusOnEllipse(theta) ;
//	std::cout << (std::abs(r*r - squareDist2D(proj,center)) < POINT_TOLERANCE) << std::endl ;
//	Ellipse etest(this) ;
//	Point q = getTangentDirection(theta) ;
//	std::cout << (q * (test - proj)) << std::endl ;
//	std::cout << minorradius << ";" << majorradius << ";" << theta << r << std::endl ;
//	proj = center + (proj - center) * r / sqrt(squareDist2D(proj, center)) ;
//	Point t1(center) ;
//	Point t2(p->x,p->y) ;
//	Point t3(proj) ;
//	std::cout << r << " ; " << sqrt(squareDist2D(t1,t2)) << " ; " << sqrt(squareDist2D(t1,t3)) << std::endl ; 
	
	p->x = proj.x ;
	p->y = proj.y ;

//	std::cout << sqrt(squareDist2D(proj,center)) << std::endl ;
	
/*	Segment seg(center, *p) ;

	Point proj = center - majoraxis * (seg.vector()*majoraxis) * (majorradius/seg.vector().norm()) - minoraxis * (seg.vector()*minoraxis) * (minorradius/seg.vector().norm()) ;

	
	p->x = center.x + r * majoraxis.x ;
	p->y = center.y + r * majoraxis.y ;*/
	
	return ;
	
}

double Ellipse::getRadiusOnEllipse(double theta) const
{
	double r = (sqradius / sqrt(pow(minorradius * cos(theta), 2) + pow(majorradius * sin(theta), 2))) ;
	return r ; 
}

double Ellipse::getRadiusOnEllipseFromFocus(double theta, bool dir) const
{
	int f = 1 ;
	if(!dir)
		f = -1 ;
	return getParameter() / (1 + f * getExcentricity() * cos(theta)) ;
}


std::vector<Point> Ellipse::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;
	
	double angle = 2 * M_PI / num_points ;
/*	double factor = 0 ;
	double r = 0 ;
	for(size_t i = 0 ; i < num_points ; i++)
	{
		r = getRadiusOnEllipseFromFocus(angle * i, true) ;
		factor = 1 / (r*r)  ;
	}

	factor = 2 * M_PI / factor ;
	angle = 0 ;
	Point p(center) ;

	for(size_t i = 0 ; i < num_points ; i++)
	{
		r = getRadiusOnEllipseFromFocus(angle, true) ;
		p  = getFocus(true) + getMajorAxis() * cos(angle) * r + getMinorAxis() * sin(angle) * r ;
		ret.push_back(p) ;
		angle += factor / (r*r) ;
	}

	return ret ;*/

	double thisangle = angle ;
	double lastangle = 0. ;
	double redfactor = 0.8 ; // factor for angle decrease < 1
	Point minoraxis = getMinorAxis() ;

//	Point firstfocus = getFocus(true) ;
//	Point secondfocus = getFocus(false) ;
//	Point bp(firstfocus) ;
	
	Point thispoint = center + majoraxis * (majorradius * cos(angle)) + minoraxis * (minorradius * sin(angle)) ;
	Point lastpoint = center + majoraxis * majorradius ;
	Point lastlastpoint = center + majoraxis * (majorradius * cos(-angle)) + minoraxis * (minorradius * sin(-angle)) ;
	double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

//	ret.push_back(lastpoint) ;
	int n_iter = 0 ;

	for (size_t i = 0 ; i< num_points + 1 ; i++)
	{
		n_iter = 0 ;
		while(((criteria > angle) || (criteria < angle * redfactor)) && n_iter < 20)
		{
			if(criteria > angle)
				{ thisangle = lastangle + (thisangle - lastangle) * redfactor ; }
			else
				{ thisangle = lastangle + (thisangle - lastangle) / redfactor ; }
			thispoint = center + majoraxis * (majorradius * cos(thisangle)) + minoraxis * (minorradius * sin(thisangle)) ;
			criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
			n_iter++ ;
		}

		ret.push_back(thispoint) ;

		lastlastpoint = lastpoint ;
		lastpoint = thispoint ;
		
		lastangle = thisangle ;
		thisangle += angle ;
		
		thispoint = center + majoraxis * (majorradius * cos(thisangle)) + minoraxis * (minorradius * sin(thisangle)) ;
	}
	
	return ret ;
}

double Ellipse::area() const
{
	return M_PI * majorradius * minorradius ;
}

double Ellipse::getSquareRadius() const
{
	return sqradius ;
}

Point Ellipse::toSmallCircle(Point p, Point origin, bool b) const
{
	Point tsc(p - origin) ;
	double tscmajor = tsc * majoraxis * minorradius / majorradius ;
	double tscminor = tsc * getMinorAxis() ;
	if(!b)
	  tscmajor = tsc * majoraxis * majorradius / minorradius ;
//	p.print() ;
//	(center + majoraxis * tscmajor + getMinorAxis() * tscminor).print() ;
	return origin + majoraxis * tscmajor + getMinorAxis() * tscminor ;
}

Point Ellipse::getTangentDirection(double theta) const
{
	double alpha = this->getAxisAngle() ;
	
	double tangentx = - majorradius * sin(theta) ;
	double tangenty = minorradius * cos(theta) ;

	Point tangent(tangentx * cos(alpha) + tangenty * sin(alpha), - tangentx * sin(alpha) + tangenty * cos(alpha)) ;

	return tangent ;
}

bool Ellipse::in(const Point &p) const 
{ 
//	std::cout << "in Ellipse::in()" << std::endl ;
	double theta = (p - center).angle() - majoraxis.angle() ;
	double dist = getRadiusOnEllipse(theta) ;
//	std::cout << dist << std::endl ;
//	std::cout << squareDist2D(p, center) << std::endl ;
	return (squareDist2D(p, center) < dist * dist) ; 
}

std::vector<Point> Ellipse::getSampleBoundingPointsOnArc(size_t num_points, double alpha, double beta) const 
{
	std::vector<Point> ret ;

	double angle = (beta - alpha) / (num_points) ;
	double thisangle = alpha + angle ;
	double lastangle = alpha ;
	double redfactor = 0.8 ; // factor for angle decrease < 1
	Point minoraxis = getMinorAxis() ;
	
	Point thispoint = center + majoraxis * (majorradius * cos(thisangle)) + minoraxis * (minorradius * sin(thisangle)) ;
	Point lastpoint = center + majoraxis * (majorradius * cos(lastangle)) + minoraxis * (minorradius * sin(lastangle)) ;
	Point lastlastpoint = center + majoraxis * (majorradius * cos(lastangle-angle)) + minoraxis * (minorradius * sin(lastangle-angle)) ;
	double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

//	ret.push_back(lastpoint) ;

	for (size_t i = 1 ; i< num_points + 1 ; i++)
	{
		while((criteria > angle) || (criteria < angle * redfactor))
		{
			if(criteria > angle)
				{ thisangle = lastangle + (thisangle - lastangle) * redfactor ; }
			else
				{ thisangle = lastangle + (thisangle - lastangle) / redfactor ; }
			thispoint = center + majoraxis * (majorradius * cos(thisangle)) + minoraxis * (minorradius * sin(thisangle)) ;
			criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
		}

		ret.push_back(thispoint) ;

		lastlastpoint = lastpoint ;
		lastpoint = thispoint ;
		
		lastangle = thisangle ;
		thisangle += angle ;
		
		thispoint = center + majoraxis * (majorradius * cos(thisangle)) + minoraxis * (minorradius * sin(thisangle)) ;
	}
	
	return ret ;


}

void Ellipse::sampleBoundingSurface (size_t num_points)
{
	std::vector<Point> bound = this->getSamplingBoundingPoints(num_points) ;

	getBoundingPoints().resize(num_points) ;

	for (size_t i = 0 ; i < num_points ; i++)
	{		
		getBoundingPoints()[i] = new Point(bound[i]) ;
	}
	
}

void Ellipse::sampleSurface (size_t num_points)
{
	if(boundingPoints.size() == 0)
//		this->sampleBoundingSurface(13*num_points*pow(getMajorRadius()/getMinorRadius(),0.666666)/8) ;
		this->sampleBoundingSurface(num_points * 2 / 3) ;
	sampled = true ;

	size_t ring = 1 + num_points / (3 * M_PI / 2) ;
	if(getMinorRadius() / getMajorRadius() < 0.7071 || ring==1)
	{	
		ring++ ;
//		std::cout << "add more points inside" << std::endl ;
	}

//	if(ring < 2)
//		ring = 2 ;
/*	std::vector<double> vangle(boundingPoints.size()) ;
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
		vangle[i] = (boundingPoints[i])->angle() - majoraxis.angle() ;
	
	std::vector<Point*> temp ;
	
	for (size_t j = 0 ; j < num_points ; ++j)
	{
		double factor = (double) (j+1) / (double) num_points ;
		double r = getRadiusOnEllipse(vangle[j]) * factor ;
		temp.push_back(new Point(center + majoraxis * r * cos(vangle[j] + majoraxis.angle()) + getMinorAxis() * r * sin(vangle[j] + majoraxis.angle()))) ;
	}*/

	std::vector<Point*> temp ;
	double newb = getMinorRadius() ;
	double newa = getMajorRadius() ;
	std::vector<double> newalist ;
	newalist.push_back(getMajorRadius()) ;

        int factor = 2 ;
        int jfactor = 0 ;
        bool jchange = false ;
//	if(getMinorRadius() / getMajorRadius() < 0.7071)// || ring==1)
//	{	
//		factor = 1 ;
//		std::cout << "add more points inside" << std::endl ;
//	}

	newb = getMinorRadius() * (ring) / (ring + 1) ;
	newa = getMinorRadius() * (ring) / (ring + 1) + (getMajorRadius() - getMinorRadius()) * (ring - 1) / (ring) ;

	if(newb/newa < 0.7071)
		ring = ring + 1 ;


	for(size_t j = 0 ; j < ring ; j++)
	{
		newb = getMinorRadius() * (ring - j) / (ring + 1) ;
		newa = getMinorRadius() * (ring - j) / (ring + 1) + (getMajorRadius() - getMinorRadius()) * (ring - j - 1) / (ring) ;
		newalist.push_back(newa) ;
                for(size_t i = 0 ; i < getBoundingPoints().size() / (factor * (jfactor + 1)) ; i++)
		{
			temp.push_back(new Point(center + 
                                                getMajorAxis() * ((getBoundingPoint(i * (factor * (jfactor + 1))) - center) * getMajorAxis()) * newa / getMajorRadius() +
                                                getMinorAxis() * ((getBoundingPoint(i * (factor * (jfactor + 1))) - center) * getMinorAxis()) * newb / getMinorRadius())) ;
		}
                if(jchange)
                {
                    jfactor++ ;
		}
                jchange != jchange ;
		if(j > 4)
			jchange = true ;
		
		factor = 2 ;
	}
/*	double r = sqrt(majorradius * minorradius) ;
	double x = majorradius * sqrt(r*r - minorradius*minorradius) / sqrt(majorradius*majorradius - minorradius*minorradius) ;
	double y = minorradius * sqrt(majorradius*majorradius - r*r) / sqrt(majorradius*majorradius - minorradius*minorradius) ;

	for(size_t j = 0 ; j < num_points / 4 ; j++)
	{
		if(j%2 == 0)
			temp.push_back(new Point(center + majoraxis * majorradius * j / (num_points / 4))) ;
		if(j%2 == 1)
			temp.push_back(new Point(center - majoraxis * majorradius * j / (num_points / 4))) ;
		if(j%2 == 0)
			temp.push_back(new Point(center + getMinorAxis() * minorradius * j / (num_points / 4))) ;
		if(j%2 == 1)
			temp.push_back(new Point(center - getMinorAxis() * minorradius * j / (num_points / 4))) ;
		if(j%2 == 0)
			temp.push_back(new Point(center + majoraxis * x * j / (num_points / 4) + getMinorAxis() * y * j / (num_points / 4))) ;
		if(j%2 == 1)
			temp.push_back(new Point(center + majoraxis * x * j / (num_points / 4) - getMinorAxis() * y * j / (num_points / 4))) ;
		if(j%2 == 0)
			temp.push_back(new Point(center - majoraxis * x * j / (num_points / 4) + getMinorAxis() * y * j / (num_points / 4))) ;
		if(j%2 == 1)
			temp.push_back(new Point(center - majoraxis * x * j / (num_points / 4) - getMinorAxis() * y * j / (num_points / 4))) ;
	}*/
	
	int toadd = 1 ;
	if(getMinorRadius() / getMajorRadius() < 0.5)// || ring==1)
	{	
		toadd = 1 + 4 * (std::min((int) newalist.size(),4) - 1) ;
//		std::cout << "add more points inside" << std::endl ;
	}
	inPoints.resize(temp.size() + toadd) ;
	inPoints[0] = new Point(center) ;
	if(toadd > 1)
	{
		for(size_t j = 0 ; j < newalist.size() - 1 ; j++)
		{
			inPoints[1 + j] = new Point(center + getMajorAxis() * (newalist[j] + newalist[j + 1] * 2) / 3) ;
			inPoints[newalist.size() + j] = new Point(center + getMajorAxis() * (newalist[j] * 2 + newalist[j + 1]) / 3) ;
			inPoints[newalist.size() * 2 + j - 1] = new Point(center - getMajorAxis() * (newalist[j] + newalist[j + 1] * 2) / 3) ;
			inPoints[newalist.size() * 3 + j - 2] = new Point(center - getMajorAxis() * (newalist[j] * 2 + newalist[j + 1]) / 3) ;
		}
	}
	for(size_t i = 0 ; i < temp.size() ; i++)
		inPoints[i+toadd] = temp[i] ;
	//std::copy(temp.begin(), temp.end(),&inPoints[1]) ;
	//std::cout << "we have " << num_points << " sample points" << std::endl ;
}

double Ellipse::getRadius() const
{
	return majorradius ;
}

std::vector<Point> Ellipse::getBoundingBox() const
{
	std::vector<Point> bbox(4) ;
	Point minoraxis = this->getMinorAxis() ;
	
	bbox[0] = center + majoraxis * majorradius + minoraxis * minorradius ;
	bbox[1] = center + majoraxis * majorradius - minoraxis * minorradius ;
	bbox[2] = center - majoraxis * majorradius - minoraxis * minorradius ;
	bbox[3] = center - majoraxis * majorradius + minoraxis * minorradius ;

	return bbox ;
}

double Ellipse::getMinorRadius() const
{
	return minorradius ;
}

Point Ellipse::getMajorAxis() const
{
	return majoraxis ;
}





