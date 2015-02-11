// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_2D.h"
#include "../polynomial/vm_base.h"
#include "../utilities/random.h"

using namespace Amie ;

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


	if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
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

	radius = sqrt((squareDist2D(p1, circumCenter)+ squareDist2D(p0, circumCenter))*.5);
	sqradius = radius*radius ;

}

void Triangle::setCenter(const Point & newCenter)
{
	Geometry::setCenter(newCenter);
	computeCircumCenter();
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

	if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}

	computeCircumCenter() ;
	computeCenter() ;

	if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
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

	if((p0->getZ() == p1->getZ())  && (p0->getZ() == p2->getZ()) &&  (p0->getZ() == 0))
	{
		if(!isTrigoOriented())
		{
			std::swap(boundingPoints[1], boundingPoints[2]) ;
		}
	}

	computeCircumCenter() ;
	computeCenter() ;

	radius = (squareDist2D(*p1, circumCenter) + squareDist2D(*p0, circumCenter) + squareDist2D(*p2, circumCenter)+squareDist2D(*p3, circumCenter))/4.;

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

std::vector<Point> OrientedRectangle::getBoundingBox() const
{
	double maxx = boundingPoints[0]->getX() ;
	double minx = boundingPoints[0]->getX() ;
	double maxy = boundingPoints[0]->getY() ;
	double miny = boundingPoints[0]->getY() ;

	for(size_t i = boundingPoints.size()/4 ; i < boundingPoints.size() ; i += boundingPoints.size()/4 )
	{
		double x = boundingPoints[i]->getX() ;
		double y = boundingPoints[i]->getY() ;
		if(x > maxx)
			maxx = x ;
		if(x < minx)
			minx = x ;
		if(y > maxy)
			maxy = y ;
		if(y < miny)
			miny = y ;
	}

	std::vector<Point> ret ;
	ret.push_back(Point(minx,miny)) ;
	ret.push_back(Point(maxx,miny)) ;
	ret.push_back(Point(maxx,maxy)) ;
	ret.push_back(Point(minx,maxy)) ;

	return  ret ;
}

void OrientedRectangle::computeCircumCenter()
{
	if (fabs(boundingPoints[1]->getY()-boundingPoints[0]->getY()) < 20*POINT_TOLERANCE_2D)
	{
		double m2 = - (boundingPoints[2]->getX()-boundingPoints[1]->getX()) / (boundingPoints[2]->getY()-boundingPoints[1]->getY());
		double mx2 = (boundingPoints[1]->getX() + boundingPoints[2]->getX()) / 2.0;
		double my2 = (boundingPoints[1]->getY() + boundingPoints[2]->getY()) / 2.0;
		double xc = (boundingPoints[1]->getX() + boundingPoints[0]->getX()) / 2.0;
		double yc = fma(m2, (xc - mx2), my2);

		circumCenter = Point(xc, yc) ;
	}
	else if (fabs(boundingPoints[2]->getY()-boundingPoints[1]->getY()) < 20*POINT_TOLERANCE_2D)
	{
		double m1 = - (boundingPoints[1]->getX()-boundingPoints[0]->getX()) / (boundingPoints[1]->getY()-boundingPoints[0]->getY());
		double mx1 = (boundingPoints[0]->getX() + boundingPoints[1]->getX()) / 2.0;
		double my1 = (boundingPoints[0]->getY() + boundingPoints[1]->getY()) / 2.0;
		double xc = (boundingPoints[2]->getX() + boundingPoints[1]->getX()) / 2.0;
		double yc = fma(m1, (xc - mx1), my1);

		circumCenter = Point(xc, yc) ;
	}
	else
	{
		double m1 = - (boundingPoints[1]->getX()-boundingPoints[0]->getX()) / (boundingPoints[1]->getY()-boundingPoints[0]->getY());
		double m2 = - (boundingPoints[2]->getX()-boundingPoints[1]->getX()) / (boundingPoints[2]->getY()-boundingPoints[1]->getY());
		double mx1 = (boundingPoints[0]->getX() + boundingPoints[1]->getX()) / 2.0;
		double mx2 = (boundingPoints[1]->getX() + boundingPoints[2]->getX()) / 2.0;
		double my1 = (boundingPoints[0]->getY() + boundingPoints[1]->getY()) / 2.0;
		double my2 = (boundingPoints[1]->getY() + boundingPoints[2]->getY()) / 2.0;
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

	Point c0 = getBoundingPoint(0) ;
	Point c1 = getBoundingPoint(boundingPoints.size()/4) ;
	Point c2 = getBoundingPoint(boundingPoints.size()*2/4) ;
	Point c3 = getBoundingPoint(boundingPoints.size()*3/4) ;

	Segment s0(c0, c1) ;
	Segment s1(c1, c2) ;
	Segment s2(c2, c3) ;
	Segment s3(c3, c0) ;

	Point p0(*p) ;
	Point p1(*p) ;
	Point p2(*p) ;
	Point p3(*p) ;

	p0 = s0.project(p0) ;
	p1 = s1.project(p1) ;
	p2 = s2.project(p2) ;
	p3 = s3.project(p3) ;

	double d0 = squareDist3D(p0, *p) ;
	double d1 = squareDist3D(p1, *p) ;
	double d2 = squareDist3D(p2, *p) ;
	double d3 = squareDist3D(p3, *p) ;

	int i = 0 ;
	double d = d0 ;
	if(d1 < d)
	{
		d = d1 ;
		i = 1 ;
	}
	if(d2 < d)
	{
		d = d2 ;
		i = 2 ;
	}
	if(d3 < d)
	{
		d = d3 ;
		i = 3 ;
	}

	switch(i)
	{
	case 0:
		p->getX() = p0.getX() ;
		p->getY() = p0.getY() ;
		break ;
	case 1:
		p->getX() = p1.getX() ;
		p->getY() = p1.getY() ;
		break ;
	case 2:
		p->getX() = p2.getX() ;
		p->getY() = p2.getY() ;
		break ;
	case 3:
		p->getX() = p3.getX() ;
		p->getY() = p3.getY() ;
		break ;
	}

	return ;

}


bool OrientedRectangle::in(const Point &p) const
{

	Point c0 = getBoundingPoint(0) ;
	Point c1 = getBoundingPoint(boundingPoints.size()/4) ;
	Point c2 = getBoundingPoint(boundingPoints.size()*2/4) ;
	Point c3 = getBoundingPoint(boundingPoints.size()*3/4) ;

	Segment s0(c0, c1) ;
	Segment s1(c1, c2) ;
	Segment s2(c2, c3) ;
	Segment s3(c3, c0) ;

	Segment t(getCenter(), p) ;

	return !(t.intersects(s0) || t.intersects(s1) || t.intersects(s2) || t.intersects(s3)) ;

	bool in = false ;

	for (size_t i = 0, j  =  boundingPoints.size()-1; i <  boundingPoints.size(); j = i++)
	{
		if (
			(((boundingPoints[i]->getY() <= p.getY() )
			  && (p.getY()<boundingPoints[j]->getY()))
			 || ((boundingPoints[j]->getY() <= p.getY())
				 && (p.getY()<boundingPoints[i]->getY())))
			&& (p.getX() < (boundingPoints[j]->getX() - boundingPoints[i]->getX()) * (p.getY() - boundingPoints[i]->getY()) / (boundingPoints[j]->getY() - boundingPoints[i]->getY()) + boundingPoints[i]->getX()))
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
	num_points = num_points + 4 - num_points%4 ;

	Point * v0(boundingPoints[0]) ;
	Point * v1(boundingPoints[boundingPoints.size()/4]) ;
	Point * v2(boundingPoints[2*boundingPoints.size()/4]) ;
	Point * v3(boundingPoints[3*boundingPoints.size()/4]) ;

	for(size_t i = 1 ; i < boundingPoints.size() ; i++)
	{
		if(i != boundingPoints.size()/4 && i != 2*boundingPoints.size()/4 && i != 3*boundingPoints.size()/4)
			delete boundingPoints[i] ;
	}
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
//	size_t n = 2*num_points ;
	this->sampleBoundingSurface(num_points) ;

	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;

	std::vector<Point> newInPoints ;


	size_t numberOfPointsAlongX = boundingPoints.size()/4 ;
	size_t numberOfPointsAlongY = boundingPoints.size()/4 ;
	size_t numberOfBoundingPoints = boundingPoints.size() ;

	Point c = getBoundingPoint(boundingPoints.size()/4)-getBoundingPoint(0) ;

	UniformDistribution uniform(-0.1/numberOfPointsAlongY,+0.1/numberOfPointsAlongY) ;
	for(size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
		for(size_t j = 1 ; j < numberOfPointsAlongY ; j++)
		{
			double dj = (double) j / (double) numberOfPointsAlongY ;
			dj += uniform.draw() ;
			double di = uniform.draw() ;
			Point p1(*boundingPoints[i]) ;
			Point p2(*boundingPoints[numberOfBoundingPoints-numberOfPointsAlongX-i]) ;
			p1 = p1*dj + p2*(1.-dj) + c*di ;
			if(in(p1))
				newInPoints.push_back(p1) ;
		}
	}


	inPoints.resize(newInPoints.size()) ;
	for(size_t i = 0 ; i < inPoints.size() ; i++)
		inPoints[i] = new Point(newInPoints[i]) ;




}

Triangle::Triangle( Point *p0,  Point *p1,  Point *p2): ConvexGeometry(3)
{
	gType = TRIANGLE ;

	assert(this->size() == 3) ;


	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;


	if((p0->getZ() == p1->getZ())  && (p0->getZ() == p2->getZ()) &&  (p0->getZ() == 0))
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
	box.push_back(getCircumCenter()+Point(-getRadius(), getRadius())) ;
	box.push_back(getCircumCenter()+Point(getRadius(), getRadius())) ;
	box.push_back(getCircumCenter()+Point(getRadius(), -getRadius())) ;
	box.push_back(getCircumCenter()+Point(-getRadius(), -getRadius())) ;

	return box ;
}

void Triangle::computeCircumCenter()
{

	if (std::abs(getBoundingPoint(1).getY()-getBoundingPoint(0).getY()) < 100.*POINT_TOLERANCE_2D)
	{
		double m2  =  (getBoundingPoint(1).getX()-getBoundingPoint(2).getX() ) / (getBoundingPoint(2).getY()-getBoundingPoint(1).getY());
		double mx2 = (getBoundingPoint(1).getX() + getBoundingPoint(2).getX()) ;
		double my2 = (getBoundingPoint(1).getY() + getBoundingPoint(2).getY()) ;
		double xc  = (getBoundingPoint(1).getX() + getBoundingPoint(0).getX()) ;
		double yc  = m2 * (xc - mx2) + my2;

		circumCenter.set(xc/2., yc/2.) ;
	}
	else if (std::abs(getBoundingPoint(2).getY()-getBoundingPoint(1).getY()) < 100.*POINT_TOLERANCE_2D)
	{
		double m1  =  (getBoundingPoint(0).getX() - getBoundingPoint(1).getX() ) / (getBoundingPoint(1).getY()-getBoundingPoint(0).getY());
		double mx1 = (getBoundingPoint(0).getX() + getBoundingPoint(1).getX()) ;
		double my1 = (getBoundingPoint(0).getY() + getBoundingPoint(1).getY()) ;
		double xc  = (getBoundingPoint(2).getX() + getBoundingPoint(1).getX()) ;
		double yc  = m1 * (xc - mx1) + my1;

		circumCenter.set(xc/2., yc/2.) ;
	}
	else
	{
		double m1  = (getBoundingPoint(0).getX()-getBoundingPoint(1).getX()) / (getBoundingPoint(1).getY()-getBoundingPoint(0).getY());
		double m2  = (getBoundingPoint(1).getX()-getBoundingPoint(2).getX()) / (getBoundingPoint(2).getY()-getBoundingPoint(1).getY());
		double mx1 = (getBoundingPoint(0).getX() + getBoundingPoint(1).getX()) ;
		double mx2 = (getBoundingPoint(1).getX() + getBoundingPoint(2).getX()) ;
		double my1 = (getBoundingPoint(0).getY() + getBoundingPoint(1).getY()) ;
		double my2 = (getBoundingPoint(1).getY() + getBoundingPoint(2).getY()) ;
		double xc  = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		double yc  = m1 * (xc - mx1) + my1;

		circumCenter.set(xc/2., yc/2.) ;
	}
}

bool Triangle::inCircumCircle(const Point & p) const
{
	if(p.getX() > circumCenter.getX()+1.01*radius)
		return false ;
	if(p.getX() < circumCenter.getX()-1.01*radius)
		return false ;
	if(p.getY() > circumCenter.getY()+1.01*radius)
		return false ;
	if(p.getY() < circumCenter.getY()-1.01*radius)
		return false ;

	if(squareDist2D(circumCenter, p) < .99*sqradius)
		return true ;

	double delta = POINT_TOLERANCE_2D ;
	Point a(p) ;
	a.getX() += delta ;
	a.getY() += delta ;
	Point c(p) ;
	c.getX() += delta ;
	c.getY() -= delta ;
	Point e(p) ;
	e.getX() -= delta ;
	e.getY() += delta ;
	Point g(p) ;
	g.getX() -= delta ;
	g.getY() -= delta ;

	return  squareDist2D(circumCenter, a) < sqradius
			&&  squareDist2D(circumCenter, c) < sqradius
			&&  squareDist2D(circumCenter, e) < sqradius
			&&  squareDist2D(circumCenter, g) < sqradius;
	double x = circumCenter.getX() -p.getX() ;
	double y = circumCenter.getY() -p.getY() ;
	return  fma(x, x, y*y)< sqradius*(1. - 100.*POINT_TOLERANCE_2D)  ;
}

bool Triangle::inCircumCircle(const Point *p) const
{
	if(p->getX() > circumCenter.getX()+1.01*radius)
		return false ;
	if(p->getX() < circumCenter.getX()-1.01*radius)
		return false ;
	if(p->getY() > circumCenter.getY()+1.01*radius)
		return false ;
	if(p->getY() < circumCenter.getY()-1.01*radius)
		return false ;

	if(squareDist2D(circumCenter, *p) < .99*sqradius)
		return true ;

	double delta = POINT_TOLERANCE_2D ;
	Point a(*p) ;
	a.getX() += delta ;
	a.getY() += delta ;
	Point c(*p) ;
	c.getX() += delta ;
	c.getY() -= delta ;
	Point e(*p) ;
	e.getX() -= delta ;
	e.getY() += delta ;
	Point g(*p) ;
	g.getX() -= delta ;
	g.getY() -= delta ;

	return  squareDist2D(circumCenter, a) < sqradius
			&&  squareDist2D(circumCenter, c) < sqradius
			&&  squareDist2D(circumCenter, e) < sqradius
			&&  squareDist2D(circumCenter, g) < sqradius;
	double x = circumCenter.getX() -p->getX() ;
	double y = circumCenter.getY() -p->getY() ;
	return  fma(x, x, y*y) < sqradius*(1. - 100.*POINT_TOLERANCE_2D)  ;
}


double Triangle::area() const
{
// 	const DelaunayTriangle * tri = dynamic_cast<const DelaunayTriangle *>(this) ;
// 	if(tri)
// 	{
// 		Segment s0(*(tri->first), *(tri->second)) ;
// 		Segment s1(*(tri->first), *(tri->third)) ;
// 		return 0.5*std::abs((s0.vector()^s1.vector()).getZ()) ;
// 	}

// 	assert(this->boundingPoints.size() == 3) ;
	int pointsInTimePlane = this->boundingPoints.size()/timePlanes() ;
// 	if(getBoundingPoint(0).getT() != 0)
// 	{
// 		pointsInTimePlane = 0 ;
// 		double init = getBoundingPoint(0).getT() ;
// 		int counter = 0 ;
// 		while(std::abs(getBoundingPoint(counter++).getT()-init) < POINT_TOLERANCE_2D)
// 			pointsInTimePlane++ ;
// 	}
	Segment s0(getBoundingPoint(0), getBoundingPoint(pointsInTimePlane/3)) ;
	Segment s1(getBoundingPoint(0), getBoundingPoint(2*pointsInTimePlane/3)) ;

	return 0.5*std::abs((s0.vector()^s1.vector()).getZ()) ;
}


void Triangle::project(Point * p) const
{
	Segment s(getCenter(), *p) ;
	if(dist(*p, getCircumCenter()) < POINT_TOLERANCE_2D)
		return ;

	Segment s0(getBoundingPoint(0), getBoundingPoint(getBoundingPoints().size()/(3*timePlanes()))) ;
	Segment s1(getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes()))) ;
	Segment s2( getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(0)) ;
	std::map<double, Point> proj ;
	Point p0 = s0.project(*p) ;
	Point p1 = s1.project(*p) ;
	Point p2 = s2.project(*p) ;
	proj[dist(p0, *p)] = p0 ;
	proj[dist(p1, *p)] = p1 ;
	proj[dist(p2, *p)] = p2 ;
	*p = proj.begin()->second ;
	return ;

//
// 	for(size_t i = 0 ; i < pts.size() ; i++)
// 	{
// 		Segment seg(pts[i], pts[(i+1)%pts.size()]) ;
// 		if(s.intersects(seg))
// 		{
// 			Point t = s.intersection(seg);
// 			p->getX() = t.getX() ;
// 			p->getY() = t.getY() ;
// 			return ;
// 		}
// 	}
//
// 	double r = getRadius() ;
// 	Point reach = (*p - getCenter()) ;
// 	Point trans = getCenter() + reach * (2.*r/reach.norm()) ;
// 	Segment sec(getCenter(), trans) ;
//
// 	for(size_t i = 0 ; i < pts.size() ; i++)
// 	{
// 		Segment seg(pts[i], pts[(i+1)%pts.size()]) ;
// 		if(sec.intersects(seg))
// 		{
// 			Point t = sec.intersection(seg);
// 			p->getX() = t.getX() ;
// 			p->getY() = t.getY() ;
// 			return ;
// 		}
// 	}

}

bool Triangle::in(const Point &p) const
{

// 	p.print() ;
// 	bool isAPoint = false ;
// 	for (int i = 0; i <  getBoundingPoints().size(); i++)
// 	{
// 		if(squareDist2D( p, getBoundingPoint(i))  < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D)
// 		{
// 			return true ;
// 		}
// 	}
//
// 	Point proj(p) ; project(&proj) ;
// 	bool isOnSurface = squareDist2D(p, proj) < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D ;
// 	if(isOnSurface)
// 		return true ;
//
// 	Segment s(p, getCenter()) ;
//
//
// 				bool ret = false ;
// 			std::vector<Point> pts ;
// 			std::multimap<double, Point> pt ;
// 			for(size_t i = 0 ; i < getBoundingPoints().size() ;  i++)
// 			{
// 				pt.insert(std::make_pair(
// 				std::abs(
// 				squareDist2D(getCircumCenter(), getBoundingPoint(i))-getRadius()*getRadius()), getBoundingPoint(i)));
// 			}
// 			std::multimap<double, Point>::const_iterator ptend = pt.begin() ;
// 			ptend++ ; ptend++ ; ptend++ ;
//
// 			for(std::multimap<double, Point>::const_iterator i = pt.begin() ; i != ptend ; ++i )
// 				pts.push_back(i->second);
//
// 			if(s.on(pts[0]) || s.on(pts[1]) || s.on(pts[2]))
// 				return true ;
//
// 			Segment sa(pts[0],pts[1]) ;
// 			Segment sb(pts[1],pts[2]) ;
// 			Segment sc(pts[2],pts[0]) ;
//
//
// 			return !sa.intersects(*this) || sb.intersects(*this) || sc.intersects(*this) ;
//
//
// 	return !s.intersects(this) || isAPoint || isOnSurface;
//
// 		bool isAPoint = false ;
	for (int i = 0; i <  getBoundingPoints().size(); i++)
	{
		if(p == getBoundingPoint(i) || squareDist2D( p, getBoundingPoint(i))  < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D)
		{
			return true ;
		}
	}

	Point proj(p) ;
	project(&proj) ;
	bool isOnSurface = squareDist2D(p, proj) < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D ;
	if(isOnSurface)
		return true ;


	proj = p ;
	proj.getT() = getBoundingPoint(0).getT() ;

	Point c = getCenter() ;
	c.getT() = proj.getT() ;

	Segment s(proj, c) ;

	size_t npts = getBoundingPoints().size()/timePlanes() ;

	if(s.on(getBoundingPoint(0)) || s.on(getBoundingPoint(npts/3)) || s.on(getBoundingPoint(npts*2/3)))
		return false ;

	Segment sa(getBoundingPoint(0),getBoundingPoint(npts/3)) ;
	Segment sb(getBoundingPoint(npts/3),getBoundingPoint(npts*2/3)) ;
	Segment sc(getBoundingPoint(npts*2/3),getBoundingPoint(0)) ;

	return !(sa.intersects(s) || sb.intersects(s) || sc.intersects(s)) ;

}


std::vector<Point> Triangle::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;
	assert(num_points%3 == 0) ;

	int n = getBoundingPoints().size() ;

	Point v0 = *boundingPoints[0] ;
	Point v1 = *boundingPoints[n/3] ;
	Point v2 = *boundingPoints[2*n/3] ;

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
	Point * v1 = &getBoundingPoint(getBoundingPoints().size()/3) ;
	Point * v2 = &getBoundingPoint(2*getBoundingPoints().size()/3) ;

	/*	for(size_t i = 1 ; i < num_points/3 ; i++)
		{
			if(i != getBoundingPoints().size()/3 && i != 2*getBoundingPoints().size()/3)
				delete boundingPoints[i] ;
		}*/


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

	sampleBoundingSurface(num_points*3) ;

	std::vector<Point> newPoints ;

	size_t end_i = boundingPoints.size()/3 ;
	Point p1 = getBoundingPoint( 0 ) ;
	Point p2 = getBoundingPoint( end_i )-p1 ;
	Point p3 = getBoundingPoint( end_i*2 )-p1 ;
	double d = (p2.norm()+p3.norm())*0.01 ;

/*	for(double i = 1./num_points ; i < 1 ; i += 1./num_points)
	{
		for(double j = 1./num_points ; j < i ; j += 1./num_points)
		{
			Point test = p1 + p2*(i+xrand) + p3*(j+yrand) ;
			if(in(test))
				newPoints.push_back(test) ;
		}
	}*/




	for(int i = 0 ; i < ((int)num_points-1) ; i+=1)
	{
		for(int j = 0 ; j < i+1 ; j+=1)
		{
			double fact = (double)(j+1)/(double)(i+1) ;
			double xrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
			double yrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
			if(this->in(getBoundingPoint(i+1)*(1.-fact) + getBoundingPoint(end_i-1-i)*fact+Point(xrand, yrand)))
				newPoints.push_back(getBoundingPoint(i+1)*(1.-fact) + getBoundingPoint(end_i-1-i)*fact+Point(xrand, yrand)) ;
		}
	}
	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;

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
	topLeft = Point(originX-0.5*x, originY+0.5*y) ;
	boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(originX+0.5*x, originY+0.5*y) ;
	boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(originX+0.5*x, originY-0.5*y) ;
	boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft =  Point(originX-0.5*x, originY-0.5*y) ;
	boundingPoints[3] = new Point(bottomLeft) ;
}

Rectangle::Rectangle(double x, double y, const Point &center) :  ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE ;
	this->center = center ;
	topLeft = Point(center.getX()-0.5*x, center.getY()+0.5*y) ;
	boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(center.getX()+0.5*x, center.getY()+0.5*y) ;
	boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(center.getX()+0.5*x, center.getY()-0.5*y) ;
	boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft =  Point(center.getX()-0.5*x, center.getY()-0.5*y) ;
	boundingPoints[3] = new Point(bottomLeft) ;
}

Rectangle::Rectangle() :  ConvexGeometry(4), size_y(2), size_x(2)
{
	gType = RECTANGLE ;
	this->center = Point(0,0) ;
	topLeft = Point(-1, 1) ;
	boundingPoints[0] = new Point(topLeft) ;
	topRight = Point(1,1) ;
	boundingPoints[1] = new Point(topRight) ;
	bottomRight = Point(1, -1) ;
	boundingPoints[2] = new Point(bottomRight) ;
	bottomLeft = Point(-1, -1) ;
	boundingPoints[3] = new Point(bottomLeft) ;


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
	if(p.getX() < getCenter().getX() - 0.5*width())
		return false ;
	if(p.getX()  > getCenter().getX() + 0.5*width())
		return false ;
	if(p.getY() > getCenter().getY() + 0.5*height())
		return false ;
	if(p.getY()  < getCenter().getY() - 0.5*height())
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
		ret.push_back(Point(center.getX()-0.5*size_x, center.getY() + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy)) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
		double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongX-1)
			randx = 0 ;
		ret.push_back(Point( center.getX()-0.5*size_x+i*distanceBetweenPointsAlongX+ randx,
							 getCenter().getY()-0.5*size_y));
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongY-1)
			randy = 0 ;
		ret.push_back(Point(center.getX()+0.5*size_x,
							center.getY()-0.5*size_y+i*distanceBetweenPointsAlongY+ randy));
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
		ret.push_back(Point(
						  center.getX() + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
						  center.getY() + 0.5*size_y)) ;
	}

	return ret ;
}

void Rectangle::sampleBoundingSurface(size_t num_points)
{
	//	assert(num_points%4 == 0) ;
	double perimeter = 2*(size_x+size_y) ;

	double distanceBetweenPointsx = std::min(perimeter/num_points, size_x) ;
	double distanceBetweenPointsy = std::min(perimeter/num_points, size_y) ;
	double dy = distanceBetweenPointsy ;
	double dx = distanceBetweenPointsx ;
	distanceBetweenPointsy = std::min(dx*1.5, dy) ;
	distanceBetweenPointsx = std::min(dy*1.5, dx) ;

	this->numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x/distanceBetweenPointsx) + 1);
	double distanceBetweenPointsAlongX = size_x/(this->numberOfPointsAlongX-1) ;

	this->numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y/distanceBetweenPointsy) + 1);
	double distanceBetweenPointsAlongY = size_y/(this->numberOfPointsAlongY-1) ;

	num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;

	for(size_t i = 0 ; i < boundingPoints.size(); i++)
		delete boundingPoints[i] ;

	boundingPoints.resize(num_points) ;

	for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
	{
// 		double randx=((2.*rand()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == 0 || i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[i] = new Point(center.getX()-0.5*size_x, center.getY() + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy) ;
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
	{
		double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongX-1)
			randx = 0 ;
		boundingPoints[numberOfPointsAlongY+i-1] = new Point( center.getX()-0.5*size_x + i*distanceBetweenPointsAlongX+ randx,
				getCenter().getY()-0.5*size_y);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
	{
		double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		if(i == numberOfPointsAlongY-1)
			randy = 0 ;
		boundingPoints[numberOfPointsAlongX+numberOfPointsAlongY+i-2] = new Point(center.getX()+0.5*size_x,
				center.getY()-0.5*size_y+i*distanceBetweenPointsAlongY+ randy);
	}
	for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
	{
		double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
		assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
		boundingPoints[2*numberOfPointsAlongY+numberOfPointsAlongX+i-3] = new Point(
			center.getX() + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
			center.getY() + 0.5*size_y) ;
	}
}

void Rectangle::sampleSurface(size_t num_points)
{

	if(std::max(size_x/size_y, size_y/size_x) < 10)
	{
		sampleBoundingSurface((size_t)round((double)num_points*2.*std::max(size_x/size_y, size_y/size_x)/(M_PI))) ;
	}
	else if(std::max(size_x/size_y, size_y/size_x) < 60)
	{
		sampleBoundingSurface((size_t)round((double)num_points*0.5*std::max(size_x/size_y, size_y/size_x)/(M_PI))) ;
	}
	else
	{
		sampleBoundingSurface((size_t)round((double)num_points*0.2*std::max(size_x/size_y, size_y/size_x)/(M_PI))) ;
	}

	size_t nip = static_cast<size_t>((numberOfPointsAlongX-2)*(numberOfPointsAlongY-2)) ;

	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;

	std::vector<Point *> newInPoints ;

	double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;
	double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;

	if(nip > 0)
	{
		for(size_t j = 0 ; j < numberOfPointsAlongY-2 ; j++)
		{
			for(size_t i = 0 ; i < numberOfPointsAlongX-2+(j+1)%2 ; i++)
			{

				double randx= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.25*distanceBetweenPointsAlongX ;
				double randy= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.25*distanceBetweenPointsAlongY ;

				newInPoints.push_back( new Point(center.getX() - 0.5*size_x + (double)(i+0.66)*distanceBetweenPointsAlongX+(double)((j)%2)*distanceBetweenPointsAlongX*.5-(double)((j+1)%2)*distanceBetweenPointsAlongX*.15+randx,
												 center.getY() - 0.5*size_y + (double)(j+1)*distanceBetweenPointsAlongY+ randy)) ;
			}
		}
	}

	inPoints.resize(newInPoints.size()) ;
	for(size_t i = 0 ; i < inPoints.size() ; i++)
		inPoints[i] = newInPoints[i] ;
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

Circle::Circle(double r, const Point & center)
{
	gType = CIRCLE ;
	this->center = center ;
	this->radius = r ;
	this->sqradius = r*r ;
}


void Circle::setRadius(double newr)
{
	double ratio = newr/(radius) ;


	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		getBoundingPoint(i).getX() = (getBoundingPoint(i).getX() - center.getX())*ratio + center.getX() ;
		getBoundingPoint(i).getY() = (getBoundingPoint(i).getY() - center.getY())*ratio + center.getY() ;
	}

	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		getInPoint(i).getX() = (getInPoint(i).getX() - center.getX())*ratio + center.getX() ;
		getInPoint(i).getY() = (getInPoint(i).getY() - center.getY())*ratio + center.getY() ;
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
	if(squareDist2D(p, &getCenter() ) < POINT_TOLERANCE_2D*POINT_TOLERANCE_2D)
	{
		p->getX() +=getRadius() ;
		return ;
	}


	Line l(*p, *p-getCenter()) ;

	std::vector<Point> inter = l.intersection(this) ;
// 	if(inter.empty() || inter.size() == 1)
// 	{
// 		p->print() ;
// 		getCenter().print() ;
// 	}
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
		ret.push_back(Point(getRadius()*cos((double)i*angle) + getCenter().getX(), getRadius()*sin((double)i*angle) + getCenter().getY()));
	}

	return ret ;
}

std::vector<Point> Circle::getSamplingBoundingPointsOnArc(size_t num_points, const Point & start, const Point & finish) const
{
	std::vector<Point> ret ;
	if(num_points == 0)
		return ret ;
	Point init(start) ;
	project(&init) ;
	init -= getCenter() ;
	Point fin(finish) ;
	project(&fin) ;
	fin -= getCenter() ;
	double angle = (init.angle() -fin.angle())/num_points ;


	for (double i = 0 ; i< num_points ; i++)
	{
		Point newPoint(init.getX()*cos(i*angle)+init.getY()*sin(i*angle), -init.getX()*sin(i*angle)+init.getY()*cos(i*angle)) ;
		newPoint+= getCenter() ;
		ret.push_back(newPoint);
	}
	return ret ;
}

void Circle::sampleBoundingSurface(size_t num_points)
{
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
		delete boundingPoints[i] ;

	boundingPoints.resize(num_points) ;
	double angle = 2.*M_PI/ (num_points) ;

	for (size_t i = 0 ; i < num_points ; i++)
	{
		double randa= 0;//((2.*(double)rand()/(RAND_MAX+1.0))-1.)*0.15*(M_PI/num_points) ;
		boundingPoints[i] = new Point(getRadius()*cos((double)i*angle) + getCenter().getX(), getRadius()*sin((double)i*angle+randa) + getCenter().getY());
// 		std::cout << "x = " << boundingPoints[i]->getX()() << ", y = " << boundingPoints[i]->getY() << std::endl ;
	}
}

void Circle::sampleSurface(size_t num_points)
{
	if(!sampled)
	{
		num_points = std::max(num_points, (size_t)6) ;
		sampleBoundingSurface(num_points*3/2) ;
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
				temp.push_back(new Point((r+randr)*cos((double)(j+0.5*(i))*angle+randa+offset) + getCenter().getX(), (r+randr)*sin((double)(j+0.5*(i))*angle+randa) + getCenter().getY()));
			}

			num_points = (size_t)(/*std::max(*/(double)num_points_start*(r/getRadius())/*, (double)8)*/) ;

			angle = 2.*M_PI/ (num_points) ;

			offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;

			if(num_points < 5)
				break ;
		}
		for(size_t i = 0 ; i < inPoints.size() ; i++)
			delete inPoints[i] ;

		inPoints.resize(temp.size() + 1) ;
		inPoints[0] = new Point(center) ;
		std::copy(temp.begin(), temp.end(),&inPoints[1]) ;

		//	for(size_t i = 0 ; i < inPoints.size() ; i++)
		//		inPoints[i]->print() ;
		//std::cout << "we have " << num_points << " sample points" << std::endl ;
	}
}

bool Circle::in(const Point & v) const
{
	if(v.getX() < getCenter().getX()-getRadius()*1.01)
		return false ;
	if(v.getX() > getCenter().getX()+getRadius()*1.01)
		return false ;
	if(v.getY() < getCenter().getY()-getRadius()*1.01)
		return false ;
	if(v.getY() > getCenter().getY()+getRadius()*1.01)
		return false ;

	return squareDist2D(v, getCenter()) < sqradius ;
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
	sampleBoundingSurface(num_points) ;
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
			temp.push_back(new Point((r+randr)*cos((double)(j+0.5*(i))*angle+randa+offset) + getCenter().getX(), (r+randr)*sin((double)(j+0.5*(i))*angle+randa) + getCenter().getY()));
		}

		num_points = (size_t)((double)num_points_start*(r/getRadius())) ;

		angle = 2.*M_PI/ (num_points) ;

		offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;
	}

	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;

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
		getBoundingPoint(i).getX() = (getBoundingPoint(i).getX() - center.getX())*ratio + center.getX() ;
		getBoundingPoint(i).getY() = (getBoundingPoint(i).getY() - center.getY())*ratio + center.getY() ;
	}

	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		getInPoint(i).getX() = (getInPoint(i).getX() - center.getX())*ratio + center.getX() ;
		getInPoint(i).getY() = (getInPoint(i).getY() - center.getY())*ratio + center.getY() ;
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

SegmentedLine::SegmentedLine(const Point & p0, const Point & p1) : NonConvexGeometry(0)
{
	gType = SEGMENTED_LINE ;

	this->center = (p0+p1)*.5 ;

	this->boundingPoints.resize(2) ;
	boundingPoints[0] = new Point(p0) ;
	boundingPoints[1] = new Point(p1) ;
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

	p->getX() = projections.begin()->second.getX() ;
	p->getY() = projections.begin()->second.getY() ;
}

void SegmentedLine::sampleSurface(size_t num_points)
{

}

bool SegmentedLine::in(const Point & v) const
{
	return false ;
	for (size_t i = 0 ; i < boundingPoints.size() ; i++)
		if(*boundingPoints[i] == v)
			return true ;

	return false ;
}


Ellipse::Ellipse(Point center, Point a, Point b) : majorAxis(a), minorAxis(b)
{
	gType = ELLIPSE ;
	this->center = center ;
}

Ellipse::Ellipse(Circle c)
{
	gType = ELLIPSE ;
	this->center = c.getCenter() ;
	majorAxis = Point(c.getRadius(),0.) ;
	minorAxis = Point(0.,c.getRadius()) ;
}

Ellipse::Ellipse(Point center, Point a, double b) : majorAxis(a)
{
	gType = ELLIPSE ;
	double b_ = std::min(std::abs(b), 1.) ;
	this->minorAxis = Point(-a.getY()*b_, a.getX()*b_) ;
	this->center = center ;
}


Ellipse::Ellipse(Point center, double a, double b)
{
	gType = ELLIPSE ;
	this->center = center ;
	double dir_x = (double) rand() / (double) RAND_MAX ;
	double dir_y = (double) rand() / (double) RAND_MAX ;
	double a_ = std::max(std::abs(a),std::abs(b)) ;
	double b_ = std::min(std::abs(a),std::abs(b)) ;
	this->majorAxis = Point(a_*dir_x,a_*dir_y) ;
	this->minorAxis = Point(-b_*dir_y,b_*dir_x) ;
}


const Point Ellipse::getFocus(bool dir) const
{
	int dirsign = 1 ;
	if(dir)
	{
		dirsign = - 1 ;
	}
	Point focus(center + majorAxis * (getExcentricity() * dirsign)) ;
	return focus ;
}

const std::pair<Point, Point> Ellipse::getBothFocus() const
{
	std::pair<Point, Point> focus(getFocus(true),getFocus(false)) ;
	return focus ;
}

void Ellipse::computeCenter()
{
}


const Point Ellipse::getPointOnEllipse(double theta) const
{
	Point p(center) ;
	p += getMajorAxis() * std::cos(theta) ;
	p += getMinorAxis() * std::sin(theta) ;
	return p ;
}


Point Ellipse::project(const Point &p) const
{
	if((p-center).norm() < POINT_TOLERANCE_2D*100.)
	{
		return p+getMajorAxis() ;
	} else {
		Point proj(p.getX(), p.getY()) ;
		this->project(&proj) ;
		return proj ;
	}
}

void Ellipse::project(Point * p) const
{
	Point test(p->getX(),p->getY()) ;
	if(test == center)
	{
//		p->getX() += getMinorAxis().getX() ;
//		p->getY() += getMinorAxis().getY() ;
		p->print() ;
		return ;
	}
	else
	{
		double alpha = majorAxis.angle() ;
		Point prot((*p-center).getX()*cos(-alpha)-(*p-center).getY()*sin(-alpha),
				   +(*p-center).getX()*sin(-alpha)+(*p-center).getY()*cos(-alpha)) ;

		Line majL(this->getCenter(), this->getMajorAxis()) ;
		Line minL(this->getCenter(), this->getMinorAxis()) ;

		if(minL.on(test))
		{
			// case point is colinear to minor axis
			Point A = center + getMinorAxis() ;
			Point B = center - getMinorAxis() ;
			Point pa(p->getX()-A.getX(),p->getY()-A.getY())  ;
			Point pb(p->getX()-B.getX(),p->getY()-B.getY())  ;
			if(pa.norm() < pb.norm())
			{
				p->getX() = A.getX() ;
				p->getY() = A.getY() ;
				return ;
			}
			else
			{
				p->getX() = B.getX() ;
				p->getY() = B.getY() ;
				return ;
			}
		}

		if(majL.on(test))
		{
			// case point is colinear to major axis
			Point C = center + getMajorAxis() ;
			Point D = center - getMajorAxis() ;
			Point pc(p->getX()-C.getX(),p->getY()-C.getY())  ;
			Point pd(p->getX()-D.getX(),p->getY()-D.getY())  ;
			if(pc.norm() < pd.norm())
			{
				p->getX() = C.getX() ;
				p->getY() = C.getY() ;
				return ;
			}
			else
			{
				p->getX() = D.getX() ;
				p->getY() = D.getY() ;
				return ;
			}
		}

		// else, we need brute force

		double tmin = 0. ;
		double tmax = M_PI * 2. ;
		int iter = 0 ;
		while(std::abs(tmin-tmax) > POINT_TOLERANCE_2D && iter < 16)
		{
			iter++ ;
			double dtheta = (tmax - tmin) / 32. ;
			Vector allDistances(33) ;
			for(size_t i = 0 ; i < allDistances.size() ; i++)
				allDistances[i] = dist( test, getPointOnEllipse( tmin+dtheta*i )) ;

			double j = 0 ;
			for(size_t i = 0 ; i < allDistances.size() ; i++)
			{
				if(allDistances[i] == allDistances.min())
				{
					j = (double) i ;
					break ;
				}
			}
			tmin = tmin+dtheta*(j-1) ;
			tmax = tmin+dtheta*(j+1) ;
		}

		double found = (tmin+tmax)/2. ;		
//		std::cout << found << "\t" ;

		p->getX() = getPointOnEllipse(found).getX() ;
		p->getY() = getPointOnEllipse(found).getY() ;

	}

}

std::vector<Point> Ellipse::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> ret ;

	double angle = 2. * M_PI / num_points ;

	double thisangle = angle ;
	double lastangle = 0. ;
	double redfactor = 0.8 ; // factor for angle decrease < 1

	ret.push_back(getPointOnEllipse(0.)) ;


	Point thispoint = getPointOnEllipse(angle) ;
	Point lastpoint = getPointOnEllipse(lastangle) ;
	Point lastlastpoint = getPointOnEllipse(-angle) ;
	double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

	int n_iter = 0 ;



	bool accepted = true ;
	while (accepted)
	{
		n_iter = 0 ;
		while(((criteria > angle) || (criteria < angle * redfactor)) && n_iter < 20)
		{
			if(criteria > angle)
			{
				thisangle = lastangle + (thisangle - lastangle) * redfactor ;
			}
			else
			{
				n_iter = 20 ;
			} //thisangle = lastangle + (thisangle - lastangle) / redfactor ; }
			thispoint = getPointOnEllipse(thisangle) ;
			criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
			n_iter++ ;
		}

//		std::cout << thisangle << "-" ;
		if(thisangle < 2.*M_PI)
			ret.push_back(thispoint) ;
		else
			accepted = false ;

		lastlastpoint = lastpoint ;
		lastpoint = thispoint ;

		lastangle = thisangle ;
		thisangle += angle ;

		thispoint = getPointOnEllipse(thisangle) ;
	}

	if(ret.size() > 1)
	{
		
		double div = (ret[0]-ret[1]).norm() ;

		while(ret.size() > 1 && (ret[ret.size()-1]-ret[0]).norm() < 0.15*div)
		{
			ret.pop_back() ;
		}
	}

	return ret ;

}

Function Ellipse::getEllipseFormFunction() const
{
	double alpha = majorAxis.angle() ;
	Function x("x") ;
	Function y("y") ;
	Function x_((x-center.getX())*cos(alpha)-(y-center.getY())*sin(-alpha)) ;
	Function y_((x-center.getX())*sin(-alpha)+(y-center.getY())*cos(alpha)) ;
	double mm0 = (getMajorRadius()*getMajorRadius()) ;
	double mm1 = (getMinorRadius()*getMinorRadius()) ;
	Function f0 = x_*x_/mm0 ;
	Function f1 = y_*y_/mm1 ;
	Function f2 = f0+f1-1. ;
	return f2 ;
}

bool Ellipse::in(const Point &p) const
{
	return VirtualMachine().eval(getEllipseFormFunction(),p) < POINT_TOLERANCE_2D ;
}

std::vector<Point> Ellipse::getSampleBoundingPointsOnArc(size_t num_points, double alpha, double beta) const
{
	std::vector<Point> ret ;

	double angle = (beta-alpha) * M_PI / num_points ;

	double thisangle = alpha + angle ;
	double lastangle = alpha ;
	double redfactor = 0.8 ; // factor for angle decrease < 1


	Point thispoint = getPointOnEllipse(angle) ;
	Point lastpoint = getPointOnEllipse(lastangle) ;
	Point lastlastpoint = getPointOnEllipse(alpha-angle) ;
	double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

	int n_iter = 0 ;

	for (size_t i = 0 ; i< num_points + 1 ; i++)
	{
		n_iter = 0 ;
		while(((criteria > angle) || (criteria < angle * redfactor)) && n_iter < 20)
		{
			if(criteria > angle)
			{
				thisangle = lastangle + (thisangle - lastangle) * redfactor ;
			}
			else
			{
				thisangle = lastangle + (thisangle - lastangle) / redfactor ;
			}
			thispoint = getPointOnEllipse(thisangle) ;
			criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
			n_iter++ ;
		}
		thispoint = getPointOnEllipse(thisangle + RandomNumber().uniform(-angle*0.05, angle*0.05)) ;

		ret.push_back(thispoint) ;

		lastlastpoint = lastpoint ;
		lastpoint = thispoint ;

		lastangle = thisangle ;
		thisangle += angle ;

		thispoint = getPointOnEllipse(thisangle) ;
	}

	double div = (ret[0]-ret[1]).norm() ;

	while((ret[ret.size()-1]-ret[0]).norm() < 0.15*div)
	{
		ret.pop_back() ;
	}

	return ret ;

}

void Ellipse::sampleBoundingSurface (size_t num_points)
{
	if( getMinorAxis().norm() < POINT_TOLERANCE_2D || getMajorAxis().norm() < POINT_TOLERANCE_2D )
		return ;
	
	if(num_points < 3)
		return ;

	std::vector<Point> bound = getSamplingBoundingPoints(num_points) ;

	for(size_t i = 0 ; i <getBoundingPoints().size() ; i++)
		delete boundingPoints[i];

	getBoundingPoints().resize(bound.size()) ;

	for (size_t i = 0 ; i < bound.size() ; i++)
	{
		getBoundingPoints()[i] = new Point(bound[i]) ;
	}

}

void Ellipse::sampleSurface (size_t num_points)
{
	if( getMinorAxis().norm() < POINT_TOLERANCE_2D || getMajorAxis().norm() < POINT_TOLERANCE_2D )
		return ;
	
	if(num_points < 3)
		num_points = 4 ;

	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;

	inPoints.resize(1) ;
	inPoints[0] = new Point(center) ;

	size_t n = num_points*getMajorRadius()/getMinorRadius() ;

	sampleBoundingSurface(n) ;
	sampled = true ;

	double dist = (getBoundingPoint(0)-getBoundingPoint(1)).norm() ;

	size_t ring = (num_points*3 / 5) - 1 ;

	if(getMinorRadius() / getMajorRadius() < 0.5)
		ring-- ;

	if(ring > 1)
	{
		std::vector<double> newa ;
		std::vector<double> newb ;
		std::vector<size_t> newn ;

		double minn = (double) num_points ;
		size_t inc = 1 ;
		if(ring > 4)
			inc++ ;

		for(size_t i = 0 ; i < ring ; i++)
		{
			newa.push_back(getMajorRadius() * (ring - i) / (ring + 1) ) ;
			newb.push_back(getMinorRadius() * (ring - i) / (ring + 1) ) ;
			newn.push_back(  minn + (size_t) ((double) num_points) * ((double) (ring-i)*(ring-i) / (double) ((ring+1)*(ring+1)))) ;
		}
		std::vector<Point> toadd ;
		int rm = 0 ;

		for(size_t i = 0 ; i < newa.size() ; i+=inc)
		{
			Ellipse elln(center,getMajorAxis()*newa[i]/getMajorRadius(),getMinorAxis()*newb[i]/getMinorRadius()) ;

//			int factor = 2 + i/2 ;
			size_t nell = newn[i] ;
			if(nell < minn) { nell = minn ; }

			std::vector<Point> pn = elln.getSamplingBoundingPoints(nell) ;
			for(size_t j = 0 ; j < pn.size() ; j++)
			{
				pn[j] += Point(RandomNumber().uniform(dist*0.1), RandomNumber().uniform(dist*0.1)) ;
				bool alone = in(pn[j]) ;
				if(!alone)
				{
					pn[j].print() ;
					std::cout << "not inside" << std::endl ;
				}

				size_t k = 0 ;
				while(alone && k < getBoundingPoints().size())
				{
					if((pn[j]-getBoundingPoint(k)).norm() < dist)
					{
						alone = false ;
					}
					k++ ;
				}
				
				
				k = 0 ;
				while(alone && k < toadd.size())
				{
					if((pn[j]-toadd[k]).norm() < dist)
					{
						alone = false ;
					}
					k++ ;
				}

				if(alone)
					toadd.push_back(pn[j]) ;
				else
					rm++ ;
			}
		}
		PointArray inTemp(inPoints) ;

		inPoints.resize(inPoints.size()+toadd.size()) ;

		for(size_t j = 0 ; j < inTemp.size() ; j++)
			inPoints[j] = inTemp[j] ;
		for(size_t j = 0 ; j < toadd.size() ; j++)
		{
			inPoints[j+inTemp.size()] = new Point(toadd[j]) ;
		}
	}
	else
	{
		if(ring == 1)
		{
			Ellipse elln(center, getMajorAxis() * 0.5, getMinorAxis() * 0.6666666) ;
			std::vector<Point> pn = elln.getSamplingBoundingPoints(num_points / 3) ;
			PointArray inTemp(inPoints) ;

			inPoints.resize(inPoints.size()+pn.size()) ;

			for(size_t j = 0 ; j < inTemp.size() ; j++)
				inPoints[j] = inTemp[j] ;
			for(size_t j = 0 ; j < pn.size() ; j++)
			{
				inPoints[j+inTemp.size()] = new Point(pn[j]) ;
			}
		}
	}
}

Ellipse Ellipse::getEllipseInLocalCoordinates() const
{
	Point c(center) ;
	Point a(getMajorRadius(), 0.) ;
	Point b(0.,getMinorRadius()) ;

	return Ellipse(c,a,b) ;

}

Point Ellipse::toLocalCoordinates(const Point & p) const
{
	double alpha = majorAxis.angle() ;
	Point prot((p-center).getX()*cos(-alpha)-(p-center).getY()*sin(-alpha),
			   +(p-center).getX()*sin(-alpha)+(p-center).getY()*cos(-alpha)) ;
	return prot ;
}




std::vector<Point> Ellipse::getBoundingBox() const
{
	std::vector<Point> bbox(4) ;

	Point A = center + majorAxis + minorAxis ;
	Point B = center + majorAxis - minorAxis ;
	Point C = center - majorAxis - minorAxis ;
	Point D = center - majorAxis + minorAxis ;

	double minx = std::min(A.getX(), std::min(B.getX(), std::min(C.getX(), D.getX()))) ;
	double maxx = std::max(A.getX(), std::max(B.getX(), std::max(C.getX(), D.getX()))) ;
	double miny = std::min(A.getY(), std::min(B.getY(), std::min(C.getY(), D.getY()))) ;
	double maxy = std::max(A.getY(), std::max(B.getY(), std::max(C.getY(), D.getY()))) ;

	bbox[0] = Point(minx, maxy) ;
	bbox[1] = Point(maxx, maxy) ;
	bbox[2] = Point(maxx, miny) ;
	bbox[3] = Point(minx, miny) ;

	return  bbox ;



	double a = getMajorRadius() ;
	bbox[0] = center + Point(a,a) ;
	bbox[1] = center + Point(a,-a) ;
	bbox[2] = center + Point(-a,-a) ;
	bbox[3] = center + Point(-a,a) ;

	return bbox ;

//	double a = getMajorRadius() ;
	double b = getMinorRadius() ;
	double alpha = getMajorAxis().angle() ;

	double thetax1 = std::acos(std::sqrt((a*a*std::cos(alpha)*std::cos(alpha)) / (a*a*std::cos(alpha)*std::cos(alpha) + b*b*std::sin(alpha)*std::sin(alpha)))) ;
	double thetax2 = std::acos(- std::sqrt((a*a*std::cos(alpha)*std::cos(alpha)) / (a*a*std::cos(alpha)*std::cos(alpha) + b*b*std::sin(alpha)*std::sin(alpha)))) ;
	double thetay1 = std::asin(std::sqrt((b*b*std::cos(alpha)*std::cos(alpha)) / (b*b*std::cos(alpha)*std::cos(alpha) + a*a*std::sin(alpha)*std::sin(alpha)))) ;
	double thetay2 = std::asin(- std::sqrt((b*b*std::cos(alpha)*std::cos(alpha)) / (b*b*std::cos(alpha)*std::cos(alpha) + a*a*std::sin(alpha)*std::sin(alpha)))) ;

	double x1 = a*std::cos(alpha)*std::cos(thetax1) - b*std::sin(alpha)*std::sin(thetax1) ;
	double x2 = a*std::cos(alpha)*std::cos(thetax2) - b*std::sin(alpha)*std::sin(thetax2) ;

	double xmin = std::min(x1,x2) ;
	double xmax = std::max(x1,x2) ;

	double y1 = a*std::sin(alpha)*std::cos(thetay1) + b*std::cos(alpha)*std::sin(thetay1) ;
	double y2 = a*std::sin(alpha)*std::cos(thetay2) + b*std::cos(alpha)*std::sin(thetay2) ;

	double ymin = std::min(y1,y2) ;
	double ymax = std::max(y1,y2) ;

	bbox[0] = center + Point(xmax,ymax) ;
	bbox[1] = center + Point(xmax,ymin) ;
	bbox[2] = center + Point(xmin,ymin) ;
	bbox[3] = center + Point(xmin,ymax) ;

	return bbox ;

}


void Polygon::computeCenter()
{
    
        for(size_t i = 0 ; i < originalPoints.size() ; i++) 
        {
            center += originalPoints[i] ;
        }
        center /= originalPoints.size() ;
}


Polygon::Polygon(const std::valarray<Point *> & points) : NonConvexGeometry(0), originalPoints(points.size())
{
    gType = POLYGON ;

    boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
        boundingPoints[i] = points[i] ;
        originalPoints[i] = *(points[i]) ;
    }
    
    computeCenter();
}

Polygon::~Polygon() { };

void Polygon::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> newPoints = getSamplingBoundingPoints(num_points) ;

    boundingPoints.resize(newPoints.size());
    
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        boundingPoints[i] = new Point(newPoints[i]) ;
    
}

std::vector<Point> Polygon::getSamplingBoundingPoints(size_t num_points) const
{
    double perimeter = 0 ;
    std::vector<Point> ret ;
    
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        perimeter += dist(originalPoints[i], originalPoints[inext]) ;
    }
    
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        double fraction = dist(originalPoints[i], originalPoints[inext])/perimeter ;
        int numPointsOnSegment = std::max(round(fraction*num_points), 2.) ;
        
        ret.push_back(originalPoints[i]);
        for(int j = 1 ; j < numPointsOnSegment-1 ; j++)
        {
            ret.push_back(originalPoints[i]*(double)j/(numPointsOnSegment-1) + originalPoints[inext]*(1.-(double)j/(numPointsOnSegment-1)) );
        }
    }
    return ret ;
}

double Polygon::getPerimeter() const 
{
    double perimeter = 0 ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        perimeter += dist(originalPoints[i], originalPoints[inext]) ;
    }
    return perimeter ;
}

void Polygon::sampleSurface(size_t num_points)
{
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;
    
    num_points *=2 ;
    sampleBoundingSurface(num_points*3);
    std::vector<Segment> segments ;
    double perimeter = 0 ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        segments.push_back(Segment(originalPoints[i], originalPoints[inext]));
        perimeter += dist(originalPoints[i], originalPoints[inext]) ;
    }
    
    double delta = perimeter/num_points ;
    
    
    std::vector<Point> newPoints ;
    double rad = getRadius() ;
    int iteration = 0 ;
    do
    {
        iteration++ ;
        std::vector<Segment> newSegments ;
        for(size_t i = 0 ; i < segments.size() ; i++)
        {
            Point vector(segments[i].normal(center)) ;
     
            vector /= vector.norm() ;
            vector *= 0.75*delta ;
            if(!in(segments[i].midPoint()+vector*iteration))
                vector *= -1 ;
                
            segments[i].set(segments[i].first()+vector , segments[i].second()+vector) ;   
            Point vs = segments[i].second()- segments[i].midPoint() ;
            Point vf = segments[i].first()- segments[i].midPoint() ;
            segments[i].set(segments[i].first()+vf/vf.norm()*2.*delta , segments[i].second()+vs/vs.norm()*2.*delta) ; 
        }
        
        for(size_t i = 0 ; i < segments.size() ; i++)
        {
            std::multimap<double,Point> candidatesf ;
            std::multimap<double,Point> candidatess ;
            for(size_t j = 0 ; j < segments.size() ; j++)
            {
                if(i == j)
                    continue ;
                
                if(segments[i].intersects(segments[j]))
                {
                    Point inter = segments[i].intersection(segments[j]) ;
                    candidatess.insert(std::make_pair(dist(segments[i].second(),inter), inter)) ;
                    candidatesf.insert(std::make_pair(dist(segments[i].first(),inter), inter)) ;
                }
            }
            
            if(candidatess.size() > 1)
            {
                Segment newSeg(candidatesf.begin()->second, candidatess.begin()->second) ;
                double n = newSeg.norm() ;
                if(n > delta)
                {
                    newSegments.push_back(newSeg);
                }
            }
        }
        
        segments.clear();
        if(newSegments.empty())
        {
            newPoints.push_back(center);
            break ;
        }
        int start = newPoints.size() ;
        Segment * testSegment = &newSegments[0] ;
        Segment * startSegment = testSegment ;
        newPoints.push_back(testSegment->second());
        int tries  = 0 ;
        bool found = true ;
        while(tries < newSegments.size() && found)
        {
           found = false ;
            for(size_t i = 0 ; i < newSegments.size() ; i++)
            {
                if(&newSegments[i] == testSegment)
                    continue ;
                
                if(dist(newPoints.back(), newSegments[i].first()) < POINT_TOLERANCE_2D)
                {
                    newPoints.push_back(newSegments[i].second());
                    testSegment = &newSegments[i] ;
                    found = true ;
                    break ;
                }
                if(dist(newPoints.back(), newSegments[i].second()) < POINT_TOLERANCE_2D)
                {
                    newPoints.push_back(newSegments[i].first());
                    testSegment = &newSegments[i] ;
                    found = true ;
                    break ;
                }
            }
            
            if(startSegment == testSegment)
                break ;
            if(!found)
            {
                tries++ ;
                Segment * testSegment = &newSegments[tries] ;
                Segment * startSegment = testSegment ;
//                 newPoints.erase(newPoints.begin()+start, newPoints.end()) ;
                newPoints.push_back(testSegment->second());
            }
 
        }
        

        
        double currentperimeter = 0 ;
        for(size_t i = start ; i < newPoints.size() ; i++ )
        {
            int inext = i+1 ;
            if(inext >= newPoints.size())
                inext = start ;
            segments.push_back(Segment(newPoints[i], newPoints[inext]));
            currentperimeter +=segments.back().norm() ;
        }
        
        if(currentperimeter < delta)
        {
            newPoints.push_back(center);
            break ;
        }
        
        for(size_t i = 0 ; i < segments.size() ; i++ )
        {
            double fraction = segments[i].norm()/currentperimeter ;
            int numPointsOnSegment = round(fraction*num_points*currentperimeter/perimeter) ;
  
            for(int j = 0 ; j < numPointsOnSegment-1 ; j++)
            {
                newPoints.push_back(segments[i].first()*(double)j/(numPointsOnSegment-1) + segments[i].second()*(1.-(double)j/(numPointsOnSegment-1)) );
            }
        }
        
    }while(!segments.empty()) ;
    
    inPoints.resize(newPoints.size());
    for(size_t i = 0 ; i < newPoints.size() ; i++ )
        inPoints[i] = new Point(newPoints[i]) ;
    
}


bool Polygon::in(const Point & v) const
{
    Point out = center + Point(0, getRadius()*1.5) ;
    Segment s(out, v) ;
    
    int interCount = 0 ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        interCount += Segment(originalPoints[i], originalPoints[inext]).intersects(s);
    }
    return interCount%2 ;
}

double Polygon::area() const 
{ 
    std::vector<Point> bb = getBoundingBox() ;
    double maxx = bb[1].getX() ;
    double minx = bb[0].getX() ;
    double miny = bb[1].getY() ;
    double maxy = bb[2].getY() ;
    std::default_random_engine generator;
    std::uniform_real_distribution< double > distributionx(minx, maxx);
    std::uniform_real_distribution< double > distributiony(miny, maxy);
    
    double incount = 0 ;
    for(size_t i = 0 ; i < 1024  ;i++)
    {
        Point test(distributionx(generator), distributiony(generator)) ;
        if(in(test))
            incount++ ;
    }
    return (maxx-minx)*(maxy-miny)*(incount/1024.) ;
    
}

double Polygon::volume() const { return 0 ;} 

void Polygon::project(Point * init) const
{
    std::map<double, Point> potential ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        Segment s(originalPoints[i], originalPoints[inext]);
        Point p = s.project(*init);
        potential[dist(p, *init)] = p ;
    }
    init->set(potential.begin()->second) ;
}

double Polygon::getRadius() const
{
    double r = dist(center, originalPoints[0]) ;
    for(const auto & p : originalPoints)
        r = std::max(r, dist(p, center)) ;
    
    return r ;
}

SpaceDimensionality Polygon::spaceDimensions() const
{
    return SPACE_TWO_DIMENSIONAL ;
}

std::vector<Point> Polygon::getBoundingBox() const 
{ 
    double maxx = originalPoints[0].getX() ;
    double minx = originalPoints[0].getX() ;
    double maxy = originalPoints[0].getY() ;
    double miny = originalPoints[0].getY() ;
    for(size_t i = 1 ; i < originalPoints.size() ; i++ )
    {
        maxx = std::max(maxx,originalPoints[i].getX()) ;
        minx = std::min(minx,originalPoints[i].getX()) ;
        maxy = std::max(maxy,originalPoints[i].getY()) ;
        miny = std::min(miny,originalPoints[i].getY()) ;
    }
    
    return {Point(minx, miny), Point(maxx, miny), Point(maxx, maxy), Point(minx,maxy)} ;
    
}


