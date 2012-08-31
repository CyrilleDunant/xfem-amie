// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_3D.h"
#include "geometry_2D.h"
#include "../utilities/itoa.h"
#include "../utilities/random.h"
#include <fstream>
#include <iomanip>

using namespace Mu ;



RegularOctahedron::RegularOctahedron(double s, double x, double y, double z) : ConvexGeometry(6), length(s)
{
	this->gType = REGULAR_OCTAHEDRON ;
	getCenter().set(x, y, z) ;
	double h = s/std::sqrt(2.) ;
	//the four points of the basis are
	boundingPoints[0] = new Point(x-h, y/*-s*.5*/, z) ;
	boundingPoints[1] = new Point(x/*-s*.5*/, y+h, z) ;
	boundingPoints[2] = new Point(x+h, y/*+s*.5*/, z) ;
	boundingPoints[3] = new Point(x/*+s*.5*/, y-h, z) ;
	//the two points are
	boundingPoints[4] = new Point(x, y, z+h) ;
	boundingPoints[5] = new Point(x, y, z-h) ;

}
	
RegularOctahedron::RegularOctahedron(double s, Point c) : ConvexGeometry(6), length(s)
{
	this->gType = REGULAR_OCTAHEDRON ;
	getCenter().set(c) ;
	double h = s/std::sqrt(2.) ;
	//the four points of the basis are
	boundingPoints[0] = new Point(c.x-h, c.y/*-s*.5*/, c.z) ;
	boundingPoints[1] = new Point(c.x/*-s*.5*/, c.y+h, c.z) ;
	boundingPoints[2] = new Point(c.x+h, c.y/*+s*.5*/, c.z) ;
	boundingPoints[3] = new Point(c.x/*+s*.5*/, c.y-h, c.z) ;
	//the two points are
	boundingPoints[4] = new Point(c.x, c.y, c.z+h) ;
	boundingPoints[5] = new Point(c.x, c.y, c.z-h) ;
}

	
void RegularOctahedron::sampleBoundingSurface(size_t num_points)
{
	std::vector<Point> samplingPoints = getSamplingBoundingPoints( num_points) ;
	
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		delete boundingPoints[i] ;
	}

	random_shuffle(samplingPoints.begin(), samplingPoints.end()) ;
	
	this->boundingPoints.resize(samplingPoints.size()) ;
	for(size_t i = 0 ; i < samplingPoints.size() ; i++)
		boundingPoints[i] = new Point(samplingPoints[i]) ;
}

std::vector<Point> RegularOctahedron::getSamplingBoundingPoints(size_t num_points) const
{
	size_t pointsOnEquator = (size_t)round(pow(num_points, 0.33333333)*6.) ;
	Rectangle(length, length, getCenter().x, getCenter().y) ;
	std::vector<Point> samplingPoints ;
	std::vector<Point> newPoints = Rectangle(length, length, getCenter().x, getCenter().y).getSamplingBoundingPoints(pointsOnEquator) ;
	size_t realPointsOnEquator = newPoints.size() ;

	Matrix rot(3,3) ;
	rot[0][0] = cos(M_PI*.25) ;rot[0][1] = sin(M_PI*.25) ;
	rot[1][0] = -sin(M_PI*.25) ;rot[1][1] = cos(M_PI*.25) ;
	rot[2][2] = 1 ;
	for(size_t i = 0 ; i < newPoints.size() ;i++)
	{
		newPoints[i] -= Point(center.x, center.y) ;
		newPoints[i] *= rot ;
		newPoints[i] += Point(center.x, center.y) ;
		newPoints[i].z = center.z ; 
	}
	
	samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;

	int iterations = 1 ;
	while(true)
	{
		double factor = (double)(realPointsOnEquator-2.*iterations)/(double)realPointsOnEquator ;
		if(factor < 0)
			break ;
		iterations++ ;
	}
	double totalIterations = iterations ;

	iterations = 1 ;
	
	double sq2 = std::sqrt(2.) ;
	while(newPoints.size() >= 4)
	{
		double factor = (totalIterations-iterations)/totalIterations ;
		if(factor <= 0)
			break ;
		newPoints = Rectangle(length*factor, length*factor, getCenter().x, getCenter().y).getSamplingBoundingPoints(realPointsOnEquator-2.*iterations) ;
		if(newPoints.size() >= 4)
		{
			for(size_t i = 0 ; i < newPoints.size() ;i++)
			{
				newPoints[i] -= Point(center.x, center.y) ;
				newPoints[i] *= rot ;
				newPoints[i] += Point(center.x, center.y) ;
				newPoints[i].z = center.z + (1.-factor)*length/sq2 ; 
			}
			samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;
		
			for(size_t i = 0 ; i < newPoints.size() ;i++)
			{
				newPoints[i].z = center.z - (1.-factor)*length/sq2; 
			}
			samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;
		
			iterations++ ;
		}
		else
			break ;
	}
	samplingPoints.push_back(*boundingPoints[4]) ;
	samplingPoints.push_back(*boundingPoints[5]) ;
	
	return samplingPoints ;
}

void RegularOctahedron::sampleSurface(size_t num_points)
{

	sampleBoundingSurface(num_points) ;
	
	int realPointsOnEquator = Rectangle(length, length, getCenter().x, getCenter().y).getSamplingBoundingPoints(round(sqrt(num_points))).size() ;
	int iterations = 1 ;
	while(true)
	{
		double factor = (double)(realPointsOnEquator-1.*iterations)/(double)realPointsOnEquator ;
		if(factor < 0)
			break ;
		iterations++ ;
	}
	
	size_t numberOfShells = round(iterations) ;
	iterations = 1 ;
	std::vector<Point> samplingPoints ;

	for(size_t i = 0 ; i < numberOfShells ; i++)
	{
		double factor = (double)(numberOfShells-1.*iterations)/(double)numberOfShells ;
		if(factor <= 0)
			break ;
		std::vector<Point> newPoints = RegularOctahedron(length*factor, getCenter()).getSamplingBoundingPoints(num_points*factor) ;


		if(newPoints.size() <= 6)
			break ;
		samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;
		
		iterations++ ;
	}
	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;
		
	random_shuffle(samplingPoints.begin(), samplingPoints.end()) ;
	
	inPoints.resize(samplingPoints.size()) ;
	for(size_t i = 0 ; i < samplingPoints.size() ; i++)
		inPoints[i] = new Point(samplingPoints[i]) ;
	
}

bool RegularOctahedron::in(const Point &p) const
{

	Point local = p-center ;
	double test = std::abs(local.x) + std::abs(local.y) + std::abs(local.z) ;
	return (test-length/std::sqrt(2.)) < POINT_TOLERANCE_3D ;

// 	return true ;
	Matrix rot(3,3) ;
	rot[0][0] = cos(-M_PI*.25) ;rot[0][1] = sin(-M_PI*.25) ;
	rot[1][0] = -sin(-M_PI*.25) ;rot[1][1] = cos(-M_PI*.25) ;
	rot[2][2] = 1 ;

	Point v = p ;//*rot ;
	v -=center ;
	v*=rot ;
	v+=center ;
	if(abs(v.x-center.x) > length)
		return false ;
	if(abs(v.y-center.y) > length)
		return false ;
	Point v_(v.x, v.z) ;
	if(!Triangle(Point(getCenter().x-.5*length, center.z), Point(getCenter().x+.5*length, getCenter().z), Point(getCenter().x, getCenter().z+0.70711*length)).in(v_) && !Triangle(Point(getCenter().x-.5*length, getCenter().z), Point(getCenter().x+.5*length, getCenter().z), Point(getCenter().x, getCenter().z-0.70711*length)).in(v_))
		return false ;
	
	v_.set(v.y, v.z) ;
	if(!Triangle(Point(getCenter().y-.5*length, getCenter().z), Point(getCenter().y+.5*length, getCenter().z), Point(getCenter().y, getCenter().z+0.70711*length)).in(v_) && !Triangle(Point(getCenter().y-.5*length, getCenter().z), Point(getCenter().y+.5*length, getCenter().z), Point(getCenter().y, getCenter().z-0.70711*length)).in(v_))
		return false ;
	
	return true;
	
}
	
double RegularOctahedron::area() const
{
	return 2.*sqrt(3.)*length*length ;
}

double RegularOctahedron::volume() const
{
	return sqrt(2.)/3.*length*length*length ;
}
		
void RegularOctahedron::project(Point * p) const
{
}

void RegularOctahedron::computeCenter()
{
}

double RegularOctahedron::getRadius() const 
{
	return .5*sqrt(2)*length ;
}
	
	
std::vector<Point> RegularOctahedron::getBoundingBox() const
{
	std::vector<Point> ret ;
	ret.push_back(Point(getCenter().x+0.5*length, getCenter().y+0.5*length, getCenter().z+0.5*length)) ;
	ret.push_back(Point(getCenter().x+0.5*length, getCenter().y+0.5*length, getCenter().z-0.5*length)) ;
	ret.push_back(Point(getCenter().x+0.5*length, getCenter().y-0.5*length, getCenter().z+0.5*length)) ;
	ret.push_back(Point(getCenter().x+0.5*length, getCenter().y-0.5*length, getCenter().z-0.5*length)) ;
	ret.push_back(Point(getCenter().x-0.5*length, getCenter().y+0.5*length, getCenter().z+0.5*length)) ;
	ret.push_back(Point(getCenter().x-0.5*length, getCenter().y+0.5*length, getCenter().z-0.5*length)) ;
	ret.push_back(Point(getCenter().x-0.5*length, getCenter().y-0.5*length, getCenter().z+0.5*length)) ;
	ret.push_back(Point(getCenter().x-0.5*length, getCenter().y-0.5*length, getCenter().z-0.5*length)) ;
	return ret ;
}
	


Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;
	boundingPoints[3] = p3 ;
	
	if(volume() < 0 )
	{
		boundingPoints[0] = p1; 
		boundingPoints[1] = p0;
	}
	
	computeCircumCenter() ;
	Vector r(4) ;
	r[0] = squareDist3D(*p0, circumCenter) ;
	r[1] = squareDist3D(*p1, circumCenter) ;
	r[2] = squareDist3D(*p2, circumCenter) ;
	r[3] = squareDist3D(*p3, circumCenter) ;
	std::sort(&r[0], &r[4]) ;
	sqradius = r[0] ;
	radius = sqrt(sqradius);
	
	assert(volume() >0 ) ;

	computeCenter() ;	
}

Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7): ConvexGeometry(8)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;
	boundingPoints[3] = p3 ;
	boundingPoints[4] = p4 ;
	boundingPoints[5] = p5 ;
	boundingPoints[6] = p6 ;
	boundingPoints[7] = p7 ;
	
	if(this->volume() < 0 )
	{
		boundingPoints[0] = p1; 
		boundingPoints[1] = p0;
	}
	computeCircumCenter() ;
	Vector r(4) ;
	r[0] = squareDist3D(*p0, circumCenter) ;
	r[1] = squareDist3D(*p1, circumCenter) ;
	r[2] = squareDist3D(*p2, circumCenter) ;
	r[3] = squareDist3D(*p3, circumCenter) ;
	std::sort(&r[0], &r[4]) ;
	sqradius = r[0] ;
	radius = sqrt(sqradius);
	assert(this->volume() >0 ) ;
	computeCenter() ;	
}

Tetrahedron::Tetrahedron(const Point &p0, const Point &p1, const Point &p2, const Point &p3): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[0] = new Point(p0) ;
	boundingPoints[1] = new Point(p1) ;
	boundingPoints[2] = new Point(p2) ;
	boundingPoints[3] = new Point(p3) ;
	if(this->volume() < 0 )
	{
		std::swap(boundingPoints[0], boundingPoints[1]) ;
	}
	computeCircumCenter() ;
	Vector r(4) ;
	r[0] = squareDist3D(p0, circumCenter) ;
	r[1] = squareDist3D(p1, circumCenter) ;
	r[2] = squareDist3D(p2, circumCenter) ;
	r[3] = squareDist3D(p3, circumCenter) ;
	std::sort(&r[0], &r[4]) ;
	sqradius = r[0] ;
	radius = sqrt(sqradius);
	computeCenter() ;	
}

Tetrahedron::Tetrahedron(): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[2] = new Point(1,0,0) ;
	boundingPoints[3] = new Point(0,1,0) ;
	boundingPoints[0] = new Point(0,0,1) ;
	boundingPoints[1] = new Point(0,0,0) ;
	if(this->volume() < 0 )
	{
		std::swap(boundingPoints[0], boundingPoints[1]) ;
	}
	computeCircumCenter() ;
	Vector r(4) ;
	r[0] = squareDist3D(*boundingPoints[0], circumCenter) ;
	r[1] = squareDist3D(*boundingPoints[1], circumCenter) ;
	r[2] = squareDist3D(*boundingPoints[2], circumCenter) ;
	r[3] = squareDist3D(*boundingPoints[3], circumCenter) ;
	std::sort(&r[0], &r[4]) ;
	sqradius = r[0] ;
	radius = sqrt(sqradius);
	computeCenter() ;	
}

XMLTree * Tetrahedron::toXML()
{
	XMLTree * tet = new XMLTree("tetrahedron") ;
	tet->addChild(boundingPoints[0]->toXML()) ;
	tet->addChild(boundingPoints[1]->toXML()) ;
	tet->addChild(boundingPoints[2]->toXML()) ;
	tet->addChild(boundingPoints[3]->toXML()) ;
	return tet ;	
}

std::vector<Point> Tetrahedron::getSamplingBoundingPoints(size_t num_points) const
{
	return std::vector<Point>() ;
}

void Tetrahedron::sampleBoundingSurface(size_t num_points)
{
	//! \todo make it do something
	return ;
}

void Tetrahedron::computeCenter()
{
	this->center = Point(0,0,0) ;
	for(size_t i = 0  ; i < boundingPoints.size() ; i++)
	{
		this->center += (*boundingPoints[i])/boundingPoints.size() ;
	}
}

double Tetrahedron::getRadius() const
{
	return radius ;
}

std::vector<Point> Tetrahedron::getBoundingBox() const
{
	
	double min_x_0 = 0, min_y_0 = 0, max_x_0 = 0, max_y_0 = 0, max_z_0 = 0, min_z_0 = 0;
	min_y_0 = getBoundingPoint(0).y ;
	max_y_0 = getBoundingPoint(0).y ;
	min_x_0 = getBoundingPoint(0).x ;
	max_x_0 = getBoundingPoint(0).x ;
	min_z_0 = getBoundingPoint(0).z ;
	max_z_0 = getBoundingPoint(0).z ;
	
	for(size_t k  =  1 ; k <  getBoundingPoints().size() ; k++)
	{
		if(getBoundingPoint(k).y < min_y_0)
			min_y_0 = getBoundingPoint(k).y ;
		if(getBoundingPoint(k).y > max_y_0)
			max_y_0 = getBoundingPoint(k).y ;
		
		if(getBoundingPoint(k).x < min_x_0)
			min_x_0 = getBoundingPoint(k).x ;
		if(getBoundingPoint(k).x > max_x_0)
			max_x_0 = getBoundingPoint(k).x ;
		
		if(getBoundingPoint(k).z < min_z_0)
			min_z_0 = getBoundingPoint(k).z ;
		if(getBoundingPoint(k).z > max_z_0)
			max_z_0 = getBoundingPoint(k).z ;
	}
				
	
	double r = getRadius() ;
	std::vector<Point> ret ;
	ret.push_back(Point(max_x_0, max_y_0, max_z_0)) ;
	ret.push_back(Point(max_x_0, max_y_0, min_z_0)) ;
	ret.push_back(Point(max_x_0, min_y_0, max_z_0)) ;
	ret.push_back(Point(max_x_0, min_y_0, min_z_0)) ;
	ret.push_back(Point(min_x_0, max_y_0, max_z_0)) ;
	ret.push_back(Point(min_x_0, max_y_0, min_z_0)) ;
	ret.push_back(Point(min_x_0, min_y_0, max_z_0)) ;
	ret.push_back(Point(min_x_0, min_y_0, min_z_0)) ;
	return ret ;
}

const Mu::Point& Tetrahedron::getCircumCenter() const
{
	return this->circumCenter ;
}

void Tetrahedron::sampleSurface(size_t num_points)
{
	//! \todo make it do something
	return ;
}

bool Tetrahedron::in(const Point & p) const
{
	if(!inCircumSphere(p))
		return false ;
	
	Point v_ ;
	
	if(boundingPoints.size() == 4)
	{
		Matrix S(4,4) ;
		S[0][0] = getBoundingPoint(2).x;
		S[0][1] = getBoundingPoint(3).x; 
		S[0][2] = getBoundingPoint(0).x ; 
		S[0][3] = getBoundingPoint(1).x ;  
		
		S[1][0] = getBoundingPoint(2).y ;
		S[1][1] = getBoundingPoint(3).y;
		S[1][2] = getBoundingPoint(0).y ; 
		S[1][3] = getBoundingPoint(1).y ;  

		S[2][0] = getBoundingPoint(2).z ;
		S[2][1] = getBoundingPoint(3).z;
		S[2][2] = getBoundingPoint(0).z ; 
		S[2][3] = getBoundingPoint(1).z ;  
		
		S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]= 1;
		
		Vector v(4) ; 
		v[0] = p.x ;
		v[1] = p.y ;
		v[2] = p.z ;
		v[3] = 1 ;

		Vector coeff = inverse4x4Matrix(S) * v ;
		
		v_ = Point(coeff[0],coeff[1],coeff[2], p.t); 
	}
	else
	{
		Matrix S(4,4) ;
		S[0][0] = getBoundingPoint(4).x;
		S[0][1] = getBoundingPoint(6).x; 
		S[0][2] = getBoundingPoint(0).x ; 
		S[0][3] = getBoundingPoint(2).x ;  
		
		S[1][0] = getBoundingPoint(4).y ;
		S[1][1] = getBoundingPoint(6).y;
		S[1][2] = getBoundingPoint(0).y ; 
		S[1][3] = getBoundingPoint(2).y ;  

		S[2][0] = getBoundingPoint(4).z ;
		S[2][1] = getBoundingPoint(6).z;
		S[2][2] = getBoundingPoint(0).z ; 
		S[2][3] = getBoundingPoint(2).z ;  
		
		S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]= 1;
		
		Vector v(4) ; 
		v[0] = p.x ;
		v[1] = p.y ;
		v[2] = p.z ;
		v[3] = 1 ;

		Vector coeff = inverse4x4Matrix(S) * v ;

		v_ = Point(coeff[0],coeff[1],coeff[2],p.t ); 
	}
	if(v_.x < -POINT_TOLERANCE_3D)
		return false ;
	if(v_.y < -POINT_TOLERANCE_3D)
		return false ;
	if(v_.z < -POINT_TOLERANCE_3D)
		return false ;
	if(v_.x+v_.y+v_.z > 1+3.*POINT_TOLERANCE_3D)
		return false ;
	return true ;
	
	
// 	if(dynamic_cast<const TetrahedralElement *>(this))
// 	{
// 		Point v_ = static_cast<const TetrahedralElement *>(this)->inLocalCoordinates(v) ;
// 		
// 		
// 		
// 		
// 		if(v.x < 0)
// 			return false ;
// 		if(v.y < 0)
// 			return false ;
// 		if(v.z < 0)
// 			return false ;
// 		if(v.x+v.y+v.z > 1)
// 			return false ;
// 		return true ;
// 	}
// 	Point  pg = getCenter() ;//(getBoundingPoint(0)+getBoundingPoint(1)+getBoundingPoint(2)+getBoundingPoint(3))/4;
// 	if(squareDist3D(pg, v) <POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
// 		return true ;
// 	
// 	if(getBoundingPoints().size() == 4)
// 	{
// 		TriPoint t0(&getBoundingPoint(0),&getBoundingPoint(1),&getBoundingPoint(2)) ;
// 		TriPoint t1(&getBoundingPoint(0),&getBoundingPoint(1),&getBoundingPoint(3)) ;
// 		TriPoint t2(&getBoundingPoint(0),&getBoundingPoint(2),&getBoundingPoint(3)) ;
// 		TriPoint t3(&getBoundingPoint(1),&getBoundingPoint(2),&getBoundingPoint(3)) ;
// 		Segment s(pg,v) ;
// 
// 		return !(s.intersects(t0) || s.intersects(t1) || s.intersects(t2) || s.intersects(t3)) || 
// 		(isCoplanar(v, getBoundingPoint(0), getBoundingPoint(1), getBoundingPoint(2))) || 
// 		(isCoplanar(v, getBoundingPoint(0), getBoundingPoint(1), getBoundingPoint(3))) ||
// 		(isCoplanar(v, getBoundingPoint(0), getBoundingPoint(2), getBoundingPoint(3))) ||
// 		(isCoplanar(v, getBoundingPoint(1), getBoundingPoint(2), getBoundingPoint(3))) ;
// 	}
// 	else
// 	{
// 		std::multimap<double, const Point *> pts ;
// 		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
// 		{
// 			pts.insert(std::make_pair(std::abs(getRadius()*getRadius()-squareDist3D(getCircumCenter(), getBoundingPoint(i))), & getBoundingPoint(i))) ;
// 		}
// 		auto p = pts.begin() ;
// 		const Point * a = p->second ; ++p ;
// 		const Point * b = p->second ; ++p ;
// 		const Point * c = p->second ; ++p ;
// 		const Point * d = p->second ;
// 		
// 		TriPoint t0(a,b,c) ;
// 		TriPoint t1(a,b,d) ;
// 		TriPoint t2(a,c,d) ;
// 		TriPoint t3(b,c,d) ;
// 		Segment s(pg,v) ;
// 
// 		
// 		return !(s.intersects(t0) || s.intersects(t1) || s.intersects(t2) || s.intersects(t3)) || 
// 		(isCoplanar(v, *a, *b, *c)) || 
// 		(isCoplanar(v, *a, *b, *d)) ||
// 		(isCoplanar(v, *a, *c, *d)) ||
// 		(isCoplanar(v, *b, *c, *d)) ;
// 	}
	
// 	double alpha;
// 	alpha =  ((getBoundingPoint(0))*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(1)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(0))))/((v-pg)*((getBoundingPoint(0))^(getBoundingPoint(1))));
// 	if (alpha>=-1 && alpha<=1) 
// 		return true;
// 	
// 	alpha = ((getBoundingPoint(1))*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(1))))/((v-pg)*((getBoundingPoint(1))^(getBoundingPoint(2))));
// 	if (alpha>=-1 && alpha<=1) 
// 		return true;
// 	
// 	alpha =  ((getBoundingPoint(2))*((getBoundingPoint(3))^(getBoundingPoint(0)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(0)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(2))))/((v-pg)*((getBoundingPoint(2))^(getBoundingPoint(3))));
// 	if (alpha>=-1 && alpha<=1) 
// 		return true;
// 	
// 	alpha =  ((getBoundingPoint(0))*((getBoundingPoint(1))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(1)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(0))))/((v-pg)*((getBoundingPoint(0))^(getBoundingPoint(1))));
// 	if (alpha>=-1 && alpha<=1) 
// 		return true;
// 	
// 	return false ;
}

double Tetrahedron::area() const
{
	if(this->getBoundingPoints().size() == 4)
	{
	Segment s0(getBoundingPoint(1), getBoundingPoint(0)) ;
	Segment s1(getBoundingPoint(1), getBoundingPoint(2)) ;
	Segment s2(getBoundingPoint(1), getBoundingPoint(3)) ;
	
	return 0.5*(((s0.vector()^s1.vector())).norm()+
	             ((s0.vector()^s2.vector())).norm()+
	             ((s1.vector()^s2.vector())).norm()+
	             (((s1.vector()-s0.vector())^s2.vector())-s0.vector()).norm());
	}
	else
	{
		Segment s0(getBoundingPoint(2), getBoundingPoint(0)) ;
		Segment s1(getBoundingPoint(2), getBoundingPoint(4)) ;
		Segment s2(getBoundingPoint(2), getBoundingPoint(6)) ;
		return 0.5*(((s0.vector()^s1.vector())).norm()+
	             ((s0.vector()^s2.vector())).norm()+
	             ((s1.vector()^s2.vector())).norm()+
	             (((s1.vector()-s0.vector())^s2.vector())-s0.vector()).norm());
	}
}

double Tetrahedron::volume() const
{
	if(this->getBoundingPoints().size() == 4 || timePlanes() > 1)
	{
		Segment s0(getBoundingPoint(1), getBoundingPoint(0)) ;
		Segment s1(getBoundingPoint(1), getBoundingPoint(2)) ;
		Segment s2(getBoundingPoint(1), getBoundingPoint(3)) ;
		
		return ((s1.vector())^(s0.vector()))*(s2.vector())/6. ;
	}
	else
	{
		Segment s0(getBoundingPoint(2), getBoundingPoint(0)) ;
		Segment s1(getBoundingPoint(2), getBoundingPoint(4)) ;
		Segment s2(getBoundingPoint(2), getBoundingPoint(6)) ;
		
		return ((s1.vector())^(s0.vector()))*(s2.vector())/6. ;
	}
}

void Tetrahedron::project(Point * p) const
{
//! \todo make it do something
	
	return ;
}

void Tetrahedron::computeCircumCenter()
{
	//Pick three points from a face
	Point * A = boundingPoints[0] ;
	Point * B = boundingPoints[1] ;
	Point * C = boundingPoints[2] ;
	//Pick last point from a face
	Point * A_ = boundingPoints[3] ;
	Point BC = (*C-*B) ;
	
	Point BC_mid = (*B+*C)*.5 ;
	
	Point BA = (*A-*B) ;
	Point BA_mid = (*B+*A)*.5 ;
	
	Point BA_ = (*A_-*B) ;
	Point BA__mid = (*B+*A_)*.5 ;
	
	Matrix planes(3,3) ;
	planes[0][0] = BC.x ; planes[0][1] = BC.y ; planes[0][2] = BC.z ; 
	planes[1][0] = BA.x ; planes[1][1] = BA.y ; planes[1][2] = BA.z ; 
	planes[2][0] = BA_.x ; planes[2][1] = BA_.y ; planes[2][2] = BA_.z ; 
	
	Vector coord(3) ;
	coord[0] = - BC.x*BC_mid.x  - BC.y*BC_mid.y  - BC.z*BC_mid.z ; 
	coord[1] = - BA.x*BA_mid.x  - BA.y*BA_mid.y  - BA.z*BA_mid.z ; 
	coord[2] = - BA_.x*BA__mid.x  - BA_.y*BA__mid.y  - BA_.z*BA__mid.z ; 
	invert3x3Matrix( planes) ;
	coord= planes*coord ;
	
	circumCenter = Point(-coord[0], -coord[1], -coord[2]) ;
	  
	
}

bool Tetrahedron::inCircumSphere(const Point & p) const
{
	return inCircumSphere(&p) ;
}

bool Tetrahedron::inCircumSphere(const Point *p) const
{
	if(p->x > circumCenter.x+1.01*radius)
		return false ;
	if(p->x < circumCenter.x-1.01*radius)
		return false ;
	if(p->y > circumCenter.y+1.01*radius)
		return false ;
	if(p->y < circumCenter.y-1.01*radius)
		return false ;
	if(p->z > circumCenter.z+1.01*radius)
		return false ;
	if(p->z < circumCenter.z-1.01*radius)
		return false ;
	
	double d = dist(&circumCenter, p) ;

	return  d-radius < POINT_TOLERANCE_3D ;
}

Hexahedron::Hexahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7)
{
	gType = HEXAHEDRON;
	
	std::vector<Point *> pts ;
	pts.push_back(p0) ;
	pts.push_back(p1) ;
	pts.push_back(p2) ;
	pts.push_back(p3) ;
	pts.push_back(p4) ;
	pts.push_back(p5) ;
	pts.push_back(p6) ;
	pts.push_back(p7) ;
	
// 	std::stable_sort(pts.begin(), pts.end(), PointLessThan()) ;
	
	this->boundingPoints.resize(8) ;
	for(size_t i = 0 ;  i < 8 ; i++)
	{
		boundingPoints[i] = pts[i] ;
	}
	
	this->center = Point(*p0 + *p1 + *p2 + *p3 + *p4+ *p5+ *p6+ *p7)*0.125;
	this->size_x = std::abs(boundingPoints[7]->x - boundingPoints[0]->x);
	this->size_y = std::abs(boundingPoints[7]->y - boundingPoints[0]->y);
	this->size_z = std::abs(boundingPoints[7]->z - boundingPoints[0]->z);
	
}

XMLTree * Hexahedron::toXML()
{
	XMLTree * hex = new XMLTree("hexahedron") ;
	hex->addChild(boundingPoints[0]->toXML()) ;
	hex->addChild(boundingPoints[1]->toXML()) ;
	hex->addChild(boundingPoints[2]->toXML()) ;
	hex->addChild(boundingPoints[3]->toXML()) ;
	hex->addChild(boundingPoints[4]->toXML()) ;
	hex->addChild(boundingPoints[5]->toXML()) ;
	hex->addChild(boundingPoints[6]->toXML()) ;
	hex->addChild(boundingPoints[7]->toXML()) ;
	return hex ;	
}


double Hexahedron::getXSize() const 
{
	return this->size_x ;
}

double Hexahedron::getYSize() const
{
	return this->size_y ;
}

double Hexahedron::getZSize() const
{
	return this->size_z ;
}

Hexahedron::Hexahedron(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5, Point p6, Point p7)
{
	gType =HEXAHEDRON  ;
	std::vector<Point> pts ;
	pts.push_back(p0) ;
	pts.push_back(p1) ;
	pts.push_back(p2) ;
	pts.push_back(p3) ;
	pts.push_back(p4) ;
	pts.push_back(p5) ;
	pts.push_back(p6) ;
	pts.push_back(p7) ;
	
	std::stable_sort(pts.begin(), pts.end()) ;
	
	this->boundingPoints.resize(8) ;
	
	for(size_t i = 0 ;  i < 8 ; i++)
	{

// 		pts[i].print() ;std::cout << std::endl ;
		boundingPoints[i] = new Point(pts[i]) ;
	}

	
	this->center =  Point(p0 + p1 + p2 + p3 + p4+ p5+ p6+ p7)*0.125;
	
	size_x = std::abs(boundingPoints[6]->x - boundingPoints[0]->x);
	size_y = std::abs(boundingPoints[6]->y - boundingPoints[0]->y);
	size_z = std::abs(boundingPoints[6]->z - boundingPoints[0]->z);
	
}

Hexahedron::Hexahedron(double x, double y, double z, double originX, double originY, double originZ) {
	gType =HEXAHEDRON  ;
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = new Point(originX -0.5*x, originY-0.5*y, originZ-0.5*z) ;
	boundingPoints[1] = new Point(originX -0.5*x, originY-0.5*y, originZ+0.5*z) ;
	boundingPoints[2] = new Point(originX -0.5*x, originY+0.5*y, originZ-0.5*z) ;
	boundingPoints[3] = new Point(originX -0.5*x, originY+0.5*y, originZ+0.5*z) ;
	boundingPoints[4] = new Point(originX +0.5*x, originY-0.5*y, originZ-0.5*z) ;
	boundingPoints[5] = new Point(originX +0.5*x, originY-0.5*y, originZ+0.5*z) ;
	boundingPoints[6] = new Point(originX +0.5*x, originY+0.5*y, originZ-0.5*z) ;
	boundingPoints[7] = new Point(originX +0.5*x, originY+0.5*y, originZ+0.5*z) ;
	
	this->center = Point(originX, originY, originZ ) ;
	
	size_x = x ;
	size_y = y ;
	size_z = z ;
}

Hexahedron::Hexahedron(double x, double y, double z, const Point & c)
{
	gType =HEXAHEDRON  ;
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = new Point(c.x -0.5*x, c.y-0.5*y, c.z-0.5*z) ;
	boundingPoints[1] = new Point(c.x -0.5*x, c.y-0.5*y, c.z+0.5*z) ;
	boundingPoints[2] = new Point(c.x -0.5*x, c.y+0.5*y, c.z-0.5*z) ;
	boundingPoints[3] = new Point(c.x -0.5*x, c.y+0.5*y, c.z+0.5*z) ;
	boundingPoints[4] = new Point(c.x +0.5*x, c.y-0.5*y, c.z-0.5*z) ;
	boundingPoints[5] = new Point(c.x +0.5*x, c.y-0.5*y, c.z+0.5*z) ;
	boundingPoints[6] = new Point(c.x +0.5*x, c.y+0.5*y, c.z-0.5*z) ;
	boundingPoints[7] = new Point(c.x +0.5*x, c.y+0.5*y, c.z+0.5*z) ;
	
	this->center = c ;
	
	size_x = x ;
	size_y = y ;
	size_z = z ;
}

Hexahedron::Hexahedron()
{
	gType = HEXAHEDRON ;
	
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = new Point(-1,-1,-1) ;
	boundingPoints[1] = new Point(-1,-1,1) ;
	boundingPoints[2] = new Point(-1,1,-1) ;
	boundingPoints[3] = new Point(-1,1,1) ;
	boundingPoints[4] = new Point(1,-1,-1) ;
	boundingPoints[5] = new Point(1,-1,1) ;
	boundingPoints[6] = new Point(1,1,-1) ;
	boundingPoints[7] = new Point(1,1,1) ;
	
	this->center =  Point(*boundingPoints[0] + 
	*boundingPoints[1] +
	*boundingPoints[2] +
	*boundingPoints[3] +
	*boundingPoints[4] +
	*boundingPoints[5] +
	*boundingPoints[6] +
	*boundingPoints[7] )*0.125;
	
	size_x = 2 ;
	size_y = 2 ;
	size_z = 2 ;
}

std::vector<Point> Hexahedron::getSamplingBoundingPoints(size_t num_points) const
{
	std::vector<Point> points ;
	
	Point point000(*boundingPoints[0]) ;
	
	Point point111(*boundingPoints[7]) ;
	
	size_t numPointsPerDirection = std::max(static_cast<size_t>(ceil(pow((double)num_points, 1./3.))),(size_t)2) ;
	
	std::vector<double> ds ;
	
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;
		
		ds.push_back(v) ;
	}
	
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point000.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [000] [010] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			Point candidate(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z) ;
			bool in = false ;
			for(size_t k = 0 ; k < points.size() ; k++)
			{
				if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					in = true ;
					break ;
				}
			}
			if(!in)
				points.push_back(candidate) ;
		}
	}
	
	//face [000] [001] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{	
			Point candidate(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j])) ;
			bool in = false ;
			for(size_t k = 0 ; k < points.size() ; k++)
			{
				if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					in = true ;
					break ;
				}
			}
			if(!in)
				points.push_back(candidate) ;
		}
	}
	
	//face [001] [011] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			Point candidate(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z) ;
			bool in = false ;
			for(size_t k = 0 ; k < points.size() ; k++)
			{
				if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					in = true ;
					break ;
				}
			}
			if(!in)
				points.push_back(candidate) ;
		}
	}
	
	//face [100] [110] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			Point candidate(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j])) ;
			bool in = false ;
			for(size_t k = 0 ; k < points.size() ; k++)
			{
				if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					in = true ;
					break ;
				}
			}
			if(!in)
				points.push_back(candidate) ;
		}
	}
	
	//face [010] [011] [110]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			Point candidate(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j])) ;
			bool in = false ;
			for(size_t k = 0 ; k < points.size() ; k++)
			{
				if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					in = true ;
					break ;
				}
			}
			if(!in)
				points.push_back(candidate) ;
		}
	}
		
	return points ;
}

void Hexahedron::sampleBoundingSurface(size_t num_points)
{
	std::vector<Point> points = getSamplingBoundingPoints(4*num_points);
	for (size_t i = 0 ; i < boundingPoints.size() ; i++)
		delete boundingPoints[i] ;

	boundingPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		boundingPoints[i] = new Point(points[i]) ;
	}
	
}

void Hexahedron::computeCenter()
{
	for(size_t i = 0  ; i < boundingPoints.size() ; i++)
	{
		this->center += (*boundingPoints[i])/boundingPoints.size() ;
	}
}

double Hexahedron::getRadius() const
{
	return dist(getCenter(), getCenter()+Point(size_x, size_y, size_z)*.5) ;
}

void Hexahedron::sampleSurface(size_t num_points)
{
	Point point000_(*boundingPoints[0]) ;
	
	Point point111_(*boundingPoints[7]) ;
	
	size_t factor = 4 ;

	sampleBoundingSurface(num_points*factor) ;
	
	int numPointsPerDirection = std::max(static_cast<size_t>(ceil(pow((double)num_points*factor, 1./3.))),(size_t)2) -2;
	int numPointsPerDirectionOnBoundary = numPointsPerDirection+2 ;
	std::vector<Point> points ;
	
	while(numPointsPerDirection > 0)
	{
		Point point000(point000_ - center) ;
		point000 *= (double)(numPointsPerDirection-1)/(double)(numPointsPerDirectionOnBoundary-1) ;
		point000 += center ;
		
		Point point111(point111_ - center) ;
		point111 *= (double)(numPointsPerDirection-1)/(double)(numPointsPerDirectionOnBoundary-1) ;
		point111 += center ;
		
		std::vector<double> ds ;
		ds.push_back(0) ;
		for(int i = 1 ; i < numPointsPerDirection-1 ; i++)
		{
			double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;
			
			ds.push_back(v) ;
		}
		ds.push_back(1) ;

		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point000.x, 
				                       point111.y*ds[i]+point000.y*(1.-ds[i]), 
				                       point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [000] [010] [100]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z)) ;
			}
		}
		
	//face [000] [001] [100]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{	
				points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [001] [011] [101]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z)) ;
			}
		}
		
	//face [100] [110] [101]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [010] [011] [110]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
		numPointsPerDirection -= 2 ;
		
	}
	std::sort(points.begin(), points.end()) ;
	auto e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	for (size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;
	
	inPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		Point randomize(.2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX,
		                .2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX,
		                .2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX
		               ) ;
		inPoints[i] = new Point(points[i]+randomize) ;
	}
	
}

bool Hexahedron::in(const Point & v) const
{
	return v.x >= (center.x - size_x*.5 - POINT_TOLERANCE_3D) &&
		   v.x <= (center.x + size_x*.5 + POINT_TOLERANCE_3D) &&
		   v.y >= (center.y - size_y*.5 - POINT_TOLERANCE_3D) &&
		   v.y <= (center.y + size_y*.5 + POINT_TOLERANCE_3D) &&
		   v.z >= (center.z - size_z*.5 - POINT_TOLERANCE_3D) &&
		   v.z <= (center.z + size_z*.5 + POINT_TOLERANCE_3D) ;
}

double Hexahedron::area() const
{
	return 2.*size_x*size_y + 2.*size_x*size_z + 2.*size_y*size_z ;
}

double Hexahedron::volume() const
{
	return size_x*size_y*size_z ;
}

std::vector<Point> Hexahedron::getBoundingBox() const
{
	std::vector<Point> ret ;
	ret.push_back(Point(center.x+0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z)) ;
	return ret ;
}

void Hexahedron::project(Point * p) const
{
	std::vector<Point> bbox = getBoundingBox() ;
	
	Plane p0(bbox[0], bbox[1], bbox[2]) ;
	Plane p1(bbox[4], bbox[5], bbox[6]) ;
	Plane p2(bbox[0], bbox[1], bbox[4]) ;
	Plane p3(bbox[2], bbox[3], bbox[6]) ;
	Plane p4(bbox[0], bbox[2], bbox[4]) ;
	Plane p5(bbox[1], bbox[3], bbox[5]) ;
	Line l0(bbox[0], bbox[1]-bbox[0]) ;
	Line l1(bbox[2], bbox[2]-bbox[3]) ;
	Line l2(bbox[4], bbox[4]-bbox[5]) ;
	Line l3(bbox[6], bbox[6]-bbox[7]) ;
	Line l4(bbox[0], bbox[2]-bbox[0]) ;
	Line l5(bbox[1], bbox[1]-bbox[3]) ;
	Line l6(bbox[4], bbox[4]-bbox[6]) ;
	Line l7(bbox[5], bbox[5]-bbox[7]) ;
	Line l8(bbox[0], bbox[0]-bbox[4]) ;
	Line l9(bbox[1], bbox[1]-bbox[5]) ;
	Line l10(bbox[2], bbox[2]-bbox[6]) ;
	Line l11(bbox[3], bbox[3]-bbox[7]) ;
	
	std::map<double, Point> targets ;
	targets[dist(p0.projection(*p), *p)] = p0.projection(*p) ; 
	targets[dist(p1.projection(*p), *p)] = p1.projection(*p) ;
	targets[dist(p2.projection(*p), *p)] = p2.projection(*p) ; 
	targets[dist(p3.projection(*p), *p)] = p3.projection(*p) ; 
	targets[dist(p4.projection(*p), *p)] = p4.projection(*p) ; 
	targets[dist(p5.projection(*p), *p)] = p5.projection(*p) ; 
	targets[dist(l0.projection(*p), *p)] = l0.projection(*p) ; 
	targets[dist(l1.projection(*p), *p)] = l1.projection(*p) ; 
	targets[dist(l2.projection(*p), *p)] = l2.projection(*p) ; 
	targets[dist(l3.projection(*p), *p)] = l3.projection(*p) ; 
	targets[dist(l4.projection(*p), *p)] = l4.projection(*p) ; 
	targets[dist(l5.projection(*p), *p)] = l5.projection(*p) ; 
	targets[dist(l6.projection(*p), *p)] = l6.projection(*p) ; 
	targets[dist(l7.projection(*p), *p)] = l7.projection(*p) ; 
	targets[dist(l8.projection(*p), *p)] = l8.projection(*p) ; 
	targets[dist(l9.projection(*p), *p)] = l9.projection(*p) ; 
	targets[dist(l10.projection(*p), *p)] = l10.projection(*p) ; 
	targets[dist(l11.projection(*p), *p)] = l11.projection(*p) ; 
	for(size_t i = 0 ; i < 8 ; i++)
		targets[dist(*p, bbox[i])] = bbox[i] ;
	
	for(auto i = targets.begin() ; i!=targets.end() ; ++i)
	{
		if(in(i->second))
		{
			*p = i->second ;
			return ;
		}
	}
	*p = targets.begin()->second ;
}

Sphere::Sphere(double r,double x, double y, double z ) : radius(r), sqradius(r*r)
{
	center = Point(x,y,z) ;
	gType = SPHERE ;
}
	
Sphere::Sphere(double r,const Point * p0 ) : radius(r) , sqradius(r*r)
{ 
	center = Point(*p0) ;
	gType = SPHERE ;
}

Sphere::Sphere(double r,const Point p0 ) :radius(r) , sqradius(r*r)
{ 
	center = Point(p0) ; 
	gType = SPHERE ;
}

Sphere::Sphere(XMLTree * xml)
{
	gType = SPHERE ;
	this->center = Point(0,0,0) ;
	this->radius = 1 ; 

	if(xml->match("sphere"))
	{
		this->center = Point(xml->getChild(0)->getChild(0)) ;
		this->radius = xml->getChild(1)->buildDouble().second ;
	}

	this->sqradius = this->radius * this->radius ;

}


XMLTree * Sphere::toXML()
{
	XMLTree * sph = new XMLTree("sphere") ;
	XMLTree * c = new XMLTree("center") ;
	c->addChild(this->getCenter().toXML()) ;
	sph->addChild(c) ;
	sph->addChild(new XMLTree("radius",radius)) ;
	return sph ;	
}

void Sphere::dumpSampleBoundingPoints(size_t n, size_t iter)
{
	Sphere sph(1.,0.,0.,0.) ;
	std::vector<Point> p = sph.getSamplingPointsOnSphere(n,1., iter) ;
	std::valarray<double> d(p.size()*3) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		d[3*i+0] = p[i].x ;
		d[3*i+1] = p[i].y ;
		d[3*i+2] = p[i].z ;
	}
	
	std::string file = "../geometry/sphere_mesh_" ;
	file += itoa(n,10) ;
	std::fstream out ;
	out.open(file.c_str(), std::ios::out|std::ios::binary) ;
//	out << p.size() << std::endl ;
	for(size_t i = 0 ; i < d.size() ; i++)
	{
//		unsigned char * buf = (unsigned char *) & d[i] ;
//		for(int j = 0 ; j < sizeof(double) ; j++)
//			out << buf[j] ;
		out << std::setprecision(16) << d[i] ;
	}
	p[0].print() ;
	out.close() ;
}

std::vector<Point> Sphere::importStandardBoundingPoints(size_t n)
{
	std::vector<Point> ret ;
	std::string file = "../geometry/sphere_mesh_" ;
	file += itoa(n,10) ;
	std::fstream in ;
	in.open(file.c_str(), std::ios::in|std::ios::binary) ;
	double dx, dy, dz ;
	for(size_t i = 0 ; i < n ; i++)
	{
		in >> dx ;
		in >> dy ;
		in >> dz ;
		ret.push_back(Point(dx,dy,dz)) ;
	}
	return ret ;
}


std::vector<Point> Sphere::getBoundingBox() const
{
	double r = getRadius() ;
	std::vector<Point> ret ;
	ret.push_back(Point(center.x+0.5*r, center.y+0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y+0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y-0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y-0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y+0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y+0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y-0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y-0.5*r, center.z-0.5*r)) ;
	return ret ;
}

std::vector<Point> Sphere::getSamplingPointsOnSphere(size_t num_points, double r, size_t iter, size_t threshold) const
{
	
	size_t import = 128 ;
	while(import < num_points)
		import *= 2 ;


	double ns = ((double) import-(double)num_points)/(double) (import/2) ;
	ns *= import ;
//	ns = std::max(ns,(double)import/4) ;

//	std::cerr << iter << "-" << ns << "-" << import << "-" << num_points << std::endl ;
	
	if(iter > ns)
	{
		std::vector<Point> standard = getStandardSamplingBoundingPointsOnSphere(num_points) ;
		for(size_t i = 0 ; i < standard.size() ; i++)
		{
			standard[i] -= getCenter() ;
			standard[i] *= r ;
			standard[i] += getCenter() ;
		}
		
		return standard ;
	}
	
	//equation of a sphere = 
	//x = r sin(theta) cos(phi)
	//y = r sin(theta) sin(phi)
	//z = r cos(theta)
	
	std::vector<Point> points ;
	if(r  < POINT_TOLERANCE_3D)
		return points ;
	
// 		first we sample a cube.

	Point point000(-r + center.x, 
	               -r + center.y, 
	               -r + center.z) ;
	Point point111(r + center.x, 
	               r + center.y, 
	               r + center.z) ;
	
	
	size_t numPointsPerDirection = static_cast<size_t>(pow((double)0.25*num_points, 1./3.)) ;
	//face [000] [010] [001]
	
	std::vector<double> ds ;
	
	if(numPointsPerDirection > 1)
	{
		for(size_t i = 0 ; i < numPointsPerDirection ; i++)
		{
			double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;
			ds.push_back(v) ;
		}
	}
	else
	{
		ds.push_back(1) ;
	}
	
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point000.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [000] [010] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z)) ;
		}
	}
	
	//face [000] [001] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{	
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [001] [011] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z)) ;
		}
	}
	
	//face [100] [110] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [010] [011] [110]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	

	
	for(size_t i = points.size() ; i < num_points ; i++)
	{
		double rx = ((double)rand()/(RAND_MAX))*2.-1. ;
		double ry = ((double)rand()/(RAND_MAX))*2.-1. ;
		double rz = ((double)rand()/(RAND_MAX))*2.-1. ;
		points.push_back(Point(rx, ry, rz)+center) ;
		project(&points[i],r) ;
	}
	
	bool haveDuplicates = true ;
	while(haveDuplicates)
	{
		haveDuplicates = false ;
		for(size_t i  = 0 ; i < points.size() ; i++)
		{
			for(size_t j  = i+1 ; j < points.size() ; j++)
			{
				if(squareDist3D(points[i], points[j])< 128*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					haveDuplicates = true ;
					points.erase(points.begin()+j) ;
					break ;
				}
			}
			
			if(haveDuplicates)
				break ;
		}
	}
	
	smooth(points, r, iter*10) ;
	
	return points ;
}

void Sphere::smooth(std::vector<Point> & points,double r, size_t iter) const
{
	if(points.empty())
		return ;
	std::valarray<Point> speeds(points.size()) ;
//	std::cout << r << std::endl ;
	Point vec ;
	double error = 2. ;
	double last_error = 1. ;
	int count = 0 ;
	double derr = 1. ;
	for(size_t i = 0 ; /*(i < iter) &&*/ std::abs(error-last_error)/last_error > POINT_TOLERANCE_3D*POINT_TOLERANCE_3D*points.size()*points.size() && (count == 0); i++)
	{
		derr = std::abs(error-last_error) ;
		if(iter && i%iter == 0)
			std::cout << derr << std::endl ;
		last_error = error ;
		error = 0. ;
		for(size_t j = 0 ; j < points.size() ; j++)
		{
			for(size_t k = j+1 ; k < points.size() ; k++)
			{
				if(squareDist3D( points[j], points[k]) > 512.*POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
				{
					vec.set(points[j].x-points[k].x,points[j].y-points[k].y,points[j].z-points[k].z) ;
					double n = vec.sqNorm() ;
					error += n ;
					vec *= (sqradius/n) ;

					speeds[j] += vec ;
					speeds[k] -= vec ;
				}
				else
				{
					count++ ;
					Point vec(r,0,0) ;

					speeds[j] += vec ;
					speeds[k] -= vec ;
				}
			}
		}
		
		
		for(size_t j = 0 ; j < points.size() ; j++)
		{
			points[j] += speeds[j];
			project(&points[j],r) ;
		}
		if(last_error < 1e-12)
			break ;
		
		speeds = Point() ;
	}
/* 	std::cout << error << std::endl ;
				std::cout << std::abs(error-last_error)/last_error << std::endl ;*/
	
// 	for(size_t j = 0 ; j < points.size() ; j++)
// 	{
// 			project(&points[j],r) ;
// 	}
}

std::vector<Point> Sphere::getStandardSamplingBoundingPointsOnSphere(size_t n) const
{
	size_t import = 128 ;
	while(import < n)
		import *= 2 ;
		
	if(import > 16384)
		import = 16384 ;
	
//	std::cout << n << "-" << import << std::endl ;
	
	std::vector<Point> p = Sphere::importStandardBoundingPoints(import) ;
	std::random_shuffle(p.begin(), p.end());
	if(n <= p.size())
		p.erase(p.begin()+n, p.end()) ;
	
//	p[0].print() ;

	Point c = getCenter() ;
	for(size_t i = 0 ; i < p.size() ; i++)
		p[i] += c ;

	double ns = ((double) import-n)/(double) (import-import/2) ;
	ns *= import ;
//	ns = std::max(ns,(double)import/4) ;
//	std::cerr << ns << "here" << std::endl ;
//	std::cout << p.size() << "-" << n << std::endl ; 

	smooth(p,1., int(ns*2)) ;
//	p[0].print() ;
//	std::cerr << "smoothed" << std::endl ;

	return p ;
}

std::vector<Point> Sphere::getSamplingBoundingPoints(size_t num_points) const
{
	return getSamplingPointsOnSphere(num_points, radius) ;
}

void Sphere::sampleBoundingSurface(size_t num_points)
{	
	std::vector<Point> points = getSamplingPointsOnSphere(num_points, radius) ;
	for(size_t i = 0 ; i < boundingPoints.size() ; i++)
		delete boundingPoints[i] ;
	
	boundingPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
// 		std::cout << points[i].x << ";" << points[i].y << ";" << points[i].z << std::endl ;
		boundingPoints[i] = new Point(points[i]) ;
	}
	
	
}

void Sphere::sampleSurface(size_t num_points) 
{
	if(num_points < 4)
		return ;

	sampleBoundingSurface(num_points) ;
	
	std::vector<Point> points ;

	size_t numPointsOnSurface = boundingPoints.size() ;
	size_t numberOfRadii = static_cast<size_t>(round(sqrt(numPointsOnSurface/6))) ;
	
	for(size_t i = 0 ; i < numberOfRadii ; i++)
	{
// 		if(i == 0)
// 			numPointsOnSurface = 7*numPointsOnSurface/8 ;
		if(i == 0)
			numPointsOnSurface = 6*numPointsOnSurface/7 ;
		if(i == 1)
			numPointsOnSurface = 5*numPointsOnSurface/6 ;
		if(i == 2)
			numPointsOnSurface = 4*numPointsOnSurface/5 ;
		if(i == 3)
			numPointsOnSurface = 3*numPointsOnSurface/4 ;
		if(i == 4)
			numPointsOnSurface = 2*numPointsOnSurface/3 ;
		if(i == 5)
			numPointsOnSurface = 1*numPointsOnSurface/2 ;
		
		double r = radius*((double)(numberOfRadii-i-1)/(double)numberOfRadii) ;
		num_points = (size_t)((double)numPointsOnSurface*((r*r)/(radius*radius)));
		if(num_points < 8 )
			break ;
		
		std::vector<Point> newPoints = getSamplingPointsOnSphere(num_points, r) ;
		
		for(size_t j = 0 ; j < newPoints.size() ; j++)
		{
// 			double r = ((double)rand()/(double)(RAND_MAX+1)*2. -1.)*.0005*(sqrt(radius)/(double)(numberOfRadii+1)) ;
// 			Point dr = (newPoints[j] - center)*r ;
			points.push_back(newPoints[j]/*+ dr*/) ;
		}
	}
	points.push_back(this->center) ;
	for(size_t i = 0 ; i < inPoints.size() ; i++)
		delete inPoints[i] ;
	
	inPoints.resize(points.size()) ;
	
	for(size_t i = 0 ; i < points.size() ; i++)
	{
// 		std::cout << points[i].x << "; " << points[i].y << ";" << points[i].z << std::endl ;
		inPoints[i] = new Point(points[i]) ;
	}
}

bool Sphere::in(const Point & v) const 
{ 
	if(v.x < center.x-1.0001*getRadius())
		return false ;
	if(v.x > center.x+1.0001*getRadius())
		return false ;
	if(v.y < center.y-1.0001*getRadius())
		return false ;
	if(v.y > center.y+1.0001*getRadius())
		return false ;
	if(v.z < center.z-1.0001*getRadius())
		return false ;
	if(v.z > center.z+1.0001*getRadius())
		return false ;
	return squareDist3D(v, getCenter()) < getRadius()*getRadius() + 2.*getRadius()*POINT_TOLERANCE_3D + POINT_TOLERANCE_3D*POINT_TOLERANCE_3D;
}

double Sphere::area() const
{
	return 4.*M_PI*getRadius()*getRadius() ;
}

double Sphere::volume() const
{
	return 1.333333333333333333*M_PI*sqradius*radius ;
}

void Sphere::project(Point * p) const
{
	project(p, getRadius()) ;
}

void Sphere::project(Point * p, double r) const
{
	//x = r sin(theta) cos(phi)
	//y = r sin(theta) sin(phi)
	//z = r cos(theta)
	int id = p->id ;
	if(squareDist3D(*p, getCenter() ) < POINT_TOLERANCE_3D*POINT_TOLERANCE_3D)
	{
		p->x =+ r ;
		return ;
	}
	
	Point p_prime = *p-getCenter() ;
	double n = p_prime.norm() ;
	if(std::abs(n - r) < POINT_TOLERANCE_3D)
		return ;
	p_prime *= r/n ;
	*p = getCenter()+p_prime ;
	p->id = id ;
	return ;

}

void Sphere::computeCenter() { } ;

double Sphere::getRadius() const 
{
	return this->radius ;
}

void Sphere::setRadius(double newr)
{

	double ratio = newr/radius ;
	
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
	this->radius = newr ;
}

void TriangulatedSurface::addPoint(Point * p)
{
	const Point * nearest = boundary[0] ;
	for(size_t i = 0 ; i < boundary.size() ; i++)
	{
		
	}
}

TriangulatedSurface::TriangulatedSurface(const Point * p0, const Point * p1, const Point * p2)
{
	mesh.push_back(TriPoint(p0, p1, p2)) ;
	boundary.push_back(p0) ;
	boundary.push_back(p1) ;
	boundary.push_back(p2) ;
	
	center = (*p0 + *p1 + *p2)/3 ;
}

TriangulatedSurface::TriangulatedSurface()
{
	mesh.push_back(TriPoint(new Point(0,0,1), new Point(1,0,0), new Point(0,0,0))) ;
	boundary.push_back(mesh[0].point[0]) ;
	boundary.push_back(mesh[0].point[1]) ;
	boundary.push_back(mesh[0].point[2]) ;
	center = Point(1./3., 0, 1./3.) ;
}

double TriangulatedSurface::area() const
{
	double ret = 0 ;
	for(size_t i = 0 ; i < mesh.size() ;i++)
	{
		ret += mesh[i].area() ;
	}
	
	return ret ;
}

void TriangulatedSurface::project(Point * p) const
{
	
}

void TriangulatedSurface::computeCenter()
{
	
}

double TriangulatedSurface::getRadius() const {return 0 ;}

std::vector<Point> TriangulatedSurface::getBoundingBox() const {return std::vector<Point>(0) ;}



	
