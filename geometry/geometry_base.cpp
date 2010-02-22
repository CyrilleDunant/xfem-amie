// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009 (added: ellipses)
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_base.h"
#include <limits>
#include <iomanip>

#include "../mesher/delaunay.h"
#include "../polynomial/vm_function_base.h"

using namespace Mu ;

Point::Point() : id(-1)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_setzero_pd() ;
	veczt = _mm_setzero_pd() ;
	#else
	x = 0 ; y = 0 ; z = 0 ; t = 0 ;
	#endif
}

#ifdef HAVE_SSE3
Point::Point(const Point & p) : vecxy(p.vecxy), veczt(p.veczt), id(p.id)
{
}
#else
Point::Point(const Point & p) 
{
	x = p.x ; y = p.y ; z = p.z ; t = p.t ; id = p.id ;
}
#endif

double Point::angle() const
{
	return atan2(y,x) ;
}

Point::Point(double x_, double y_) : id(-1)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_setr_pd(x_, y_) ;
	veczt =_mm_setzero_pd() ;
	#else
	x= x_ ; y = y_ ; z = 0 ; t = 0 ;
	#endif
}

Point::Point(double x_, double y_, double z_): id(-1)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_setr_pd(x_, y_) ;
	veczt = _mm_set_sd(z_) ;
	#else
		x= x_ ; y = y_ ; z = z_ ; t = 0 ;
	#endif
}

Point::Point(double x_, double y_, double z_, double t_): id(-1)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_setr_pd(x_, y_) ;
	veczt = _mm_setr_pd(z_, t_) ;
	#else
	x= x_ ; y = y_ ; z = z_ ; t = t_ ;
	#endif
}

void Point::print() const
{
	std::cout << " ( id = " << id << std::flush ;
	std::cout << " ; "<< std::setprecision(14) << x << std::flush ;
	std::cout << "; " << std::setprecision(14)<< y << std::flush ;
	std::cout << "; " << std::setprecision(14)<< z << std::flush ;
	std::cout << "; " << std::setprecision(14)<< t << ") " << std::endl;
}

double Point::norm() const
{
#ifdef HAVE_SSE4
	vecdouble r0 ;
	r0.vec = _mm_dp_pd(vecxy, vecxy, 61) ;
	r0.vec += _mm_dp_pd(veczt, veczt, 62) ;
	return sqrt(r0.val[0]+ r0.val[1]);
#elif defined HAVE_SSE3
	vecdouble rzt ;
	rzt.vec = _mm_add_pd(_mm_mul_pd(veczt, veczt), _mm_mul_pd(vecxy, vecxy)) ;
	return sqrt(rzt.val[0]+ rzt.val[1]);
#else 
	return sqrt(x*x+y*y+z*z+t*t) ;
#endif
}

double Point::sqNorm() const
{
#ifdef HAVE_SSE4
	vecdouble r0 ;
	r0.vec = _mm_dp_pd(vecxy, vecxy, 61) ;
	r0.vec += _mm_dp_pd(veczt, veczt, 62) ;	
	return r0.val[0]+ r0.val[1];
#elif HAVE_SSE3
	vecdouble rzt ;
	rzt.vec = _mm_add_pd(_mm_mul_pd(veczt, veczt), _mm_mul_pd(vecxy, vecxy)) ;
	return rzt.val[0]+ rzt.val[1];
#else 
	return x*x+y*y+z*z+t*t ;
#endif
}

void Point::setX(double v) 
{ 
	x = v ;
}

void Point::setY(double v) 
{ 
	y = v ;
}

void Point::setZ(double v) 
{ 
	z = v ;
}

void Point::setT(double v) 
{ 
	t = v ;
}

void Point::set(const Point & p)
{
	#ifdef HAVE_SSE3
	vecxy = p.vecxy ;
	veczt = p.veczt ;
	#else
	x = p.x ; y = p.y ; z = p.z ; t = p.t ;
	#endif
}

void Point::set(const Point * p)
{
	#ifdef HAVE_SSE3
	vecxy = p->vecxy ;
	veczt = p->veczt ;
	#else
	x = p->x ; y = p->y ; z = p->z ; t = p->t ;
	#endif
}


void Point::set(double v, double vv)
{
	x = v ; 
	y = vv ;
}

void Point::set(double v, double vv, double vvv)
{
	x = v ; 
	y = vv ;
	z = vvv ;
}

void Point::set(double v, double vv, double vvv, double vvvv)
{
	x = v ; 
	y = vv ;
	z = vvv ;
	t = vvvv ;
}

bool Point::operator==(const Point &p) const
{
	
	if(std::abs(p.x-x) > 2.*POINT_TOLERANCE)
		return false ;
	if(std::abs(p.y-y) > 2.*POINT_TOLERANCE)
		return false ;
	if(std::abs(p.z-z) > 2.*POINT_TOLERANCE)
		return false ;
	if(std::abs(p.t-t) > 2.*POINT_TOLERANCE)
		return false ;
	
	double delta = POINT_TOLERANCE ;
	Point a(p) ; a.x += delta ; a.y += delta ; a.z += delta ;
	Point b(p) ; b.x += delta ; b.y += delta; b.z -= delta ;
	Point c(p) ; c.x += delta ; c.y -= delta; c.z += delta ;
	Point d(p) ; d.x += delta ; d.y -= delta; d.z -= delta ;
	Point e(p) ; e.x -= delta ; e.y += delta; e.z += delta ;
	Point f(p) ; f.x -= delta ; f.y += delta; f.z -= delta ;
	Point g(p) ; g.x -= delta ; g.y -= delta; g.z += delta ;
	Point h(p) ; h.x -= delta ; h.y -= delta; h.z -= delta ;

	return squareDist( &p, this) < POINT_TOLERANCE*POINT_TOLERANCE 
		|| squareDist( &a, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &b, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &c, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &d, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &e, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &f, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &g, this) < POINT_TOLERANCE*POINT_TOLERANCE
		|| squareDist( &h, this) < POINT_TOLERANCE*POINT_TOLERANCE ;
}

bool Point::operator!=(const Point & p) const
{
	return !(*this == p) ;
}

Point Point::operator-(const Point &p) const
{
	Point ret((*this)) ;
	#ifdef HAVE_SSE3
	ret.vecxy = _mm_sub_pd(ret.vecxy, p.vecxy) ;
	ret.veczt = _mm_sub_pd(ret.veczt, p.veczt) ;
	ret.id = std::max(id, p.id) ;
	#else 
	ret.x -= p.x ;
	ret.y -= p.y ;
	ret.z -= p.z ;
	ret.t -= p.t ;
	#endif
	return ret ;
}

Point Point::operator-(const Vector &p) const
{
	Point ret((*this)) ;
	ret.x -= p[0] ; 
	ret.y -= p[1] ; 
	if(p.size() > 2)
		ret.z -= p[2] ; 
	if(p.size() > 3)
		ret.t -= p[3] ; 
	ret.id = id ;
	return ret ; 
}

Point Point::operator+(const Point &p) const
{
	Point ret((*this)) ;
	#ifdef HAVE_SSE3
	ret.vecxy = _mm_add_pd(ret.vecxy, p.vecxy) ;
	ret.veczt = _mm_add_pd(ret.veczt, p.veczt) ;
	#else 
	ret.x += p.x ;
	ret.y += p.y ;
	ret.z += p.z ;
	ret.t += p.t ;
	#endif
	ret.id = std::max(id, p.id) ;
	return ret ;
}

void Point::operator+=(const Point &p)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_add_pd(vecxy, p.vecxy) ;
	veczt = _mm_add_pd(veczt, p.veczt) ;
	#else 
	x += p.x ;
	y += p.y ;
	z += p.z ;
	t += p.t ;
	#endif
}

void Point::operator-=(const Point &p)
{
	#ifdef HAVE_SSE3
	vecxy = _mm_sub_pd(vecxy, p.vecxy) ;
	veczt = _mm_sub_pd(veczt, p.veczt) ;
	#else 
	x -= p.x ;
	y -= p.y ;
	z -= p.z ;
	t -= p.t ;
	#endif
}

Point Point::operator+(const Vector &p) const
{
	Point ret((*this)) ;
	ret.x += p[0] ; 
	ret.y += p[1] ; 
	if(p.size() > 2)
		ret.z += p[2] ; 
	if(p.size() > 3)
		ret.t += p[3] ; 
	ret.id = id ;
	return ret ; 
}

Point Point::operator/(const double p) const 
{
	Point ret((*this)) ;
	double inv = 1./p ;
	
	#ifdef HAVE_SSE3
	__m128d temp = _mm_load1_pd(&inv) ;
	ret.vecxy = _mm_mul_pd(ret.vecxy, temp) ;
	ret.veczt = _mm_mul_pd(ret.veczt, temp) ;
	#else 
	ret.x *= inv ;
	ret.y *= inv ;
	ret.z *= inv ;
	ret.t *= inv ;
	#endif
	ret.id = id ;
	return ret ; 
}

Geometry::~Geometry() 
{

	for(size_t i = 0 ; i < this->inPoints.size() ; i++)
	{
		delete inPoints[i] ;
		inPoints[i] = NULL ;
	}
}

bool Point::operator <(const Point &p) const 
{
// 	if(p == *this)
// 		return false ;
// 	
	if(x < p.x)
		return true ;
	else if(x > p.x)
		return false ;
	
	if(y < p.y)
		return true ;
	else if(y > p.y)
		return false ;
	
	if(z < p.z)
		return true ;
	else if (z > p.z)
		return false ;
	
	if(t < p.t)
		return true ;
	
	return false ;
}

bool Point::operator >(const Point &p) const 
{
	if(p == *this)
		return false ;
	
	double tol = POINT_TOLERANCE ;
	return (y > p.y ) 
		|| (( std::abs(y - p.y) < tol) 
		    && (x > p.x)) 
		|| (( std::abs(y - p.y) < tol) 
		    && ( std::abs(x - p.x) < tol) 
		    && (z> p.z)) 
		||(( std::abs(y - p.y) < tol) 
		   && ( std::abs(x - p.x) < tol) 
		   && ( std::abs(z - p.z) < tol) 
		   && (z> p.z));
}

Point Point::operator*(const double p)  const 
{
	Point ret((*this)) ;
	#ifdef HAVE_SSE3
	__m128d temp = _mm_load1_pd(&p) ;

	ret.vecxy = _mm_mul_pd(ret.vecxy, temp) ;
	ret.veczt = _mm_mul_pd(ret.veczt, temp) ;
	#else
	ret.x *= p ;
	ret.y *= p ;
	ret.z *= p ;
	ret.t *= p ;
	#endif
	ret.id = id;
	return ret ; 
}

double Point::operator*(const Point &p) const
{
#ifdef HAVE_SSE4
	vecdouble r ;
	r.vec = _mm_dp_pd(p.vecxy, vecxy, 61) ;
	r.vec += _mm_dp_pd(p.veczt, veczt, 62) ;
	return r.val[0] + r.val[1];
#elif defined HAVE_SSE3
	vecdouble r ;
	r.vec = _mm_add_pd(_mm_mul_pd(p.vecxy, vecxy), _mm_mul_pd(p.veczt, veczt)) ;
	return r.val[0] + r.val[1];
#endif
	return p.x*x+p.y*y+p.z*z+p.t*t ;
	
}



double Point::operator*(const Vector &p) const
{
	double ret = x*p[0] + y*p[1] ;
	if(p.size() > 2)
		ret+=z*p[2] ;
	if(p.size() > 3)
		ret+=t*p[3] ;
	return ret ; 
}


Point Point::operator^(const Point &p) const
{
	Point ret ;

	ret.x = y*p.z - z*p.y ; //fma(y,p.z,  -z*p.y) ;
	ret.y = z*p.x - x*p.z ;//fma(z,p.x , -x*p.z) ;
	ret.z = x*p.y - y*p.x ; //fma(x,p.y , -y*p.x) ;
	
	ret.id = std::max(id, p.id) ;
	return ret ;
}

Point Point::operator^(const Vector &p) const
{
	Point ret ;
	ret.x = y*p[2] - z*p[1] ;
	ret.y = z*p[0] - x*p[2] ;
	ret.z = x*p[1] - y*p[0] ;
	ret.id = id;
	return ret ;
}


PointSet::PointSet() : boundingPoints(0)
{
	this->chullEndPos = 0;
}


std::vector<Point> Geometry::getBoundingBox() const
{
	return std::vector<Point>() ;
}

const std::valarray<Point *> & ConvexGeometry::getBoundingPoints() const
{ 
	return boundingPoints ; 
}

std::valarray<Point *> & ConvexGeometry::getBoundingPoints() 
{ 
	return boundingPoints ; 
}

const Point & ConvexGeometry::getBoundingPoint(size_t i) const
{ 
	return *boundingPoints[i] ; 
}

Point & ConvexGeometry::getBoundingPoint(size_t i) 
{ 
	return *boundingPoints[i] ; 
}

const Point & ConvexGeometry::getPoint(size_t i) const
{
	if (i < inPoints.size())
	{
		return *inPoints[i] ;
	}
	
	return *boundingPoints[i-inPoints.size()] ;
}

Point & ConvexGeometry::getPoint(size_t i)
{
	if (i < inPoints.size())
	{
		return *inPoints[i] ;
	}
	
	return *boundingPoints[i-inPoints.size()] ;
}

std::valarray<Point *> & NonConvexGeometry::getBoundingPoints() 
{ 
	return boundingPoints ; 
}

const std::valarray<Point *> & NonConvexGeometry::getBoundingPoints() const
{ 
	return boundingPoints ; 
}

const Point & NonConvexGeometry::getBoundingPoint(size_t i) const
{ 
	return *boundingPoints[i] ; 
}

Point & NonConvexGeometry::getBoundingPoint(size_t i) 
{ 
	return *boundingPoints[i] ; 
}

const Point & Geometry::getInPoint(size_t i) const
{ 
	return *inPoints[i] ; 
}


void Geometry::setCenter(const Point & newCenter)
{
	Point delta = newCenter-getCenter();
	getCenter() = newCenter ;
	
	for(size_t  i = 0 ; i < getInPoints().size() ;i++)
		getInPoint(i) += delta ;
	for(size_t  i = 0 ; i < getBoundingPoints().size() ;i++)
		getBoundingPoint(i) += delta ;
}

Point & Geometry::getInPoint(size_t i) 
{ 
	return *inPoints[i] ; 
}

void NonConvexGeometry::setBoundingPoint(size_t i, Point * p)
{
	boundingPoints[i] = p ;
}

void NonConvexGeometry::setBoundingPoints(const PointArray & nb)
{
	boundingPoints.resize(nb.size()) ;
	boundingPoints = nb ;
}

void ConvexGeometry::setBoundingPoint(size_t i, Point * p)
{
	boundingPoints[i] = p ;
}

void ConvexGeometry::setBoundingPoints(const PointArray & nb)
{
	boundingPoints.resize(nb.size()) ;
	boundingPoints = nb ;
}

double & Point::operator[](size_t i)
{
	return (*(&x+i)) ;
}
double Point::operator[](size_t i) const
{
	return (*(&x+i)) ;
}


Plane::Plane(const Point & origin, const Point & vector) : p(origin), v(vector)
{
}

bool Plane::intersects(const Line &l) const
{
	double d = p*v ;
	double numerator = d - l.origin()*v ;
	double denominator = l.vector()*v ;

	// line is in the plane
	if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
		return true ;
	
	// the intersection exists and is unique
	if(std::abs(denominator) > POINT_TOLERANCE)
		return true ;
	
	return false ;
}

bool Plane::intersects(const Segment &l) const
{
	double d = p*v ;
	double numerator = d - l.second()*v ;
	double denominator = l.vector()*v ;
	
	// segment is in the plane
	if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
		return true ;
	
	// the intersection exists and is unique
	if(std::abs(denominator) > POINT_TOLERANCE)
	{
		double t = numerator/denominator ;
		
		return t > 0 && t < 1 ;
	}
	
	return false ;
}

bool Plane::intersects(const Geometry *g) const
{
	switch (g->getGeometryType())
	{
	case SPHERE:
		{
			double n = v.norm() ;
			double d = v*p ;
			d /= n ;
			
			return d >= g->getRadius() ;
			
		}
	case HEXAHEDRON:
		{
			std::vector<Point> bbox = g->getBoundingBox() ;
			double maxx =  bbox[0].x ;
			double minx =  bbox[0].x ;
			double maxy =  bbox[0].y ;
			double miny =  bbox[0].y ;
			double maxz =  bbox[0].z ;
			double minz =  bbox[0].z ;
			
			for(size_t i = 1 ; i < bbox.size() ; i++)
			{
				if(bbox[i].x > maxx)
					maxx = bbox[i].x ;
				if(bbox[i].x < minx)
					minx = bbox[i].x ;
				if(bbox[i].y > maxy)
					maxy = bbox[i].y ;
				if(bbox[i].y < miny)
					miny = bbox[i].y ;
				if(bbox[i].z > maxz)
					maxz = bbox[i].z ;
				if(bbox[i].z < minz)
					minz = bbox[i].z ;
			}
			
			Point corner1 (minx, miny, minz) ;
			Point corner2 (minx, miny, maxz) ;
			Point corner3 (minx, maxy, minz) ;
			Point corner4 (minx, maxy, maxz) ;
			Point corner5 (maxx, miny, minz) ;
			Point corner6 (maxx, miny, maxz) ;
			Point corner7 (maxx, maxy, minz) ;
			Point corner8 (maxx, maxy, maxz) ;
			
			Segment side1(corner1, corner2) ;
			if(intersects(side1))
				return true ;
			Segment side2(corner1, corner3) ;
			if(intersects(side2))
				return true ;
			Segment side3(corner1, corner5) ;
			if(intersects(side3))
				return true ;
			Segment side4(corner2, corner4) ;
			if(intersects(side4))
				return true ;
			Segment side5(corner2, corner6) ;
			if(intersects(side5))
				return true ;
			Segment side6(corner3, corner4) ;
			if(intersects(side6))
				return true ;
			Segment side7(corner3, corner7) ;
			if(intersects(side7))
				return true ;
			Segment side8(corner4, corner8) ;
			if(intersects(side8))
				return true ;
			Segment side9(corner5, corner6) ;
			if(intersects(side9))
				return true ;
			Segment side10(corner5, corner7) ;
			if(intersects(side10))
				return true ;
			Segment side11(corner6, corner8) ;
			if(intersects(side11))
				return true ;
			Segment side12(corner7, corner8) ;
			if(intersects(side12))
				return true ;

			return false ;
		}
	case TETRAHEDRON:
		{
			Point corner1 (g->getBoundingPoint(0)) ;
			Point corner2 (g->getBoundingPoint(1)) ;
			Point corner3 (g->getBoundingPoint(2)) ;
			Point corner4 (g->getBoundingPoint(3)) ;
			
			Segment side1(corner1, corner2) ;
			if(intersects(side1))
				return true ;
			Segment side2(corner1, corner3) ;
			if(intersects(side2))
				return true ;
			Segment side3(corner2, corner4) ;
			if(intersects(side3))
				return true ;
			Segment side4(corner3, corner4) ;
			if(intersects(side4))
				return true ;
			
			return false ;
		}
	default:
		{
			std::cout << "intersection type not implemented" << std::endl ;
			return false ;
		}
	}
}

bool Plane::intersects(const Plane & plane) const
{
	return v/v.norm() != plane.vector()/plane.vector().norm() ;
}

bool Plane::on(const Point &point) const
{
	double d = v*p ;
	double dtest = point*v ;
	return std::abs(d-dtest) < POINT_TOLERANCE ;
}

std::vector<Point> Plane::intersection(const Geometry * g) const 
{
	switch (g->getGeometryType())
	{
	case SPHERE:
		{
			Point vec = projection(g->getCenter())-g->getCenter() ;
			double l = dist(g->getCenter(),projection(g->getCenter())) ;
			Point centerOfIntersection = g->getCenter() + vec/vec.norm()*l ;
			double radiusOfIntersection = sqrt(g->getRadius()*g->getRadius() - l*l) ;
			OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
			
			size_t num_points = std::max((int)round(7.*sqrt(g->getBoundingPoints().size())*radiusOfIntersection/g->getRadius()), 8 ) ;
			C.sampleSurface(num_points) ;
			
			return C.getSamplingBoundingPoints(num_points) ;
			
		}
	case HEXAHEDRON:
		{
			std::vector<Point> ret ;
			std::vector<Point> bbox = g->getBoundingBox() ;
			double maxx =  bbox[0].x ;
			double minx =  bbox[0].x ;
			double maxy =  bbox[0].y ;
			double miny =  bbox[0].y ;
			double maxz =  bbox[0].z ;
			double minz =  bbox[0].z ;
			
			for(size_t i = 1 ; i < bbox.size() ; i++)
			{
				if(bbox[i].x > maxx)
					maxx = bbox[i].x ;
				if(bbox[i].x < minx)
					minx = bbox[i].x ;
				if(bbox[i].y > maxy)
					maxy = bbox[i].y ;
				if(bbox[i].y < miny)
					miny = bbox[i].y ;
				if(bbox[i].z > maxz)
					maxz = bbox[i].z ;
				if(bbox[i].z < minz)
					minz = bbox[i].z ;
			}
			
			Point corner1 (minx, miny, minz) ;
			Point corner2 (minx, miny, maxz) ;
			Point corner3 (minx, maxy, minz) ;
			Point corner4 (minx, maxy, maxz) ;
			Point corner5 (maxx, miny, minz) ;
			Point corner6 (maxx, miny, maxz) ;
			Point corner7 (maxx, maxy, minz) ;
			Point corner8 (maxx, maxy, maxz) ;
			
			Segment side1(corner1, corner2) ;
			if(intersects(side1))
				ret.push_back(intersection(side1)) ;
			Segment side2(corner1, corner3) ;
			if(intersects(side2))
				ret.push_back(intersection(side2)) ;
			Segment side3(corner1, corner5) ;
			if(intersects(side3))
				ret.push_back(intersection(side3)) ;
			Segment side4(corner2, corner4) ;
			if(intersects(side4))
				ret.push_back(intersection(side4)) ;
			Segment side5(corner2, corner6) ;
			if(intersects(side5))
				ret.push_back(intersection(side5)) ;
			Segment side6(corner3, corner4) ;
			if(intersects(side6))
				ret.push_back(intersection(side6)) ;
			Segment side7(corner3, corner7) ;
			if(intersects(side7))
				ret.push_back(intersection(side7)) ;
			Segment side8(corner4, corner8) ;
			if(intersects(side8))
				ret.push_back(intersection(side8)) ;
			Segment side9(corner5, corner6) ;
			if(intersects(side9))
				ret.push_back(intersection(side9)) ;
			Segment side10(corner5, corner7) ;
			if(intersects(side10))
				ret.push_back(intersection(side10)) ;
			Segment side11(corner6, corner8) ;
			if(intersects(side11))
				ret.push_back(intersection(side11)) ;
			Segment side12(corner7, corner8) ;
			if(intersects(side12))
				ret.push_back(intersection(side12)) ;
			
			return ret ;
		}
	case TETRAHEDRON:
		{
			std::vector<Point> ret ;
			Point corner1 (g->getBoundingPoint(0)) ;
			Point corner2 (g->getBoundingPoint(1)) ;
			Point corner3 (g->getBoundingPoint(2)) ;
			Point corner4 (g->getBoundingPoint(3)) ;
			
			Segment side1(corner1, corner2) ;
			if(intersects(side1))
				ret.push_back(intersection(side1)) ;
			Segment side2(corner1, corner3) ;
			if(intersects(side2))
				ret.push_back(intersection(side2)) ;
			Segment side3(corner2, corner4) ;
			if(intersects(side3))
				ret.push_back(intersection(side3)) ;
			Segment side4(corner3, corner4) ;
			if(intersects(side4))
				ret.push_back(intersection(side4)) ;
			
			return ret ;
		}
	default:
		{
			std::cout << "intersection type not implemented" << std::endl ;
			return std::vector<Point>() ;
		}
	}
}

Point Plane::intersection(const Line &l) const
{
	double d = p*v ;
	double numerator = d - l.origin()*v ;
	double denominator = l.vector()*v ;
	
	// segment is in the plane
	if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
		return p ;
	
	//there is no intersection, but we need to return something
	return p ;
	
	// the intersection exists and is unique
	double t = numerator/denominator ;
	return l.origin()+l.vector()*t ;

}

Point Plane::intersection(const Segment &l) const
{
	double d = p*v ;
	double numerator = d - l.second()*v ;
	double denominator = l.vector()*v ;
	
	// segment is in the plane
	if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
		return l.midPoint() ;
	
	// the intersection exists and is unique
	if(std::abs(denominator) > POINT_TOLERANCE)
	{
		double t = numerator/denominator ;
		
		return l.second()+l.vector()*t ;
	}
	
	return p ;
}

Line Plane::intersection(const Plane & plane) const 
{
	Point vec = v^plane.vector() ;
	double d1 = p*v ;
	double d2 = plane.vector()*plane.origin() ;
	if(!((v.x != 0  && plane.vector().y != 0) || (v.y != 0  && plane.vector().x != 0)))
	{
		Matrix m(2, 2) ;
		m[0][0] = v.x ; m[0][1] = v.y ;
		m[1][0] = plane.vector().x ; m[1][1] = plane.vector().y ;
		
		invert2x2Matrix(m) ;
		
		Vector d (2) ; d[0] = d1 ; d[1] = d2 ;
		
		Vector xy = m*d ;
		
		return Line(Point(xy[0], xy[1], 0), vec) ;
	}
	else if(!((v.x != 0  && plane.vector().z != 0) || (v.z != 0  && plane.vector().x != 0)))
	{
		Matrix m(2, 2) ;
		m[0][0] = v.x ; m[0][1] = v.z ;
		m[1][0] = plane.vector().x ; m[1][1] = plane.vector().z ;
		
		invert2x2Matrix(m) ;
		
		Vector d (2) ; d[0] = d1 ; d[1] = d2 ;
		
		Vector xz = m*d ;
		
		return Line(Point(xz[0], 0, xz[1]), vec) ;
	}
	else if(!((v.y != 0  && plane.vector().z != 0) || (v.z != 0  && plane.vector().y != 0)))
	{
		Matrix m(2, 2) ;
		m[0][0] = v.y ; m[0][1] = v.z ;
		m[1][0] = plane.vector().x ; m[1][1] = plane.vector().z ;
		
		invert2x2Matrix(m) ;
		
		Vector d (2) ; d[0] = d1 ; d[1] = d2 ;
		
		Vector yz = m*d ;
		
		return Line(Point(0, yz[0], yz[1]), vec) ;
	}
	
	return Line(Point(0, 0, 0), vec) ;
}

const Point & Plane::vector() const 
{
	return v ;
}

const Point & Plane::origin() const 
{
	return p ;
}

Point Plane::projection(const Point &p ) const 
{
	return intersection(Line(p, v)) ;
}



const std::valarray<Point* > & Geometry::getInPoints() const
{
	return inPoints ;
}

std::valarray<Point* > & Geometry::getInPoints() 
{
	return inPoints ;
}

bool Geometry::intersects(const Geometry *g) const
{
	switch(getGeometryType())
	{
	case TRIANGLE :
		{
			int numpoints = getBoundingPoints().size() ;
			Segment s0(getBoundingPoint(0), getBoundingPoint(numpoints/3)) ;
			Segment s1(getBoundingPoint(numpoints/3), getBoundingPoint(2*numpoints/3)) ;
			Segment s2(getBoundingPoint(0),getBoundingPoint(2*numpoints/3)) ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				int gnumpoints = g->getBoundingPoints().size() ;
				Segment t0(g->getBoundingPoint(0), g->getBoundingPoint(gnumpoints/3)) ;
				Segment t1(g->getBoundingPoint(gnumpoints/3), g->getBoundingPoint(2*gnumpoints/3)) ;
				Segment t2(g->getBoundingPoint(0), g->getBoundingPoint(2*gnumpoints/3)) ;
				
				return s0.intersects(t0) || s0.intersects(t1) || s0.intersects(t2)
					|| s1.intersects(t0) || s1.intersects(t1) || s1.intersects(t2) 
					|| s2.intersects(t0) || s2.intersects(t1) || s2.intersects(t2) ;
				
			}
			
			return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g) ;
		}
	case RECTANGLE:
		{
			std::vector<Point> box = this->getBoundingBox() ;
			Segment s0(box[0], box[1]) ;
			Segment s1(box[1], box[2]) ;
			Segment s2(box[2], box[3]) ;
			Segment s3(box[3], box[0]) ;

			return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g)  || s3.intersects(g);
		}
	case SEGMENTED_LINE:
		{
			std::vector<Segment> segs ;
			for(size_t i = 0 ; i < this->getBoundingPoints().size()-1 ; i++)
			{
				segs.push_back(Segment(this->getBoundingPoint(i), this->getBoundingPoint(i+1))) ;
			}
			bool intersects = false ;
			for(size_t i = 0 ; i < segs.size() ; i++)
			{
				intersects = intersects || segs[i].intersects(g) ;
			}
			return intersects ;
		}
	case CIRCLE:
		{
			
			if(g->getGeometryType() == CIRCLE)
			{
				double birad = getRadius()+g->getRadius() ;
				if(getCenter().x + birad < g->getCenter().x)
					return false ;
				if(getCenter().x - birad > g->getCenter().x)
					return false ;
				if(getCenter().y + birad < g->getCenter().y)
					return false ;
				if(getCenter().y - birad > g->getCenter().y)
					return false ;

				return squareDist2D( getCenter(), g->getCenter()) < birad*birad ;
			}
			
			std::vector<Segment> segs ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
			}
			if(g->getGeometryType() == RECTANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
			}
			if(g->getGeometryType() == SEGMENTED_LINE)
			{
				for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
				{
					segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
				}
			}
			
			bool intersects = false ;
			for(size_t i = 0 ; i < segs.size() ; i++)
			{
				intersects = intersects || segs[i].intersects(this) ;
			}
			return intersects ;
		}
	case ELLIPSE:
		{
//			std::cout << "ellipse" << std::endl ;
			if(g->getRadius() < this->getRadius())
			{  
//				std::cout << g->intersects(this) << std::endl ;
				return g->intersects(this) ;
			}

			std::vector<Segment> segs ;
			bool isInSegments = false ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				isInSegments = true ;
			}
			if(g->getGeometryType() == RECTANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
				isInSegments = true ;
			}
			if(g->getGeometryType() == SEGMENTED_LINE)
			{
				for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
				{
					segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
				}
				isInSegments = true ;
			}
			
			if(isInSegments)
			{
				bool intersects = false ;
				for(size_t i = 0 ; i < segs.size() ; i++)
				{
					intersects = intersects || segs[i].intersects(this) ;
				}
				return intersects ;
			}

			if(g->getGeometryType() == ELLIPSE && (g->in(this->getCenter()) || this->in(g->getCenter())))
				return true ;

/*			if(g->getGeometryType() == ELLIPSE)
			{
				Circle thiscircle(this->getRadius(), this->getCenter()) ;	
				Circle othercircle(g->getRadius(), g->getCenter()) ;
				Geometry * gcircle = &othercircle ;
				return thiscircle.intersects(gcircle) ;
			}*/

			// this intersection is not complete and must be refined... later...
			std::vector<Point> box = this->getBoundingBox() ;
			Segment s0(box[0], box[1]) ;
			Segment s1(box[1], box[2]) ;
			Segment s2(box[2], box[3]) ;
			Segment s3(box[3], box[0]) ;
			
			return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g)  || s3.intersects(g) ;
/*			Circle largecircle(this->getRadius(), this->getCenter()) ;
			return largecircle.intersects(g) ;*/

		}
	case SPHERE:
		{
			if(g->getGeometryType() == SPHERE)
			{
			  if(in(g->getCenter()) || g->in(getCenter()))
				return true ;
				
				double birad = getRadius()+g->getRadius() ;
// 				if(((getCenter().x + birad < g->getCenter().x) ||
// 					return false ;
// 				if(getCenter().x - birad > g->getCenter().x))
// 					return false ;
// 				if((getCenter().y + birad < g->getCenter().y) ||
// 					return false ;
// 				if(getCenter().y - birad > g->getCenter().y))
// 					return false ;
// 				if((getCenter().z + birad < g->getCenter().z) ||
// 					return false ;
// 				if(getCenter().z - birad > g->getCenter().z)))
// 					return false ;

				return dist(getCenter(), g->getCenter()) < birad ;

			}
			if(g->getGeometryType() == HEXAHEDRON)
			{
				std::vector<Point> bbox = g->getBoundingBox() ;
				if(g->in(getCenter()))
					return (getCenter().x+getRadius() < bbox[7].x) ||
				         (getCenter().x-getRadius() > bbox[0].x)    ||
				         (getCenter().y+getRadius() < bbox[7].y)    ||
				         (getCenter().y-getRadius() > bbox[0].y)    ||
				         (getCenter().z+getRadius() < bbox[7].z)    ||
					     (getCenter().z-getRadius() > bbox[0].z)     ;

				if(!g->in(getCenter()))
					return (getCenter().x+getRadius() > bbox[0].x) ||
					  (getCenter().x-getRadius() < bbox[7].x)    ||
					  (getCenter().y+getRadius() > bbox[0].y)    ||
					  (getCenter().y-getRadius() < bbox[7].y)    ||
					  (getCenter().z+getRadius() > bbox[0].z)    ||
					  (getCenter().z-getRadius() < bbox[7].z) ;
			}
			if(g->getGeometryType() == TETRAHEDRON)
			{
				int incount = 0 ;
				int outcount = 0 ;
				for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
				{
					if(in(g->getBoundingPoint(i)))
						incount++ ;
					else
						outcount++ ;
				}
				
				if(incount && outcount)
					return true ;
				
				Segment s(g->getCenter(), getCenter()) ;
				TriPoint t0(&g->getBoundingPoint(0), &g->getBoundingPoint(1), &g->getBoundingPoint(2)) ;
				TriPoint t1(&g->getBoundingPoint(0), &g->getBoundingPoint(1), &g->getBoundingPoint(3)) ;
				TriPoint t2(&g->getBoundingPoint(0), &g->getBoundingPoint(2), &g->getBoundingPoint(3)) ;
				TriPoint t3(&g->getBoundingPoint(1), &g->getBoundingPoint(2), &g->getBoundingPoint(3)) ;
				if(s.intersects(t0))
					return in(s.intersection(t0)[0]) ;
				if(s.intersects(t1))
					return in(s.intersection(t1)[0]) ;
				if(s.intersects(t2))
					return in(s.intersection(t2)[0]) ;
				if(s.intersects(t3))
					return in(s.intersection(t3)[0]) ;

				return false ;
				
			}
			return false ;
			
		}
	case HEXAHEDRON:
		{
			if(g->getGeometryType() == SPHERE)
			{
				std::vector<Point> bbox = getBoundingBox() ;
				return (((g->getCenter().x+g->getRadius() < bbox[7].x) ||
				         (g->getCenter().x-g->getRadius() > bbox[0].x) ||
				         (g->getCenter().y+g->getRadius() < bbox[7].y) ||
				         (g->getCenter().y-g->getRadius() > bbox[0].y) ||
				         (g->getCenter().z+g->getRadius() < bbox[7].z) ||
				         (g->getCenter().z-g->getRadius() > bbox[0].z)    )&&in(g->getCenter())) ||
					(((g->getCenter().x+g->getRadius() > bbox[0].x) ||
					  (g->getCenter().x-g->getRadius() < bbox[7].x)    ||
					  (g->getCenter().y+g->getRadius() > bbox[0].y)    ||
					  (g->getCenter().y-g->getRadius() < bbox[7].y)    ||
					  (g->getCenter().z+g->getRadius() > bbox[0].z)    ||
					  (g->getCenter().z-g->getRadius() > bbox[7].z)    )&& !in(g->getCenter()))
					;
			}
			if(g->getGeometryType() == HEXAHEDRON)
			{
				size_t count = 0 ;
				size_t count_alt = 0 ;
				
				std::vector<Point> myBB = getBoundingBox() ;
				std::vector<Point> gBB = g->getBoundingBox() ;
				
				for(size_t i = 0 ; i < 8 ; i++)
				{
					if(in(gBB[i]))
						count++ ;
					if(g->in(myBB[i]))
						count_alt++ ;
				}
				
				return count != 0 && count !=8 && count_alt != 0 && count_alt !=8;

			}
			return false ;
		}
	case TETRAHEDRON:
		{
			if(g->getGeometryType() == SPHERE)
			{
				Segment s(g->getCenter(), getCenter()) ;
				TriPoint t0(&getBoundingPoint(0), &getBoundingPoint(1), &getBoundingPoint(2)) ;
				TriPoint t1(&getBoundingPoint(0), &getBoundingPoint(1), &getBoundingPoint(3)) ;
				TriPoint t2(&getBoundingPoint(0), &getBoundingPoint(2), &getBoundingPoint(3)) ;
				TriPoint t3(&getBoundingPoint(1), &getBoundingPoint(2), &getBoundingPoint(3)) ;
				if(s.intersects(t0))
					return g->in(s.intersection(t0)[0]) ;
				if(s.intersects(t1))
					return g->in(s.intersection(t1)[0]) ;
				if(s.intersects(t2))
					return g->in(s.intersection(t2)[0]) ;
				if(s.intersects(t3))
					return g->in(s.intersection(t3)[0]) ;

				return false ;
			}
		}
	default:
		{
			return false ;
		}
	}
}

std::vector<Point> Geometry::intersection(const Geometry * g) const
{
	std::vector<Point> ret ;

	switch(getGeometryType())
	{
	case TRIANGLE :
		{
			Segment s0(this->getBoundingPoint(0), 
			           this->getBoundingPoint(this->getBoundingPoints().size()/3)) ;
			Segment s1(this->getBoundingPoint(this->getBoundingPoints().size()/3),
			           this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
			Segment s2(this->getBoundingPoint(0), 
			           this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
			std::vector<Point> intersection = s0.intersection(g) ;
			ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			intersection = s1.intersection(g) ;
			ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			intersection = s2.intersection(g) ;
			ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			return ret ;
		}
	case RECTANGLE:
		{
			std::vector<Point> box = this->getBoundingBox() ;
			Segment s0(box[0], box[1]) ;
			Segment s1(box[1], box[2]) ; 
			Segment s2(box[2], box[3]) ; 
			Segment s3(box[3], box[0]) ;
			if(g->getGeometryType() == RECTANGLE)
			{
				std::vector<Point> intersection = s0.intersection(g) ;
				std::vector<Point> it = s1.intersection(g) ;
				intersection.insert(intersection.end(), it.begin(), it.end()) ;
				it = s2.intersection(g) ;
				intersection.insert(intersection.end(), it.begin(), it.end()) ;
				it = s3.intersection(g) ;
				intersection.insert(intersection.end(), it.begin(), it.end()) ;
				
				std::sort(intersection.begin(), intersection.end()) ;
				std::vector<Point>:: iterator e = std::unique(intersection.begin(), intersection.end()) ;
				intersection.erase(e, intersection.end()) ;
				return intersection ;
			}
			std::vector<Point> intersection = s0.intersection(g) ;
			
			double perimetre = s0.norm()+s1.norm()+s2.norm()+s3.norm() ;
			
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)4*g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}

			}
			if(!intersection.empty())
				ret.push_back(intersection.back()) ;

			intersection = s1.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)4*g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection.back()) ;

			intersection = s2.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)4*g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection.back()) ;

			intersection = s3.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)4*g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection.back()) ;

			std::stable_sort(ret.begin(), ret.end()) ;
			std::vector<Point>:: iterator e = std::unique(ret.begin(), ret.end()) ;
			ret.erase(e, ret.end()) ;
			
			return ret ;
		}
	case SEGMENTED_LINE:
		{
			std::vector<Segment> segs ;
			for(size_t i = 0 ; i < this->getBoundingPoints().size()-1 ; i++)
			{
				segs.push_back(Segment(this->getBoundingPoint(i), this->getBoundingPoint(i+1))) ;
			}
			for(size_t i = 0 ; i < segs.size() ; i++)
			{
				std::vector<Point> intersection = segs[i].intersection(g) ;
				ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			}
			return ret ;
		}
	case CIRCLE:
		{
			
			if(g->getGeometryType() == CIRCLE)
			{
				double r = getRadius() ;
				double R = g->getRadius() ;


				//first quadratic equation
				double a = 1. ;
				double b = 1. - 2.*R*R ;
				double c = R*R*R*R - r*r ;
				
				double delta = b*b - 4.*a*c ;
				
				if(delta < 0)
					return ret ;
				
				double y_squared_0 = (-b + sqrt(delta))/(2.*a) ;
				if(delta == 0)
				{
					if(y_squared_0 < 0)
						return ret ;
					
					
					
					if((r*r-y_squared_0) >= 0)
					{
						double y = sqrt(y_squared_0) ;
						double x = sqrt(r*r-y_squared_0) ;
						Point A(x,y) ;
						if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						   std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
							ret.push_back(A) ;
						
						Point B(-x,y) ;
						if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						   std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
							ret.push_back(B) ;
						
						Point C(-x,-y) ;
						if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						   std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
							ret.push_back(C) ;
						
						Point D(x,-y) ;
						if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						   std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
							ret.push_back(D) ;
						

					}
					
					return ret ;
				}

				double y_squared_1 = (-b - sqrt(delta))/(2.*a) ;

				if(y_squared_0 >= 0)
				{
					if((r*r-y_squared_0) >= 0)
					{
						double x = sqrt(r*r-y_squared_0) ;
						double y = sqrt(y_squared_0) ;
						
						Point A(x,y) ;
						if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						)
							ret.push_back(A) ;
						
						Point B(-x,y) ;
						if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						)
							ret.push_back(B) ;
						
						Point C(-x,-y) ;
						if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						)
							ret.push_back(C) ;
						
						Point D(x,-y) ;
						if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
						std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE 
						)
							ret.push_back(D) ;
					}
				}
				if(y_squared_1 >= 0)
				{
					if((r*r-y_squared_1) >= 0)
					{
						
						double x = sqrt(r*r-y_squared_1) ;
						double y = sqrt(y_squared_1) ;
						
						Point A(x,y) ;
						if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE &&
						std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE 
						)
							ret.push_back(A) ;
						
						Point B(-x,y) ;
						if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE &&
						std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE 
						)
							ret.push_back(B) ;
						
						Point C(-x,-y) ;
						if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE &&
						std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE 
						)
							ret.push_back(C) ;
						
						Point D(x,-y) ;
						if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE &&
						std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE 
						)
							ret.push_back(D) ;
					}
				}
				
				return ret ;
					

			}
			
			std::vector<Segment> segs ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3), g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
			}
			if(g->getGeometryType() == RECTANGLE)
			{
				return g->intersection(this) ;
			}
			if(g->getGeometryType() == SEGMENTED_LINE)
			{
				for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
				{
					segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
				}
			}
			
			for(size_t i = 0 ; i < segs.size() ; i++)
			{
				std::vector<Point>intersection = segs[i].intersection(this) ;
				ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			}
			return ret ;
		}
	case ELLIPSE:
		{

			std::vector<Segment> segs ;
			bool isInSegments = false ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				segs.push_back(Segment(g->getBoundingPoint(0),
				                       g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
				isInSegments = true ;
			}
			if(g->getGeometryType() == RECTANGLE)
			{
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
				segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
				isInSegments = true ;
			}
			if(g->getGeometryType() == SEGMENTED_LINE)
			{
				for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
				{
					segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
				}
				isInSegments = true ;
			}
			
			if(isInSegments)
			{
				for(size_t i = 0 ; i < segs.size() ; i++)
				{
//					std::cout << "in Geometry::intersection" << std::endl ;
//					segs[i].first().print() ;
//					segs[i].second().print() ;
//					std::cout << "in Segment::intersection" << std::endl ;
					std::vector<Point> vtemp = segs[i].intersection(this) ;
					for(size_t j = 0 ; j < vtemp.size() ; j++)
						{
//							vtemp[j].print() ;
							ret.push_back(vtemp[j]) ;
						}
				}
				return ret ;
			}

//			std::vector<Segment> segs ;
			for(size_t i = 0 ; i < getBoundingPoints().size() - 1 ; i++)
			{
//				getBoundingPoint(i).print() ;
				segs.push_back(Segment(getBoundingPoint(i),getBoundingPoint(i+1))) ;
			}
			segs.push_back(Segment(getBoundingPoint(getBoundingPoints().size()-1),getBoundingPoint(0))) ;
			
			for(size_t i = 0 ; i < segs.size() ; i++)
			{
				std::vector<Point>intersection = segs[i].intersection(g) ;
//				segs[i].first().print() ;
//				segs[i].second().print() ;
				for(size_t j = 0 ; j < intersection.size() ; j++)
				{
//					intersection[j].print() ;
					ret.push_back(intersection[j]) ;
				}
//				ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
			}
			return ret ;

		}
	case SPHERE:
		{
			if(g->getGeometryType() == SPHERE)
			{
				double dc = dist(getCenter(), g->getCenter()) ;
				double r0 = getRadius() ;
				double r1 = g->getRadius() ;
				
				if (dc < 1e-8)
					return ret ;
				
				double cosAngleInSelf = ((r0*r0+dc*dc-r1*r1)/(2.*r0*dc)) ;
				
				double d = cosAngleInSelf*r0/*0.5*getRadius() +0.5*g->getRadius() - dc*/;
				
				Point v = g->getCenter()-getCenter() ;

				Point centerOfIntersection = getCenter() + v*(d/dc) ;
				double radiusOfIntersection = sqrt(r0*r0 - d*d) ;
				
				OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
				
				size_t num_points = (size_t)(2*sqrt(g->getBoundingPoints().size())) ;
				
				if(num_points < 3 )
					return ret ;
				
				C.sampleBoundingSurface(num_points) ;
				
				for(size_t i = 0 ;  i < num_points ; i++)
				{
					ret.push_back(C.getBoundingPoint(i)) ;
				}

				return ret ;
				
			}
			if(g->getGeometryType() == HEXAHEDRON)
			{
				std::vector<Point> bbox = g->getBoundingBox() ;
				double maxx =  bbox[0].x ;
				double minx =  bbox[0].x ;
				double maxy =  bbox[0].y ;
				double miny =  bbox[0].y ;
				double maxz =  bbox[0].z ;
				double minz =  bbox[0].z ;
				
				for(size_t i = 1 ; i < bbox.size() ; i++)
				{
					if(bbox[i].x > maxx)
						maxx = bbox[i].x ;
					if(bbox[i].x < minx)
						minx = bbox[i].x ;
					if(bbox[i].y > maxy)
						maxy = bbox[i].y ;
					if(bbox[i].y < miny)
						miny = bbox[i].y ;
					if(bbox[i].z > maxz)
						maxz = bbox[i].z ;
					if(bbox[i].z < minz)
						minz = bbox[i].z ;
				}
				
				
				if(std::abs(center.x - minx) <= getRadius())
				{
					double d = std::abs(center.x - minx) ;
					Point v(-1,0,0) ;
					Point centerOfIntersection(minx, center.y, center.z) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;					
						else if(C.getBoundingPoint(i).y > maxy)
							ret.push_back(Point(C.getBoundingPoint(i).x, maxy, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).y < miny)
							ret.push_back(Point(C.getBoundingPoint(i).x, miny,C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).z > maxz)
							ret.push_back(Point(C.getBoundingPoint(i).x, C.getBoundingPoint(i).y, maxz)) ;
						else if(C.getBoundingPoint(i).z < minz)
							ret.push_back(Point(C.getBoundingPoint(i).x, C.getBoundingPoint(i).y, minz)) ;
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
				if(std::abs(center.x - maxx) < getRadius())
				{
					double d = std::abs(center.x - maxx) ;
					Point v(1,0,0) ;
					Point centerOfIntersection(maxx, center.y, center.z) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ; ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
	
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
				
				if(std::abs(center.y - miny) < getRadius())
				{
					double d = std::abs(center.y - miny) ;
					Point v(0,-1,0) ;
					Point centerOfIntersection(center.x, miny, center.z) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
					C.sampleSurface(num_points) ;
					
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
						else if(C.getBoundingPoint(i).x > maxx)
							ret.push_back(Point(maxx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).x < minx)
							ret.push_back(Point(minx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).z > maxz)
							ret.push_back(Point(C.getBoundingPoint(i).x, C.getBoundingPoint(i).y, maxz)) ;
						else if(C.getBoundingPoint(i).z < minz)
							ret.push_back(Point(C.getBoundingPoint(i).x, C.getBoundingPoint(i).y, minz)) ;
					}
					
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
				if(std::abs(center.y - maxy) < getRadius())
				{
					double d = std::abs(center.y - maxy) ;
					Point v(0,1,0) ;
					Point centerOfIntersection(center.x, maxy, center.z) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
						else if(C.getBoundingPoint(i).x > maxx)
							ret.push_back(Point(maxx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).x < minx)
							ret.push_back(Point(minx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).z > maxz)
							ret.push_back(Point(C.getBoundingPoint(i).x,  C.getBoundingPoint(i).y,maxz)) ;
						else if(C.getBoundingPoint(i).z < minz)
							ret.push_back(Point(C.getBoundingPoint(i).x,  C.getBoundingPoint(i).y,minz)) ;
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
				
				if(std::abs(center.z - minz) < getRadius())
				{
					double d = std::abs(center.z - minz) ;
					Point v(0,0,-1) ;
					Point centerOfIntersection(center.x, center.y, minz) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
						else if(C.getBoundingPoint(i).x > maxx)
							ret.push_back(Point(maxx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).x < minx)
							ret.push_back(Point(minx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).y > maxy)
							ret.push_back(Point(C.getBoundingPoint(i).x, maxy, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).y < miny)
							ret.push_back(Point(C.getBoundingPoint(i).x, miny, C.getBoundingPoint(i).z)) ;
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
				if(std::abs(center.z - maxz) < getRadius())
				{
					double d = std::abs(center.z - maxz) ;
					Point v(0,0,1) ;
					Point centerOfIntersection(center.x, center.y, maxz) ;
					double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
					OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;
					
					size_t num_points = (size_t)round(7.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
						else if(C.getBoundingPoint(i).x > maxx)
							ret.push_back(Point(maxx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).x < minx)
							ret.push_back(Point(minx,C.getBoundingPoint(i).y, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).y > maxy)
							ret.push_back(Point(C.getBoundingPoint(i).x, maxy, C.getBoundingPoint(i).z)) ;
						else if(C.getBoundingPoint(i).y < miny)
							ret.push_back(Point(C.getBoundingPoint(i).x, miny, C.getBoundingPoint(i).z)) ;
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
			}
			if(g->getGeometryType() == TETRAHEDRON)
			{
				Segment s0(g->getBoundingPoint(0), g->getBoundingPoint(1)) ;
				Segment s1(g->getBoundingPoint(0), g->getBoundingPoint(2)) ;
				Segment s2(g->getBoundingPoint(0), g->getBoundingPoint(3)) ;
				Segment s3(g->getBoundingPoint(1), g->getBoundingPoint(2)) ;
				Segment s4(g->getBoundingPoint(1), g->getBoundingPoint(3)) ;
				Segment s5(g->getBoundingPoint(2), g->getBoundingPoint(3)) ;
				std::vector<Segment *> intersectingSegments ;
				if(s0.intersects(this))
				{
					std::vector<Point> inter = s0.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				if(s1.intersects(this))
				{
					std::vector<Point> inter = s1.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				if(s2.intersects(this))
				{
					std::vector<Point> inter = s2.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				if(s3.intersects(this))
				{
					std::vector<Point> inter = s3.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				if(s4.intersects(this))
				{
					std::vector<Point> inter = s4.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				if(s5.intersects(this))
				{
					std::vector<Point> inter = s5.intersection(this) ;
					ret.insert(ret.end(), inter.begin(), inter.end()) ;
				}
				
				TriPoint t0(&g->getBoundingPoint(0), &g->getBoundingPoint(1), &g->getBoundingPoint(2)) ;
				TriPoint t1(&g->getBoundingPoint(0), &g->getBoundingPoint(1), &g->getBoundingPoint(3)) ;
				TriPoint t2(&g->getBoundingPoint(0), &g->getBoundingPoint(2), &g->getBoundingPoint(3)) ;
				TriPoint t3(&g->getBoundingPoint(1), &g->getBoundingPoint(2), &g->getBoundingPoint(3)) ;
				
				Plane p0(*t0.point[0], t0.normal) ;
				Plane p1(*t1.point[0], t1.normal) ;
				Plane p2(*t2.point[0], t2.normal) ;
				Plane p3(*t3.point[0], t3.normal) ;
				
// 				if(p0.intersects(this))
// 				{
// 					std::vector<Point> inter = p0.intersection(this) ;
// 					for(size_t i = 0 ; i < inter.size() ; i++)
// 					{
// 						if(t0.in(inter[i]))
// 							ret.push_back(inter[i]) ;
// 					}
// 				}
// 				if(p1.intersects(this))
// 				{
// 					std::vector<Point> inter = p1.intersection(this) ;
// 					for(size_t i = 0 ; i < inter.size() ; i++)
// 					{
// 						if(t1.in(inter[i]))
// 							ret.push_back(inter[i]) ;
// 					}
// 				}
// 				if(p2.intersects(this))
// 				{
// 					std::vector<Point> inter = p2.intersection(this) ;
// 					for(size_t i = 0 ; i < inter.size() ; i++)
// 					{
// 						if(t2.in(inter[i]))
// 							ret.push_back(inter[i]) ;
// 					}
// 				}
// 				if(p3.intersects(this))
// 				{
// 					std::vector<Point> inter = p3.intersection(this) ;
// 					for(size_t i = 0 ; i < inter.size() ; i++)
// 					{
// 						if(t3.in(inter[i]))
// 							ret.push_back(inter[i]) ;
// 					}
// 				}
			}
			return ret ;
		}
	case TETRAHEDRON:
		{
			return g->intersection(this) ;
		}
	case HEXAHEDRON:
		{

			if(g->getGeometryType() == SPHERE)
				return g->intersection(this) ;
		}
	default:
		{
			return ret ;
		}
	}
	
}


ConvexGeometry::ConvexGeometry(): ConvexPolygon((size_t)0)
{
}

ConvexGeometry::ConvexGeometry(size_t s) : ConvexPolygon(s)
{
}

/*bool ConvexGeometry::intersects(const Geometry * g) const
{
	std::cout << "convex" << std::endl ;
	Point inbetween = (this->getCenter() + g->getCenter()) * 0.5 ;
	Point thispoint(this->getCenter()) ;
	double tempdist = 1. / (inbetween - this->getCenter()).norm() ;
	double maxrad = 1.1 * std::max(this->getRadius(), g->getRadius()) ;
	thispoint.x = thispoint.x - (inbetween - this->getCenter()).y * tempdist * maxrad * 2 + inbetween.x ;
	thispoint.y = thispoint.y + (inbetween - this->getCenter()).x * tempdist * maxrad * 2+ inbetween.y ;
	Point nextpoint(thispoint) ;
	g->project(&nextpoint) ;
	if(this->in(nextpoint))
		return true ;
	double lastdist = 0. ;
	double thisdist = squareDist2D(thispoint,nextpoint) ;
	bool isthis = true ;
	while(std::abs(lastdist - thisdist) > POINT_TOLERANCE)
	{
//		if(isthis)
//			this->project(&nextpoint) ;
//		else
//			g->project(&nextpoint) ;
		lastdist = thisdist ;
		thispoint = nextpoint ;
		if(isthis)
		{
			g->project(&nextpoint) ;
			if(this->in(nextpoint))
				return true ;
		}
		else
		{
			this->project(&nextpoint) ;
			if(g->in(nextpoint))
				return true ;
		}
		thisdist = squareDist2D(thispoint,nextpoint) ;
		isthis = !isthis ;
	}
	return (thisdist < POINT_TOLERANCE) ;
}*/

const Point & Geometry::getCenter() const
{
	return this->center;
}

size_t & Geometry::timePlanes()
{
	return this->time_planes ;
}

size_t Geometry::timePlanes() const
{
	return this->time_planes ;
}

Point & Geometry::getCenter() 
{
	return this->center;
}

PointSet::PointSet(size_t npoints) 
{
	this->boundingPoints.resize(npoints) ; 
	for(size_t i = 0 ; i < npoints ;i++)
		boundingPoints[i] = NULL ;
	this->chullEndPos = 0;
} ;

size_t ConvexGeometry::size() const
{
	return 	boundingPoints.size() + inPoints.size() ;
}

size_t NonConvexGeometry::size() const
{
	return 	boundingPoints.size() + inPoints.size() ;
}

double PointSet::x(size_t i) 
{
	return boundingPoints[i]->x ; 
}

double PointSet::y(size_t i) 
{ 
	return boundingPoints[i]->y ; 
}

double PointSet::z(size_t i) 
{ 
	return boundingPoints[i]->z ; 
}

void PointSet::setX(size_t i, double vv) 
{ 
	boundingPoints[i]->Point::setX(vv) ; 
}

void PointSet::setY(size_t i, double vv) 
{ 
	boundingPoints[i]->Point::setY(vv) ; 
}

void PointSet::setZ(size_t i, double vv) 
{ 
	boundingPoints[i]->Point::setZ(vv) ; 
}

void PointSet::set(size_t i, Point * p) 
{ 
	delete boundingPoints[i] ;
	boundingPoints[i] = p; 
}

void PointSet::set(size_t i, double x, double y) 
{ 
	int id = -1 ;
	if(boundingPoints[i] != NULL)
		id = boundingPoints[i]->id ;
	
	delete boundingPoints[i] ;
	
	boundingPoints[i] = new Point(x,y) ; 
	boundingPoints[i]->id = id ;
}

void PointSet::set(size_t i, double x, double y, double z) 
{ 
	int id = -1 ;
	if(boundingPoints[i] != NULL)
		id = boundingPoints[i]->id ;
	
	delete boundingPoints[i] ;
	
	boundingPoints[i] = new Point(x,y,z) ; 
	boundingPoints[i]->id = id ;
}

Point * PointSet::operator[](size_t i) 
{ 
	return boundingPoints[i] ; 
}

Point * PointSet::operator[](size_t i) const
{ 
	return boundingPoints[i] ; 
}

Point * PointSet::getPoint(size_t i) const
{ 
	return boundingPoints[i] ; 
}

Point * PointSet::getPoint(size_t i)
{ 
	return boundingPoints[i] ; 
}

PointSet::const_iterator PointSet::begin() const
{
	return &boundingPoints[0] ; 
}

PointSet::const_iterator PointSet::end() const
{
	return  &boundingPoints[boundingPoints.size()] ;
}

PointSet::iterator PointSet::begin()
{
	return &boundingPoints[0] ; 
}

PointSet::iterator PointSet::end()
{
	return  &boundingPoints[boundingPoints.size()] ;
}

ConvexPolygon * PointSet::convexHull() const
{
		//!have we allready done that ?
	
	if(this->chullEndPos != 0 )
	{
		ConvexPolygon *hull = new ConvexPolygon(chullEndPos) ;
		std::copy(begin(), end(), hull->begin()) ;
		return hull ;
	}
	
	ConvexPolygon *ret = new ConvexPolygon(this) ;
// 	chullEndPos = ret->size() ;
	
	return ret ;
	
}


/*bool PointSet::in(const Point & p)  const 
{
	ConvexPolygon * hull = convexHull() ;
	bool ret = hull->in(p) ;
	delete hull ;
	return ret ;
}*/

size_t PointSet::size() const
{ 
	return  boundingPoints.size() ;
}

Point PointSet::computeCenter() const 
{
	Point ret;
	for (size_t i = 0 ; i < boundingPoints.size() ; i++)
	{
		ret.x += boundingPoints[i]->x/boundingPoints.size() ;
		ret.y += boundingPoints[i]->y/boundingPoints.size() ;
	}
	return ret ;
}

void PointSet::removePoint(size_t index)
{
	std::valarray<Point *> n(size()-1) ;
	//std::copy((*this)[0], (*this)[index-1], (*n)[0]) ;
	for (size_t i = 0 ; i < index ; i++)
	{
		n[i] = boundingPoints[i] ;
	}
	for (size_t i = index+1 ; i < size() ; i++)
	{
		n[i] = boundingPoints[i] ;
	}
	
	boundingPoints.resize(n.size()) ;
	boundingPoints = n ;
}



Nurb::Nurb(std::vector<double> knot, size_t degree)
{
	this->knot = knot ;
	this->degree = degree;
}

Nurb::Nurb(std::vector<Point> controlPoint, std::vector<double> weight,std::vector<double> knot, size_t degree)
{
	this->controlPoint = controlPoint ;
	this->weight = weight ;
	this->knot = knot ;
	this->degree = degree ;
}

Nurb::Nurb(std::vector<Point> controlPoint, std::vector<double> weight, std::vector<double> knot)
{
	this->controlPoint = controlPoint ;
	this->weight = weight ;
	this->knot = knot ;
	
}

int Nurb::computeDegree(const std::vector<double> &vec,const std::vector<Point> &Pv)
{
	return vec.size() - Pv.size() - 1;
}

std::vector<double> Nurb::transVector(std::vector<double> &vec)
{
	double diff = vec.back() - vec.front();
	for(size_t i = 0 ; i != vec.size(); ++i)
	{
		vec[i] = ( vec[i] - vec.front() ) / diff ;
	}
	return vec;
}

double Nurb::getBasis(double u, int i, int k, const std::vector<double> &vec)
{
	if (k == 0)
	{
		if (vec[i] <= u && u <= vec[i+1])
		{
		return 1;
		}
	return 0;
	}
	double denom1 = vec[i+k]-vec[i];
	double denom2 = vec[i+k+1]-vec[i+1];
	double part1 = 0; double part2 = 0;
	if (denom1 > 0)
	{
	part1 = ((u-vec[i])/denom1)*getBasis(u,i,k-1,vec);
	}
	if (denom2 > 0)
	{
	part2 = ((vec[i+k+1]-u)/denom2)*getBasis(u,i+1,k-1,vec);
	}
	return part1 + part2;
}

Point Nurb::getNurbPoint(double u, double r)
{
	sort(knot.begin(),knot.end());
	transVector(knot);
	int k = computeDegree(knot,controlPoint);
	Point outPointNom;
	double outDenom = 0;
	for(size_t i = 0 ; i != controlPoint.size() ; ++i)
	{
		double koef = getBasis(u,i,k,knot);
		if(koef > r)
		{
			outPointNom += controlPoint[i]*weight[i]*koef ;
			outDenom += weight[i]*koef ;
		}
	}
	(outPointNom / outDenom).print();
	return outPointNom / outDenom;
}

std::vector<Point> Nurb::sampleNurb(size_t n)
{	
 	std::vector<Point>  nurbPoint;
 	double step = 1./n;
 	for(double u = 0; u != 1; u += step) 
 	{
 		nurbPoint.push_back(getNurbPoint(u));
 	}
	return nurbPoint;
}

void WeightedPoint::print() const 
{
	std::cout << " ( id = " << id << " ; "<< x << ", " << y << ", "<< z << ", "<< w <<") " << std::endl ;
}

WeightedPoint WeightedPoint::operator*(const double & p)  const
{
	return WeightedPoint(x*p, y*p, z*p, w*p); 
}

WeightedPoint WeightedPoint::operator+(const WeightedPoint &wp)  const
{
	return WeightedPoint(x + wp.x, y + wp.y,z + wp.z, w + wp.w); 
}

WeightedPoint WeightedPoint::operator/(const double p) const 
{
	return WeightedPoint(x/p, y/p, z/p, w/p);
}

// WeightedPoint::WeightedPoint() 
// {
// 	bas = Point();
// 	wei = 1;
// }
// 
// WeightedPoint::WeightedPoint(Point P, double w) 
// {
// 	bas = Point(P.x,P.y,P.z);
// 	wei = w;
// }
// 
// void WeightedPoint::print() const 
// {
// 	std::cout << " ( id = " << bas.id << " ; "<< bas.x << ", " << bas.y << ", "<< bas.z << ", "<< wei <<") " << std::endl ;
// }
// 
// WeightedPoint WeightedPoint::operator+(const WeightedPoint &wP) const
// {
// 	
// 	WeightedPoint ret((*this)) ;
// 	ret.bas += wP.bas ; 
// 	ret.wei += wP.wei ; 
// 	return ret ;
// }
// 
// WeightedPoint WeightedPoint::operator*(const double p)  const
// {
// 	WeightedPoint ret((*this)) ;
// 	ret.bas *= p ;
// 	ret.wei *= p ; 
// 	return ret ; 
// }

// Point Nurb::pointOnNurb(double u)
// {
// 	size_t i = 0;
// 
// 	while (u >= knot[i]) 
// 		++i;
// 
// 	size_t leftBound = i - 1;
// 	
// 	size_t numberOfCycles, knotMultiplicity  = 0;
// 
// 	if ( u == knot[leftBound] )
// 	{
// 		knotMultiplicity = count(knot.begin(),knot.end(),u);
// 		numberOfCycles = degree - knotMultiplicity;
// 	}
// 	else numberOfCycles = degree;
// 
// 	size_t r = leftBound - knotMultiplicity + 1;
// 
// 	std::valarray<WeightedPoint> originalWPoint(controlPoint.size());
// 	for (size_t i = 0; i != originalWPoint.size(); ++i)
// 	{	
// 		originalWPoint[i] = WeightedPoint(controlPoint[i]*weight[i],weight[i]);
// 	}
// 
//  	double coeff = 0;
// 	std::valarray<WeightedPoint> generatedWPoint(originalWPoint.size());
// 	for (size_t j = 1; j <= numberOfCycles; ++j)
// 	{
// 		for (size_t k = leftBound - degree + j; k != r ; ++k)
// 		{
// 			coeff = ( u - knot[k] ) / ( knot[k+1+degree-j] - knot[k] );	
// 			generatedWPoint[k] = originalWPoint[k-1] * ( 1 - coeff) + originalWPoint[k] *coeff;
// 		}
// 		
// 
// 		for (size_t p = 0; p != originalWPoint.size(); ++p)
// 		{
// 			originalWPoint[p]=generatedWPoint[p];
// 		}
// 	}
// 	
// 	Point resultPoint = originalWPoint[r - 1].getPoint() / originalWPoint[r - 1].getWeight();
// 	resultPoint.print();
//  	return resultPoint;
// }

Point Nurb::pointOnNurb(double u)
{
	size_t i = 0;

	while (u >= knot[i]) 
		++i;

	size_t leftBound = i - 1;
	
	size_t numberOfCycles, knotMultiplicity  = 0;

	if ( u == knot[leftBound] )
	{
		knotMultiplicity = count(knot.begin(),knot.end(),u);
		numberOfCycles = degree - knotMultiplicity;
	}
	else numberOfCycles = degree;

	size_t r = leftBound - knotMultiplicity + 1;

	std::valarray<WeightedPoint> originalWPoint(controlPoint.size());
	for (size_t i = 0; i != originalWPoint.size(); ++i)
	{	
		originalWPoint[i] = WeightedPoint(controlPoint[i]*weight[i],weight[i]);
	}

 	double coeff = 0;
	std::valarray<WeightedPoint> generatedWPoint(originalWPoint.size());
// 	std::valarray<WeightedPoint> forSubdivision(2*numberOfCycles-1);
	for (size_t j = 1; j <= numberOfCycles; ++j)
	{
		for (size_t k = leftBound - degree + j; k != r ; ++k)
		{
			coeff = ( u - knot[k] ) / ( knot[k+1+degree-j] - knot[k] );	
			generatedWPoint[k] = originalWPoint[k-1] * ( 1 - coeff) + originalWPoint[k] *coeff;
		}
		

		for (size_t k = 0; k != originalWPoint.size(); ++k)
		{
			originalWPoint[k]=generatedWPoint[k];
		}
	}
	
	Point resultPoint = originalWPoint[r - 1] / originalWPoint[r - 1].w;
	resultPoint.print();
 	return resultPoint;
}

// Nurb::Point intersection(const Line *l) const
// {
// double epsilon = 0;
// while (epsilon > 1e-08) 
// {
/*
}*/

Geometry::Geometry(): inPoints(0),gType(NULL_GEOMETRY)
{
	sampled = false ;
	time_planes = 1 ;
}

Geometry::Geometry(size_t numPoints):inPoints(numPoints), gType(NULL_GEOMETRY)
{
	time_planes = 1 ;
}

GeometryType Geometry::getGeometryType() const
{
	return gType ;
}

// ConvexPolygon::ConvexPolygon()
// {
// 	std::cout << "calling default ConvexPolygon constructor" << std::endl ;
// }

ConvexPolygon::ConvexPolygon(size_t npoints) : PointSet(npoints)
{
}

ConvexPolygon::ConvexPolygon(const PointSet * po) : PointSet(po->size())
{
	std::vector<Point*> *points = new std::vector<Point*>(po->size()) ;
	std::copy(po->begin(), po->end(), points->begin()) ;
	
	Point *pivot( (*points)[0] ) ;
	for(size_t i = 1 ; i < points->size() ; i++)
	{
		if (pivot->y < (*points)[i]->y)
			pivot = (*points)[i] ;
	}
	
	//! we build a map ordered by the angle to the pivot.
	
	std::map< double, Point *>  pointSet ;
	
	for(size_t i = 0 ; i < points->size() ; i++)
	{
		if ((*(*points)[i]) != (*pivot))
		{
			double angle = atan2(pivot->y-(*points)[i]->y, pivot->x-(*points)[i]->x) ;
			
			if(pointSet.find(angle) != pointSet.end())
			{
				if(squareDist((*pointSet[angle]), (*pivot)) < squareDist((*(*points)[i]), (*pivot)))
				{
					pointSet[angle] = (*points)[i] ;
				}
			}
			else
				pointSet[angle] = (*points)[i] ;
		}
	}
	
	pointSet[-10] = pivot ;
	
	//!we build the convex hull : we are on the convex hull if and only if we only turned left.
	
	std::vector<Point *> temphull ;
	std::vector<Point *> orderedset ;
	
	for(std::map<double, Point*>::const_iterator i = pointSet.begin() ; i!= pointSet.end() ; ++i)
	{
		orderedset.push_back(i->second) ;
	}
	
	temphull.push_back(orderedset[0] ) ;
	temphull.push_back(orderedset[1] ) ;
	
	for(std::vector< Point *>::iterator i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
	{
		for(size_t j = 0 ; j < temphull.size() ; j++)
			
		//! this is a usual cross product of two vectors...
			if(  ((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			     ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y - temphull[temphull.size()-2]->y ) > 
			     POINT_TOLERANCE )
			{
				temphull.push_back(*i) ;
			}	
		else
		{
			while( !(((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			         ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y -temphull[temphull.size()-2]->y ) > 
			         POINT_TOLERANCE))
			{
				temphull.pop_back();
			}
			temphull.push_back(*i) ;
		}
		
		
	}
	
	delete points ;
	std::copy(temphull.begin(), temphull.end(),this->begin()) ;
}


/*bool ConvexPolygon::in(const Point & p) const 
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
	
}*/


bool ConvexPolygon::isTrigoOriented()  const
{
	
	for (size_t i = 0 ;  i <  boundingPoints.size()-1; i++)
	{
		Point v_0 = *boundingPoints[(i+1)%boundingPoints.size()] 
			- *boundingPoints[(i)%boundingPoints.size()] ;
		Point v_1 = *boundingPoints[(i+2)%boundingPoints.size()] 
			- *boundingPoints[(i+1)%boundingPoints.size()] ;
		if((v_0^v_1).z < 0 )
			return false ;
	}
	
	return true ;
	
}

NonConvexGeometry::NonConvexGeometry() : PointSet(1)
{
	orderedSet.resize(1) ;
	orderedSet[0] = inPoints[0] ;
}

NonConvexGeometry::NonConvexGeometry(size_t numPoints) : PointSet(numPoints)
{
	orderedSet.resize(numPoints) ;
	for(size_t i = 0 ; i < numPoints ; i++)
		orderedSet[i] = boundingPoints[i] ;
}

NonConvexGeometry::NonConvexGeometry(const PointArray & p)
{
	boundingPoints.resize(p.size()) ;
	std::copy(&p[0], &p[p.size()], begin()) ;//&boundingPoints[0]) ;
}

const Point & NonConvexGeometry::getPoint(size_t i) const
{
	if (i < inPoints.size())
		return *inPoints[i] ;
	return *boundingPoints[i-inPoints.size()] ;
}

Point & NonConvexGeometry::getPoint(size_t i) 
{
	if (i < inPoints.size())
		return *inPoints[i] ;
	return *boundingPoints[i-inPoints.size()] ;
}

Line::Line(const Point & origin, const Point &vector)
{
	p = origin ;
	v = vector ;
}

Line::Line(const Segment & base)
{
	p = base.first() ;
	v = base.vector() ;
}

Line::Line(const Line & l, const Point & through)
{
	p = through ;
	v = l.vector() ;
}

Line::Line()
{
	p = Point(0,0) ;
	v = Point(1,1) ;
}

bool Line::intersects(const Line &l) const
{
	return v.x * l.vector().y - v.y * l.vector().x != 0 ;
}

bool Line::intersects(const Segment &s) const
{
	
	if (std::abs(-v.y * s.vector().x + v.x * s.vector().y) <= std::numeric_limits<double>::epsilon())
		return false ;
	
	Matrix m(2,2) ;
	Vector vv(2) ;
	
	m[0][0] = s.vector().x ; m[0][1] = -v.x ;
	m[1][0] = s.vector().y ; m[1][1] = -v.y ;
	
	vv[0] = p.x-s.second().x  ; vv[1] = p.y-s.second().y   ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * vv ;
	
	return fac[0] < 1 && fac[0] > 0 ;
}

bool Line::on(const Point &m) const
{
	return isAligned(m, p, (p+v)) ;
}

const Point & Line::vector() const
{
	return v ;
}

const Point & Line::origin() const
{
	return p ;
}

Point Line::intersection(const Line &l) const
{
	Matrix m(2,2) ;
	Vector vec(2) ;
	
	m[0][0] = v.x ; m[0][1] = -l.vector().x ;
	m[1][0] = v.y ; m[1][1] = -l.vector().y ;
	
	vec[0] = l.origin().x - p.x ; vec[1] = l.origin().y - p.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * vec ;

	return p + v*fac[0];
	
// 	double t = 0;
// 	if(v.x != 0 && v.y != 0)
// 		t = ((p.y - l->origin()->y) + v.y/v.x * (p.x - l->origin()->x)) / (v.x/v.y - l->vector()->y) ;
// 	else if (v.x == 0)
// 		t = ( p.x - l->origin()->x ) / l->vector()->x ;
// 	else if (v.y == 0)
// 		t = ( p.y - l->origin()->y ) / l->vector()->y ;
// 	
// 	return (*l->origin()) + (*l->vector())*t ;
}

Point Line::intersection(const Segment &s) const
{
	Matrix m(2,2) ;
	Vector vec(2) ;
	
	m[0][0] = v.x ; m[0][1] = -s.vector().x ;
	m[1][0] = v.y ; m[1][1] = -s.vector().y ;
	
	vec[0] = s.first().x - p.x ; vec[1] = s.first().y - p.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * vec ;
	
	return p + v*fac[0];
}


bool Line::intersects(const TriPoint &g) const
{
	if(isCoplanar(g.point[1],g.point[0], g.point[2],&p))
		return true ;

	Matrix mat(3,3) ;
	
	mat[0][0] = v.x; mat[0][1] = g.point[1]->x-g.point[0]->x; mat[0][2] = g.point[2]->x-g.point[0]->x; 
	mat[1][0] = v.y; mat[1][1] = g.point[1]->y-g.point[0]->y; mat[1][2] = g.point[2]->y-g.point[0]->y; 
	mat[2][0] = v.z; mat[2][1] = g.point[1]->z-g.point[0]->z; mat[2][2] = g.point[2]->z-g.point[0]->z; 
	return abs(det(mat)) > POINT_TOLERANCE ;


}

std::vector<Point> Line::intersection(const TriPoint &s) const
{
	Point u = *s.point[1]-*s.point[0] ;
	Point v = *s.point[2]-*s.point[0] ;
	double a = -(s.normal*(p-*s.point[0])) ;
	double b = v*s.normal ;
	
	if(b < POINT_TOLERANCE)
	{
		if(abs(a) > POINT_TOLERANCE)
			return std::vector<Point>(0) ;
		
		Triangle t(*s.point[0], *s.point[1], *s.point[2]) ;
		return this->intersection(&t) ;
	}
	
	double r = a/b ;
	
	std::vector<Point> ret ;
	ret.push_back( p+v*r) ;
	return ret ;
}

bool Line::intersects(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
			return squareDist(projection(g->getCenter()), g->getCenter()) < g->getRadius()*g->getRadius() ;
		}
	case ELLIPSE:
		{
			double a = g->getRadius() ;
			double b = dynamic_cast<const Ellipse*>(g)->getMinorRadius() ;
			double c = this->vector() * dynamic_cast<const Ellipse*>(g)->getMinorAxis() / (this->vector().norm() * dynamic_cast<const Ellipse*>(g)->getMinorAxis().norm()) ;
			Line L(g->getCenter(),dynamic_cast<const Ellipse*>(g)->getMinorAxis()) ;
			Point I = this->intersection(L) ;
			double d = (g->getCenter() - I).norm() ;

			if(abs(c) == 1)
			{
				double coordx = (this->origin() - g->getCenter()) * dynamic_cast<const Ellipse*>(g)->getMajorAxis() ;
				// case: Line is too far away				
				if(coordx > g->getRadius())
					return false ; 
				else
					return true ;
				}

			c = c / sqrt(1 - c*c) ;
			double A = (1 / (a*a) + (c*c) / (b*b)) ;
			double B = 2*c*d/(b*b) ;
			double C = (d*d) / (b*b) - 1 ;

			double delta = B * B - 4 * A * C ;
			return !(delta < 0) ;
		}
	case TRIANGLE:
		{
			bool ret = false ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case RECTANGLE:
		{

			bool ret = false ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case CONVEX_POLYGON:
		{

			bool ret = false ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case SEGMENTED_LINE:
		{

			bool ret = false ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case SPHERE:
		{
			Point center(g->getCenter()) ;
			double uc = v*g->getCenter() ;
			double delta = uc*uc - v.sqNorm()*g->getCenter().sqNorm() - g->getRadius()*g->getRadius() ;
			return delta >= 0 ;
		}
	default:
		{
			std::cout << "case not solved : " << g->getGeometryType() << std::endl ;
		return false ;
		}
	}
}


std::vector<Point> Line::intersection(const Geometry * g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
			double a = v.sqNorm() ;
			double b = p.x*v.x + p.y*v.y ;
			double c = p.sqNorm()-g->getRadius()*g->getRadius() ;
			double delta = b*b - 4*a*c ;
			
			if(delta == 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p+v*(-b/(2.*a))) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p+v*((-b + sqrt(delta))/(2.*a))) ;
				ret.push_back(p+v*((-b - sqrt(delta))/(2.*a))) ;
				return ret ;
			}

			return std::vector<Point>(0) ;
		}
	case ELLIPSE:
		{
			double a = dynamic_cast<const Ellipse*>(g)->getMajorRadius() ;
			double b = dynamic_cast<const Ellipse*>(g)->getMinorRadius() ;
			double c = this->vector() * dynamic_cast<const Ellipse*>(g)->getMinorAxis() / (this->vector().norm() * dynamic_cast<const Ellipse*>(g)->getMinorAxis().norm()) ;
			Line L(g->getCenter(),dynamic_cast<const Ellipse*>(g)->getMinorAxis()) ;
			Point I = this->intersection(L) ;
			I = g->getCenter() - I ;
			double d = I.norm() ;
			if((I ^ (dynamic_cast<const Ellipse*>(g)->getMajorAxis())).z < 0)
				d = - d ;

			double coordx = 0 ;
			double coordy = 0 ;

			// case: Line is parallel to minor axis
			if(abs(c) == 1)
			{
				std::vector<Point> ret ;
				coordx = (this->origin() - g->getCenter()) * dynamic_cast<const Ellipse*>(g)->getMajorAxis() ;
				// case: Line is too far away				
				if(coordx > g->getRadius())
					return ret ; 
				// case: Line is tangent
				if(coordx == g->getRadius())
					{ ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx) ; }
				// any other case: get two points on ellipse
				else
				{
					coordy = sqrt(b * b * ( 1 - (coordx * coordx) / (a * a))) ;
					ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx + dynamic_cast<const Ellipse*>(g)->getMinorAxis() * coordy) ;
					ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx - dynamic_cast<const Ellipse*>(g)->getMinorAxis() * coordy) ;
				}
				return ret ;
			}

//			std::cout << c <<  ";" << d << std::endl ;
			c = c / sqrt(1 - c*c) ;
			if((this->vector() * dynamic_cast<const Ellipse*>(g)->getMajorAxis()) < 0)
				c = -c ;
//			std::cout << c <<  ";" << d << std::endl ;

			double A = (1 / (a*a) + (c*c) / (b*b)) ;
			double B = 2*c*d/(b*b) ;
			double C = (d*d) / (b*b) - 1 ;

			double delta = B * B - 4 * A * C ;

			if(delta == 0)
			{
				std::vector<Point> ret ;
				coordx = - B / (2 * A) ;
				coordy = c * coordx + d ;
				ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx + dynamic_cast<const Ellipse*>(g)->getMinorAxis() * coordy) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				coordx = - (B + sqrt(delta)) / (2 * A) ;
				coordy = c * coordx + d ;
				ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx + dynamic_cast<const Ellipse*>(g)->getMinorAxis() * coordy) ;
				coordx = - (B - sqrt(delta)) / (2 * A) ;
				coordy = c * coordx + d ;
				ret.push_back(g->getCenter() + dynamic_cast<const Ellipse*>(g)->getMajorAxis() * coordx + dynamic_cast<const Ellipse*>(g)->getMinorAxis() * coordy) ;
				return ret ;
			}
		}
	case TRIANGLE:
		{

			std::vector<Point> ret ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				if(s.intersects(*this))
				{
					Point inter =  s.intersection(*this) ;
					ret.push_back(inter) ;
				}
			}
			
			return ret ;
		}
	case RECTANGLE:
		{
			
			std::vector<Point> ret ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				if(s.intersects(*this))
					ret.push_back(intersection(s)) ;
			}
			
			return ret ;
		}
	case CONVEX_POLYGON:
		{
			
			std::vector<Point> ret ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				if(s.intersects(*this))
					ret.push_back(intersection(s)) ;
			}
			
			return ret ;
		}
	case SEGMENTED_LINE:
		{
			
			std::vector<Point> ret ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
				if(s.intersects(*this))
					ret.push_back(intersection(s)) ;
			}
			
			return ret ;
		}
	case SPHERE:
		{
			Point center(g->getCenter()) ;
			double uc = v*g->getCenter() ;
			double delta = uc*uc - v.sqNorm()*g->getCenter().sqNorm() - g->getRadius()*g->getRadius() ;
			if(delta == 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p + v * uc) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				ret.push_back(p + v * (uc + sqrt(delta))) ;
				ret.push_back(p + v * (uc - sqrt(delta))) ;
				return ret ;
			}
		}
	default:
	{
		std::cout << "line-geometry type " << g->getGeometryType() << " not implemented" << std::endl ;
		return std::vector<Point>() ;
	}
	}
}

Point Line::projection(const Point &m ) const
{
// 	m.print() ;
	Point d = v/v.norm() ;
	Point w = m-p ;
	Point dir = d*(w*d) ;
	return p+dir ;
// 	double nu = (-v.x-v.y-v.z+v.x*m.x+v.y*m.y+v.z*m.z)/(v.x*v.x+v.y*v.y+v.z*v.z) ;
// 	Point can0 = p+v*nu ;
// 	Point can1 = p-v*nu ;
// // 	can0.print() ;
// // 	can1.print() ;
// 	if(squareDist3D(m, can0) < squareDist3D(m, can1))
// 		return can0 ;
// 	
// 	return can1 ;
}

Point Segment::normal() const
{
	if(this->norm() > 1e-8)
		return Point(this->f.y-this->s.y, this->s.x-this->f.x)/this->norm() ;
	
	return Point(0, 0) ;
}

Point Segment::normal(const Point & inside) const
{
	double nx = vec.y ;
	double ny = -vec.x ;
	double dx = mid.x - inside.x ;
	double dy = mid.y - inside.y ;
	Point n(nx, ny) ;
	
	if((nx*dx+ny*dy) > 0)
		return n/n.norm() ;
	else
		return n/(-n.norm()) ;
}

Segment::Segment(const Point & p0, const Point & p1)
{
	f = p0 ;
	s = p1 ;
	mid = (p0+p1)*0.5;
	vec = f-s ;
}

double Segment::norm() const
{
	return sqrt((f.y-s.y)*(f.y-s.y)+(f.x-s.x)*(f.x-s.x)) ;
}

std::vector<std::pair<Point, double> > Segment::getGaussPoints() const
{
	std::vector< std::pair<Point, double> > gp ;
	Point a = f*0.788675134594813+ s*(1.-0.788675134594813) ;
	Point b = s*0.788675134594813+ f*(1.-0.788675134594813) ;
	gp.push_back(std::pair<Point, double>(a, 1.0)) ;
	gp.push_back(std::pair<Point, double>(b, 1.0)) ;
	return gp ;
}

Segment::Segment()
{
	f = Point(0,0) ;
	s = Point(0,0) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}
Segment::~Segment()
{
}


void Segment::print() const 
{
	std::cout << "[ (" << f.x << ", " << f.y << ") ; (" << s.x << ", " << s.y << ") ]" << std::endl ;
}

bool Segment::intersects(const Line & l) const
{
	if (std::abs(-vec.x*l.vector().y+l.vector().x*vec.y) <= std::numeric_limits<double>::epsilon())
		return false ;
// 	if (-vec.x * l.vector().y + vec.y * l.vector().x == 0)
// 		return false ;
	
	Matrix m(2,2) ;
	Vector vv(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -l.vector().x ;
	m[1][0] = vec.y ; m[1][1] = -l.vector().y ;
	
	vv[0] = l.origin().x - s.x ; vv[1] = l.origin().y - s.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * vv ;
	return fac[0] < 1 && fac[0] > 0 ;
	
}

bool Segment::intersects(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
// 			return !intersection(g).empty() ;
// 			print() ;
			Line l(f, vec) ;
			Point proj = l.projection(g->getCenter()) ;
			if(g->in(f) && g->in(s))
				return false ;
			if(!g->in(f) && !g->in(s) && g->in(proj) && on(proj))
				return true ;
			if((g->in(f) && !g->in(s)) || (!g->in(f) && g->in(s)))
				return true ;
			
			return false ;
// 			
			double a = vec.sqNorm() ;
			double b = -(s.x-g->getCenter().x)*2.*vec.x - (s.y-g->getCenter().y)*2.*vec.y ;
			double c = (s.x-g->getCenter().x)*(s.x-g->getCenter().x) + (s.y-g->getCenter().y)*(s.y-g->getCenter().y)-g->getRadius()*g->getRadius() ;
// 			double delta = b*b - 4*a*c ;
			return b*b - 4.*a*c >= 0;
		}
	case ELLIPSE:
		{
//			vec.print() ;
			Line l(f,s-f) ;
			std::vector<Point> in = l.intersection(g) ;
//			std::cout << in.size() << std::endl ;
			if(in.size() == 0)
				return false ;
			for(size_t i = 0 ; i < in.size() ; i++)
			{
//				std::cout << isAligned(f,s,in[i]) <<  isAligned(f,in[i],s)  <<  isAligned(in[i],f,s)  << std::endl ;
//				in[i].print() ;
//				std::cout << this->on(in[i]) << std::endl ;
				if(this->on(in[i]))
					return true ;
			}
			return false ;
		}
	case TRIANGLE:
		{
			bool ret = false ;	
			
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;

				ret = ret || s.intersects(*this) /*|| isAligned(g->getBoundingPoint(i), &f, &this->s)*/;
			}
			
			return ret ;
		}
	case RECTANGLE:
		{	
			std::vector<Point> box = g->getBoundingBox() ;
			Segment s0(box[0], box[1]) ;
			Segment s1(box[1], box[2]) ;
			Segment s2(box[2], box[3]) ;
			Segment s3(box[3], box[0]) ;
			
			return intersects(s0) || intersects(s1) || intersects(s2) || intersects(s3) ;
		}
	case SEGMENTED_LINE :
		{
			bool ret = false ;	
			
			for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case CONVEX_POLYGON:
		{
			
			bool ret = false ;	
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				ret = ret || s.intersects(*this) ;
			}
			
			return ret ;
		}
	case SPHERE:
		{
			if(g->in(s) && g->in(f))
				return false ;
			if(g->in(s) || g->in(f))
				return true ;
			Point center(g->getCenter()) ;
			Point v = s-f ;
			double uc = v*g->getCenter() ;
			double delta = uc*uc - v.sqNorm()*g->getCenter().sqNorm() - g->getRadius()*g->getRadius() ;
			return delta >= 0 ;
		}
	default:
		return false ;
	}
}

bool Segment::intersects(const TriPoint *g) const
{
	if(isCoplanar(f, *g->point[0],*g->point[1],*g->point[2]))
		return g->in(f) ;
	if(isCoplanar(s, *g->point[0],*g->point[1],*g->point[2]))
		return g->in(s) ;
	Vector vec(3) ;
	vec[0] = f.x - g->point[0]->x ;
	vec[1] = f.y - g->point[0]->y ;
	vec[2] = f.z - g->point[0]->z ;

	Matrix mat(3,3) ;
	mat[0][0] = f.x-s.x; mat[0][1] = g->point[1]->x-g->point[0]->x; mat[0][2] = g->point[2]->x-g->point[0]->x; 
	mat[1][0] = f.y-s.y; mat[1][1] = g->point[1]->y-g->point[0]->y; mat[1][2] = g->point[2]->y-g->point[0]->y; 
	mat[2][0] = f.z-s.z; mat[2][1] = g->point[1]->z-g->point[0]->z; mat[2][2] = g->point[2]->z-g->point[0]->z; 
	Vector tuv = inverse3x3Matrix(mat)*vec ;

	return tuv.max() <=1 && tuv.min() >= 0 &&tuv[1] +tuv[2] <= 1  ;

// 	Plane p(*g->point[0], g->normal) ;
// 	Line l(f, vec) ;
// 	if(p.intersects(l))
// 	{
// 		Point i = p.intersection(l) ;
// 		Point v0(*g->point[0] - f) ;
// 		Point v1(*g->point[1] - f) ;
// 		Point n = v0^v1 ;
// 		n /= n.norm() ;
// 		if((i*n-f*n) >= 0)
// 			return true ;
// 	}
// 	return false ;
}

std::vector<Point> Segment::intersection(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case TRIANGLE:
		{
			std::vector<Point> ret ;
			
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				if(s.intersects(*this))
					ret.push_back( s.intersection(*this)) ;
			}
			std::sort(ret.begin(), ret.end()) ;
			std::vector<Point>::iterator e = std::unique(ret.begin(), ret.end()) ;
			ret.erase(e, ret.end()) ;
			return ret ;
		}
	case RECTANGLE:
		{
			std::vector<Point> ret ;
			
			std::vector<Point> bbox = g->getBoundingBox() ;
			Segment s0(bbox[0], bbox[1]) ; 
			
			if(s0.intersects(*this))
				ret.push_back( s0.intersection(*this)) ;
			
			Segment s1(bbox[1], bbox[2]) ; 
			
			if(s1.intersects(*this))
				ret.push_back( s1.intersection(*this)) ;
			
			Segment s2( bbox[2],  bbox[3]) ; 
			
			if(s2.intersects(*this))
				ret.push_back( s2.intersection(*this)) ;
			
			Segment s3(bbox[3],bbox[0]) ; 
			
			if(s3.intersects(*this))
				ret.push_back( s3.intersection(*this)) ;
			
			std::sort(ret.begin(), ret.end()) ;
			std::vector<Point>::iterator e = std::unique(ret.begin(), ret.end()) ;
			ret.erase(e, ret.end()) ;
			return ret ;
		}
	case CIRCLE:
		{
			if(g->in(f) && g->in(s))
				return std::vector<Point>(0) ;
			
			if((g->in(f) && !g->in(s)) || (!g->in(f) && g->in(s)))
			{
				Line l(f, vec) ;
				Point proj = l.projection(g->getCenter()) ;
				double d = sqrt(g->getRadius()*g->getRadius() - squareDist2D(proj, g->getCenter())) ;
				Point unitVector = vec/vec.norm() ;
				Point candidate = proj + unitVector*d ;
				if(on(candidate))
				{
					std::vector<Point> ret ;
					ret.push_back(candidate) ;
					return ret ;
				}
				else
				{
					std::vector<Point> ret ;
					ret.push_back(proj - unitVector*d) ;
					return ret ;
				}
				
			}
			else
			{
				std::vector<Point> ret ;
				Line l(f, vec) ;
				Point proj = l.projection(g->getCenter()) ;
				double dd = g->getRadius()*g->getRadius() - squareDist2D(proj, g->getCenter()) ;
				if(dd < 0)
					return ret ;
				
				double d = sqrt(dd) ;
				Point unitVector = vec/vec.norm() ;
				Point candidateA = proj + unitVector*d ;
				Point pa(candidateA) ; g->project(&pa) ;
				Point candidateB = proj - unitVector*d ;
				Point pb(candidateB) ; g->project(&pb) ;
// 				if(dist(candidateA, pa) < POINT_TOLERANCE)
// 					ret.push_back(candidateA) ;
// 				if(dist(candidateB, pb) < POINT_TOLERANCE)
// 					ret.push_back(candidateB) ;
				if(on(candidateA))
					ret.push_back(candidateA) ;
				if(on(candidateB))
					ret.push_back(candidateB) ;
				
				return ret ;
			}
		}
	case ELLIPSE:
		{
//			std::cout << "in Segment::intersection" << std::endl ;
//			f.print() ;
//			s.print() ;
//			vec.print() ;
//			(s-f).print() ;
//			std::cout << " " << std::endl ;
			Line l(f,s-f) ;
//			l.origin().print() ;
//			l.vector().print() ;
			std::vector<Point> ret = l.intersection(g) ;
//			if(ret.size() > 0)
//				{f.print() ; s.print() ;}
			for(size_t i = ret.size() ; i > 0 ; i--)
			{
//				ret[i-1].print() ;
				if(!(this->on(ret[i-1])))
					ret.erase(ret.begin()+i-1) ;
			}
			return ret ;

		}
	case SEGMENTED_LINE:
		{
			std::vector<Point> ret ;
			for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
			{
				Segment test(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
				if(test.intersects(*this))
					ret.push_back(test.intersection(*this)) ;
			}
			
			return ret ;
		}
	case SPHERE:
		{
			double a = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
			double b = (f.x-g->getCenter().x)*2.*vec.x + (f.y-g->getCenter().y)*2.*vec.y + (f.z-g->getCenter().z)*2.*vec.z;
			double c = (f.x-g->getCenter().x)*(f.x-g->getCenter().x) + (f.y-g->getCenter().y)*(f.y-g->getCenter().y)+ (f.z-g->getCenter().z)*(f.z-g->getCenter().z)-g->getRadius()*g->getRadius() ;
			double delta = b*b - 4.*a*c ;
			
			if(delta == 0)
			{
				std::vector<Point> ret ;
				Point A(f+vec*(-b/(2.*a))) ;
				if(on(A))
					ret.push_back(A) ;
				return ret ;
			}
			else if (delta > 0)
			{
				std::vector<Point> ret ;
				Point A(f+vec*(-b + sqrt(delta))/(2.*a)) ;
				if(on(A))
					ret.push_back(A) ;
				Point B(f+vec*(-b - sqrt(delta))/(2.*a)) ;
				if(on(B))
					ret.push_back(B) ;
				return ret ;
			}
			else
			{
				return std::vector<Point>(0) ;
			}
		}
		
	default:
		return std::vector<Point>(0) ;
	}
}


bool Segment::intersects(const Segment & l) const
{
	if (isAligned(l.first(), s, f) && isAligned(l.second(), s, f))
	{
		return l.on(f) || l.on(s) || on(l.first()) || on(l.second()) ;
	}
	
	if(l.on(f) || l.on(s) || on(l.first()) || on(l.second()))
		return true ;
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	if(std::abs(vec.x) < POINT_TOLERANCE && std::abs(l.vector().x) < POINT_TOLERANCE)
	{
		return false ;
	}	

	if(std::abs(vec.y) < POINT_TOLERANCE && std::abs(l.vector().y) < POINT_TOLERANCE)
	{
		return false ;
	}


	m[0][0] = -vec.x ; m[0][1] = l.vector().x ;
	m[1][0] = -vec.y ; m[1][1] = l.vector().y ;

	v[0] = -l.first().x + f.x ; v[1] = -l.first().y + f.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * v ;
	
	Point intersect = f + vec*fac[0];
// 	intersect.print() ;
	return on(intersect) && l.on(intersect) ;
	
}

Point Segment::project(const Point & p) const
{
	if(isAligned(f, s, p))
	{
		double d = squareDist3D(f, s) ;
		if( squareDist3D(f, p) < d && squareDist2D(s, p) < d )
			return p ;
		else if(squareDist3D(f, p) < squareDist2D(s, p))
			return f ;

		return s ;
	}
	Line l(s, Point(vec.x, vec.y)) ;
	
	Point candidate = l.projection(p) ;
	if(on(candidate))
	{
		return candidate ;
	}


	if(squareDist3D(p, f) > squareDist3D(p, s))
		return s ;

	return f ;
}

bool Segment::intersects(const Point & a, const Point & b) const
{
	return intersects(Segment(a,b));
}

bool Segment::on(const Point &p) const
{
// 	std::cout << "plyf" << std::flush ;	
// 		f.print() ;
// 		p.print() ;
// 		s.print() ;
	if(!isAligned(p, f, s))
		return false ;
// 	
// 	std::cout << "plaf" << std::flush ;	

	if(std::abs(vec.x) > 100.*POINT_TOLERANCE)
	{
// 		std::cout << "plof " << std::flush ;	

		if(s.x < f.x)
			return (s.x <= p.x) && (p.x <= f.x) ;
		return (f.x <= p.x) && (p.x <= s.x) ;
	}
	else if(std::abs(vec.y) > 100.*POINT_TOLERANCE)
	{
// 		std::cout << "pluf" << std::endl ;
		if(s.y < f.y)
			return (s.y <= p.y) && (p.y <= f.y) ;
		return (f.y <= p.y) && (p.y <= s.y) ;
	}

	return false ;
	
}

void Segment::setFirst(const Point & p) 
{
	f = p ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

void Segment::setFirst(double x, double y)
{
	f.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

void Segment::setSecond(const Point & p)
{
	s = p ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

void Segment::setSecond(double x, double y)
{
	s.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

void Segment::set(const Point & p0, const Point & p1)
{
	f.set(p0) ;
	s.set(p1) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

void Segment::set(double x0, double y0, double x1, double y1)
{
	f.set(x0, y0) ;
	s.set(x1, y1) ;
	mid = f*0.5 + s*0.5;
	vec = f-s ;
}

const Point & Segment::first() const
{
	return f ;
}

const Point & Segment::second() const
{
	return s ;
}

const Point & Segment::midPoint() const
{
	return mid ;
}

const Point & Segment::vector() const
{
	return vec ;
}


Point Segment::intersection(const Line & l) const
{
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -l.vector().x ;
	m[1][0] = vec.y ; m[1][1] = -l.vector().y ;
	
	v[0] = l.origin().x - f.x ; v[1] = l.origin().y - f.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * v ;
	
	return f + vec*fac[0];
}

Point TriPoint::projection(const Point & p) const
{
	Plane plane(*point[0], normal) ;
	
	Point planProj = plane.projection(p) ;
	
	if(in(planProj))
		return planProj ;
	
	Segment s0(*point[0], *point[1]) ;
	Segment s1(*point[0], *point[2]) ;
	Segment s2(*point[1], *point[2]) ;
	
	Point sproj0 = s0.project(planProj) ;
	double d0 = squareDist3D(p, sproj0) ;
	Point sproj1 = s1.project(planProj) ;
	double d1 = squareDist3D(p, sproj1) ;
	Point sproj2 = s2.project(planProj) ;
	double d2 = squareDist3D(p, sproj2) ;
	
	if(d0 <= d1 && d0 <= d2)
		return sproj0 ;
	if(d1 <= d0 && d1 <= d2)
		return sproj1 ;
	
	return sproj2 ;
	
}

bool TriPoint::in(const Point & p) const
{
	Point u = *point[1]-*point[0] ;
	Point v = *point[2]-*point[0] ;
	
	double uu = u*u ;
	double uv = u*v ;
	double vv = v*v ;
	Point w = p - *point[0] ;
	double wu = w*u ;
	double wv = w*v ;
	double d = uv*uv-uu*vv ;
	
	double q = (uv*wv - vv*wu) / d;
	
	if(q < 0. || q > 1.)
		return false ;
	
	double t = (uv * wu - uu * wv) / d;
	
	if (t < 0. || (q + t) > 1.) 
		return false;
	
	return true;
}

bool Segment::intersects(const TriPoint &g) const
{
	if(isCoplanar(f, *g.point[0],*g.point[1],*g.point[2]))
		return g.in(f) ;
	if(isCoplanar(s, *g.point[0],*g.point[1],*g.point[2]))
		return g.in(s) ;
	Vector vec(3) ;
	vec[0] = f.x - g.point[0]->x ;
	vec[1] = f.y - g.point[0]->y ;
	vec[2] = f.z - g.point[0]->z ;

	Matrix mat(3,3) ;
	mat[0][0] = f.x-s.x; mat[0][1] = g.point[1]->x-g.point[0]->x; mat[0][2] = g.point[2]->x-g.point[0]->x; 
	mat[1][0] = f.y-s.y; mat[1][1] = g.point[1]->y-g.point[0]->y; mat[1][2] = g.point[2]->y-g.point[0]->y; 
	mat[2][0] = f.z-s.z; mat[2][1] = g.point[1]->z-g.point[0]->z; mat[2][2] = g.point[2]->z-g.point[0]->z; 
	Vector tuv = inverse3x3Matrix(mat)*vec ;

	return tuv.max() <=1 && tuv.min() >= 0 &&tuv[1] +tuv[2] <= 1  ;
/*
	Point u = *s.point[1]-*s.point[0] ;
	Point v = *s.point[2]-*s.point[0] ;
	double a = -(s.normal*(f-*s.point[0])) ;
	double b = v*s.normal ;
	if(b < POINT_TOLERANCE)
	{
		if(abs(a) > POINT_TOLERANCE)
			return false ;
		
		Triangle t(*s.point[0], *s.point[1], *s.point[2]) ;
		return this->intersects(&t) ;
	}
	
	double r = a/b ;
	
	if(r < 0 || r > 1)
		return false ;
	
	Point planeIntersection = f+v*r ;
	
	double uu = u*u ;
	double uv = u*v ;
	double vv = v*v ;
	Point w = planeIntersection - *s.point[0] ;
	double wu = w*u ;
	double wv = w*v ;
	double d = uv*uv-uu*vv ;
	
	double q = (uv*wv - vv*wu) / d;
	
	if(q < 0. || q > 1.)
		return false ;
	
	double t = (uv * wu - uu * wv) / d;
	
	if (t < 0. || (q + t) > 1.) 
		return false;
	
	return true;*/
}

Point Segment::intersection(const Segment &l) const                                                                                                                 
{                                                                                                                                                                   
	if (isAligned(l.first(), s, f) && isAligned(l.second(), s, f))                                                                                               
	{                                                                                                                                                            
	if(on(l.first()) && on(l.second())) ;                                                                                                                
		return l.midPoint() ;                                                                                                                        
	if(on(l.first()))                                                                                                                                    
		return l.first() ;                                                                                                                           
	if(on(l.second()))                                                                                                                                   
		return l.second() ;                                                                                                                          
	if(l.on(f) && l.on(s)) ;                                                                                                                             
		return midPoint() ;                                                                                                                          
	if(l.on(f))                                                                                                                                          
		return f ;                                                                                                                                   
	if(l.on(s))                                                                                                                                          
		return s ;                                                                                                                                   
	}                                                                                                                                                            

	Matrix m(2,2) ;                                                                                                                                              
	Vector v(2) ;                                                                                                                                                
                                                                                                                                                              
	m[0][0] = vec.x ; m[0][1] = -l.vector().x ;                                                                                                                  
	m[1][0] = vec.y ; m[1][1] = -l.vector().y ;                                                                                                                  
                                                                                                                                                              
	v[0] = l.first().x - f.x ; v[1] = l.first().y - f.y ;                                                                                                        
                                                                                                                                                              
	invert2x2Matrix(m) ;                                                                                                                                         
                                                                                                                                                              
	Vector fac = m * v ;                                                                                                                                         
	return f + vec*fac[0];
}

std::vector<Point> Segment::intersection(const TriPoint &g) const
{
	std::vector<Point> ret ;
	if(isCoplanar(f, *g.point[0],*g.point[1],*g.point[2]))
		if( g.in(f))
		{
			ret.push_back(f) ;
		}

	if(isCoplanar(s, *g.point[0],*g.point[1],*g.point[2]))
		if( g.in(s))
		{
			ret.push_back(s) ;
		}

	Vector vec(3) ;
	vec[0] = f.x - g.point[0]->x ;
	vec[1] = f.y - g.point[0]->y ;
	vec[2] = f.z - g.point[0]->z ;

	Matrix mat(3,3) ;
	mat[0][0] = f.x-s.x; mat[0][1] = g.point[1]->x-g.point[0]->x; mat[0][2] = g.point[2]->x-g.point[0]->x; 
	mat[1][0] = f.y-s.y; mat[1][1] = g.point[1]->y-g.point[0]->y; mat[1][2] = g.point[2]->y-g.point[0]->y; 
	mat[2][0] = f.z-s.z; mat[2][1] = g.point[1]->z-g.point[0]->z; mat[2][2] = g.point[2]->z-g.point[0]->z; 
	Vector tuv = inverse3x3Matrix(mat)*vec ;
	ret.push_back(f+ (s-f)*tuv[0]) ;
	return ret ;
}

bool isInTriangle(const Point & test, const Point&  p0, const Point & p1, const Point  &p2) 
{
	return isOnTheSameSide( test, p0, p1, p2) && isOnTheSameSide(test, p1, p0, p2) && isOnTheSameSide(test, p2, p1, p2) ;
}

bool isOnTheSameSide(const Point & test, const Point & witness, const Point & f0, const Point & f1) 
{
	Point frontier(f1-f0) ;
	Point yes(witness-f0) ;
	Point perhaps(test-f0) ;
	return (frontier^yes).z*(frontier^perhaps).z > -POINT_TOLERANCE ;
}

bool isOnTheSameSide(const Point * test, const Point *witness, const Point *f0, const Point *f1) 
{
	Point frontier(*f1-*f0) ;
	Point yes(*witness-*f0) ;
	Point perhaps(*test-*f0) ;
	return (frontier^yes).z*(frontier^perhaps).z > -POINT_TOLERANCE ;
}

bool isOnTheSameSide(const Point & test, const Point & witness, const Point & f0, const Point & f1, const Point & f2) 
{
	Point f2test(f2-test) ;
	Point f1test(f1-test) ;
	Point f0test(f0-test) ;
	Point f2witness(f2-witness) ;
	Point f1witness(f1-witness) ;
	Point f0witness(f0-witness) ;
	return ((f2test.x*(f1test.y*f0test.z - f0test.y*f1test.z)-f2test.y*(f1test.x*f0test.z - f0test.x*f1test.z)+f2test.z*(f1test.x*f0test.y - f0test.x*f1test.y))*(f2witness.x*(f1witness.y*f0witness.z - f0witness.y*f1witness.z)-f2witness.y*(f1witness.x*f0witness.z - f0witness.x*f1witness.z)+f2witness.z*(f1witness.x*f0witness.y - f0witness.x*f1witness.y)) > 0) ;
}

bool isOnTheSameSide(const Point * test, const Point * witness, const Point * f0, const Point * f1, const Point * f2) 
{
	
	return (((f2->x-test->x)*((f1->y-test->y)*(f0->z-test->z) - (f0->y-test->y)*(f1->z-test->z))-(f2->y-test->y)*((f1->x-test->x)*(f0->z-test->z) - (f0->x-test->x)*(f1->z-test->z))+(f2->z-test->z)*((f1->x-test->x)*(f0->y-test->y) - (f0->x-test->x)*(f1->y-test->y)))*((f2->x-witness->x)*((f1->y-witness->y)*(f0->z-witness->z) - (f0->y-witness->y)*(f1->z-witness->z))-(f2->y-witness->y)*((f1->x-witness->x)*(f0->z-witness->z) - (f0->x-witness->x)*(f1->z-witness->z))+(f2->z-witness->z)*((f1->x-witness->x)*(f0->y-witness->y) - (f0->x-witness->x)*(f1->y-witness->y))) > 0) ;
}

double dist(const Point & v1, const Point & v2)
{
#ifdef HAVE_SSE4
		__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	r.vec += _mm_dp_pd(temp, temp, 62) ;
	return sqrt(r.val[0]+ r.val[1] );
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return sqrt(rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1]);
#else 
	double x = v1.x-v2.x ;
	double y = v1.y-v2.y ;
	double z = v1.z-v2.z ;
	double t = v1.t-v2.t ;
	return sqrt(x*x+y*y+z*z+t*t) ;
#endif
}

double dist(const Point * v1, const Point * v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	r.vec += _mm_dp_pd(temp, temp, 62) ;
	return sqrt(r.val[0]+ r.val[1] );
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return sqrt(rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1]);
#else 
	double x = v1->x-v2->x ;
	double y = v1->y-v2->y ;
	double z = v1->z-v2->z ;
	double t = v1->t-v2->t ;
	return sqrt(x*x+y*y+z*z+t*t) ;
#endif
}



double squareDist(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	r.vec += _mm_dp_pd(temp, temp, 62) ;
	return r.val[0]+ r.val[1] ;
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else 
	double x = v1.x-v2.x ;
	double y = v1.y-v2.y ;
	double z = v1.z-v2.z ;
	double t = v1.t-v2.t ;
	return x*x+y*y+z*z+t*t ;
#endif
}

double squareDist(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	r.vec += _mm_dp_pd(temp, temp, 62) ;
	return r.val[0]+ r.val[1];
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else 
	double x = v1->x-v2->x ;
	double y = v1->y-v2->y ;
	double z = v1->z-v2->z ;
	double t = v1->t-v2->t ;
	return x*x+y*y+z*z+t*t ;
#endif
}

double squareDist2D(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	return r.val[0] ;
#elif defined HAVE_SSE3
	vecdouble rzt ;
	rzt.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	return rzt.val[0]+ rzt.val[1];
#else 
	double x = v1.x-v2.x ;
	double y = v1.y-v2.y ;
	return x*x+y*y ;
#endif
}

double squareDist2D(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r ;
	temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	r.vec = _mm_dp_pd(temp, temp, 61) ;
	return r.val[0] ;
#elif defined HAVE_SSE3
	vecdouble rzt ;
	rzt.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	return rzt.val[0]+ rzt.val[1] ;
#else 
	double x = v1->x-v2->x ;
	double y = v1->y-v2->y ;
	return x*x+y*y ;
#endif
}

double squareDist3D(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r0 ;
	temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
	r0.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	r0.vec += _mm_dp_pd(temp, temp, 62) ;
	return r0.val[0]+ r0.val[1] ;
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else 
	double x = v1.x-v2.x ;
	double y = v1.y-v2.y ;
	double z = v1.z-v2.z ;
	return x*x+y*y+z*z ;
#endif
}



double squareDist3D(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
	__m128d temp ;
	vecdouble r0 ;
	temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
	r0.vec = _mm_dp_pd(temp, temp, 61) ;
	temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	r0.vec += _mm_dp_pd(temp, temp, 62) ;
	return r0.val[0]+ r0.val[1] ;
#elif defined HAVE_SSE3
	vecdouble rzt ;
	vecdouble rxy ;
	rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
	rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
	rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
	return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else 
	double x = v1->x-v2->x ;
	double y = v1->y-v2->y ;
	double z = v1->z-v2->z ;
	return x*x+y*y+z*z ;
#endif
}


ConvexPolygon* convexHull(const std::vector<Point *> * points)
{
	Point *pivot = (*points)[0]  ;
	for(size_t i = 1 ; i < points->size() ; i++)
	{
		if (pivot->y < (*points)[i]->y)
			pivot = (*points)[i] ;
	}
	std::cout << "pivot = " << pivot->x << ", " << pivot->y << std::endl ;
	
	//!then we build a map ordered by the angle to the pivot.
	
	std::map< double, Point* >  pointSet ;
	
	for(size_t i = 0 ; i < points->size() ; i++)
	{
		if ((*(*points)[i]) != (*pivot))
		{
			double angle = atan2(pivot->y-(*points)[i]->y, pivot->x-(*points)[i]->x) ;
			
			if(pointSet.find(angle) != pointSet.end())
			{
				std::cout << "already there..." << std::endl ;
				if(squareDist(pointSet[angle], pivot) < squareDist((*points)[i], pivot))
				{
					pointSet[angle] = (*points)[i] ;
				}
			}
			else
				pointSet[angle] = (*points)[i] ;
		}
	}
	
	pointSet[-10] = pivot ;
	
	//!we build the convex hull : we are on the convex hull if and only if we only turned left.
	
	std::vector<Point*> temphull ;
	std::vector<Point*> orderedset ;
	
	for(std::map<double, Point*>::const_iterator i = pointSet.begin() ; i!= pointSet.end() ; ++i)
	{
		std::cout << "point " <<  i->second->x << ", " <<  i->second->y << std::endl ;
		orderedset.push_back(i->second) ;
	}
	
	temphull.push_back(orderedset[0] ) ;
	temphull.push_back(orderedset[1] ) ;
	
	for(std::vector< Point *>::iterator i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
	{
		for(size_t j = 0 ; j < temphull.size() ; j++)
			std::cout << "(" << temphull[j]->x << ", " << temphull[j]->y << ")" << std::endl ;
		
		//! this is a usual cross product of two vectors...
		if(  ((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
		     ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y - temphull[temphull.size()-2]->y ) > 
		     1e-6 )
		{
			temphull.push_back(*i) ;
			std::cout << "new point in hull = " <<  (*i)->x << ", " <<  (*i)->y << std::endl ;
		}	
		else
		{
			while( !(((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			         ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y -temphull[temphull.size()-2]->y ) > 
			         1e-6))
			{
				std::cout << "out of the hull = " <<  (*temphull.rbegin())->x << ", " <<  (*temphull.rbegin())->y << std::endl ;
				temphull.pop_back();
			}
			temphull.push_back(*i) ;
			std::cout << "new point in hull = " <<  (*i)->x << ", " <<  (*i)->y << std::endl ;
		}
		
		
	}
	
	ConvexPolygon *hull = new ConvexPolygon(temphull.size()) ;
	
	std::copy(temphull.begin(), temphull.end(),hull->begin()) ;
	
	for(size_t i = 0 ; i < hull->size() ; i++)
		std::cout << "(" << (*hull)[i]->x << ", " << (*hull)[i]->y << ")" << std::endl ;
	
	return hull ;
} 

OrientableCircle::OrientableCircle(double r,double x, double y, double z, Point n )
{
	gType = ORIENTABLE_CIRCLE ;
	this->center = Point(x,y,z) ;
	this->normal = n/ sqrt(n.x*n.x+n.y*n.y+n.z*n.z) ;
	this->radius = r ;
}

OrientableCircle::OrientableCircle(double r,const Point * p0, Point normal ) 
{
	gType = ORIENTABLE_CIRCLE ;
	this->center = *p0 ;
	this->normal = normal/normal.norm() ;
	this->radius = r ;
}

OrientableCircle::OrientableCircle(double r,const Point p0, Point normal ) 
{
	gType = ORIENTABLE_CIRCLE ;
	this->center = p0 ;
	this->normal = normal/normal.norm() ;
	this->radius = r ;
}

OrientableCircle::OrientableCircle()
{
	gType = ORIENTABLE_CIRCLE ;
	this->center = Point() ;
	this->normal = Point(0,0,1) ;
	this->radius = 1 ;
}


std::vector<Point> OrientableCircle::getSamplingBoundingPoints(size_t num_points) const
{
	Vector start(3) ; start[0] = -normal.y*radius ;  start[1] = normal.x*radius ;  start[2] = 0 ; 
	
	if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]))
	{
		start[0] = 0 ;  start[1] = normal.z*radius ;  start[2] = -normal.y*radius ; 
	}
	
	if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
	{
		start[0] = normal.z*radius ;  start[1] = 0 ;  start[2] = -normal.x*radius ; 
	}
	
	
	double t = (1./(double)num_points)*M_PI ;
	
	double q0 = cos(t) ; 
	double st = sin(t) ;
	double q1 = st * normal.x ;
	double q2 = st * normal.y ;
	double q3 = st * normal.z ;
	
	Matrix R(3,3) ;
	
	R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3; R[0][1] = 2.*(q1*q2 - q0*q3) ; R[0][2] = 2.*(q1*q3 + q0*q2) ;
	R[1][0] = 2.*(q2*q1 + q0*q3) ; R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ; R[1][2] = 2.*(q2*q3 - q0*q1) ;
	R[2][0] = 2.*(q3*q1 - q0*q2) ; R[2][1] = 2.*(q3*q2 + q0*q1) ; R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;
	
	std::vector<Point> ret ;
	for(size_t i = 0 ; i < num_points ; i++)
	{
		start=R*start ;
		ret.push_back(Point( start[0] + center.x,
		                     start[1] + center.y, 
		                     start[2] + center.z)) ;
	}
	
	return ret ;
}

void OrientableCircle::sampleBoundingSurface(size_t num_points)
{
// 	std::cout << "sample OC" << std::endl ;
	
	Vector start(3) ; start[0] = -normal.y*radius ;  start[1] = normal.x*radius ;  start[2] = 0 ; 
	
	if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]))
	{
		start[0] = 0 ;  start[1] = normal.z*radius ;  start[2] = -normal.y*radius ; 
	}
	
	if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
	{
		start[0] = normal.z*radius ;  start[1] = 0 ;  start[2] = -normal.x*radius ; 
	}
	
	
	double t = (1./(double)num_points)*M_PI ;
		
	double q0 = cos(t) ; 
	double st = sin(t) ;
	double q1 = st * normal.x ;
	double q2 = st * normal.y ;
	double q3 = st * normal.z ;
		
	Matrix R(3,3) ;
		
	R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3; R[0][1] = 2.*(q1*q2 - q0*q3) ; R[0][2] = 2.*(q1*q3 + q0*q2) ;
	R[1][0] = 2.*(q2*q1 + q0*q3) ; R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ; R[1][2] = 2.*(q2*q3 - q0*q1) ;
	R[2][0] = 2.*(q3*q1 - q0*q2) ; R[2][1] = 2.*(q3*q2 + q0*q1) ; R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;
	
	this->boundingPoints.resize(num_points) ;
	for(size_t i = 0 ; i < num_points ; i++)
	{
		start=R*start ;
		boundingPoints[i] = new Point( start[0] + center.x,
		                           start[1] + center.y, 
		                           start[2] + center.z) ;
	}
	
}

void OrientableCircle::sampleSurface(size_t num_points)
{
	sampleBoundingSurface(num_points) ;
	
	std::vector<Point> toAdd ;
	
	double dr = 4.*M_PI*radius/num_points ;
	double rad = radius-dr ;
	
	while(rad > POINT_TOLERANCE)
	{
		
		size_t nPoints = (size_t)round((double)num_points*(rad/radius)) ;
		if(nPoints < 3)
			break ;
		
		Vector start(3) ; start[0] = -normal.y*rad ;  start[1] = normal.x*rad ;  start[2] = 0 ; 
		
		if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
		{
			start[0] = 0 ;  start[1] = normal.z*rad ;  start[2] = -normal.y*rad ; 
		}
		
		if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
		{
			start[0] = normal.z*rad ;  start[1] = 0 ;  start[2] = -normal.x*rad ; 
		}
		
		double t = (1./(double)nPoints)*M_PI ;
		double st = sin(t) ;
		double q0 = cos(t) ; 
		double q1 = st * normal.x ;
		double q2 = st * normal.y ;
		double q3 = st * normal.z ;
		
		Matrix R(3,3) ;
		
		R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3; R[0][1] = 2.*(q1*q2 - q0*q3) ; R[0][2] = 2.*(q1*q3 + q0*q2) ;
		R[1][0] = 2.*(q2*q1 + q0*q3) ; R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ; R[1][2] = 2.*(q2*q3 - q0*q1) ;
		R[2][0] = 2.*(q3*q1 - q0*q2) ; R[2][1] = 2.*(q3*q2 + q0*q1) ; R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;
		
		for(size_t i = 0 ; i < nPoints ; i++)
		{
			start=R*start ;
			toAdd.push_back( Point( start[0] + center.x,
			                        start[1] + center.y, 
			                        start[2] + center.z)) ;
		}
		rad -= dr ;
	}
	toAdd.push_back(center) ;
	this->inPoints.resize(toAdd.size()) ;
	for(size_t i = 0 ; i < toAdd.size() ;i++)
	{
		inPoints[i] = new Point(toAdd[i]) ;
	}
}

bool OrientableCircle::in(const Point & v) const
{
	return false ;
}

double OrientableCircle::area() const
{
	return radius*radius*M_PI ;
}

double OrientableCircle::volume() const
{
	return 0 ;
}

void OrientableCircle::project(Point * p) const
{
	return ;
}

void OrientableCircle::computeCenter()
{
	return ;
}

double OrientableCircle::getRadius() const
{
	return radius ;
}

bool isCoplanar(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2)  
{
	return isCoplanar(*test, *f0, *f1, *f2) ;
} ;

double signedAlignement(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1)
{
	Point a(f1) ; a -= test ;
	Point b(f0) ; b -= test ;
	return (a^b).z ;
// 	return (f1.x-test.x)*(f0.y-test.y) - (f0.x-test.x)*(f1.y-test.y);
}

bool isAligned(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1) 
{

	Point centre = (test+f0+f1)*.3333333333333333333333333 ;
	Point f0_(f0-centre) ;
	Point f1_(f1-centre) ;
	Point test_(test-centre) ;
	double scale = 100.*sqrt(std::max(std::max(f0_.sqNorm(), f1_.sqNorm()), test_.sqNorm())) ;
	f0_ /=scale ;
	f1_ /=scale ;
	test_ /=scale ;
	if(test == f1 || test == f0)
		return true ;
	if(std::abs((f0 - test) * (f1 - test)) < POINT_TOLERANCE)
		return false ;

	if (std::abs(signedAlignement(test_, f0_, f1_)) >= 2.*POINT_TOLERANCE)
		return false ;
	
	double mdist = sqrt(std::max(squareDist2D(f0_, f1_), std::max(squareDist2D(f0_, test_), squareDist2D(f1_, test_)))) ;

	double delta = .25*POINT_TOLERANCE*mdist ;

//	Point a(test_) ; a.x += delta ;
	Point b(test_) ; b.x += delta ; b.y += delta;
//	Point c(test_) ; c.y += delta;
	Point d(test_) ; d.x -= delta ; d.y += delta; 
//	Point e(test_) ; e.x -= delta ;
	Point f(test_) ; f.x -= delta ; f.y -= delta; 
//	Point g(test_) ; g.y -= delta;
	Point h(test_) ; h.x += delta ; h.y -= delta; 
	

	if(std::abs(signedAlignement(b, f0_, f1_)) >= 2.*POINT_TOLERANCE)
		return false ;

	if(std::abs(signedAlignement(d, f0_, f1_)) >= 2.*POINT_TOLERANCE)
		return false ;
	if(std::abs(signedAlignement(f, f0_, f1_)) >= 2.*POINT_TOLERANCE)
		return false ;
	if(std::abs(signedAlignement(h, f0_, f1_))  >= 2.*POINT_TOLERANCE)
		return false ;
	return true ;

} 

bool isAligned(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1)  
{
	return isAligned(*test, *f0, *f1) ;
} ;

bool isCoplanar(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  
{

	if(test == f1)
		return true ;
	if(test == f0)
		return true ;
	if(test == f2)
		return true ;
	
	Mu::Point A(f1) ; A -= f0 ;
	Mu::Point B(f2) ; B -= f1 ;
	Mu::Point C(f2) ; C -= test ;

	double c0 = signedCoplanarity(test, f0, f1, f2) ;
	double c02 = c0*c0 ;
	if(c02 > 8.*std::numeric_limits<double>::epsilon()*A.sqNorm())
		return false ;
	if(c02 > 8.*std::numeric_limits<double>::epsilon()*B.sqNorm())
		return false ;
	if(c02 > 8.*std::numeric_limits<double>::epsilon()*C.sqNorm())
		return false ;
	

	Point normal = A^B ;
	normal /= normal.norm()*4.*sqrt(std::numeric_limits<double>::epsilon()) ;
	
	Point a(test) ; a += normal ;
	Point b(test) ; b -= normal ;

	double c1 = signedCoplanarity(a, f0, f1, f2) ;
	double c2 = signedCoplanarity(b, f0, f1, f2) ;

// 	if(c0 > 0)
// 		return c1 < 0 || c2 < 0 ;
// 	if(c0 < 0)
// 		return c1 > 0 || c2 > 0 ;
// 	
// 	return true ;
// 	
	bool positive = c0 > 0 || c1 > 0 || c2 > 0 ;
	bool negative = c0 < 0 || c1 < 0 || c2 < 0 ;
	return  positive && negative ;
} ;

double coplanarity(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2)  
{
	return coplanarity(*test, *f0, *f1, *f2) ;
} ;

double coplanarity(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  
{

	Mu::Point A(f0) ; A -= f1 ;
	Mu::Point B(f2) ; B -= f1 ; 
	Mu::Point C(f2) ; C -= test ; 

	return  std::abs(triProduct(A, B, C))  ;
} ;

double signedCoplanarity(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  
{

	Mu::Point A(f0) ; A -= f1 ;
	Mu::Point B(f2) ; B -= f1 ; 
	Mu::Point C(f2) ; C -= test ; 

	return  triProduct(A, B, C)  ;
} ;

double signedCoplanarity(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2)  
{
	return signedCoplanarity(*test, *f0, *f1, *f2) ;
} ;

double triProduct(const Mu::Point &A, const Mu::Point &B, const Mu::Point &C)
{
	Point temp(A^B) ;
	return temp*C ;
	return (A.y*B.z - A.z*B.y)*C.x + (A.z*B.x - A.x*B.z)*C.y + (A.x*B.y - A.y*B.x)*C.z ;
}
