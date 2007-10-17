// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_base.h"

using namespace Mu ;

Point::Point() 
{
	this->x = 0 ; 
	this->y = 0 ;
	this->z = 0 ;
	this->t = 0 ;
	this->id = -1 ;
}

Point::Point(const Point & p) 
{
	this->x = p.x ; 
	this->y = p.y ;
	this->z = p.z ;
	this->t = p.t ;
	this->id = p.id ;
}

double Point::angle() const
{
	assert(z==0) ;
	return atan2(y,x) ;
}

Point::Point(double x, double y)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x ; 
	this->y = y ;
	this->z = 0 ;
	this->t = 0 ;
	this->id = -1 ;
}

Point::Point(double x, double y, double z)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x ; 
	this->y = y ;
	this->z = z ;
	this->t = 0 ;
	this->id = -1 ;
}

Point::Point(double x, double y, double z, double t)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x ; 
	this->y = y ;
	this->z = z ;
	this->t = t ;
	this->id = -1 ;
}

void Point::print() const
{
	std::cout << " ( id = " << id << " ; "<< round(x*1000.)/1000. << "; " << round(y*1000.)/1000. << "; "<< round(z*1000.)/1000.<< "; "<< round(t*1000.)/1000. << ") " << std::endl;
}

double Point::norm() const
{
	return sqrt(fma(x, x, fma(y, y, z*z))) ;
	return sqrt(x*x + y*y + z*z) ;
}

double Point::sqNorm() const
{
	return fma(x, x, fma(y, y, z*z)) ;
	return x*x + y*y + z*z;
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
	x = p.x ; 
	y = p.y;
	z = p.z;
	t = p.t;
}

void Point::set(const Point * p)
{
	x = p->x ; 
	y = p->y;
	z = p->z;
	t = p->t;
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

	return dist( &p, this) < POINT_TOLERANCE ;
/*
	return  std::abs(x-p.x) < 1e-8 && 
		std::abs(y-p.y) < 1e-8 &&  std::abs(z-p.z) < 1e-8;*/
}

bool Point::operator!=(const Point & p) const
{
// 	return squareDist( p, *this) >= 1e-8 ;
	
	return   dist(this, &p) >  POINT_TOLERANCE ;
}

Point Point::operator-(const Point &p) const
{
	Point ret((*this)) ;
	ret.x -= p.x; 
	ret.y -= p.y ; 
	ret.z -= p.z ; 
	ret.t -= p.t ; 
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
	return ret ; 
}

Point Point::operator+(const Point &p) const
{
	Point ret((*this)) ;
	ret.x += p.x; 
	ret.y += p.y ; 
	ret.z += p.z ; 
	ret.t += p.t ; 
	return ret ;
}

void Point::operator+=(const Point &p)
{
	x += p.x; 
	y += p.y ; 
	z += p.z ; 
	t += p.t ; 
}

void Point::operator-=(const Point &p)
{
	x -= p.x; 
	y -= p.y ; 
	z -= p.z ; 
	t -= p.t ; 
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
	return ret ; 
}

Point Point::operator/(const double p) const 
{
	Point ret((*this)) ;
	ret.x /= p ; 
	ret.y /= p ; 
	ret.z /= p ; 
	ret.t /= p ; 
	return ret ; 
}


bool Point::operator <(const Point &p) const 
{
	if(p == *this)
		return false ;
	
	double tol = POINT_TOLERANCE ;
	return (y < p.y ) 
		|| (( std::abs(y - p.y) < tol) 
		    && (x < p.x)) 
		|| (( std::abs(y - p.y) < tol) 
		    && ( std::abs(x - p.x) < tol) 
		    && (z < p.z)) 
		|| (( std::abs(y - p.y) < tol) 
		    && ( std::abs(x - p.x) < tol) 
		    && ( std::abs(z - p.z) < tol) 
		    && (t < p.t));
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
	ret.x *= p ; 
	ret.y *= p ; 
	ret.z *= p ; 
	ret.t *= p ; 
	return ret ; 
}


double Point::operator*(const Point &p) const
{
	return fma(x,p.x,fma(y,p.y, fma(z, p.z, t*p.t))) ;
	return x*p.x + y*p.y + z*p.z+ t*p.t;
}

double Point::operator*(const Vector &p) const
{
	double ret = x*p[0] + y*p[1] ;
	if(p.size() >2)
		ret+=z*p[2] ;
	if(p.size() >3)
		ret+=t*p[3] ;
	return ret ; 
}

Point Point::operator^(const Point &p) const
{
	Point ret ;
	ret.x = fma(y,p.z,  -z*p.y) ;
	ret.y = fma(z,p.x , -x*p.z) ;
	ret.z = fma(x,p.y , -y*p.x) ;
	
	return ret ;
}

Point Point::operator^(const Vector &p) const
{
	Point ret ;
	ret.x = y*p[2] - z*p[1] ;
	ret.y = z*p[0] - x*p[2] ;
	ret.z = x*p[1] - y*p[0] ;
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
			Segment s0(this->getBoundingPoint(0), this->getBoundingPoint(this->getBoundingPoints().size()/3)) ;
			Segment s1(this->getBoundingPoint(this->getBoundingPoints().size()/3), this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
			Segment s2(this->getBoundingPoint(0), this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
			
			if(g->getGeometryType() == TRIANGLE)
			{
				Segment t0(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/3)) ;
				Segment t1(g->getBoundingPoint(g->getBoundingPoints().size()/3), g->getBoundingPoint(2*g->getBoundingPoints().size()/3)) ;
				Segment t2(g->getBoundingPoint(0), g->getBoundingPoint(2*g->getBoundingPoints().size()/3)) ;
				
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
				return squareDist( this->getCenter(), g->getCenter()) < (this->getRadius() + g->getRadius())*(this->getRadius() + g->getRadius()) ;
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
	case SPHERE:
		{
			if(g->getGeometryType() == SPHERE)
			{
				if(dist(getCenter(), g->getCenter()) < 1e-8)
					return false ;
				if( dist(getCenter(), g->getCenter()) > getRadius()+g->getRadius()) 
					return false;
				if(dist(getCenter(), g->getCenter())+std::min(getRadius(),g->getRadius()) < std::max(getRadius(),g->getRadius()))
					return false ;
				
				return true ;
			}
			if(g->getGeometryType() == HEXAHEDRON)
			{
				std::vector<Point> bbox = g->getBoundingBox() ;
				
				if(g->in(getCenter()))
					return (getCenter().x+getRadius() > bbox[0].x) ||
				         (getCenter().x-getRadius() < bbox[7].x)    ||
				         (getCenter().y+getRadius() > bbox[0].y)    ||
				         (getCenter().y-getRadius() < bbox[7].y)    ||
				         (getCenter().z+getRadius() > bbox[0].z)    ||
					     (getCenter().z-getRadius() < bbox[7].z)     ;

				if(!g->in(getCenter()))
					return (getCenter().x+getRadius() < bbox[0].x) ||
					  (getCenter().x-getRadius() > bbox[7].x)    ||
					  (getCenter().y+getRadius() < bbox[0].y)    ||
					  (getCenter().y-getRadius() > bbox[7].y)    ||
					  (getCenter().z+getRadius() < bbox[0].z)    ||
					  (getCenter().z-getRadius() > bbox[7].z) ;
			}
			return false ;
			
		}
	case HEXAHEDRON:
		{
			if(g->getGeometryType() == SPHERE)
			{
				
				
// 				return true ;
				
				return (((g->getCenter().x+g->getRadius() > getBoundingBox()[0].x) ||
				         (g->getCenter().x-g->getRadius() < getBoundingBox()[7].x)    ||
				         (g->getCenter().y+g->getRadius() > getBoundingBox()[0].y)    ||
				         (g->getCenter().y-g->getRadius() < getBoundingBox()[7].y)    ||
				         (g->getCenter().z+g->getRadius() > getBoundingBox()[0].z)    ||
				         (g->getCenter().z-g->getRadius() < getBoundingBox()[7].z)    )&&in(g->getCenter())) ||
					(((g->getCenter().x+g->getRadius() < getBoundingBox()[0].x) ||
					  (g->getCenter().x-g->getRadius() > getBoundingBox()[7].x)    ||
					  (g->getCenter().y+g->getRadius() < getBoundingBox()[0].y)    ||
					  (g->getCenter().y-g->getRadius() > getBoundingBox()[7].y)    ||
					  (g->getCenter().z+g->getRadius() < getBoundingBox()[0].z)    ||
					  (g->getCenter().z-g->getRadius() > getBoundingBox()[7].z)    )&& !in(g->getCenter()))
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
			Segment s0(this->getBoundingPoint(0), this->getBoundingPoint(this->getBoundingPoints().size()/3)) ;
			Segment s1(this->getBoundingPoint(this->getBoundingPoints().size()/3), this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
			Segment s2(this->getBoundingPoint(0), this->getBoundingPoint(2*this->getBoundingPoints().size()/3)) ;
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
			std::vector<Point> intersection = s0.intersection(g) ;
			
			double perimetre = s0.norm()+s1.norm()+s2.norm()+s3.norm() ;
			
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
				
// 				ret.push_back(intersection[i]*3./5.+intersection[i+1]*2./5.) ;
// 				ret.push_back(intersection[i]*2./5.+intersection[i+1]*3./5.) ;
// 				ret.push_back(intersection[i]*1./5.+intersection[i+1]*4./5.) ;
			}
			if(!intersection.empty())
				ret.push_back(intersection[intersection.size()-1]) ;

			intersection = s1.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection[intersection.size()-1]) ;

			intersection = s2.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection[intersection.size()-1]) ;

			intersection = s3.intersection(g) ;
			for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
			{
				ret.push_back(intersection[i]) ;
				double l = dist(intersection[i], intersection[i+1]) ;
				size_t jmax = (size_t)round((double)g->getBoundingPoints().size()*(l/perimetre)) ;
				for(size_t j = 1 ; j < jmax ; j++ )
				{
					ret.push_back(intersection[i]*(double)(jmax-j)/(double)(jmax)+intersection[i+1]*(double)(j)/(double)(jmax)) ;
				}
			}
			if(!intersection.empty())
				ret.push_back(intersection[intersection.size()-1]) ;

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
				double a = 1 ;
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
					
					
					
					if(r*r-y_squared_0 >= 0)
					{
						double y = sqrt(y_squared_0) ;
						double x = sqrt(r*r-y_squared_0) ;
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
					
					return ret ;
				}

				double y_squared_1 = (-b - sqrt(delta))/(2.*a) ;

				if(y_squared_0 >= 0)
				{
					if(r*r-y_squared_0 >= 0)
					{
						double x = sqrt(r*r-y_squared_0) ;
						double y = sqrt(y_squared_0) ;
						
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
				if(y_squared_1 >= 0)
				{
					if(r*r-y_squared_1 >= 0)
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
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
						if(g->in(C.getInPoint(i)) && in(C.getBoundingPoint(i)))
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ; ;
					C.sampleSurface(num_points) ;
					for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
					{
						if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getBoundingPoint(i)) ;
	
					}
					for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
					{
						if(g->in(C.getInPoint(i)) && in(C.getBoundingPoint(i)))
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
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
					
					size_t num_points = (size_t)round(3.*sqrt(getBoundingPoints().size())*radiusOfIntersection/getRadius()) ;
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
						if(g->in(C.getInPoint(i)) && in(C.getBoundingPoint(i)))
							ret.push_back(C.getInPoint(i)) ;
					}
				}
			}
			return ret ;
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

const Point & Geometry::getCenter() const
{
	return this->center;
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


bool PointSet::in(const Point & p)  const 
{
	ConvexPolygon * hull = convexHull() ;
	bool ret = hull->in(p) ;
	delete hull ;
	return ret ;
}

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
}

Geometry::Geometry(size_t numPoints):inPoints(numPoints), gType(NULL_GEOMETRY)
{
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
			     std::numeric_limits<double>::epsilon() )
			{
				temphull.push_back(*i) ;
			}	
		else
		{
			while( !(((*i)->y - (*temphull.rbegin())->y)*((*temphull.rbegin())->x - temphull[temphull.size()-2]->x ) - 
			         ((*i)->x - (*temphull.rbegin())->x)*((*temphull.rbegin())->y -temphull[temphull.size()-2]->y ) > 
			         std::numeric_limits<double>::epsilon()))
			{
				temphull.pop_back();
			}
			temphull.push_back(*i) ;
		}
		
		
	}
	
	delete points ;
	std::copy(temphull.begin(), temphull.end(),this->begin()) ;
}


bool ConvexPolygon::in(const Point & p) const 
{
	
	bool in = false ;
	
	for (size_t i = 0, j  =  boundingPoints.size()-1; i <  boundingPoints.size(); j = i++)
	{
		if ((((boundingPoints[i]->y <= p.y ) && (p.y<boundingPoints[j]->y)) || ((boundingPoints[j]->y <= p.y) && (p.y<boundingPoints[i]->y))) &&
		    (p.x < (boundingPoints[j]->x - boundingPoints[i]->x) * (p.y - boundingPoints[i]->y) / (boundingPoints[j]->y - boundingPoints[i]->y) + boundingPoints[i]->x))
			in = !in;
	}
	
	return in ;
	
}


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

NonConvexGeometry::NonConvexGeometry(const std::valarray<Point *> & p)
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
	
	if (v.x * s.vector().y - v.y * s.vector().x == 0)
		return false ;
	
	Matrix m(2,2) ;
	Vector vv(2) ;
	
	m[0][0] = v.x ; m[0][1] = -s.vector().x ;
	m[1][0] = v.y ; m[1][1] = -s.vector().y ;
	
	vv[0] = s.first().x - p.x ; vv[1] = s.first().y - p.y ; 
	
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
// 	double t = 0;
// 	if(v.x != 0 && v.y != 0)
// 		t = ((p.y - s->first()->y) + v.y/v.x * (p.x - s->first()->x)) / (v.x/v.y - s->vector()->y) ;
// 	else if (v.x == 0)
// 		t = ( p.x - s->first()->x ) / s->vector()->x ;
// 	else if (v.y == 0)
// 		t = ( p.y - s->first()->y ) / s->vector()->y ;
// 	
// 	return (*s->first()) + (*s->vector())*t ;
	
	Matrix m(2,2) ;
	Vector vec(2) ;
	
	m[0][0] = v.x ; m[0][1] = -s.vector().x ;
	m[1][0] = v.y ; m[1][1] = -s.vector().y ;
	
	vec[0] = s.first().x - p.x ; vec[1] = s.first().y - p.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * vec ;
	
	return p + v*fac[0];
}


bool Line::intersects(const Geometry *g) const
{
	switch(g->getGeometryType())
	{
	case CIRCLE:
		{
			return squareDist(projection(g->getCenter()), g->getCenter()) < g->getRadius()*g->getRadius() ;
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
	case TRIANGLE:
		{

			std::vector<Point> ret ;
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				Point inter =  s.intersection(*this) ;
				ret.push_back(inter) ;
			}
			
			return ret ;
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
	
	Point normalVector(this->v.y, -this->v.x) ;
	
	Line perp(m, normalVector) ;
	return perp.intersection(*this) ;
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
		return n/(-1.*n.norm()) ;
}

Segment::Segment(const Point & p0, const Point & p1)
{
	f = p0 ;
	s = p1 ;
	mid = p0*0.5 + p1*0.5;
	vec = p1-p0 ;
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
	
	if (vec.x * l.vector().y - vec.y * l.vector().x == 0)
		return false ;
	
	Matrix m(2,2) ;
	Vector vv(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -l.vector().x ;
	m[1][0] = vec.y ; m[1][1] = -l.vector().y ;
	
	vv[0] = l.origin().x - f.x ; vv[1] = l.origin().y - f.y ; 
	
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
			return !intersection(g).empty() ;
// 			Line l(f, vec) ;
// 			
// 			return squareDist(l.projection(g->getCenter()), *g->getCenter()) < g->getRadius()*g->getRadius() ;
// 			
			double a = vec.sqNorm() ;
			double b = (f.x-g->getCenter().x)*2.*vec.x + (f.y-g->getCenter().y)*2.*vec.y ;
			double c = (f.x-g->getCenter().x)*(f.x-g->getCenter().x) + (f.y-g->getCenter().y)*(f.y-g->getCenter().y)-g->getRadius()*g->getRadius() ;
// 			double delta = b*b - 4*a*c ;
			return b*b - 4.*a*c >= 0;
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
			bool ret = false ;	
			
			for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
			{
				Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
				ret = ret || s.intersects(*this) /*|| isAligned(g->getBoundingPoint(i), &f, &this->s)*/;
			}
			
			return ret ;
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
	default:
		return false ;
	}
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
			std::vector<Point>::iterator e = std::unique(ret.begin(), ret.end()) ;
			ret.erase(e, ret.end()) ;
			return ret ;
		}
	case RECTANGLE:
		{
			std::vector<Point> ret ;
			
			Segment s0(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4)) ;
			
			if(s0.intersects(*this))
				ret.push_back( s0.intersection(*this)) ;
			
			Segment s1(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4)) ;
			
			if(s1.intersects(*this))
				ret.push_back( s1.intersection(*this)) ;
			
			Segment s2(g->getBoundingPoint(2*g->getBoundingPoints().size()/4), g->getBoundingPoint(3*g->getBoundingPoints().size()/4)) ;
			
			if(s2.intersects(*this))
				ret.push_back( s2.intersection(*this)) ;
			
			Segment s3(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(0)) ;
			
			if(s3.intersects(*this))
				ret.push_back( s3.intersection(*this)) ;
			
			return ret ;
		}
	case CIRCLE:
		{
			double a = vec.x*vec.x + vec.y*vec.y ;
			double b = (f.x-g->getCenter().x)*2.*vec.x + (f.y-g->getCenter().y)*2.*vec.y ;
			double c = (f.x-g->getCenter().x)*(f.x-g->getCenter().x) + (f.y-g->getCenter().y)*(f.y-g->getCenter().y)-g->getRadius()*g->getRadius() ;
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


bool Segment::intersects(const Segment & s) const
{
	if (std::abs(vec.x * s.vector().y - vec.y * s.vector().x) < POINT_TOLERANCE)
	{
		return isAligned(s.first(), s.second(), mid) ;
// 		return false ;
	}
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -s.vector().x ;
	m[1][0] = vec.y ; m[1][1] = -s.vector().y ;
	
	v[0] = s.first().x - f.x ; v[1] = s.first().y - f.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * v ;
	
	return fac[0] < 1.+POINT_TOLERANCE && fac[0] > -POINT_TOLERANCE && fac[1] < 1.+POINT_TOLERANCE && fac[1] > -POINT_TOLERANCE;
	
}

bool Segment::intersects(const Point & a, const Point & b) const
{
// 	if (std::abs(vec.x * (b->y - a->y) - vec.y * (a->x - b->x)) 
// 	{
// // 		return isAligned(a, b, &mid) ;
// 		return false ;
// 	}
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = (a.x - b.x) ;
	m[1][0] = vec.y ; m[1][1] = (a.y - b.y) ;
	
	if(std::abs(det(m)) < POINT_TOLERANCE)
		return false ;
	
	v[0] = a.x - f.x ; v[1] = a.y - f.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * v ;
	
	return fac[0] < 1.+POINT_TOLERANCE && fac[0] > -POINT_TOLERANCE && fac[1] < 1.+POINT_TOLERANCE && fac[1] > -POINT_TOLERANCE;
	
}

bool Segment::on(const Point &p) const
{
	return std::abs((vec.norm()- (f-p).norm() - (p-s).norm())) < POINT_TOLERANCE ;
}

void Segment::setFirst(const Point & p) 
{
	f = p ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setFirst(double x, double y)
{
	f.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setSecond(const Point & p)
{
	s = p ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::setSecond(double x, double y)
{
	s.set(x, y) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::set(const Point & p0, const Point & p1)
{
	f.set(p0) ;
	s.set(p1) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
}

void Segment::set(double x0, double y0, double x1, double y1)
{
	f.set(x0, y0) ;
	s.set(x1, y1) ;
	mid = f*0.5 + s*0.5;
	vec = s-f ;
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

Point Segment::intersection(const Segment &s) const
{
	
	Matrix m(2,2) ;
	Vector v(2) ;
	
	m[0][0] = vec.x ; m[0][1] = -s.vector().x ;
	m[1][0] = vec.y ; m[1][1] = -s.vector().y ;
	
	v[0] = s.first().x - f.x ; v[1] = s.first().y - f.y ; 
	
	invert2x2Matrix(m) ;
	
	Vector fac = m * v ;
	
	return f + vec*fac[0];
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
	return (frontier^yes).z*(frontier^perhaps).z > -1e-12 ;
}

bool isOnTheSameSide(const Point * test, const Point *witness, const Point *f0, const Point *f1) 
{
	Point frontier(*f1-*f0) ;
	Point yes(*witness-*f0) ;
	Point perhaps(*test-*f0) ;
	return (frontier^yes).z*(frontier^perhaps).z > -1e-12 ;
}

bool isOnTheSameSide(const Point & test, const Point & witness, const Point & f0, const Point & f1, const Point & f2) 
{
	
	return (((f2.x-test.x)*((f1.y-test.y)*(f0.z-test.z) - (f0.y-test.y)*(f1.z-test.z))-(f2.y-test.y)*((f1.x-test.x)*(f0.z-test.z) - (f0.x-test.x)*(f1.z-test.z))+(f2.z-test.z)*((f1.x-test.x)*(f0.y-test.y) - (f0.x-test.x)*(f1.y-test.y)))*((f2.x-witness.x)*((f1.y-witness.y)*(f0.z-witness.z) - (f0.y-witness.y)*(f1.z-witness.z))-(f2.y-witness.y)*((f1.x-witness.x)*(f0.z-witness.z) - (f0.x-witness.x)*(f1.z-witness.z))+(f2.z-witness.z)*((f1.x-witness.x)*(f0.y-witness.y) - (f0.x-witness.x)*(f1.y-witness.y))) > 0) ;
}

bool isOnTheSameSide(const Point * test, const Point * witness, const Point * f0, const Point * f1, const Point * f2) 
{
	
	return (((f2->x-test->x)*((f1->y-test->y)*(f0->z-test->z) - (f0->y-test->y)*(f1->z-test->z))-(f2->y-test->y)*((f1->x-test->x)*(f0->z-test->z) - (f0->x-test->x)*(f1->z-test->z))+(f2->z-test->z)*((f1->x-test->x)*(f0->y-test->y) - (f0->x-test->x)*(f1->y-test->y)))*((f2->x-witness->x)*((f1->y-witness->y)*(f0->z-witness->z) - (f0->y-witness->y)*(f1->z-witness->z))-(f2->y-witness->y)*((f1->x-witness->x)*(f0->z-witness->z) - (f0->x-witness->x)*(f1->z-witness->z))+(f2->z-witness->z)*((f1->x-witness->x)*(f0->y-witness->y) - (f0->x-witness->x)*(f1->y-witness->y))) > 0) ;
}

double dist(const Point & v1, const Point & v2)
{
	double x = v2.x-v1.x ;
	double y = v2.y-v1.y ;
	double z = v2.z-v1.z ;
	double t = v2.t-v1.t ;
	return sqrt(fma(x,x, fma(y,y, fma(z, z, t*t)))) ;
	return sqrt (x*x+y*y+z*z +t*t) ;
}

double dist(const Point * v1, const Point * v2)
{
	double x = v2->x-v1->x ;
	double y = v2->y-v1->y ;
	double z = v2->z-v1->z ;
	double t = v2->t-v1->t ;
	return sqrt(fma(x,x, fma(y,y, fma(z, z, t*t)))) ;
	return sqrt (x*x+y*y+z*z+t*t) ;
}


double squareDist(const  Point &v1, const Point & v2)
{
	double x = v2.x-v1.x ;
	double y = v2.y-v1.y ;
	double z = v2.z-v1.z ;
	double t = v2.t-v1.t ;
	return fma(x,x, fma(y,y, fma(z, z, t*t))) ;
	return (v2.x-v1.x)*(v2.x-v1.x)+(v2.y-v1.y)*(v2.y-v1.y) + (v2.z-v1.z)*(v2.z-v1.z)+ (v2.t-v1.t)*(v2.t-v1.t);
}

double squareDist(const Point *v1, const Point *v2)
{
	double x = v2->x-v1->x ;
	double y = v2->y-v1->y ;
	double z = v2->z-v1->z ;
	double t = v2->t-v1->t ;
	return fma(x,x, fma(y,y, fma(z, z, t*t))) ;
	return (v2->x-v1->x)*(v2->x-v1->x)+(v2->y-v1->y)*(v2->y-v1->y) + (v2->z-v1->z)*(v2->z-v1->z)+ (v2->t-v1->t)*(v2->t-v1->t);
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

