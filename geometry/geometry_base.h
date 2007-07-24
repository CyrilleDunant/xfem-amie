// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_BASE_H_
#define __GEOMETRY_BASE_H_
#include "../matrixops.h"

#include <map>

#define GEO_DERIVED_OBJECT( __geo_type__)    virtual const std::valarray<Point *> & getBoundingPoints() const \
{                                                          \
	return this->__geo_type__::getBoundingPoints() ;       \
}                                                          \
virtual std::valarray<Point *> & getBoundingPoints()       \
{                                                          \
	return this->__geo_type__::getBoundingPoints() ;       \
}                                                          \
virtual GeometryType getGeometryType() const               \
{                                                          \
return this->__geo_type__::gType ;                         \
}                                                          \
virtual const Point & getBoundingPoint(size_t i) const     \
{                                                          \
return this->__geo_type__::getBoundingPoint(i) ;           \
}                                                          \
virtual Point & getBoundingPoint(size_t i)                 \
{                                                          \
return this->__geo_type__::getBoundingPoint(i) ;           \
}                                                          \
virtual std::vector<Point> getBoundingBox() const          \
{                                                          \
return this->__geo_type__::getBoundingBox() ;              \
}                                                          \
virtual void setBoundingPoint(size_t i, Point * p)         \
{                                                          \
	this->__geo_type__::setBoundingPoint(i,p) ;            \
}                                                          \
virtual void setBoundingPoints(std::valarray<Point *> * nb)\
{                                                          \
	this->__geo_type__::setBoundingPoints(nb) ;            \
}                                                          \
virtual void setInPoints(std::valarray<Point *> * nb)      \
{                                                          \
	this->__geo_type__::setInPoints(nb) ;                  \
}                                                          \
virtual void project(Point * p) const                      \
{                                                          \
this->__geo_type__::project(p) ;                           \
}                                                          \
virtual double getRadius() const                           \
{                                                          \
return this->__geo_type__::getRadius() ;                   \
}                                                          \
virtual void sampleBoundingSurface(size_t n)               \
{                                                          \
this->__geo_type__::sampleBoundingSurface(n) ;             \
}                                                          \
virtual void sampleSurface(size_t n)                       \
{                                                          \
this->__geo_type__::sampleSurface(n) ;                     \
}                                                          \
virtual SpaceDimensionality spaceDimensions() const        \
{                                                          \
return this->__geo_type__::spaceDimensions() ;             \
}                                                          \
virtual bool intersects(const Geometry * g) const          \
{                                                          \
return this->__geo_type__::intersects(g) ;                 \
}                                                          \
virtual std::vector<Point> intersection(const Geometry * g) const          \
{                                                          \
return this->__geo_type__::intersection(g) ;               \
}                                                          \
virtual bool in(const Point & p) const                     \
{                                                          \
return this->__geo_type__::in(p) ;                         \
}                                                          \
const Point & getPoint(size_t i) const                     \
{                                                          \
return this->__geo_type__::getPoint(i) ;                   \
}                                                          \
Point & getPoint(size_t i)                                 \
{                                                          \
return this->__geo_type__::getPoint(i) ;                   \
}                                                          \
virtual const Point &  getInPoint(size_t i) const          \
{                                                          \
return this->__geo_type__::getInPoint(i) ;                 \
}                                                          \
virtual Point &  getInPoint(size_t i)                      \
{                                                          \
return this->__geo_type__::getInPoint(i) ;                 \
}                                                          \
virtual const std::valarray<Mu::Point*> & getInPoints() const    \
{                                                          \
return this->__geo_type__::getInPoints() ;                 \
}                                                          \
virtual  std::valarray<Mu::Point*> & getInPoints()     \
{                                                          \
return this->__geo_type__::getInPoints() ;                 \
}                                                          \
virtual size_t size() const                                \
{                                                          \
return this->__geo_type__::size() ;                        \
}                                                          \
virtual size_t sides() const                               \
{                                                          \
return this->__geo_type__::sides() ;                       \
}                                                          \
virtual Point & getCenter()                                 \
{                                                          \
return this->__geo_type__::getCenter() ;                   \
}                                                          \
virtual const Point & getCenter() const                     \
{                                                          \
		return this->__geo_type__::getCenter() ;           \
}                                                          \
virtual double area() const                                \
{                                                          \
return this->__geo_type__::area() ;                        \
}                                                          \
virtual double volume() const                              \
{                                                          \
return this->__geo_type__::volume() ;                      \
}


namespace Mu
{

static const double POINT_TOLERANCE =  1e-12 ; // 4.*std::numeric_limits<double>::epsilon() ;

typedef enum
{
	NULL_GEOMETRY,
	CIRCLE,
	TRIANGLE,
	RECTANGLE,
	PARALLELOGRAMME,
	CONVEX_POLYGON,
	SEGMENTED_LINE,
	ORIENTABLE_CIRCLE,
	CLOSED_NURB,
	TETRAHEDRON = 16,
	HEXAHEDRON = 17,
	SPHERE = 18,
} GeometryType ;

typedef enum
{
	SPACE_ONE_DIMENSIONAL,
	SPACE_TWO_DIMENSIONAL,
	SPACE_THREE_DIMENSIONAL
} SpaceDimensionality ;


class Segment ;
class ConvexPolygon ;

struct Point
{
	double x ;
	double y ;
	double z ;
	double t ;
	int id ;
	Point();
	Point(double x, double y) ;
	Point(double x, double y, double z) ;
	Point(double x, double y, double z, double t) ;
	
	void setX(double v) ;
	void setY(double v) ;
	void setZ(double v) ;
	void setT(double v) ;

	void set(double v, double vv) ;
	void set(double v, double vv, double vvv) ;
	void set(double v, double vv, double vvv, double vvvv) ;
	void set(const Point & p) ;
	void set(const Point * p) ;
	
	
	/** Robust point comparator.
	 * 
	 * @param p \c Point to compare to.
	 * @return is within 2*std::numeric_limits\<double\>::epsilon() of p.
	 */
	bool operator==(const Point &p) const ;
	
	/** Robust point comparator.
	 * 
	 * @param p \c Point to compare to.
	 * @return is outside of 2*std::numeric_limits\<double\>::epsilon() of p.
	 */
	bool operator!=(const Point &p) const ;
	
	/** Operator defining order in the plane.
	 * 
	 * @param p \c Point to compare to 
	 * @return this is bottommost then leftmost compared to p.
	 */
	bool operator <(const Point & p) const ;
	
	/** Operator defining order in the plane.
	 * 
	 * @param p \c Point to compare to 
	 * @return this is topmost then rightmost compared to p.
	 */
	bool operator >(const Point &p) const ;
	
	Point operator-(const Point &p) const ;
	Point operator-(const Vector &p) const ;
	Point operator+(const Point &p) const ;
	Point operator+(const Vector &p) const ;
	Point operator/(const double p) const ;
	Point operator*(const double p) const ;
	double operator*(const Point &p) const ;
	double operator*(const Vector &p) const ;
	Point operator^(const Point &p) const ;
	Point operator^(const Vector &p) const ;
	
	void operator+=(const Point &p) ;
	void operator-=(const Point &p) ;
	void operator*=(const double d) {
		x *= d ; 
		y *= d ; 
		z *= d ; 
	}
	
	void operator*=(const Matrix & m){
		Vector vec(3) ;
		vec[0] = x ; vec[1] = y ; vec[2] = z ; 
		vec = vec*m ;
		
		x = vec[0] ; y = vec[1] ; z = vec[2] ;
	}
	void operator/=(const double d) {
		x /= d ; 
		y /= d ; 
		z /= d ; 
	}
	
	/** Returns the norm of the Vector. The norm is simply
	 * \f$ \sqrt{x^2 + y^2} \f$.
	 * 
	 * @return the norm of the vector, 
	 */
	double norm() const;
	double sqNorm() const;
	void print() const  ;
	double angle() const ;
	double & operator[](size_t i) ;
	double operator[](size_t i) const ;
	
	
} ;


class Geometry 
{
protected:
	
	std::valarray<Point *> * inPoints ;
	bool sampled ;
	
	Point center ;
	
	virtual void computeCenter() = 0;
	
	GeometryType gType ;
	
public:
	
	Geometry() ;
	Geometry(size_t numPoints) ;
	virtual ~Geometry() 
	{
		if(this->inPoints != NULL)
		{
			for(size_t i = 0 ; i < this->inPoints->size() ; i++)
			{
				delete (*this->inPoints)[i] ;
				(*this->inPoints)[i] = NULL ;
			}
			delete this->inPoints ;
			inPoints = NULL ;
		}
		
// 		for(size_t i = 0 ; i < boundingPoints->size() ; i++)
// 		{
// 			delete (*boundingPoints)[i] ;
// 			(*boundingPoints)[i] = NULL ;
// 		}
// 		delete this->boundingPoints ;
	}
	
	virtual const std::valarray<Point *> & getBoundingPoints() const = 0;
	virtual std::valarray<Point *> & getBoundingPoints() = 0;
	virtual const Point & getBoundingPoint(size_t i) const = 0;
	virtual Point & getBoundingPoint(size_t i)  = 0;
	virtual const std::valarray<Point *> & getInPoints() const ;
	virtual std::valarray<Point *> & getInPoints() ;
	virtual const Point & getInPoint(size_t i) const ;
	virtual Point & getInPoint(size_t i);
	virtual void setBoundingPoint(size_t i, Point * p) = 0;
	virtual void setBoundingPoints(std::valarray<Point *> * nb) = 0;
	virtual void setInPoints(std::valarray<Point *> * nb) {delete inPoints ; inPoints=nb ; }
	virtual const Point & getPoint(size_t i) const = 0;
	virtual Point & getPoint(size_t i)  = 0;
	virtual const Point & getCenter() const ;
	virtual Point & getCenter() ;
	virtual void project(Point *) const = 0;
	
	virtual GeometryType getGeometryType() const ;
	
	virtual double getRadius() const = 0;
	
	virtual void sampleBoundingSurface(size_t num_points) = 0 ;
	virtual void sampleSurface(size_t num_points) = 0 ;
	virtual bool in(const Point & p)const  = 0;
	
	virtual size_t size() const = 0 ;
	virtual double area() const = 0;
	virtual double volume() const = 0;
	virtual size_t sides() const { return 3 ; }
	
	virtual bool intersects(const Geometry *) const ;
	virtual std::vector<Point> intersection(const Geometry *) const ;
	virtual std::vector<Point> getBoundingBox() const = 0;
	
	virtual SpaceDimensionality spaceDimensions() const = 0 ;
} ;

class Line
{
protected:
	Point p ;
	Point v ;
public:
	Line(const Point & origin, const Point & vector) ;
	Line(const Segment & base) ;
	Line(const Line & l, const Point & through) ;
	Line() ;
	
	bool intersects(const Line &l) const;
	bool intersects(const Segment &s) const;
	bool intersects(const Geometry *g) const;
	bool on(const Point &p) const;
	
	std::vector<Point> intersection(const Geometry * g) const ;
	Point intersection(const Line &l) const;
	Point intersection(const Segment &l) const;
	
	const Point & vector() const ;
	const Point & origin() const ;
	
	Point projection(const Point &p ) const ;
	
};

class Segment
{
protected:
	Point f ;
	Point s ;
	Point mid ;
	Point vec ;
	
public:
	Segment(const Point & p0, const  Point & p1) ;
	Segment() ;
	virtual ~Segment() ;
	
	bool intersects(const Line &l) const;
	bool intersects(const Segment &s) const;
	bool intersects(const Geometry *g) const;
	bool on(const Point &p) const;
	
	void setFirst(const Point & p) ;
	void setFirst(double x, double y) ;
	void setSecond(const Point & p) ;
	void setSecond(double x, double y) ;
	void set(const Point & p0,const  Point & p1) ;
	void set(double x0, double y0, double x1, double y1) ;
	Point normal() const ;
	Point normal(const Point & inside) const ;
	
	void print() const ;
	
	const Point & first() const;
	const Point & second() const;
	
	/** midpoint of the segment. is recalculated if the endpoints change.*/
	const Point & midPoint() const;
	const Point & vector() const ;
	double norm() const ;
	
	
	Point intersection(const Line &l) const;
	Point intersection(const Segment &l) const;
	bool intersects(const Point & a, const Point & b) const ;
	std::vector<Point> intersection(const Geometry * g) const;
	
	std::vector<std::pair<Point, double> > getGaussPoints() const ;
} ;


class PointSet
{
protected:
	std::valarray<Point *> * boundingPoints  ;
	size_t chullEndPos ;
	
public:
	PointSet() ;
	PointSet(size_t npoints) ;
	
	virtual ~PointSet()
	{
		if(this->boundingPoints != NULL)
		{
			for(size_t i = 0 ; i < this->boundingPoints->size() ; i++)
			{
				if((*this->boundingPoints)[i] != NULL)
				{
					delete (*this->boundingPoints)[i] ;
				}
			}
			delete this->boundingPoints ;
		}
	}
	double x(size_t i);
	double y(size_t i) ;
	double z(size_t i) ;
	
	void setX(size_t i, double v);
	void setY(size_t i, double v) ;
	void setZ(size_t i, double v) ;
	void set(size_t i, Point *p);
	void set(size_t i, double x, double y) ;
	void set(size_t i, double x, double y, double z) ;
	
	Point * operator[](size_t i) ;
	Point * operator[](size_t i) const;
	Point * getPoint(size_t i) const;
	Point * getPoint(size_t i) ;
	
	virtual bool in(const Point & p) const ;
	virtual size_t size() const;
	
	void removePoint(size_t index) ;

	Point computeCenter() const;
	
	ConvexPolygon * convexHull() const;
	
	typedef Point** iterator;
	
	iterator begin() const ;
	iterator end() const ;
} ;

class WeightedPoint : public Point
{ 
		
public:	/// Hold the weight of the point
	double w;

	WeightedPoint()  {this->w = 1.0;} 

	/// WeightedPoint constructor
	WeightedPoint(double x, double y, double z, double w) : Point(x,y,z), w(w) {} ;
	WeightedPoint(const Point &p, double w) : Point(p), w(w) {} ;
	WeightedPoint(const Point *p, double w) : Point(*p), w(w) {} ;

	/// Overload operators
	WeightedPoint operator+(const WeightedPoint &wP) const;
	WeightedPoint operator*(const double &w) const;
	WeightedPoint operator/(const double p)  const;

	/// Print WeightedPoint 
	void print() const ;
	
};

class Nurb {
protected:
	/// Hold the degree for testing purpose
	size_t degree;
	
	/// Hold the control points of the Nurb
	std::vector<Point> controlPoint;

	/// Hold the weights of the control points
	std::vector<double> weight;

	/// Hold the knot vector of the Nurb
	std::vector<double> knot;

public:
	/// Nurb constructor
	Nurb(std::vector<double> knot, size_t degree);

	/// Nurb constructor
	Nurb(std::vector<Point> controlPoint, std::vector<double> weight,std::vector<double> knot, size_t degree);

	/// Nurb constructor
	Nurb(std::vector<Point> controlPoint, std::vector<double> weight, std::vector<double> knot);

	/// Compute degree of the nurb based on the knowledge of the size of knot vector and number of control points
	int computeDegree(const std::vector<double> &vec,const std::vector<Point> &Pv);
	
	/// Transform knot vector into vector (0....1)
	std::vector<double> transVector(std::vector<double> &vec);

	/// Get the basis function - recursion
	double getBasis(double u, int i, int k, const std::vector<double> &vec); 

        /// Get point on the nurb
	/**
        * @param u parameter checked against the knot vector
	* @param r stopping parameter for the recursion
        */
	Point getNurbPoint(double u, double r = 0.001);
	
	/// Get nurb
	std::vector<Point> sampleNurb(size_t n);
	/**
        * @param n wished number of nurb points
        */
	
	/// Subdivide Nurb 
	Point pointOnNurb(double u);
	
	/// Check if a Nurb intersects a line
// 	bool intersects(const Line *l) const;

	/// Find the intersection - bisection method
// 	Point intersection(const Line *l) const;
	
};

class NonConvexGeometry : public PointSet, public Geometry
{
protected:
	std::valarray<Point *> orderedSet ;
	std::vector<size_t> stopPos ;
	
public:
	NonConvexGeometry() ;
	NonConvexGeometry(size_t numPoints) ;
	NonConvexGeometry(const std::valarray<Point *> & p) ;
	virtual ~NonConvexGeometry() { } ;
	
	virtual const std::valarray<Point * > & getBoundingPoints() const ;
	virtual std::valarray<Point * > & getBoundingPoints() ;
	virtual const Point & getBoundingPoint(size_t i) const ;
	virtual Point & getBoundingPoint(size_t i) ;
	virtual const Point & getPoint(size_t i) const ;
	virtual Point & getPoint(size_t i) ;
	
	virtual void setBoundingPoint(size_t i, Point * p) ;
	virtual void setBoundingPoints(std::valarray<Point *> * nb) ;
	virtual size_t size() const ;
	virtual double area() const = 0;
	
	virtual void project(Point *) const = 0;
	
} ;


class ConvexPolygon : public PointSet
{
public:
	ConvexPolygon(size_t npoints) ;
	ConvexPolygon(const PointSet * po) ;
	virtual ~ConvexPolygon() { } ;
	
	virtual bool in(const Point & p) const;
	virtual bool isTrigoOriented()  const ;
} ;



class ConvexGeometry :  public ConvexPolygon, public Geometry
{
public:
	ConvexGeometry() ;
	ConvexGeometry(size_t numPoints) ;
	virtual ~ConvexGeometry() { } ;
	
	virtual const std::valarray<Point *> & getBoundingPoints() const ;
	virtual std::valarray<Point *> & getBoundingPoints();
	virtual const Point & getBoundingPoint(size_t i) const ;
	virtual Point & getBoundingPoint(size_t i) ;
	
	virtual void setBoundingPoint(size_t i, Point * p) ;
	virtual void setBoundingPoints(std::valarray<Point *> * nb) ;
	
	virtual void sampleBoundingSurface(size_t num_points) = 0 ;
	virtual void sampleSurface(size_t num_points) = 0 ;
	virtual const Point & getPoint(size_t i) const ;
	virtual Point & getPoint(size_t i)  ;
	virtual size_t size() const ;
	virtual double area() const = 0;
	virtual void project(Point *) const = 0;
	
} ;



class OrientableCircle : public ConvexGeometry
{
protected:
	
	double radius ;
	Point normal ;
	
public:
	
	
	OrientableCircle(double radius,double x, double y, double z, Point normal ) ;
	
	OrientableCircle(double radius,const Point * p0, Point normal ) ;
	
	OrientableCircle(double radius,const Point p0, Point normal ) ; 
	
	OrientableCircle() ; 
	
	virtual ~OrientableCircle() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point & v) const ;
	
	virtual double area() const ;
	virtual double volume() const ;
	
	virtual void project(Point * p) const;
	
	virtual void computeCenter() ;
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_THREE_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;


} ;

/** Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
inline bool isAligned(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1)  
{
	return ( std::abs((f1.x-test.x)*(f0.y-test.y) - (f0.x-test.x)*(f1.y-test.y)) < Mu::POINT_TOLERANCE) ;
} ;


inline bool isAligned(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)
{
	return ( std::abs((f2.x-test.x)*((f1.y-test.y)*(f0.z-test.z) - (f0.y-test.y)*(f1.z-test.z))-(f2.y-test.y)*((f1.x-test.x)*(f0.z-test.z) - (f0.x-test.x)*(f1.z-test.z))+(f2.z-test.z)*((f1.x-test.x)*(f0.y-test.y) - (f0.x-test.x)*(f1.y-test.y))) < Mu::POINT_TOLERANCE) ;
} ;
/** Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
inline bool isAligned(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1)  
{
	return ( std::abs((f1->x-test->x)*(f0->y-test->y) - (f0->x-test->x)*(f1->y-test->y)) < Mu::POINT_TOLERANCE) ;
} ;

inline bool isCoplanar(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2)  
{
// 	((*f0-*f1)^(*f2-*f1)).print() ; (*test-*f1).print() ;
// 	std::cout << "coplanar ? " << std::abs(((*f0-*f1)^(*f2-*f1))*(*test-*f1)) << std::endl ;
// 	return std::abs(((*f0-*f1)^(*f2-*f1))*((*f0-*f1)^(*test-*f1))-((*f0-*f1)^(*f2-*f1)).norm()*((*f0-*f1)^(*test-*f1)).norm()) < 10*std::numeric_limits<double>::epsilon() ;
	Mu::Point A (*f0-*f1) ;
	Mu::Point B (*f2-*f1) ;
	Mu::Point C (*f2-*test) ;
	
	return  std::abs((A^B)*C) < Mu::POINT_TOLERANCE ;
} ;

inline bool isCoplanar(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  
{
// 	return std::abs(((f0-f1)^(f2-f1))*((f0-f1)^(test-f1))-((f0-f1)^(f2-f1)).norm()*((f0-f1)^(test-f1)).norm()) < 10*std::numeric_limits<double>::epsilon() ;
	Mu::Point A (f0-f1) ;
	Mu::Point B (f2-f1) ; 
	Mu::Point C (f2-test) ; 
	return  std::abs((A^B)*C) < Mu::POINT_TOLERANCE ;
} ;


inline bool isAligned(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1, const Mu::Point *f2)  
{
	return ( std::abs((f2->x-test->x)*((f1->y-test->y)*(f0->z-test->z) - (f0->y-test->y)*(f1->z-test->z))-(f2->y-test->y)*((f1->x-test->x)*(f0->z-test->z) - (f0->x-test->x)*(f1->z-test->z))+(f2->z-test->z)*((f1->x-test->x)*(f0->y-test->y) - (f0->x-test->x)*(f1->y-test->y))) < Mu::POINT_TOLERANCE) ;
} ;

/** Test if a point is in a triangle defined by three points.
 * 
 * @param test point to test.
 * @param p0 vertex 0.
 * @param p1 vertex 1.
 * @param p2 vertex 2.
 * @return true if test is in.
 */
bool isInTriangle(const Mu::Point &test, const Mu::Point &p0, const Mu::Point &p1, const Mu::Point &p2)  ;

/** Test wether two points lie on the same demi-plane.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point &test, const Mu::Point &witness, const Mu::Point &f0, const Mu::Point &f1)  ;
bool isOnTheSameSide(const Mu::Point *test, const Mu::Point *witness, const Mu::Point *f0, const Mu::Point *f1)  ;
bool isOnTheSameSide(const Mu::Point &test, const Mu::Point &witness, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  ;
bool isOnTheSameSide(const Mu::Point * test, const Mu::Point * witness, const Mu::Point * f0, const Mu::Point * f1, const Mu::Point * f2)  ;
//bool isAligned(const Point test, const Point f0, const Point f1)  ;

/**Return the distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Mu::Point &v1, const Mu::Point &v2) ;

/**Return the distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Mu::Point * v1, const Mu::Point * v2) ;

/** Return the convex hull of a set of points.
 * 
 * @param points 
 * @return a convex polygon (all boundary points anti-cockwise-ordered).
 */
Mu::ConvexPolygon* convexHull(const std::vector<Mu::Point> * points) ;


struct PointLessThan
{
	bool operator()(Mu::Point * p1, Mu::Point *p2)
	{
		return *p1 < *p2 ;
	}
} ;



/**Return the sqare distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point &v1, const Mu::Point &v2) ;

/**Return the sqare distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point *v1, const Mu::Point *v2) ;


#endif
