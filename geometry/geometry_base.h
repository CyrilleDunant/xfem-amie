// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009 (added: ellipses, level set)
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_BASE_H_
#define __GEOMETRY_BASE_H_
#include "../utilities/matrixops.h"
#include "../utilities/xml.h"

#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif

#include <map>

#define GEO_DERIVED_OBJECT( __geo_type__)    virtual const PointArray & getBoundingPoints() const \
{                                                          \
	return this->__geo_type__::getBoundingPoints() ;       \
}                                                          \
protected:\
virtual void computeCenter() {this->__geo_type__::computeCenter() ;}                                  \
public:\
virtual PointArray & getBoundingPoints()       \
{                                                          \
	return this->__geo_type__::getBoundingPoints() ;       \
}                                                          \
virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points)  const   \
{                                                           \
	return this->__geo_type__::getSamplingBoundingPoints(num_points) ; \
}                                                           \
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
virtual void setBoundingPoints(const PointArray & nb)\
{                                                          \
	this->__geo_type__::setBoundingPoints(nb) ;            \
}                                                          \
virtual void setInPoints(const PointArray & nb)      \
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
virtual const PointArray & getInPoints() const    \
{                                                          \
return this->__geo_type__::getInPoints() ;                 \
}                                                          \
virtual  PointArray & getInPoints()     \
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
}                                                          \
const __geo_type__ * getPrimitive() const                         \
{                                                          \
return dynamic_cast<const __geo_type__ *>(this) ;                 \
}                                                                 \
__geo_type__ * getPrimitive()                          \
{                                                          \
return dynamic_cast<__geo_type__ *>(this) ;                 \
}                                                                 \
virtual size_t & timePlanes()                         \
{                                                          \
 return this->__geo_type__::timePlanes() ;                 \
}                                                                 \
virtual double overlapFraction(const Geometry * src) const                         \
{                                                          \
 return this->__geo_type__::overlapFraction(src) ;                 \
}                                                                 \
virtual size_t timePlanes() const                         \
{                                                          \
return this->__geo_type__::timePlanes() ;                 \
}

#define LEVEL_SET_DERIVED_OBJECT()    virtual const PointArray & getBoundingPoints() const \
{                                                          \
	return this->LevelSet::getBoundingPoints() ;       \
}                                                          \
virtual PointArray & getBoundingPoints()       \
{                                                          \
	return this->LevelSet::getBoundingPoints() ;       \
}                                                          \
virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points)  const   \
{                                                           \
	return this->LevelSet::getSamplingBoundingPoints(num_points) ; \
}                                                           \
virtual GeometryType getGeometryType() const               \
{                                                          \
return this->LevelSet::gType ;                         \
}                                                          \
virtual const Point & getBoundingPoint(size_t i) const     \
{                                                          \
return this->LevelSet::getBoundingPoint(i) ;           \
}                                                          \
virtual Point & getBoundingPoint(size_t i)                 \
{                                                          \
return this->LevelSet::getBoundingPoint(i) ;           \
}                                                          \
virtual std::vector<Point> getBoundingBox() const          \
{                                                          \
return this->LevelSet::getBoundingBox() ;              \
}                                                          \
virtual void setBoundingPoint(size_t i, Point * p)         \
{                                                          \
	this->LevelSet::setBoundingPoint(i,p) ;            \
}                                                          \
virtual void setBoundingPoints(const PointArray & nb)\
{                                                          \
	this->LevelSet::setBoundingPoints(nb) ;            \
}                                                          \
virtual void setInPoints(const PointArray & nb)      \
{                                                          \
	this->LevelSet::setInPoints(nb) ;                  \
}                                                          \
virtual void project(Point * p) const                      \
{                                                          \
this->LevelSet::project(p) ;                           \
}                                                          \
virtual double getRadius() const                           \
{                                                          \
return this->LevelSet::getRadius() ;                   \
}                                                          \
virtual void sampleBoundingSurface(size_t n)               \
{                                                          \
this->LevelSet::sampleBoundingSurface(n) ;             \
}                                                          \
virtual void sampleSurface(size_t n)                       \
{                                                          \
this->LevelSet::sampleSurface(n) ;                     \
}                                                          \
virtual SpaceDimensionality spaceDimensions() const        \
{                                                          \
return this->LevelSet::spaceDimensions() ;             \
}                                                          \
virtual bool intersects(const Geometry * g) const          \
{                                                          \
return this->LevelSet::intersects(g) ;                 \
}                                                          \
virtual std::vector<Point> intersection(const Geometry * g) const          \
{                                                          \
return this->LevelSet::intersection(g) ;               \
}                                                          \
const Point & getPoint(size_t i) const                     \
{                                                          \
return this->LevelSet::getPoint(i) ;                   \
}                                                          \
Point & getPoint(size_t i)                                 \
{                                                          \
return this->LevelSet::getPoint(i) ;                   \
}                                                          \
virtual const Point &  getInPoint(size_t i) const          \
{                                                          \
return this->LevelSet::getInPoint(i) ;                 \
}                                                          \
virtual Point &  getInPoint(size_t i)                      \
{                                                          \
return this->LevelSet::getInPoint(i) ;                 \
}                                                          \
virtual const PointArray & getInPoints() const    \
{                                                          \
return this->LevelSet::getInPoints() ;                 \
}                                                          \
virtual  PointArray & getInPoints()     \
{                                                          \
return this->LevelSet::getInPoints() ;                 \
}                                                          \
virtual size_t size() const                                \
{                                                          \
return this->LevelSet::size() ;                        \
}                                                          \
virtual size_t sides() const                               \
{                                                          \
return this->LevelSet::sides() ;                       \
}                                                          \
virtual Point & getCenter()                                 \
{                                                          \
return this->LevelSet::getCenter() ;                   \
}                                                          \
virtual const Point & getCenter() const                     \
{                                                          \
		return this->LevelSet::getCenter() ;           \
}                                                          \
virtual double area() const                                \
{                                                          \
return this->LevelSet::area() ;                        \
}                                                          \
virtual double volume() const                              \
{                                                          \
return this->LevelSet::volume() ;                      \
}                                                          \
const LevelSet * getPrimitive() const                         \
{                                                          \
return dynamic_cast<const LevelSet *>(this) ;                 \
}                                                                 \
LevelSet * getPrimitive()                          \
{                                                          \
return dynamic_cast<LevelSet *>(this) ;                 \
}                                                                 \
virtual size_t & timePlanes()                         \
{                                                          \
 return this->LevelSet::timePlanes() ;                 \
}                                                                 \
virtual size_t timePlanes() const                         \
{                                                          \
return this->LevelSet::timePlanes() ;                 \
}



namespace Mu
{

static const double POINT_TOLERANCE =  1e-10 ;//std::numeric_limits<double>::epsilon() ;

typedef enum
{
	NULL_GEOMETRY,
	CIRCLE,
	LAYERED_CIRCLE,
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
	LAYERED_SPHERE = 19,
	REGULAR_OCTAHEDRON,
	ELLIPSE,
	LEVEL_SET
} GeometryType ;

/** \brief defines the possible vertices, axes and faces of a bounding box */
typedef enum
{
  TOP,
  LEFT,
  BOTTOM,
  RIGHT,
  FRONT,
  BACK,
  TOP_LEFT,
  TOP_RIGHT,
  BOTTOM_LEFT,
  BOTTOM_RIGHT,
  FRONT_LEFT,
  FRONT_RIGHT,
  BACK_LEFT,
  BACK_RIGHT,
  FRONT_TOP,
  FRONT_BOTTOM,
  TOP_LEFT_FRONT,
  TOP_LEFT_BACK,
  BOTTOM_LEFT_FRONT,
  BOTTOM_LEFT_BACK,
  TOP_RIGHT_FRONT,
  TOP_RIGHT_BACK,
  BOTTOM_RIGHT_FRONT,
  BOTTOM_RIGHT_BACK,
	TOP_BEFORE,
  LEFT_BEFORE,
  BOTTOM_BEFORE,
  RIGHT_BEFORE,
  FRONT_BEFORE,
  BACK_BEFORE,
  TOP_LEFT_BEFORE,
  TOP_RIGHT_BEFORE,
  BOTTOM_LEFT_BEFORE,
  BOTTOM_RIGHT_BEFORE,
  FRONT_LEFT_BEFORE,
  FRONT_RIGHT_BEFORE,
  BACK_LEFT_BEFORE,
  BACK_RIGHT_BEFORE,
  FRONT_TOP_BEFORE,
  FRONT_BOTTOM_BEFORE,
  TOP_LEFT_FRONT_BEFORE,
  TOP_LEFT_BACK_BEFORE,
  BOTTOM_LEFT_FRONT_BEFORE,
  BOTTOM_LEFT_BACK_BEFORE,
  TOP_RIGHT_FRONT_BEFORE,
  TOP_RIGHT_BACK_BEFORE,
  BOTTOM_RIGHT_FRONT_BEFORE,
  BOTTOM_RIGHT_BACK_BEFORE,
  TOP_AFTER,
  LEFT_AFTER,
  BOTTOM_AFTER,
  RIGHT_AFTER,
  FRONT_AFTER,
  BACK_AFTER,
  TOP_LEFT_AFTER,
  TOP_RIGHT_AFTER,
  BOTTOM_LEFT_AFTER,
  BOTTOM_RIGHT_AFTER,
  FRONT_LEFT_AFTER,
  FRONT_RIGHT_AFTER,
  BACK_LEFT_AFTER,
  BACK_RIGHT_AFTER,
  FRONT_TOP_AFTER,
  FRONT_BOTTOM_AFTER,
  TOP_LEFT_FRONT_AFTER,
  TOP_LEFT_BACK_AFTER,
  BOTTOM_LEFT_FRONT_AFTER,
  BOTTOM_LEFT_BACK_AFTER,
  TOP_RIGHT_FRONT_AFTER,
  TOP_RIGHT_BACK_AFTER,
  BOTTOM_RIGHT_FRONT_AFTER,
  BOTTOM_RIGHT_BACK_AFTER,
} BoundingBoxPosition ;

typedef enum
{
	SPACE_ONE_DIMENSIONAL,
	SPACE_TWO_DIMENSIONAL,
	SPACE_THREE_DIMENSIONAL
} SpaceDimensionality ;


class Segment ;
class ConvexPolygon ;

/** \brief Four dimensionnal point. Storage of the coordinates makes use of the processors's vector extensions*/
struct Point
{
	#ifdef HAVE_SSE3
	union
	{
		__m128d vecxy ;
		struct
		{
	#endif
			double x ;
			double y ;
	#ifdef HAVE_SSE3
		} ;
	};
	union
	{
		__m128d veczt ;
		struct
		{
	#endif
			double z ;
			double t ;
	#ifdef HAVE_SSE3
		} ;
	};
	#endif
	/** \brief ID of the point, useful for mesh indexing*/
	int id ;
	/** \brief default constructor, coordiantes are nil, and id is -1 */
	Point();
	/** \brief initialise x and y values
	*
	* @param x
	* @param y
	*/
	Point(double x, double y) ;
	
	/** \brief initialise x, y and z values
	*
	* @param x
	* @param y
	* @param z
	*/
	Point(double x, double y, double z) ;
	/** \brief initialise x, y, z and t values
	*
	* @param x
	* @param y
	* @param z
	* @param t
	*/
	Point(double x, double y, double z, double t) ;
	
	/** \brief copy-constructor.*/
	Point(const Point & p) ;
	
	Point& operator = (const Point & p) ;

	/** \brief constructor from a XML tree item*/
	Point(XMLTree * xml) ;

	/** \brief export to XML.
	* Data structure is <point><vector> 4 x y z t </vector><id> id </id></point>
	*/
	XMLTree * toXML() ;

	
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
	
	/** \brief dot-product*/
	double operator*(const Point &p) const ;
	
	/** \brief dot-product*/
	double operator*(const Vector &p) const ;
	
	/** \brief cross-product*/
	Point operator^(const Point &p) const ;
	
	/** \brief cross-product*/
	Point operator^(const Vector &p) const ;
	
	void operator+=(const Point &p) ;
	void operator-=(const Point &p) ;
	void operator*=(const double d) {
		x *= d ; 
		y *= d ; 
		z *= d ; 
		t *= d ; 
	}
	
	void operator*=(const Matrix & m){
		if(m.numCols() == 3)
		{
			Vector vec(3) ;
			vec[0] = x ; vec[1] = y ; vec[2] = z ; 
			vec = vec*m ;
			
			x = vec[0] ; y = vec[1] ; z = vec[2] ;
		}
		else if(m.numCols() == 4)
		{
			Vector vec(4) ;
			vec[0] = x ; vec[1] = y ; vec[2] = z ; vec[3] = t ; 
			vec = vec*m ;
			
			x = vec[0] ; y = vec[1] ; z = vec[2] ;  t = vec[3] ;
		}
	}

	Point operator* (const Matrix & m) const {
		if(m.numCols() == 3)
		{
			Vector vec(3) ;
			vec[0] = x ; vec[1] = y ; vec[2] = z ; 
			vec = vec*m ;
			return Point(vec[0], vec[1], vec[2]) ;
		}
		else if(m.numCols() == 4)
		{
			Vector vec(4) ;
			vec[0] = x ; vec[1] = y ; vec[2] = z ; vec[3] = t ; 
			vec = vec*m ;
			
			return Point(vec[0], vec[1], vec[2], vec[3]) ;
		}
		
		return Point() ;
	}
	void operator/=(const double d) {
		double inv = 1./d ;
		x *= inv ; 
		y *= inv ; 
		z *= inv ; 
		t *= inv ; 
	}
	
	/** Returns the norm of the Vector. The norm is simply
	 * \f$ \sqrt{x^2 + y^2 + z^2 + t^2} \f$.
	 * 
	 * @return the norm of the vector, 
	 */
	double norm() const;
	
	/** Returns the sqare norm of the Vector. The norm is simply
	 * \f$ x^2 + y^2 + z^2 + t^2 \f$.
	 * 
	 * @return the square norm of the vector, 
	 */
	double sqNorm() const;
	
	/** \brief print point coords and id*/
	void print() const  ;
	
	/** \brief return atan2(y, x)*/
	double angle() const ;
	
	/** \brief access x, y, z or t*/
	double & operator[](size_t i) ;
	
	/** \brief access x, y, z or t*/
	double operator[](size_t i) const ;
	
	
} ;

/** \brief light class for 3D triangle-line intersection computation*/
struct TriPoint
{
	Point normal ;

	std::valarray<Point const *> point ;
	std::vector<TriPoint *> neighbour ;
	
	/** \brief constructor from three points */
	TriPoint(const Point * p0, const Point * p1, const Point * p2) : point(3) 
	{
		point[0] = p0 ;
		point[1] = p1 ;
		point[2] = p2 ;
		normal = (*p0-*p1)^(*p2-*p1) ;
	}
	
	/** \brief return area of the triangle*/
	double area() const
	{
		return .5*normal.norm() ;
	}

	/** \brief return true is the argument is in*/
	bool in(const Point & p) const ;
	
	Point projection(const Point & p) const ;

};


typedef  std::valarray<Point *> PointArray;

/** \brief Basic interface for a geometrical object.
* 
* A geometrical object is defined by a set of points, in and out discrimination and intersection and projection operators
*/
class Geometry 
{
protected:
	
	PointArray inPoints ;
	bool sampled ;
	
	Point center ;
	
	virtual void computeCenter() = 0;
	
	GeometryType gType ;
	size_t time_planes ;
	
public:
	
	Geometry() ;
	Geometry(size_t numPoints) ;
	virtual ~Geometry() ;

	virtual XMLTree * toXML() ;
	
	/** \brief Acessor get the bounding points of the geometry*/
	virtual const PointArray & getBoundingPoints() const = 0;
	
	/** \brief Acessor get the bounding points of the geometry*/
	virtual PointArray & getBoundingPoints() = 0;
	
	/** \brief Acessor get the ith bounding point of the geometry*/
	virtual const Point & getBoundingPoint(size_t i) const = 0;
	
	/** \brief Acessor get the ith bounding point of the geometry*/
	virtual Point & getBoundingPoint(size_t i)  = 0;
	
	/** \brief Acessor get the inside points of the geometry*/
	virtual const PointArray & getInPoints() const ;
	
	/** \brief Acessor get the inside points of the geometry*/
	virtual PointArray & getInPoints() ;
	
	/** \brief Acessor get the ith inside point of the geometry*/
	virtual const Point & getInPoint(size_t i) const ;
	
	/** \brief Acessor get the ith inside point of the geometry*/
	virtual Point & getInPoint(size_t i);
	
	/** \brief Acessor set the ith bounding point of the geometry*/
	virtual void setBoundingPoint(size_t i, Point * p) = 0;
	
	/** \brief Acessor set all bounding poinst of the geometry*/
	virtual void setBoundingPoints(const PointArray & nb) = 0;
	
	/** \brief Acessor set all inside poinst of the geometry*/
	virtual void setInPoints(const PointArray & nb) {inPoints.resize(nb.size()) ; inPoints=nb ; }
	
	/** \brief Acessor get the ith point of the geometry*/
	virtual const Point & getPoint(size_t i) const = 0;
	
	/** \brief Acessor get the ith point of the geometry*/
	virtual Point & getPoint(size_t i)  = 0;
	
	/** \brief Acessor get the center of the geometry*/
	virtual const Point & getCenter() const ;
	
	/** \brief Acessor get the center of the geometry*/
	virtual Point & getCenter() ;
	
	/** \brief Project the argument on the geometry*/
	virtual void project(Point *) const = 0;
	
	/** \brief Set a new center for the geometry*/
	virtual void setCenter(const Point & newCenter) ;
	
	/** \brief Return the geometry type*/
	virtual GeometryType getGeometryType() const ;
	
	/** \brief Return the circumscribing radius*/
	virtual double getRadius() const = 0;
	
	/** \brief Sample the bounding surface with a given number of sampling points. the points are stored as boundingPoints*/
	virtual void sampleBoundingSurface(size_t num_points) = 0 ;

	/** \brief Return points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const = 0 ;
	
	/** \brief Sample the bounding surface and the surface of the Geometry. The points are stored as inpoints and bounding points */
	virtual void sampleSurface(size_t num_points) = 0 ;
	
	/** \brief Return true  if the argument is in the geometry*/
	virtual bool in(const Point & p)const  = 0;
	
	/** \brief Return the total number of points stored as inpoints and boundingpoints*/
	virtual size_t size() const = 0 ;
	
	/** \brief Area of the Geometry*/
	virtual double area() const = 0;
	
	/** \brief Volume of the Geometry*/
	virtual double volume() const = 0;
	
	/** \brief Return the number of sides of the geometry. This only makes sense if the geometry is a polyhedron, or polygon*/
	virtual size_t sides() const { return 3 ; }
	
	/** \brief return true if this geometry intersects with the argument*/
	virtual bool intersects(const Geometry *) const ;
	
	/** \brief return the intersection points between this geometry and the argument*/
	virtual std::vector<Point> intersection(const Geometry *) const ;
	
	/** \brief Get the bounding box. 
   * The points are topLeft, topRight, bottomRight, bottomLeft or in 3D:
	 * center.x+0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z) ;
	 * center.x+0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z) ;
	 * center.x+0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z) ;
	 * center.x+0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z) ;
	 * center.x-0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z) ;
	 * center.x-0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z) ;
	 * center.x-0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z) ;
	 * center.x-0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z) ;
	 */
	virtual std::vector<Point> getBoundingBox() const = 0;
	
	/** \brief Return the number of space dimensions of the geometry (SPACE_TWO_DIMENSIONAL or SPACE_THREE_DIMENSIONAL)*/
	virtual SpaceDimensionality spaceDimensions() const = 0 ;

	/** \brief Return the number of time slices present in this geometry*/
	virtual size_t timePlanes() const ;
	virtual size_t & timePlanes() ;
	
	/** \brief return the fraction of the area of the current geometry also lying in the target geometry. */
	virtual double overlapFraction(const Geometry * target) const ;

} ;

/** \brief Line used for intersection computation. It also provides a projection method*/
class Line
{
protected:
	Point p ;
	Point v ;
public:
	/** \brief Constructor. The line is constructed from a point through which it passes and a direction*/
	Line(const Point & origin, const Point & vector) ;
	
	/** \brief Constructor. Build the line as the prolongation of a Segment.*/
	Line(const Segment & base) ;
	
	/** \brief Constructor. build the line as the parallel to a given line and a point to pass through*/
	Line(const Line & l, const Point & through) ;
	
	/** \brief default constructor, Lines passing through 0,0 and with direction (1, 0)*/
	Line() ;
	
	/** \brief return true if this lines intersects with the argument*/
	bool intersects(const Line &l) const;
	
	/** \brief return true if this line intersects with the argument*/
	bool intersects(const Segment &s) const;
	
	/** \brief return true if this line intersects with the argument*/
	bool intersects(const Geometry *g) const;

	/** \brief return true if this line intersects with the argument*/
	bool intersects(const TriPoint &g) const;
	
	/** \brief return true if the argument lies on this Line*/
	bool on(const Point &p) const;
	
	/** \brief Return the set of intersection Point s between this line and the argument*/
	std::vector<Point> intersection(const Geometry * g) const ;
	
	/** \brief Return the intersection point between this Line and the argument. 
	 *A Point will always be returned, so the user should first check whether the lines intersect*/
	Point intersection(const Line &l) const;
	
	/** \brief Return the intersection point between this Line and the argument. 
	 *A Point will always be returned, so the user should first check whether the line and segment intersect*/
	Point intersection(const Segment &l) const;

	/** \brief Return the intersection point between this Line and the argument. 
	 *A Point will always be returned, so the user should first check whether the line and segment intersect*/
	std::vector<Point> intersection(const TriPoint &l) const;
	
	/** \brief return the direction vector*/
	const Point & vector() const ;
	
	/** \brief return the origin point*/
	const Point & origin() const ;
	
	/** \brief Return the projection of the argument on this line*/
	Point projection(const Point &p ) const ;
	
};

/** \brief Helper class used to compute intersection between geometries in three dimensions*/
class Plane
{
protected:
	Point p ;
	Point v ;
public:
	
	/** \brief Contructor. construct a plane from a point and a normal vector.*/
	Plane(const Point & origin, const Point & vector) ;
	
	/** \brief Contructor. construct a plane from three points*/
	Plane(const Point & a, const Point & b, const Point &c) ;
	
	/** \brief Return true if this plane intersects the argument*/
	bool intersects(const Line &l) const;
	
	/** \brief Return true if this plane intersects the argument*/
	bool intersects(const Segment &s) const;
	
	/** \brief Return true if this plane intersects the argument*/
	bool intersects(const Geometry *g) const;
	
	/** \brief return true if this plane intersects thr argument. Two parallel planes will intersect.*/
	bool intersects(const Plane &) const;
	
	/** \brief Return true if the argument is on the plane*/
	bool on(const Point &p) const;
	
	/** \brief return the set of intersection points with the argument Geometry.*/
	std::vector<Point> intersection(const Geometry * g) const ;
	
	/** \brief return the intersection point with the argument.
	*
	* This function will always return a point so the user should first check for the existence of the intersection, using interscts()
	*/
	Point intersection(const Line &l) const;
	
	/** \brief return the intersection point with the argument.
	*
	* This function will always return a point so the user should first check for the existence of the intersection, using interscts()
	*/
	Point intersection(const Segment &l) const;
	
	/** \brief return the intersection Line with the argument.
	*
	* This function will always return a point so the user should first check for the existence of the intersection, using interscts()
	*/
	Line intersection(const Plane &) const ;
	
	/** \brief return the normal vector*/
	const Point & vector() const ;
	
	/** \brief return the origin*/
	const Point & origin() const ;
	
	/** \brief return the projected point from the argument on the plane*/
	Point projection(const Point &p ) const ;
	
};

/** \brief Helper class used to compute geometry-geometry intersections*/
class Segment
{
protected:
	Point f ;
	Point s ;
	Point mid ;
	Point vec ;
	
public:
	double * getFirstX() {return &f.x ;}
	double * getFirstY() {return &f.y ;}
	double * getFirstZ() {return &f.z ;}
	double * getSecondX() {return &s.x ;}
	double * getSecondY() {return &s.y ;}
	double * getSecondZ() {return &s.z ;}
public:
	
	/** \brief Constructor. construct a segment from two endpoints*/
	Segment(const Point & p0, const  Point & p1) ;
	
	Segment(const Segment & s) ;
	Segment & operator =(const Segment &s) ;
	
	/** \brief Default Constructor, return a segment (0,0), (1,0)*/
	Segment() ;
	virtual ~Segment() ;
	
	/** \brief return true if this Segment intersects the argument*/
	bool intersects(const Line &l) const;
	
	/** \brief return true if this Segment intersects the argument*/
	bool intersects(const Segment &s) const;
	
	/** \brief return true if this Segment intersects the argument*/
	bool intersects(const Geometry *g) const;
	
	/** \brief return true if this Segment intersects the argument*/
	bool intersects(const TriPoint *g) const;
	
	/** \brief return true if this Segment intersects the argument*/
	bool intersects(const Plane &g) const;
	
	/** \brief return true if the argument is on the segment*/
	bool on(const Point &p) const;
	
	/** \brief Accessor, set the first point*/
	void setFirst(const Point & p) ;
	
	/** \brief Accessor, set the first point*/
	void setFirst(double x, double y) ;
	
	/** \brief Accessor, set the second point*/
	void setSecond(const Point & p) ;
	
	/** \brief Accessor, set the second point*/
	void setSecond(double x, double y) ;
	
	/** \brief Accessor, set the endpoints */
	void set(const Point & p0,const  Point & p1) ;
	
	/** \brief Accessor, set the endpoints */
	void set(double x0, double y0, double x1, double y1) ;
	
	/** \brief Return a normal to the segment.*/
	Point normal() const ;
	
	/** \brief Return a normal to the segment. The argument gives the "inside" of the Segment*/
	Point normal(const Point & inside) const ;
	
	void print() const ;
	
	/** \brief return the first Point*/
	const Point & first() const;
	
	/** \brief return the second Point*/
	const Point & second() const;
	
	/** \brief midpoint of the segment. is recalculated if the endpoints change.*/
	const Point & midPoint() const;
	/** \brief Vector obtained by computing first - second.*/
	const Point & vector() const ;
	
	/** \brief return the norm of the Segment*/
	double norm() const ;
	
	/** \brief return the projected point to the segment. This point will be one of the endpoint, 
	* if the original point is not in the infinite band defined by the segment
	*/
	Point project(const Point & p) const ;
	
	/** \brief return the intersection between the segment and the line. This always return a point, 
	 * so the user should check whether line and segment intersect
	 */
	Point intersection(const Line &l) const;
	
		/** \brief return the intersection between the segment and the line. This always return a point, 
	 * so the user should check whether segments intersect
	 */
	Point intersection(const Segment &l) const;
	
		/** \brief segment-secment intersection, with the segment segment defined by the arguments
	 */
	bool intersects(const Point & a, const Point & b) const ;

	/** \brief segment-triangle intersection
	 */
	bool intersects(const TriPoint & a) const ;
	
	/** \brief Return the intersection between this segment and the Argument. */
	std::vector<Point> intersection(const Geometry * g) const;
	
	/** \brief segment-triangle intersection
	 */
	std::vector<Point> intersection(const TriPoint & a) const ;

	/** \brief Return two Gauss points*/
	std::vector<std::pair<Point, double> > getGaussPoints() const ;
} ;

/** \brief basic container for Geometry: the Pointset. it contains an algorithm to compute the convex hull.*/
class PointSet
{
protected:
	PointArray boundingPoints  ;
	size_t chullEndPos ;
	
public:
	
	/** \brief default constructor. creates an empty pointset*/
	PointSet() ;
	
	/** \brief Construct a pointset with a given size*/
	PointSet(size_t npoints) ;
	
	virtual ~PointSet()
	{

		for(size_t i = 0 ; i < this->boundingPoints.size() ; i++)
		{
			delete boundingPoints[i] ;
		}

	}
	
	/** \brief accessor return the x coordinate of the ith point of the pointset*/
	double x(size_t i);
	
	/** \brief accessor return the y coordinate of the ith point of the pointset*/
	double y(size_t i) ;
	
	/** \brief accessor return the z coordinate of the ith point of the pointset*/
	double z(size_t i) ;
	
	/** \brief accessor set the x coordinate of the ith point of the pointset*/
	void setX(size_t i, double v);
	
	/** \brief accessor set the y coordinate of the ith point of the pointset*/
	void setY(size_t i, double v) ;
	
	/** \brief accessor set the z coordinate of the ith point of the pointset*/
	void setZ(size_t i, double v) ;
	
	/** \brief Accessor, set the ith point of the set*/
	void set(size_t i, Point *p);
	
	/** \brief Accessor, set the x, y coordinates of the ith point*/
	void set(size_t i, double x, double y) ;
	
	/** \brief Acessor, set the x, y, z coordinates of the ith point*/
	void set(size_t i, double x, double y, double z) ;
	
	Point * operator[](size_t i) ;
	Point * operator[](size_t i) const;
	Point * getPoint(size_t i) const;
	Point * getPoint(size_t i) ;
	
//	virtual bool in(const Point & p) const ;
	virtual size_t size() const;
	
	void removePoint(size_t index) ;

	Point computeCenter() const;
	
	ConvexPolygon * convexHull() const;
	
	typedef Point** iterator;
	typedef Point* const* const_iterator;
	
	const_iterator begin() const ;
	const_iterator end() const ;
	iterator begin() ;
	iterator end() ;
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

	virtual XMLTree * toXML() {return new XMLTree("nurb") ; } ;
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
        * @param u wished number of nurb points
        */
	
	/// Subdivide Nurb 
	Point pointOnNurb(double u);
	
	/// Check if a Nurb intersects a line
// 	bool intersects(const Line *l) const;

	/// Find the intersection - bisection method
// 	Point intersection(const Line *l) const;
	
};

/** \brief Class defining the interface for a non-convex Geometry */
class NonConvexGeometry : public PointSet, public Geometry
{
protected:
	PointArray orderedSet ;
	std::vector<size_t> stopPos ;
	
public:
	
	/** \brief Default constructor, create empty Geometry*/
	NonConvexGeometry() ;
	
	/** \brief Constructor, create Geometry with @param numPoints uninitialised points*/
	NonConvexGeometry(size_t numPoints) ;
	
	/** \brief Constructor, create a Geometry from an array of Points*/
	NonConvexGeometry(const PointArray & p) ;
	virtual ~NonConvexGeometry() { } ;
	
	/** \brief Accessor, return the bounding Points. */
	virtual const PointArray & getBoundingPoints() const ;
	
	/** \brief Accessor, return the bounding Points*/
	virtual PointArray & getBoundingPoints() ;

	/** \brief Accessor, return the ith Bounding Point*/
	virtual const Point & getBoundingPoint(size_t i) const ;
	
	/** \brief Accessor, return the ith Bounding Point*/
	virtual Point & getBoundingPoint(size_t i) ;
	
	/** \brief Accessor, return the ith Point*/
	virtual const Point & getPoint(size_t i) const ;
	
	/** \brief Accessor, return the ith Point*/
	virtual Point & getPoint(size_t i) ;
	
	/** \brief Accessor, set the ith Bounding Point*/
	virtual void setBoundingPoint(size_t i, Point * p) ;
	
	/** \brief Accessor, set the Bounding Points*/
	virtual void setBoundingPoints(const PointArray & nb) ;
	
	/** \brief Return the total number of Points*/
	virtual size_t size() const ;
	
	/** \brief return the Area */
	virtual double area() const = 0;
	
	/** \brief projection operator*/
	virtual void project(Point *) const = 0;
	
} ;

/** \brief Convex Polygon, defined from a Pointset*/
class ConvexPolygon : public PointSet
{
public:
	/** \brief Constructor. Create npoints uninitialised points*/
	ConvexPolygon(size_t npoints) ;
	
	/** \brief Constructor. Construct from a pointset*/
	ConvexPolygon(const PointSet * po) ;
	virtual ~ConvexPolygon() { } ;

	virtual XMLTree * toXML() {return new XMLTree("convex polygon") ; } ;
	
	/** \brief return true if the argument is in the Polygon*/
//	virtual bool in(const Point & p) const;
	
	/** \brief return true is the points are trigonometrically oriented*/
	virtual bool isTrigoOriented()  const ;
} ;

/** \brief General convex Geometry.*/
class ConvexGeometry :  public ConvexPolygon, public Geometry
{
public:
	
	/** \brief default constructor. Create empty geometry */
	ConvexGeometry() ;
	
	/** \brief Constructor. Create npoints uninitialised points*/
	ConvexGeometry(size_t numPoints) ;
	virtual ~ConvexGeometry() { } ;
	
	/** \brief Accessor, return the bounding Points. */
	virtual const PointArray & getBoundingPoints() const ;
	
	/** \brief Accessor, return the bounding Points. */
	virtual PointArray & getBoundingPoints();
	
	/** \brief Accessor, return the ith Bounding Point*/
	virtual const Point & getBoundingPoint(size_t i) const ;
	
	/** \brief Accessor, return the ith Bounding Point*/
	virtual Point & getBoundingPoint(size_t i) ;
	
	/** \brief Accessor, set the ith Bounding Point*/
	virtual void setBoundingPoint(size_t i, Point * p) ;
	
	/** \brief Accessor, set the Bounding Points*/
	virtual void setBoundingPoints(const PointArray & nb) ;
	
	/** \brief accessor, return the ith point*/
	virtual const Point & getPoint(size_t i) const ;
	
	/** \brief accessor, return the ith point*/
	virtual Point & getPoint(size_t i)  ;
	
	/** \brief Return the total number of points*/
	virtual size_t size() const ;
	virtual double area() const = 0;
	virtual void project(Point *) const = 0;
//	virtual bool intersects(const Geometry * g) const ;
	
} ;

/** \brief Helper class. sample-able circle given a normal and a center in 3-space*/
class OrientableCircle : public ConvexGeometry
{
protected:
	
	double radius ;
	Point normal ;
	
public:
	
	/** \brief Construct a Circle from a center, a radius and a normal*/
	OrientableCircle(double radius,double x, double y, double z, Point normal ) ;
	
	/** \brief Construct a Circle from a center, a radius and a normal*/
	OrientableCircle(double radius,const Point * p0, Point normal ) ;
	
	/** \brief Construct a Circle from a center, a radius and a normal*/
	OrientableCircle(double radius,const Point p0, Point normal ) ; 
	
	/** \brief default circle, the trigonometric circle*/
	OrientableCircle() ; 
	
	virtual ~OrientableCircle() { } ;

	virtual XMLTree * toXML() ;

	/** \brief Get points sampling the circle boundary*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const;
	
	/** \brief Sample the circle boundary*/
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	/** \brief sample the circle surface*/
	virtual void sampleSurface(size_t num_points);
	
	/** \brief return true if the argument lies in the circle*/
	virtual bool in(const Point & v) const ;
	
	/** \brief return the area of the circle*/
	virtual double area() const ;
	
	/** \brief return 0*/
	virtual double volume() const ;
	
	/** \brief Project the argument on the circle boundary*/
	virtual void project(Point * p) const;
	
	/** \brief do nothing */
	virtual void computeCenter() ;
	
	/** \brief return the radius*/
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_THREE_DIMENSIONAL ;
	}
	
	/** \brief Return empty vector*/
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;



} ;

/** \brief Return the triproduct of the three arguments (A ^ B * C)*/
double signedAlignement(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1) ;

/** \brief Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
bool isAligned(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1)  ;

/** \brief Check the alignment of three points.
 * 
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
bool isAligned(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1)  ;

/** \brief return true if the four points are coplanar*/
bool isCoplanar(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2) ;

/** \brief return true if the four points are coplanar*/
bool isCoplanar(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double coplanarity(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double coplanarity(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double signedCoplanarity(const Mu::Point &test, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double signedCoplanarity(const Mu::Point *test, const Mu::Point *f0, const Mu::Point *f1,const Mu::Point *f2) ;

/**  \brief Test if a point is in a triangle defined by three points.
 * 
 * @param test point to test.
 * @param p0 vertex 0.
 * @param p1 vertex 1.
 * @param p2 vertex 2.
 * @return true if test is in.
 */
bool isInTriangle(const Mu::Point &test, const Mu::Point &p0, const Mu::Point &p1, const Mu::Point &p2)  ;

/**  \brief Test wether two points lie on the same demi-plane.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point &test, const Mu::Point &witness, const Mu::Point &f0, const Mu::Point &f1)  ;

/**  \brief Test wether two points lie on the same demi-plane.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point *test, const Mu::Point *witness, const Mu::Point *f0, const Mu::Point *f1)  ;

/**  \brief Test wether two points lie on the same demi-space.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
  * @param f1 third point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point &test, const Mu::Point &witness, const Mu::Point &f0, const Mu::Point &f1, const Mu::Point &f2)  ;

/**  \brief Test wether two points lie on the same demi-space.
 * 
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @param f1 third point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Mu::Point * test, const Mu::Point * witness, const Mu::Point * f0, const Mu::Point * f1, const Mu::Point * f2)  ;
//bool isAligned(const Point test, const Point f0, const Point f1)  ;

/** \brief Return the distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Mu::Point &v1, const Mu::Point &v2) ;

/** \brief Return the distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Mu::Point * v1, const Mu::Point * v2) ;

/** \brief Return the convex hull of a set of points.
 * 
 * @param points 
 * @return a convex polygon (all boundary points anti-cockwise-ordered).
 */
Mu::ConvexPolygon* convexHull(const std::vector<Mu::Point> * points) ;

/** \brief Functor for point comparison in STL containers */
struct PointLessThan
{
	bool operator()(Mu::Point * p1, Mu::Point *p2)
	{
		return *p1 < *p2 ;
	}
} ;

/** \brief Return the square distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point &v1, const Mu::Point &v2) ;

/** \brief Return the square distance between two points
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Mu::Point *v1, const Mu::Point *v2) ;

/** \brief Return the square distance between two points, 2D case
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist2D(const Mu::Point &v1, const Mu::Point &v2) ;

/** \brief Return the square distance between two points, 2D case
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist2D(const Mu::Point *v1, const Mu::Point *v2) ;

/** \brief Return the square distance between two points, 2D case
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist3D(const Mu::Point &v1, const Mu::Point &v2) ;

double triProduct(const Mu::Point &A, const Mu::Point &B, const Mu::Point &C) ;

/** \brief Return the square distance between two points, 2D case
 * 
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist3D(const Mu::Point *v1, const Mu::Point *v2) ;

/** \brief Functor for point comparison in STL containers */
struct PointEqTol
{
	double tol ;
	PointEqTol(double t) : tol(t) {}
	bool operator()(const Mu::Point & m, const Mu::Point & p)
	{
		return squareDist(m,p) < tol ;
	}
} ;

/** \brief Functor for point comparison in STL containers */
struct PointLess_Than_x
{
	bool operator()(const Mu::Point & p1, const Mu::Point & p2) const
	{
		return p1.x < p2.x ;
	}
} ;

/** \brief Functor for point comparison in STL containers */
struct PointLess_Than_y
{
	bool operator()(const Mu::Point & p1, const Mu::Point & p2) const
	{
		return p1.y < p2.y ;
	}
} ;


#endif
