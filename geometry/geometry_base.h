// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_BASE_H_
#define __GEOMETRY_BASE_H_
#include "../utilities/matrixops.h"

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
virtual void setCenter(const Point & newCenter)            \
{                                                          \
      this->__geo_type__::setCenter (newCenter)        ;  \
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



namespace Amie
{

static const double POINT_TOLERANCE_2D =  1e-10 ;//std::numeric_limits<double>::epsilon() ;
static const double POINT_TOLERANCE_3D =  1e-10 ;//std::numeric_limits<double>::epsilon() ;

/** \brief defines the implemented geometries */
typedef enum
{
    nullptr_GEOMETRY,
    CIRCLE,
    LAYERED_CIRCLE,
    TRIANGLE,
    RECTANGLE,
    PARALLELOGRAMME,
    CONVEX_POLYGON,
    SEGMENTED_LINE,
    POLYGON,
    ORIENTABLE_CIRCLE,
    CLOSED_NURB,
    TETRAHEDRON = 16,
    HEXAHEDRON = 17,
    SPHERE = 18,
    LAYERED_SPHERE = 19,
    REGULAR_OCTAHEDRON,
    POLYGON_PRISM,
    LOFTED_POLYGON,
    ELLIPSE,
    LEVEL_SET,
    TIME_DEPENDENT_CIRCLE,
} GeometryType ;

/** \brief defines the possible geometrical transformations  */
typedef enum
{
    ROTATE,
    SCALE,
    TRANSLATE,
} GeometricTransformationType ;

/** \brief defines the possible vertices, axes and faces of a bounding box */
typedef enum
{
    TOP,
    LEFT,
    BOTTOM,
    RIGHT,
    FRONT,
    BACK,
    BEFORE,
    NOW,
    AFTER,
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
    BOTTOM_BACK,
    TOP_BACK,
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
    BOTTOM_BACK_BEFORE,
    TOP_BACK_BEFORE,
    TOP_NOW,
    LEFT_NOW,
    BOTTOM_NOW,
    RIGHT_NOW,
    FRONT_NOW,
    BACK_NOW,
    TOP_LEFT_NOW,
    TOP_RIGHT_NOW,
    BOTTOM_LEFT_NOW,
    BOTTOM_RIGHT_NOW,
    FRONT_LEFT_NOW,
    FRONT_RIGHT_NOW,
    BACK_LEFT_NOW,
    BACK_RIGHT_NOW,
    FRONT_TOP_NOW,
    BOTTOM_BACK_NOW,
    TOP_BACK_NOW,
    FRONT_BOTTOM_NOW,
    TOP_LEFT_FRONT_NOW,
    TOP_LEFT_BACK_NOW,
    BOTTOM_LEFT_FRONT_NOW,
    BOTTOM_LEFT_BACK_NOW,
    TOP_RIGHT_FRONT_NOW,
    TOP_RIGHT_BACK_NOW,
    BOTTOM_RIGHT_FRONT_NOW,
    BOTTOM_RIGHT_BACK_NOW,
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
    BOTTOM_BACK_AFTER,
    TOP_BACK_AFTER,
} BoundingBoxPosition ;

/** \brief defines the different space dimensionnalities of the mesh.*/
typedef enum
{
    SPACE_ONE_DIMENSIONAL = 1,
    SPACE_TWO_DIMENSIONAL,
    SPACE_THREE_DIMENSIONAL
} SpaceDimensionality ;


class Segment ;
class ConvexPolygon ;

struct PtP ;
struct PtV ;

/** \brief Four dimensionnal point. Storage of the coordinates makes use of the processors's vector extensions*/
class Point
{
public :
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
public :

    /** \brief initialise x, y, z and t values
    *
    * @param x
    * @param y
    * @param z
    * @param t
    */
    Point(double x = 0, double y = 0, double z = 0, double t = 0) ;

    /** \brief copy-constructor.*/
    Point(const Point & p) ;

    /** \brief copy-constructor.*/
    Point& operator = (const Point & p) ;

    void setX(double v) ;
    void setY(double v) ;
    void setZ(double v) ;
    void setT(double v) ;
    void setId(int v) {
        id = v ;
    } ;

    const double & getT() const {
        return t ;
    } ;
    double & getT() {
        return t ;
    } ;
    const double & getZ() const {
        return z ;
    } ;
    double & getZ() {
        return z ;
    } ;
    const double & getY() const {
        return y ;
    } ;
    double & getY() {
        return y ;
    } ;
    const double & getX() const {
        return x ;
    } ;
    double & getX() {
        return x ;
    } ;

    const int & getId() const {
        return id ;
    } ;
    int & getId() {
        return id ;
    } ;

    void set(double v = 0, double vv = 0, double vvv = 0, double vvvv = 0) ;
    void set(const Point & p) ;

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
} ;


/** Robust point comparator.
    *
    * @param p \c Point to compare to.
    * @return is within 2*std::numeric_limits\<double\>::epsilon() of p.
    */
bool operator==(const Point & p_, const Point &p)  ;

/** Robust point comparator.
    *
    * @param p \c Point to compare to.
    * @return is outside of 2*std::numeric_limits\<double\>::epsilon() of p.
    */
bool operator!=(const Point & p_,const Point &p)  ;

/** Operator defining order in the plane.
    *
    * @param p \c Point to compare to
    * @return this is bottommost then leftmost compared to p.
    */
bool operator <(const Point & p_,const Point & p)  ;

/** Operator defining order in the plane.
    *
    * @param p \c Point to compare to
    * @return this is topmost then rightmost compared to p.
    */
bool operator >(const Point & p_,const Point &p)  ;

Point operator-(const Point & p_,const Point &p)  ;
Point operator-(const Point & p_,const Vector &p)  ;
Point operator+(const Point & p_,const Point &p)  ;
Point operator+(const Point & p_,const Vector &p)  ;
Point operator/(const Point & p_,const double p)  ;
Point operator*(const Point & p_,const double p)  ;
Point operator*(double v, const Point & p_)  ;

Point operator-(const Point & p_) ;

/** \brief dot-product*/
double operator*(const Point & p_,const Point &p)  ;

/** \brief dot-product*/
double operator*(const Point & p_,const Vector &p)  ;

/** \brief cross-product*/
PtP operator^(const Point & p_,const Point &p)  ;

/** \brief cross-product*/
PtV operator^(const Point & p_,const Vector &p) ;

Point & operator+=(Point &  p_ ,const Point &p) ;

Point & operator-=( Point & p_,const Point &p) ;

Point & operator*=( Point & p, const double d) ;

Point & operator*=( Point &p, const Matrix & m) ;

Point operator *( const Matrix & m, const Point &p) ;

Point operator* (const Point & p, const Matrix & m) ;

Point & operator/=( Point & p, const double d) ;


struct PtP
{
    const Point & f ;
    const Point & s ;
    PtP(const Point & f, const Point & s) :f(f), s(s) {};

    operator const Point() const
    {
        Point ret( f.getY()*s.getZ() - f.getZ()*s.getY(), f.getZ()*s.getX() - f.getX()*s.getZ(), f.getX()*s.getY() - f.getY()*s.getX()) ;
        ret.setId(std::max(f.getId(), s.getId()) ) ;
        return ret ;
    };

    double operator *(const Point & p) const
    {
        return (f.getY()*s.getZ() - f.getZ()*s.getY())*p.getX() + (f.getZ()*s.getX() - f.getX()*s.getZ())*p.getY() + (f.getX()*s.getY() - f.getY()*s.getX())*p.getZ() ;
    }

    Point operator * (const double & d) const
    {
        Point ret(d*(f.getY()*s.getZ() - f.getZ()*s.getY()), d*(f.getZ()*s.getX() - f.getX()*s.getZ()), d*(f.getX()*s.getY() - f.getY()*s.getX())) ;
        ret.setId(std::max(f.getId(), s.getId()) ) ;
        return ret ;
    }
    Point operator / (const double & d) const
    {
        Point ret((f.getY()*s.getZ() - f.getZ()*s.getY())/d, (f.getZ()*s.getX() - f.getX()*s.getZ())/d, (f.getX()*s.getY() - f.getY()*s.getX())/d) ;
        ret.setId(std::max(f.getId(), s.getId()) ) ;
        return ret ;
    }

    Point operator - (const Point & p) const
    {
        Point ret(f.getY()*s.getZ() - f.getZ()*s.getY()-p.getX(), f.getZ()*s.getX() - f.getX()*s.getZ()-p.getY(), f.getX()*s.getY() - f.getY()*s.getX() -p.getZ()) ;
        ret.setId(std::max(f.getId(), s.getId()) ) ;
        return ret ;
    }

    Point operator + (const Point & p) const
    {
        Point ret(f.getY()*s.getZ() - f.getZ()*s.getY() + p.getX(), f.getZ()*s.getX() - f.getX()*s.getZ() + p.getY(), f.getX()*s.getY() - f.getY()*s.getX() + p.getZ()) ;
        ret.setId(std::max(f.getId(), s.getId()) ) ;
        return ret ;
    }

    double norm() const
    {
        return sqrt( (f.getY()*s.getZ() - f.getZ()*s.getY()) * (f.getY()*s.getZ() - f.getZ()*s.getY()) + ( f.getZ()*s.getX() - f.getX()*s.getZ()) *  ( f.getZ()*s.getX() - f.getX()*s.getZ()) + ( f.getX()*s.getY() - f.getY()*s.getX() ) * ( f.getX()*s.getY() - f.getY()*s.getX() )) ;
    }

    double getZ() const
    {
        return f.getX()*s.getY() - f.getY()*s.getX() ;
    }
} ;

struct PtV
{
    const Point & f ;
    const Vector & s ;
    PtV(const Point & f, const Vector & s) :f(f), s(s) { };

    operator const Point() const
    {
        Point ret( f.getY()*s[2] - f.getZ()*s[1],f.getZ()*s[0] - f.getX()*s[2],f.getX()*s[1] - f.getY()*s[0]) ;
        ret.getId() = f.getId() ;
        return ret ;
    };
} ;

/** \brief light class for 3D triangle-line intersection computation*/
struct TriPoint
{
    Point normal ;
    Point center ;

    std::valarray<Point const *> point ;
    std::vector<TriPoint *> neighbour ;

    /** \brief constructor from three points */
    TriPoint(const Point * p0, const Point * p1, const Point * p2) ;

    /** \brief return area of the triangle*/
    double area() const;

    Vector normalv() const;

    Vector normalv(const Point & p) const ;

    const Point & getCenter() const {
        return center ;
    } ;
    const Point & first() const {
        return *point[0] ;
    } ;
    const Point & second() const {
        return *point[1] ;
    } ;
    const Point & third() const {
        return *point[2] ;
    } ;

    /** \brief return true is the argument is in*/
    bool in(const Point & p) const ;

    Point projection(const Point & p) const ;

    /** \brief Return Gauss points*/
    std::valarray<std::pair<Point, double> > getGaussPoints(bool timeDependent = false) const ;


};


typedef  std::valarray<Point *> PointArray;

/** \brief Basic interface for a geometrical object.
*
* A geometrical object is defined by a set of points, in and out discrimination and intersection and projection operators
*
* New method for the Geometry class must be added to the GEO_DERIVED_OBJECT defined above, otherwise AMIE will not compile!
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

    /** \brief Acessor set all bounding poinst of the geometry*/
    virtual void setBoundingPoints(const PointArray & nb) = 0;

    /** \brief Acessor set all inside poinst of the geometry*/
    virtual void setInPoints(const PointArray & nb) {
        inPoints.resize(nb.size()) ;
        inPoints=nb ;
    }

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

    /** \brief Area of the Geometry*/
    virtual double area() const = 0;

    /** \brief Volume of the Geometry*/
    virtual double volume() const = 0;

    /** \brief return true if this geometry intersects with the argument*/
    virtual bool intersects(const Geometry *) const ;

    /** \brief return the intersection points between this geometry and the argument*/
    virtual std::vector<Point> intersection(const Geometry *) const ;

    /** \brief Get the bounding box.
     * The points are topLeft, topRight, bottomRight, bottomLeft or in 3D:
     * center.getX()+0.5*size_x, center.getY()+0.5*size_y, center.getZ()+0.5*size_z) ;
     * center.getX()+0.5*size_x, center.getY()+0.5*size_y, center.getZ()-0.5*size_z) ;
     * center.getX()+0.5*size_x, center.getY()-0.5*size_y, center.getZ()+0.5*size_z) ;
     * center.getX()+0.5*size_x, center.getY()-0.5*size_y, center.getZ()-0.5*size_z) ;
     * center.getX()-0.5*size_x, center.getY()+0.5*size_y, center.getZ()+0.5*size_z) ;
     * center.getX()-0.5*size_x, center.getY()+0.5*size_y, center.getZ()-0.5*size_z) ;
     * center.getX()-0.5*size_x, center.getY()-0.5*size_y, center.getZ()+0.5*size_z) ;
     * center.getX()-0.5*size_x, center.getY()-0.5*size_y, center.getZ()-0.5*size_z) ;
     */
    virtual std::vector<Point> getBoundingBox() const = 0;

    /** \brief Return the number of space dimensions of the geometry (SPACE_TWO_DIMENSIONAL or SPACE_THREE_DIMENSIONAL)*/
    virtual SpaceDimensionality spaceDimensions() const = 0 ;

    /** \brief Return the number of time slices present in this geometry*/
    virtual size_t timePlanes() const ;

    virtual size_t & timePlanes() ;

    /** \brief return the fraction of the area of the current geometry also lying in the target geometry. */
    virtual double overlapFraction(const Geometry * target) const ;

// 	template<class GeoType> const GeoType * getPrimitive() const { return dynamic_cast<const GeoType *>(this) ; }
// 	template<class GeoType> GeoType * getPrimitive() { return dynamic_cast<GeoType *>(this) ; }


} ;


/** \brief Transforms the geometry*/
void transform(Geometry * g,  GeometricTransformationType transformation, const Point & p) ;

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

    void setOrigin( const Point & p_ ) {
        p = p_  ;
    }
    void setVector( const Point & p_ ) {
        v = p_  ;
    }

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
    double * getFirstX() {
        return &f.getX() ;
    }
    double * getFirstY() {
        return &f.getY() ;
    }
    double * getFirstZ() {
        return &f.getZ() ;
    }
    double * getSecondX() {
        return &s.getX() ;
    }
    double * getSecondY() {
        return &s.getY() ;
    }
    double * getSecondZ() {
        return &s.getZ() ;
    }
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

    /** \brief Return a normal to the segment.
     *	The vector is of norm 1.
     */
    Point normal() const ;

    Vector normalv() const
    {

        double n = norm();
        if(n > POINT_TOLERANCE_2D)
        {
            Vector ret(2) ;
            ret[0] = -vec.getY()/n ;
            ret[1] = vec.getX()/n ;
            return ret ;
        }

        return Vector(0., 2) ;
    }
    Vector normalv(const Point & p) const ;


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
    std::valarray<std::pair<Point, double> > getGaussPoints(bool timeDependent = false) const ;
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

    WeightedPoint()  {
        this->w = 1.0;
    }

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
    virtual std::vector<Point> getBoundingBox() const {
        return std::vector<Point>(0) ;
    }
} ;





/** \brief Return the triproduct of the three arguments (A ^ B * C)*/
double signedAlignement(const Amie::Point &test, const Amie::Point &f0, const Amie::Point &f1) ;

/** \brief Check the alignment of three points.
 *
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
bool isAligned(const Amie::Point &test, const Amie::Point &f0, const Amie::Point &f1)  ;

/** \brief Check the alignment of three points.
 *
 * @param test first point.
 * @param f0 second point.
 * @param f1 third point.
 * @return true if the points are aligned.
 */
bool isAligned(const Amie::Point *test, const Amie::Point *f0, const Amie::Point *f1)  ;

/** \brief return true if the four points are coplanar*/
bool isCoplanar(const Amie::Point *test, const Amie::Point *f0, const Amie::Point *f1,const Amie::Point *f2, double renorm = 1.) ;

/** \brief return true if the four points are coplanar*/
bool isCoplanar(const Amie::Point &test, const Amie::Point &f0, const Amie::Point &f1, const Amie::Point &f2, double renorm = 1.) ;

int coplanarCount( Amie::Point * const* pts, int numpoints, const Amie::Point &f0, const Amie::Point &f1, const Amie::Point &f2, double renorm = 1.) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double coplanarity(const Amie::Point &test, const Amie::Point &f0, const Amie::Point &f1, const Amie::Point &f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double coplanarity(const Amie::Point *test, const Amie::Point *f0, const Amie::Point *f1,const Amie::Point *f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double signedCoplanarity(const Amie::Point &test, const Amie::Point &f0, const Amie::Point &f1, const Amie::Point &f2) ;

/** \brief Compute a value increasing with decreasing coplanarity of the points*/
double signedCoplanarity(const Amie::Point *test, const Amie::Point *f0, const Amie::Point *f1,const Amie::Point *f2) ;

/**  \brief Test if a point is in a triangle defined by three points.
 *
 * @param test point to test.
 * @param p0 vertex 0.
 * @param p1 vertex 1.
 * @param p2 vertex 2.
 * @return true if test is in.
 */
bool isInTriangle(const Amie::Point &test, const Amie::Point &p0, const Amie::Point &p1, const Amie::Point &p2)  ;

/**  \brief Test wether two points lie on the same demi-plane.
 *
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Amie::Point &test, const Amie::Point &witness, const Amie::Point &f0, const Amie::Point &f1, double norm = 1.)  ;

/**  \brief Test wether two points lie on the same demi-plane.
 *
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Amie::Point *test, const Amie::Point *witness, const Amie::Point *f0, const Amie::Point *f1, double norm = 1.)  ;

/**  \brief Test wether two points lie on the same demi-space.
 *
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
  * @param f1 third point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Amie::Point &test, const Amie::Point &witness, const Amie::Point &f0, const Amie::Point &f1, const Amie::Point &f2, double renorm = 1.)  ;

/**  \brief Test wether two points lie on the same demi-space.
 *
 * @param test first point to test.
 * @param witness second point to test.
 * @param f0 first point defining the plane boundary.
 * @param f1 second point defining the plane boundary.
 * @param f1 third point defining the plane boundary.
 * @return true if both points are on the same side of the demi-plane.
 */
bool isOnTheSameSide(const Amie::Point * test, const Amie::Point * witness, const Amie::Point * f0, const Amie::Point * f1, const Amie::Point * f2, double renorm = 1.)  ;
//bool isAligned(const Point test, const Point f0, const Point f1)  ;

/** \brief Return the distance between two points
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Amie::Point &v1, const Amie::Point &v2) ;

/** \brief Return the distance between two points
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ \sqrt{(x_0-x_1)^2 + (y_0-y_1)^2} \f$
 */
double dist(const Amie::Point * v1, const Amie::Point * v2) ;

/** \brief Return the convex hull of a set of points.
 *
 * @param points
 * @return a convex polygon (all boundary points anti-cockwise-ordered).
 */
Amie::ConvexPolygon* convexHull(const std::vector<Amie::Point> * points) ;

/** \brief Functor for point comparison in STL containers */
struct PointLessThan
{
    bool operator()(Amie::Point * p1, Amie::Point *p2)
    {
        return *p1 < *p2 ;
    }
} ;

Matrix rotateToVector (Amie::Point * toRotate, const Amie::Point & toVector) ;
/** \brief Return the square distance between two points
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Amie::Point &v1, const Amie::Point &v2) ;

/** \brief Return the square distance between two points
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist(const Amie::Point *v1, const Amie::Point *v2) ;

/** \brief Return the square distance between two points, 2D case
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist2D(const Amie::Point &v1, const Amie::Point &v2) ;

/** \brief Return the square distance between two points, 2D case
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist2D(const Amie::Point *v1, const Amie::Point *v2) ;

/** \brief Return the square distance between two points, 2D case
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist3D(const Amie::Point &v1, const Amie::Point &v2) ;

double triProduct(const Amie::Point &A, const Amie::Point &B, const Amie::Point &C) ;

/** \brief Return the square distance between two points, 2D case
 *
 * @param v1 first point.
 * @param v2 second point.
 * @return \f$ (x_0-x_1)^2 + (y_0-y_1)^2 \f$
 */
double squareDist3D(const Amie::Point *v1, const Amie::Point *v2) ;

/** \brief Functor for point comparison in STL containers */
struct PointEqTol
{
    double tol ;
    PointEqTol(double t) : tol(t) {}
    bool operator()(const Amie::Point & m, const Amie::Point & p)
    {
        return squareDist(m,p) < tol ;
    }
} ;

/** \brief Functor for point comparison in STL containers */
struct PointLess_Than_x
{
    bool operator()(const Amie::Point & p1, const Amie::Point & p2) const
    {
        return p1.getX() < p2.getX() ;
    }
} ;

/** \brief Functor for point comparison in STL containers */
struct PointLess_Than_y
{
    bool operator()(const Amie::Point & p1, const Amie::Point & p2) const
    {
        return p1.getY() < p2.getY() ;
    }
} ;

} ;

#endif

