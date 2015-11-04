// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_3D_H_
#define __GEOMETRY_3D_H_

#include "geometry_base.h"
#include "geometry_2D.h"
#include<list>
#include<deque>

namespace Amie
{

/** \brief Tetrahedron geometry. Serves as a base for tetrahedral elements.*/
class Tetrahedron : public ConvexGeometry
{

protected:
    double cachedarea ;
    double cachedvolume ;
    double radius ;
    double sqradius ;
    Point circumCenter ;
    virtual void computeCircumCenter() ;
    virtual double computeArea()  ;
    virtual double computeVolume()  ;

public:

    /** \brief Constructor. Construct a tetrahedron from 4 points*/
    Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3) ;

    /** \brief Constructor. Construct a "quadratic" tetrahedron from six points.
     *useful as a base for elements to be directly assembled, as opposed to produced by a mesh.*/
    Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7) ;

    /** \brief Constructor, construct a tetrahedron from four points*/
    Tetrahedron(const Point & p0, const Point & p1, const Point &p2, const Point &p3) ;

    Tetrahedron() ;

    virtual ~Tetrahedron() { } ;

    /** \brief Sample The bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual void sampleBoundingSurface(double linearDensity) ;

    /** \brief Get set of points sampling the bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual std::vector<Point> getSamplingBoundingPoints(double linearDensity) const ;

    /** \brief sample the volume of the tetrahedron*/
    virtual void sampleSurface(double linearDensity);

    /** \brief return in if the argument lies in the tetrahedron, as defined by its four vertices.*/
    virtual bool in(const Point &v) const ;

    /** \brief return the surface of the tetrahedron*/
    virtual double area() const ;

    /** \brief return the volume of the tetrahedron*/
    virtual double volume() const ;

    /** \brief project the argument on the surface of the tetrahedron*/
    virtual void project(Point * p) const;

    /** \brief return true if the argument lies in the circumsphere of the tetrahedron*/
    virtual bool inCircumSphere(const Point & p) const ;

    /** \brief return true if the argument lies in the circumsphere of the tetrahedron*/
    virtual bool inCircumSphere(const Point *p) const ;

    virtual void computeCenter() ;
    void scale(double s)
    {
        circumCenter *= s ;
        radius *= s ;
        sqradius *= s*s ;
    } ;

    /** \brief return pointer to the circumcenter of the tetrahedron*/
    virtual const Point & getCircumCenter() const;

    /** \brief return the radius of the circumsphere*/
    virtual double getRadius() const ;

    virtual SpaceDimensionality spaceDimensions() const
    {
        return SPACE_THREE_DIMENSIONAL ;
    }

    virtual std::vector<Point> getBoundingBox() const ;

} ;

/** \brief Hexahedral geometry with axes aligned to the coordinate axes
*
* This geometry serves as a base for hexahedral elements
*/
class Hexahedron : public ConvexGeometry
{
protected:

    double size_y ;
    double size_x ;
    double size_z ;

public:

    /** \brief Construct Hex from six points*/
    Hexahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7) ;

    /** \brief Construct Hex from six points*/
    Hexahedron(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5, Point p6, Point p7) ;

    /** \brief construct hex from length, width and breadth and origin.*/
    Hexahedron(double x, double y, double z, double originX, double originY, double originZ) ;

    /** \brief construct hex from length, width and breadth and origin.*/
    Hexahedron(double x, double y, double z, const Point & c) ;

    Hexahedron() ;

    virtual ~Hexahedron() { } ;

    /** \brief Sample The bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual void sampleBoundingSurface(double linearDensity) ;

    /** \brief Get set of points sampling the bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual std::vector<Point> getSamplingBoundingPoints(double linearDensity) const ;

    /** \brief sample the volume of the tetrahedron*/
    virtual void sampleSurface(double linearDensity);

    /** \brief return true if the argument lies in the geometry*/
    virtual bool in(const Point &v) const ;

    /** \brief return the area*/
    virtual double area() const ;

    /** \brief the volume*/
    virtual double volume() const ;

    /** \brief project the argument on the surface of the hex*/
    virtual void project(Point * p) const;

    virtual void computeCenter() ;

    /** \brief return the circumradius*/
    virtual double getRadius() const ;

    virtual SpaceDimensionality spaceDimensions() const
    {
        return SPACE_THREE_DIMENSIONAL ;
    }

    virtual std::vector<Point> getBoundingBox() const ;

    virtual double getXSize() const final;
    virtual double getYSize() const final;
    virtual double getZSize() const final;

} ;

/** \brief Regular octahedron, with xOy aligned base geometry with axes aligned to the coordinate axes
*
* This geometry serves as a base for hexahedral elements
*/
class RegularOctahedron : public ConvexGeometry
{
protected:

    double length ;

public:

    /** \brief Construct Oct from base size and centre*/
    RegularOctahedron(double size, double x, double y, double z) ;

    /** \brief Construct Oct from base size and centre*/
    RegularOctahedron(double size, Point c) ;

    virtual ~RegularOctahedron() { } ;


    /** \brief Sample The bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual void sampleBoundingSurface(double linearDensity) ;

    /** \brief Get set of points sampling the bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual std::vector<Point> getSamplingBoundingPoints(double linearDensity) const ;

    /** \brief sample the volume of the tetrahedron*/
    virtual void sampleSurface(double linearDensity);

    /** \brief return true if the argument lies in the geometry*/
    virtual bool in(const Point &v) const ;

    /** \brief return the area*/
    virtual double area() const ;

    /** \brief the volume*/
    virtual double volume() const ;

    /** \brief project the argument on the surface of the hex*/
    virtual void project(Point * p) const;

    virtual void computeCenter() ;

    /** \brief return the circumradius*/
    virtual double getRadius() const ;

    virtual SpaceDimensionality spaceDimensions() const
    {
        return SPACE_THREE_DIMENSIONAL ;
    }

    virtual std::vector<Point> getBoundingBox() const ;

} ;

class PolygonPrism :  public NonConvexGeometry
{
protected:
    virtual void computeCenter() ;
    
    Polygon base ;
    Point axis ;
    Point origin ;
    Matrix rotationMatrix ;
public:
    PolygonPrism(const std::valarray<Point *> & points, const Point & vector, const Point & origin) ;
    
    virtual ~PolygonPrism() ;

    virtual void sampleBoundingSurface(size_t num_points) ;
    
    virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
    
    virtual void sampleSurface(size_t num_points) ;
    
    virtual bool in(const Point & v) const ;
    
    virtual double area() const ;
    
    virtual double volume() const ;

    virtual void project(Point * init) const ;
    
    virtual double getRadius() const ;
    
    virtual SpaceDimensionality spaceDimensions() const ;
    
    virtual std::vector<Point> getBoundingBox() const ;
} ;

class LoftedPolygonPrism :  public NonConvexGeometry
{
    friend class Geometry ;
protected:
    virtual void computeCenter() ;
    
    Polygon base ;
    double sweepNorm(double end = 1.) const ;
    double sweepNorm( const Point & offset, double end = 1.) const ;
    
    bool centerInXYPlane ;
    bool centerInXZPlane ;
    bool centerInYZPlane ;
    
    Point vstart ;
    Point vend ;
    
    Point tanAtCentre ;
    std::vector<Point> interpolationPoints ;
    std::pair<Point,Point> interpolatingPointAndTangent(double t) const ;
    std::pair<Point,Point> interpolatingPointAndTangent(double t, const Point & offset) const ;
    Point projectTest(const Amie::Point& v, const std::pair< Amie::Point, Amie::Point >& toCenter) const ;
    
    Matrix rotateToVector(const Amie::Point& fromvector, const Amie::Point& tovector) const ;
    Matrix rotateFromVector(const Point & fromvector,const Point & vector) const ;
    std::pair<Point,Point> projectToCenterLine(const Point & p, double * coordinate = nullptr) const ;
    std::pair<Point,Point>projectToOffsetCenterLine(const Point & p, const Point & offset) const ;
     std::vector<double> isoDistribute(int npoints) const;
     std::vector<double> isoDistribute(const Point & offset, int npoints) const;
public:
    LoftedPolygonPrism(const std::valarray<Point *> & points, const std::vector<Point> & interpolationPoints) ;
    
    virtual ~LoftedPolygonPrism() ;
    virtual void sampleBoundingSurface(size_t num_points) ; 
    virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ; 
    virtual void sampleSurface(size_t num_points) ;   
    virtual bool in(const Point & v) const ;  
    virtual double area() const ;  
    virtual double volume() const ;
    virtual void project(Point * init) const ; 
    virtual double getRadius() const ; 
    virtual SpaceDimensionality spaceDimensions() const ;  
    virtual std::vector<Point> getBoundingBox() const ;
    virtual std::pair<Point, Point> getEndNormals() const ;
} ;

/** \brief Class for the description triangulated free surfaces. Work in progress.*/
class TriangulatedSurface : public NonConvexGeometry
{
protected:
    std::vector<const Point *> boundary ;
    std::vector<TriPoint> mesh ;

public:
    void addPoint(Point *) ;

public:

    TriangulatedSurface(const Point * p0, const Point * p1, const Point * p2) ;

    TriangulatedSurface() ;

    virtual ~TriangulatedSurface() { } ;

    virtual void sampleBoundingSurface(size_t num_points) { };
    virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const {
        return std::vector<Point>() ;
    };

    virtual void sampleSurface(size_t num_points) { };
    virtual bool in(const Point &v) const {
        return false ;
    }

    virtual double area() const ;
    virtual double volume() const {
        return 0. ;
    }

    virtual void project(Point * p) const;

    virtual void computeCenter() ;
    virtual double getRadius() const ;

    virtual SpaceDimensionality spaceDimensions() const
    {
        return SPACE_THREE_DIMENSIONAL ;
    }

    virtual std::vector<Point> getBoundingBox() const ;
} ;

/** \brief Spherical geometry from center and radius*/
class Sphere : public ConvexGeometry
{
protected:

    double radius ;
    double sqradius ;

    std::vector<Point> getSamplingPointsOnSphere(double linearDensity, double radius, size_t iter = 100, size_t threshold = 500) const ;
    virtual void project(Point * p, double r) const;
    void smooth(std::vector<Point> & points, double r, size_t iter = 100) const ;

public:

    /** \brief Constructor. Build sphere from radius and center coordinates*/
    Sphere(double radius,double x, double y, double z ) ;

    /** \brief Constructor. Build sphere from radius and center*/
    Sphere(double radius,const Point * p0 ) ;

    /** \brief Constructor, build Sphere from radius and center */
    Sphere(double radius,const Point p0 ) ;

    /** \brief default constructor, build sphere of radius 1 and center (0,0,0)*/
    Sphere() ;

    virtual ~Sphere() { } ;

    /** \brief Sample The bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual void sampleBoundingSurface(double linearDensity) ;

    /** \brief Get set of points sampling the bounding surface.
     *
     * @param num_points number of points to use for the sampling.
     */
    virtual std::vector<Point> getSamplingBoundingPoints(double linearDensity) const ;

    std::vector<Point> getStandardSamplingBoundingPointsOnSphere(size_t n)	const ;

    /** \brief sample the volume of the sphere*/
    virtual void sampleSurface(double linearDensity);

    /** \brief return true if the argument lies in the sphere*/
    virtual bool in(const Point &v) const ;

    /** \brief return \f$ 4/3 \pi r^3\f$ */
    virtual double area() const ;

    /** \brief return \f$ 4 \pi r^2\f$ */
    virtual double volume() const ;

    /** \brief Project argument on the sphere surface. If the argument is the sphere center, do nothing*/
    virtual void project(Point * p) const;

    virtual void computeCenter() ;
    virtual double getRadius() const ;

    /** \brief set new radius for the sphere, if the sphere is sampled, scale the sampling points accordingly*/
    virtual void setRadius(double newr);

    virtual SpaceDimensionality spaceDimensions() const
    {
        return SPACE_THREE_DIMENSIONAL ;
    }

    virtual std::vector<Point> getBoundingBox() const ;

    static void dumpSampleBoundingPoints(double linearDensity, size_t iter = 50) ;
    static std::vector<Point> importStandardBoundingPoints(size_t n) ;

} ;

} 
#endif
