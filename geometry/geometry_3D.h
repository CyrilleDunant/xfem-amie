// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_3D_H_
#define __GEOMETRY_3D_H_

#include "geometry_base.h"
#include<list>
#include<deque>

namespace Mu
{

/** \brief Tetrahedron geometry. Serves as a base for tetrahedral elements.*/
class Tetrahedron : public ConvexGeometry
{

protected:
	double radius ;
	double sqradius ;
	Point circumCenter ;
	void computeCircumCenter() ;
	
	
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

	virtual XMLTree * toXML() ;
	
	/** \brief Sample The bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;

	/** \brief Get set of points sampling the bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
/** \brief sample the volume of the tetrahedron*/
	virtual void sampleSurface(size_t num_points);

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

	virtual XMLTree * toXML() ;

	
	/** \brief Sample The bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;

	/** \brief Get set of points sampling the bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;

/** \brief sample the volume of the tetrahedron*/
	virtual void sampleSurface(size_t num_points);

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

	virtual double getXSize() const ;
	virtual double getYSize() const ;
	virtual double getZSize() const ;
	
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

	virtual XMLTree * toXML() {return new XMLTree("regular octahedron") ; } ;

	
	/** \brief Sample The bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;

	/** \brief Get set of points sampling the bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;

/** \brief sample the volume of the tetrahedron*/
	virtual void sampleSurface(size_t num_points);

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
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const { return std::vector<Point>() ;};
	
	virtual void sampleSurface(size_t num_points) { };
	virtual bool in(const Point &v) const { return false ; }
	
	virtual double area() const ;
	virtual double volume() const { return 0. ;}
	
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

	std::vector<Point> getSamplingPointsOnSphere(size_t num_points, double radius, size_t iter = 50, size_t threshold = 500) const ;
	virtual void project(Point * p, double r) const;
	void smooth(std::vector<Point> & points, double r, size_t iter = 50) const ;
	
public:
	
	/** \brief Constructor. Build sphere from radius and center coordinates*/
	Sphere(double radius,double x, double y, double z ) ;
	
/** \brief Constructor. Build sphere from radius and center*/
	Sphere(double radius,const Point * p0 ) ;
	
/** \brief Constructor, build Sphere from radius and center */
	Sphere(double radius,const Point p0 ) ; 
	
/** \brief default constructor, build sphere of radius 1 and center (0,0,0)*/
	Sphere() ; 

	Sphere(XMLTree * xml) ;
	
	virtual ~Sphere() { } ;

	virtual XMLTree * toXML() ;
	
	/** \brief Sample The bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;

	/** \brief Get set of points sampling the bounding surface.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;

	std::vector<Point> getStandardSamplingBoundingPointsOnSphere(size_t n)	const ;

/** \brief sample the volume of the sphere*/
	virtual void sampleSurface(size_t num_points);
	
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

	static void dumpSampleBoundingPoints(size_t n, size_t iter = 50) ;
	static std::vector<Point> importStandardBoundingPoints(size_t n) ;

} ;

} ;
#endif
