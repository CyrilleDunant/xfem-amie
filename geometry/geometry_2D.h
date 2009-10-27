// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_2D_H_
#define __GEOMETRY_2D_H_

#include "geometry_base.h"

namespace Mu
{

/** \brief Triangle lying in 2-space.*/
class Triangle : public ConvexGeometry
{
protected:
	double radius ;
	double sqradius ;
	Point circumCenter ;
	
	void computeCircumCenter() ;
	virtual void computeCenter() ;
	
public:
	/** \brief Default constructor, create a (0,1) (0,0) (1,0) triangle*/
	Triangle() ;
	
	/** \brief Construct a Triangle from three points*/
	Triangle(const Point & p0, const Point & p1, const Point & p2) ;
	
	/** \brief Construct a Triangle from three points*/
	Triangle( Point *p0,  Point *p1,  Point *p2) ;
	virtual ~Triangle() { } ;
	
	/** \brief Return true if the argument lies in the Circumcenter of this triangle*/
	virtual bool inCircumCircle(const Point & p) const ;
	
	/** \brief Return true if the argument lies in the Circumcenter of this triangle*/
	virtual bool inCircumCircle(const Point *p) const ;
	
	/** \brief return reference to the circumcenter*/
	virtual const Point & getCircumCenter() const;
	
	/** \brief return circumRadius*/
	virtual double getRadius() const ;
	
	/** \brief sample the triangle borders*/
	virtual void sampleBoundingSurface(size_t num_points)  ;
	
	/** \brief sample the triangle surface*/
	virtual void sampleSurface(size_t num_points)  ;
	
	/** \brief return true if the argument is in the triangle*/
	virtual bool in(const Point & p) const ;
	
	/** \brief return 3*/
	virtual size_t sides() const { return 3 ; }
	
	/** \brief return a set of sampling points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
	/** \brief project the argument to the triangle boundary*/
	virtual void project(Point *) const ;
	
	/** \brief return area*/
	virtual double area() const ;
	
	/** \brief return 0*/
	virtual double volume() const { return 0 ;} 
	
	virtual void print() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const ;
	
} ;

/** \brief Rectangle class. Sides are aligned to the axes.*/
class Rectangle : public ConvexGeometry
{
protected:
	double size_y ;
	double size_x ;
	
	virtual void computeCenter() ;
	
	Point topLeft ;
	Point topRight ;
	Point bottomLeft ;
	Point bottomRight ;
	
	size_t numberOfPointsAlongX ;
	size_t numberOfPointsAlongY ;
	
public: 
	
	/** \brief constructor, build Rectangle from side length and center*/
	Rectangle(double x, double y, double originX, double originY) ;
	
	/** \brief constructor, build Rectangle from side length and center*/
	Rectangle(double x, double y, const Point &center) ;
	
	/** \brief Default constructor. build rectangle from (-1,-1) to (1,1)*/
	Rectangle() ;
	virtual ~Rectangle() { } ;

	/** \brief Sample the bounding surface with a given number of sampling points. the points are stored as boundingPoints*/
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	/** \brief Sample the surface with a given number of sampling points. the points are stored as inPoints*/
	virtual void sampleSurface(size_t num_points) ; 
	
	/** \brief get points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
	/** \brief return the width*/
	virtual double width() const ;
	
	/** \brief return the height*/
	virtual double height() const ;
	
	/** \brief compute the area. A = x*y*/
	virtual double area() const ;
	
	/** \brief return 0*/
	virtual double volume() const { return 0 ;} 

	/** \brief Return the circumscribing radius*/
	virtual double getRadius() const ;
	
	/** \brief Return 4*/
	virtual size_t sides() const { return 4 ; }
	
	/** \brief Return true if the argument is in*/
	virtual bool in(const Point & p) const;
	
	/** \brief Project the argument on the Recatangle Surface*/
	virtual void project(Point *) const;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	/** \brief return the boundingBox*/
	virtual std::vector<Point> getBoundingBox() const ;

} ;

/** \brief Rectangle normal to an axis in 3-space
* 
* This is used as a helper class for regular grids
*/
class OrientedRectangle : public ConvexGeometry
{
protected:
	double radius ;
	Point circumCenter ;
	
	void computeCircumCenter() ;
	virtual void computeCenter() ;
	
public:
	/** \brief Default constructor. Constructs a rectangle from -1,-1 to 1,1*/
	OrientedRectangle() ;
	
	/** \brief construct a Rectangle from four points. The points will be assumed to be planar*/
	OrientedRectangle(const Point & p0, const Point & p1, const Point&  p2, const Point & p3) ;
	
	/** \brief construct a Rectangle from four points. The points will be assumed to be planar*/
	OrientedRectangle( const Point *p0,  const Point *p1,  const Point *p2,  const Point *p3) ;
	virtual ~OrientedRectangle() { } ;
	
	/** \brief Return true if the argument is in the circumcircle*/
	virtual bool inCircumCircle(const Point p) const ;
	
	/** \brief Return true if the argument is in the circumcircle*/
	virtual bool inCircumCircle(const Point *p) const ;
	
	/** \brief Return the center of the rectangle*/
	virtual Point getCircumCenter() const;
	
	/** \brief Return the circumradius*/
	virtual double getRadius() const ;
	
	/** \brief Sample the bounding surface*/
	virtual void sampleBoundingSurface(size_t num_points)  ;
	
	/** \brief return set of points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
	/** \brief sample the surface of the rectangle*/
	virtual void sampleSurface(size_t num_points)  ;
	
	/** \brief return false*/
	virtual bool is1D() const ;
	
	/** \brief return true is the argument is in*/
	virtual bool in(const Point & p) const;
	
	/** \brief return 4*/
	virtual size_t sides() const { return 4 ; }
	
	/** \brief project to the boundary of the rectangle*/
	virtual void project(Point *) const ;
	
	/** \brief return x*y*/
	virtual double area() const ;
	
	/** \brief retunr 0*/
	virtual double volume() const { return 0 ;} 


	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;

/** \brief Circle defined from a center and a radius, in the XY plane*/
class Circle : public ConvexGeometry
{
protected:
	
	double radius ;
	double sqradius ;
	virtual void computeCenter() ;
	
public:
	/** \brief Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	Circle(double r,double originX,double originY) ;
	
	/** \brief Construct a circle. This constructor is provided for convenience and a copy of the point is made and storder in the internal points.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	Circle(double r, const Point *center) ; 
	
	/** \brief Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	Circle(double r, const Point center) ; 
	virtual ~Circle() { } ;
	
	/** \brief Sample the bounding Surface.
	 * 
	 * num_points Points are placed equidistantly on the surface. The first point is placed at \f$ \theta = 0 \f$.
	 * 
	 * @param num_points points to place on the surface.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	/** \brief return set of points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
	/** \brief return set of points sampling the bounding surface, given two angles (in radians)*/
	virtual std::vector<Point> getSamplingBoundingPointsOnArc(size_t num_points, const Point & start, const Point & finish) const ;
	
	/** \brief Sample the disc.
	 * 
	 * Points are placed in concentric circles, rotated from half the base angle each time. The spacing of the circle is so calculated as to have triangles as nearly equilateral as a regular spacing would allow.
	 * 
	 * @param num_points number of points <b>on the boundary</b>.
	 */
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point &v) const ;
	
	/** \brief Return the circle Radius.
	 * 
	 * @return the radius.
	 */
	virtual double getRadius() const ;

	/** \brief Return the circle Radius.
	 * 
	 * @return the radius.
	 */
	virtual double getSquareRadius() const ;
	
	/** \brief Calculate the area.
	 * 
	 * Using the usual \f$ \pi r^2 \f$.
	 * 
	 * @return the area.
	 */
	virtual double area() const ;
	
	/** \brief return 0*/
	virtual double volume() const { return 0 ;} 

	virtual void setRadius(double newr);
		
	/** \brief Project the point on the circle.
	 * 
	 * @param  p Point to project.
	 */
	virtual void project(Point * p) const;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const ;
	
	const Circle * getGeometry() const ;
	
	Circle * getGeometry() ;
	
} ;

/** \brief Set of concentric circles, from a center and a set of radii*/
class LayeredCircle : public Circle
{
protected:
	std::vector<double> radiuses ;
	
public:
	/** \brief Construct a circle.
	 * 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	LayeredCircle(std::vector<double> radii,double originX,double originY) ;
	
	/** \brief Construct a circle. This constructor is provided for convenience and a copy of the point is made and storder in the internal points.
	 * 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param center Center of the circle.
	 */
	LayeredCircle(std::vector<double> radii, const Point center) ; 
	
	/** \brief Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	LayeredCircle(double r, const Point center) ; 
	virtual ~LayeredCircle() { } ;
	
	/** \brief Sample the disc.
	 * 
	 * Points are placed in concentric circles, rotated from half the base angle each time. The spacing of the circle is so calculated as to have triangles as nearly equilateral as a regular spacing would allow.
	 * 
	 * @param num_points number of points <b>on the boundary</b>.
	 */
	virtual void sampleSurface(size_t num_points);

	/** \brief Set the external radius, all the other radii are scale accordingly*/
	virtual void setRadius(double newr);
	
	/** \brief add new radius. The new radius is automatically layered according to its value.*/
	virtual void addRadius(double newr);

	virtual const std::vector<double> & getRadii() const ;
	
} ;

/** \brief Segmented line. */
class SegmentedLine :  public NonConvexGeometry
{
protected:
	virtual void computeCenter() ;
	
public:
	SegmentedLine(const std::valarray<Point *> & points) ;
	virtual ~SegmentedLine() { };
	
	/** \brief Sample The bounding surface.
	 * 
	 * This function does nothing. It could, however be used to place points between the kinks.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
	/** \brief Sample The bounding surface.
	 * 
	 * This function does nothing. It could, however be used to place points between the kinks.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleSurface(size_t num_points) ;
	
	
	/** \brief Is point in ?
	 * 
	 * @param v point to test.
	 * @return true if one of the kinks or extremities.
	 */
	virtual bool in(const Point & v) const;
	
	/** \brief return 0*/
	virtual double area() const { return 0 ;}
	
	/** \brief return 0*/
	virtual double volume() const { return 0 ;} 

	
	/** \brief Return first extremity.
	 * 
	 * @return the first point of the set. This is equivalent to getBoundingPoint(0).
	 */
	virtual Point * getHead() const ;
	
	/** \brief Return second extremity.
	 * 
	 * @return the last point of the set. This is equivalent to getBoundingPoint(getBoundingPoints().size()-1).
	 */
	virtual Point * getTail() const ;
		
	virtual void project(Point *) const;
	
	/** \brief project the point on the segmented line. The projected point is the nearest point of the line to the original point*/
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_ONE_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;	

} ;

#endif
