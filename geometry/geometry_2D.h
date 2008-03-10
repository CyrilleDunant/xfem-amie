// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#ifndef __GEOMETRY_2D_H_
#define __GEOMETRY_2D_H_

#include "geometry_base.h"

namespace Mu
{

class Triangle : public ConvexGeometry
{
protected:
	double radius ;
	double sqradius ;
	Point circumCenter ;
	
	void computeCircumCenter() ;
	virtual void computeCenter() ;
	
public:
	Triangle() ;
	Triangle(const Point & p0, const Point & p1, const Point & p2) ;
	Triangle( Point *p0,  Point *p1,  Point *p2) ;
	virtual ~Triangle() { } ;
	
	virtual bool inCircumCircle(const Point & p) const ;
	virtual bool inCircumCircle(const Point *p) const ;
	virtual Point getCircumCenter() const;
	virtual double getRadius() const ;
	virtual void sampleBoundingSurface(size_t num_points)  ;
	virtual void sampleSurface(size_t num_points)  ;
	virtual bool in(const Point & p) const ;
	virtual size_t sides() const { return 3 ; }
	
	virtual void project(Point *) const ;
	
	virtual double area() const ;
	virtual double volume() const { return 0 ;} 
	
	virtual void print() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const ;
	
} ;

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
	
	
	Rectangle(double x, double y, double originX, double originY) ;
	Rectangle(double x, double y, Point &center) ;
	Rectangle() ;
	virtual ~Rectangle() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) ;
	virtual void sampleSurface(size_t num_points) ; 
	
	virtual double width() const ;
	virtual double height() const ;
	virtual double area() const ;
	virtual double volume() const { return 0 ;} 

	virtual double getRadius() const ;
	
	virtual size_t sides() const { return 4 ; }
	
	virtual bool in(const Point & p) const;
	
	virtual void project(Point *) const;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const ;

	
} ;

class Parallelogramme : public ConvexGeometry
{
protected:
	double radius ;
	Point circumCenter ;
	
	void computeCircumCenter() ;
	virtual void computeCenter() ;
	
public:
	Parallelogramme() ;
	Parallelogramme(const Point & p0, const Point & p1, const Point&  p2, const Point & p3) ;
	Parallelogramme( const Point *p0,  const Point *p1,  const Point *p2,  const Point *p3) ;
	virtual ~Parallelogramme() { } ;
	
	virtual bool inCircumCircle(const Point p) const ;
	virtual bool inCircumCircle(const Point *p) const ;
	virtual Point getCircumCenter() const;
	virtual double getRadius() const ;
	virtual void sampleBoundingSurface(size_t num_points)  ;
	virtual void sampleSurface(size_t num_points)  ;
	virtual bool is1D() const ;
	virtual bool in(const Point & p) const;
	virtual size_t sides() const { return 4 ; }
	
	virtual void project(Point *) const ;
	
	virtual double area() const ;
	virtual double volume() const { return 0 ;} 


	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;


class Circle : public ConvexGeometry
{
protected:
	
	double radius ;
	double sqradius ;
	virtual void computeCenter() ;
	
public:
	/** Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	Circle(double r,double originX,double originY) ;
	
	/** Construct a circle. This constructor is provided for convenience and a copy of the point is made and storder in the internal points.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	Circle(double r, const Point *center) ; 
	
	/** Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	Circle(double r, const Point center) ; 
	virtual ~Circle() { } ;
	
	/** Sample the bounding Surface.
	 * 
	 * num_points Points are placed equidistantly on the surface. The first point is placed at \f$ \theta = 0 \f$.
	 * 
	 * @param num_points points to place on the surface.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	
	/** Sample the disc.
	 * 
	 * Points are placed in concentric circles, rotated from half the base angle each time. The spacing of the circle is so calculated as to have triangles as nearly equilateral as a regular spacing would allow.
	 * 
	 * @param num_points number of points <b>on the boundary</b>.
	 */
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point &v) const ;
	
	/** Return the circle Radius.
	 * 
	 * @return the radius.
	 */
	virtual double getRadius() const ;
	
	/** Calculate the area.
	 * 
	 * Using the usual \f$ \pi r^2 \f$.
	 * 
	 * @return the area.
	 */
	virtual double area() const ;
	virtual double volume() const { return 0 ;} 

	virtual void setRadius(double newr);
		
	/** Project the point on the circle.
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


class LayeredCircle : public Circle
{
protected:
	
	double radius ;
	double sqradius ;
	std::vector<double> radiuses ;
	
public:
	/** Construct a circle.
	 * 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	LayeredCircle(std::vector<double> radii,double originX,double originY) ;
	
	/** Construct a circle. This constructor is provided for convenience and a copy of the point is made and storder in the internal points.
	 * 
	 * @param r radii of the circles to construct. Will be reordered.
	 * @param center Center of the circle.
	 */
	LayeredCircle(std::vector<double> radii, const Point center) ; 
	
	/** Construct a circle.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center Center of the circle.
	 */
	LayeredCircle(double r, const Point center) ; 
	virtual ~LayeredCircle() { } ;
	
	/** Sample the disc.
	 * 
	 * Points are placed in concentric circles, rotated from half the base angle each time. The spacing of the circle is so calculated as to have triangles as nearly equilateral as a regular spacing would allow.
	 * 
	 * @param num_points number of points <b>on the boundary</b>.
	 */
	virtual void sampleSurface(size_t num_points);

	virtual void setRadius(double newr);
	virtual void addRadius(double newr);

	virtual const std::vector<double> & getRadii() const ;
	
} ;

class SegmentedLine :  public NonConvexGeometry
{
protected:
	virtual void computeCenter() ;
	
public:
	SegmentedLine(const std::valarray<Point *> & points) ;
	virtual ~SegmentedLine() { };
	
	/** Sample The bounding surface.
	 * 
	 * This function does nothing. It could, however be used to place points between the kinks.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	
	/** Sample The bounding surface.
	 * 
	 * This function does nothing. It could, however be used to place points between the kinks.
	 * 
	 * @param num_points number of points to use for the sampling.
	 */
	virtual void sampleSurface(size_t num_points) ;
	
	
	/** Is point in ?
	 * 
	 * @param v point to test.
	 * @return true if one of the kinks or extremities.
	 */
	virtual bool in(const Point & v) const;
	
	virtual double area() const { return 0 ;}
	virtual double volume() const { return 0 ;} 

	
	/** Return first extremity.
	 * 
	 * @return the first point of the set. This is equivalent to getBoundingPoint(0).
	 */
	virtual Point * getHead() const ;
	
	/** Return second extremity.
	 * 
	 * @return the last point of the set. This is equivalent to getBoundingPoint(getBoundingPoints().size()-1).
	 */
	virtual Point * getTail() const ;
		
	virtual void project(Point *) const;
	
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_ONE_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
} ;	

} ;

#endif
