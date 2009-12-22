// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009 (added: ellipses)
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

/** \brief Ellipse defined from a center, two radii and the main direction (as a normalized unit-vector), in the XY plane.*/
class Ellipse : public ConvexGeometry
{
protected:
	
	double majorradius ;
	double minorradius ;
	double sqradius ;
	double excentricity ;
	Point majoraxis ;
	virtual void computeCenter() ;
	
public:
	/** \brief Construct an ellipse.
	 * 
	 * @param a major radius of the ellipse to construct.
	 * @param b minor radius of the ellipse to construct.	 
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 * @param axisX global x coordinate of the main axis direction.
	 * @param axisY global y coordinate of the main axis direction.
	 */
	Ellipse(double a, double b, double originX, double originY, double axisX, double axisY) ;
	
	/** \brief Construct an ellipse. This constructor is provided for convenience and a copy of the point is made and storder in the internal points.
	 * 
	 * @param a major radius of the ellipse to construct.
	 * @param b minor radius of the ellipse to construct.	 
	 * @param center Center of the ellipse.
	 * @param axis Axis direction of the ellipse.
	 */
	Ellipse(double a, double b, const Point *center, const Point *axis) ; 
	
	/** \brief Construct an ellipse.
	 * 
	 * @param a major radius of the ellipse to construct.
	 * @param b minor radius of the ellipse to construct.	 
	 * @param center Center of the ellipse.
	 * @param axis Axis direction of the ellipse.
	 */
	Ellipse(double a, double b, const Point center, const Point axis) ; 
	
	/** \brief Construct an ellipse with a random axis direction.
	 * 
	 * @param a major radius of the ellipse to construct.
	 * @param b minor radius of the ellipse to construct.	 
	 * @param center Center of the ellipse.
	 */
	Ellipse(double a, double b, const Point center) ; 

	/** \brief Construct an ellipse with a random axis direction.
	 * 
	 * @param a major radius of the ellipse to construct.
	 * @param b minor radius of the ellipse to construct.	 
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	Ellipse(double a, double b, double x, double y) ; 

	/** \brief Copy an ellipse.
	 * 
	 * @param e ellipse to copy.
	 */
	Ellipse(const Ellipse &e) ;

	virtual ~Ellipse() { } ;
	
	/** \brief Computes a * b and store the result in sqradius
	 *
	 */
	virtual void setSqRadius() {this->sqradius = majorradius * minorradius ;} ;

	/** \brief Computes the ellipse excentricity
	 *
	 */
	virtual void setExcentricity() {this->excentricity = sqrt(1 - ((minorradius * minorradius) / (majorradius * majorradius))) ; } ;

	/** \brief Returns the ellipse excentricity
	 *
	 * @return the ellipse excentricity
	 */
	virtual double getExcentricity() const {return excentricity ; } ;

	/** \brief Returns the ellipse polar parameter
	 *
	 * @return the ellipse polar parameter
	 */
	virtual double getParameter() const {return majorradius * (1 - excentricity*excentricity) ; } ;

	/** \brief Returns the ellipse main axis angle with the (Ox) axis
	 *
	 * @return the angle (in radians)
	 */
	virtual double getAxisAngle() const ;

	/** \brief Returns one of the ellipse focus
	 *
	 * @param dir the direction of the focus along the main axis
	 * @return the focus
	 */
	virtual Point getFocus(bool dir) const ;

	/** \brief Returns both ellipse foci
	 *
	 * @return the foci as a pair of points, First point is the focus in the same direction as the main axis. Second point is the focus in the opposite direction.
	 */
	virtual const std::pair<Point, Point> getBothFocus() ;

	/** \brief Returns the ellipse minor axis
	 *
	 * @return the minor radius as a normalized unit-vector.
	 */
	virtual Point getMinorAxis() const ;

	/** \brief Changes the main axis
	 *
	 * @param newaxis the new axis direction
	 */
	virtual void setAxis(Point newaxis) ;

	/** \brief Returns the radius from the ellipse center as a function of the angle with the main axis.
	 *
	 * @param theta angle (in radians)
	 * @return the radius
	 */
	virtual double getRadiusOnEllipse(double theta) const ;

	/** \brief Returns the ellipse major radius.
	 *
	 * @return the ellipse major radius
	 */
	virtual double radius() const ;

	/** \brief Returns a times b.
	 *
	 * @return the ellipse major radius
	 */
	virtual double getSquareRadius() const ;

	/** \brief Returns the ellipse major radius.
	 *
	 * @return the ellipse major radius
	 */
	virtual double getRadius() const ;

	/** \brief Returns the rectangle containing the ellipse
	 *
	 * @return the points of the rectangle
	 */
	virtual std::vector<Point> getBoundingBox() const ;

	/** \brief Projects a point in local coordinates, so that the ellipse is transformed into a circle with a radius equal to the ellipse minor radius. This tranformation does not change the angles. Essentially, this is an homothety along the major axis, while coordinate along the minor axis is unchanged.
	 *
	 * @param the point to project
	 * @return the point in the new coordinates
	 */
	virtual Point toSmallCircle(Point p, Point origin, bool b) const ;

	/** \brief Returns the ellipse tangent as a function of the angle with the major axis
	 *
	 * @param angle with the ellipse major radius (in radians)
	 * @return direction of the ellipse tangent
	 */
	virtual Point getTangentDirection(double theta) const ;

	/** \brief Returns bounding points on arc (defined by two angles)
	 *
	 * @param num_points number of points to sample
	 * @param alpha first angle with major axis (in radians)
	 * @param beta second angle with major axis (in radians)
	 * @return set of bounding points
	 */
	virtual std::vector<Point> getSampleBoundingPointsOnArc(size_t num_points, double alpha, double beta) const ;

	/** \brief Returns the radius from one of the ellipse focus as a function of the angle with the main axis (standard polar definition of the ellipse)
	 *
	 * @param theta angle (in radians)
	 * @param dir direction of the focus along the major axis
	 * @return the radius
	 */
	virtual double getRadiusOnEllipseFromFocus(double theta, bool dir) const ;



	/** \brief Sample the bounding Surface.
	 * 
	 * num_points Points are placed equidistantly on the surface. The first point is placed at \f$ \theta = 0 \f$.
	 * 
	 * @param num_points points to place on the surface.
	 */
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	/** \brief return set of points sampling the bounding surface*/
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const ;
	
//	/** \brief return set of points sampling the bounding surface, given two angles (in radians)*/
//	virtual std::vector<Point> getSamplingBoundingPointsOnArc(size_t num_points, const Point & start, const Point & finish) const ;
	
	/** \brief Sample the ellipse surface.
	 * 
	 * Points are placed along concentric ellipses. The ellipses are distorded so that the smallest ellipse is, in fact, a circle. Also, the number of points on each ellipse decreases with the radius. Finally, points are added on the main axis if the ellipse is too flat.
	 * 
	 * @param num_points number of points <b>on the boundary</b>.
	 */
	virtual void sampleSurface(size_t num_points) ;
	virtual bool in(const Point &v) const ;
	
	/** \brief Return the ellipse major radius.
	 * 
	 * @return the major radius.
	 */
	virtual double getMajorRadius() const { return majorradius ; } ;

	/** \brief Return the ellipse minor radius.
	 * 
	 * @return the minor radius.
	 */
	virtual double getMinorRadius() const ;

	/** \brief Return the ellipse major axis.
	 * 
	 * @return the major axis.
	 */
	virtual Point getMajorAxis() const ;

	/** \brief Calculate the area.
	 * 
	 * Using the usual \f$ \pi a b \f$.
	 * 
	 * @return the area.
	 */
	virtual double area() const ;
	
	/** \brief return 0*/
	virtual double volume() const { return 0 ;} 

	/** \brief Changes the two ellipse radii
	 * 
	 * @param newa the new major radius
	 * @param newb the new minor radius
	 */
	virtual void setRadius(double newa, double newb);
		
	/** \brief Project the point on the ellipse.
	 * 
	 * @param  p Point to project.
	 */
	Point project(Point p) const;

	virtual void project(Point * p) const /*{return ; }*/ ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_TWO_DIMENSIONAL ;
	}
		
	const Ellipse * getGeometry() const ;
	
	Ellipse * getGeometry() ;
	
} ;

} ;

#endif
