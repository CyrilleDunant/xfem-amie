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

class Tetrahedron : public ConvexGeometry
{

protected:
	double radius ;
	Point circumCenter ;
	void computeCircumCenter() ;
	
	
public:
	
	Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3) ;
	
	Tetrahedron(Point p0, Point p1, Point p2, Point p3) ; 
	
	Tetrahedron() ; 
	
	virtual ~Tetrahedron() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point &v) const ;
	
	virtual double area() const ;
	virtual double volume() const ;
		
	virtual void project(Point * p) const;
	
	virtual bool inCircumSphere(const Point p) const ;
	virtual bool inCircumSphere(const Point *p) const ;
	
	virtual void computeCenter() ;
	
	virtual const Point * getCircumCenter() const;
	
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_THREE_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const { return std::vector<Point>(0) ;}
	
} ;

class Hexahedron : public ConvexGeometry
{
protected:
	
	double size_y ;
	double size_x ;
	double size_z ;
	
public:
	
	Hexahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7) ;
	
	Hexahedron(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5, Point p6, Point p7) ; 
	
	Hexahedron(double x, double y, double z, double originX, double originY, double originZ) ;
	
	Hexahedron() ; 
	
	virtual ~Hexahedron() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point &v) const ;
	
	virtual double area() const ;
	virtual double volume() const ;
		
	virtual void project(Point * p) const;
	
	virtual void computeCenter() ;
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

struct TriPoint
{
	Point normal ;

	std::valarray<Point *> point ;
	std::vector<TriPoint *> neighbour ;
	
	TriPoint(Point * p0, Point * p1, Point * p2) : point(3) 
	{
		point[0] = p0 ;
		point[1] = p1 ;
		point[2] = p2 ;
		normal = (*p0-*p1)^(*p2-*p1) ;
	}
	
	double area() const
	{
		return .5*((*point[0]-*point[1])^(*point[2]-*point[1])).norm() ;
	}
	
};


class TriangulatedSurface : public NonConvexGeometry
{
protected:
	std::vector<Point *> boundary ;
	std::vector<TriPoint> mesh ;
	
public:
	void addPoint(Point *) ;
	
public:
	
	TriangulatedSurface(Point * p0, Point * p1, Point * p2) ;
	
	TriangulatedSurface() ;
	
	virtual ~TriangulatedSurface() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) { };
	
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


class Sphere : public ConvexGeometry
{
protected:
	
	double radius ;
	
	std::vector<Point> getSamplingPointsOnSphere(size_t num_points, double radius) const ;
	virtual void project(Point * p, double r) const;
	void smooth(std::vector<Point> & points, double r) const ;
	
public:
	
	
	Sphere(double radius,double x, double y, double z ) ;
	
	Sphere(double radius,const Point * p0 ) ;
	
	Sphere(double radius,const Point p0 ) ; 
	
	Sphere() ; 
	
	virtual ~Sphere() { } ;
	
	virtual void sampleBoundingSurface(size_t num_points) ;
	
	virtual void sampleSurface(size_t num_points);
	virtual bool in(const Point &v) const ;
	
	virtual double area() const ;
	virtual double volume() const ;
		
	virtual void project(Point * p) const;
	
	virtual void computeCenter() ;
	virtual double getRadius() const ;
	
	virtual SpaceDimensionality spaceDimensions() const
	{
		return SPACE_THREE_DIMENSIONAL ;
	}
	
	virtual std::vector<Point> getBoundingBox() const ;

} ;

} ;
#endif
