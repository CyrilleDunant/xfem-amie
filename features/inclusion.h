//
// C++ Interface: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCLUSION_H__
#define __INCLUSION_H__

#include "features.h"

namespace Mu
{

/** \brief Circular inclusion*/
class Inclusion :  virtual public Circle,   public Feature
{
public:
/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Inclusion(Feature *father, double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	Inclusion(Feature *father, double radius,  Point center) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Inclusion(double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	Inclusion(double radius, Point center) ;
	
/** \brief do nothing */
	virtual void addSamplePoints(PointSet * po ) { };

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
/** \brief Do nothing*/
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
/** \brief return false */
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;

	virtual void sample(size_t n) ;
	
} ;

/** \brief Triangular inclusion*/
class TriangularInclusion :  public Triangle,  public Feature
{
public:

/** \brief construct a triangular inclusion from three points
*
* @param father father feature
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularInclusion(Feature *father, const Point & a, const Point & b, const Point & c) ;

/** \brief construct a triangular inclusion from three points
*
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularInclusion(const Point & a, const Point & b, const Point & c) ;
	
/** \brief Do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void computeCenter()
	{
		return this->Triangle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a triangular inclusion" << std::endl ;
	}
	
/** \brief return false */
	virtual bool isVoid(const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Triangle) ;

	virtual void sample(size_t n) ;
	
} ;








} ;

#endif
