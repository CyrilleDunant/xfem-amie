//
// C++ Interface: pore
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PORE_H__
#define __PORE_H__

#include "features.h"

namespace Mu
{

/** \brief Pore class. This Feature has a VoidBehaviour by default*/
class Pore : public Feature, public Circle
{
public:

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Pore(Feature *father, double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	Pore(Feature *father, double radius, Point center) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Pore(double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	Pore(double radius, Point center) ;
	
/** \brief do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle> * dt)  ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron> * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
/** \brief Do nothing*/
	virtual Point * pointAfter(size_t i) ;
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual void print() const
	{
		std::cout << "I am a pore" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &p) const { return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;
	
	/** Sample the bounding Surface
	 * 
	 * @param n number of points for the sampling.
	 */
	virtual void sample(size_t n) ;
	
	
} ;

/** \brief Triangular pore  This Feature has a VoidBehaviour by default*/
class TriangularPore : public Feature, public Triangle
{
public:
/** \brief construct a triangular inclusion from three points
*
* @param father father feature
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularPore(Feature *father, const Point & a, const Point & b, const Point & c) ;

/** \brief construct a triangular inclusion from three points
*
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularPore(const Point & a, const Point & b, const Point & c) ;

	/** \brief Do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle> * dt)  ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron> * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void computeCenter()
	{
		return this->Triangle::computeCenter() ;
	}
	
	virtual void print() const
	{
		std::cout << "I am a pore" << std::endl ;
	}
	
	virtual bool isVoid( const Point & p) const {return false ; /*return boundary2->in(p) ;*/}
	
	
public:
	
	GEO_DERIVED_OBJECT(Triangle) ;
	
	/** Sample the bounding Surface
	 * 
	 * @param n number of points for the sampling.
	 */
	virtual void sample(size_t n) ;
	
	
} ;


} ;

#endif
