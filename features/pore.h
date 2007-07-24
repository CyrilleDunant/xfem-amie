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

class Pore : public Feature, public Circle
{
public:
	Pore(Feature *father, double radius, double x, double y) ;
	Pore(Feature *father, double radius, Point center) ;
	Pore(double radius, double x, double y) ;
	Pore(double radius, Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual void print() const
	{
		std::cout << "I am a pore" << std::endl ;
	}
	
	virtual bool isVoid( const Point &p) const {/*return dist(*p, center) < 0.90*sqrt(radius) ;*/ return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;
	
	/** Sample the bounding Surface
	 * 
	 * @param n number of points for the sampling.
	 */
	virtual void sample(size_t n) ;
	
	
} ;

class TriangularPore : public Feature, public Triangle
{
public:
	TriangularPore(Feature *father, const Point & a, const Point & b, const Point & c) ;
	TriangularPore(const Point & a, const Point & b, const Point & c) ;

	
	virtual void addSamplePoints(PointSet * po ) { };
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
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
