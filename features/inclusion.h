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

class Inclusion :  public Circle,  public Feature
{
public:
	Inclusion(Feature *father, double radius, double x, double y) ;
	Inclusion(Feature *father, double radius,  Point center) ;
	Inclusion(double radius, double x, double y) ;
	Inclusion(double radius, Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };

	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;

	virtual void sample(size_t n) ;
	
} ;


class TriangularInclusion :  public Triangle,  public Feature
{
public:
	TriangularInclusion(Feature *father, const Point & a, const Point & b, const Point & c) ;
	TriangularInclusion(const Point & a, const Point & b, const Point & c) ;
	
	virtual void addSamplePoints(PointSet * po ) { };

	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void computeCenter()
	{
		return this->Triangle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a triangular inclusion" << std::endl ;
	}
	
	virtual bool isVoid(const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Triangle) ;

	virtual void sample(size_t n) ;
	
} ;








} ;

#endif
