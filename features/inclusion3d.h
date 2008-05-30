//
// C++ Interface: inclusion3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCLUSION3D_H__
#define __INCLUSION3D_H__

#include "features.h"

namespace Mu
{

class Inclusion3D :  public Sphere,  public Feature
{
	friend class Sphere ;
public:
	Inclusion3D(Feature *father, double radius, double x, double y, double z) ;
	Inclusion3D(Feature *father, double radius,  Point center) ;
	Inclusion3D(double radius, double x, double y, double z) ;
	Inclusion3D(double radius, Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt)  ;
	
	virtual void computeCenter()
	{
		return this->Sphere::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	virtual void sample(size_t n) ;
	
} ;


class VirtualInclusion3D :  public Sphere,  public Feature
{
	friend class Sphere ;
public:
	VirtualInclusion3D(Feature *father, double radius, double x, double y, double z) ;
	VirtualInclusion3D(Feature *father, double radius,  Point center) ;
	VirtualInclusion3D(double radius, double x, double y, double z) ;
	VirtualInclusion3D(double radius, Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt)  ;
	
	virtual void computeCenter()
	{
		return this->Sphere::computeCenter() ;
	}
	
	virtual bool inBoundary(const Point *) const ;
	
	virtual bool inBoundary(const Point) const ;
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	virtual void sample(size_t n) ;
	
} ;

} ;

#endif
