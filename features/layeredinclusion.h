//
// C++ Interface: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LAYERED_INCLUSION_H__
#define __LAYERED_INCLUSION_H__

#include "features.h"

namespace Mu
{

class LayeredInclusion :  public LayeredCircle,  public Feature
{
protected:
	std::vector<Form *> layeredBehaviour ;
public:
	LayeredInclusion(Feature *father, std::vector<double> radii,double originX,double originY) ;
	LayeredInclusion(Feature *father,std::vector<double> radii, const Point center) ;
	LayeredInclusion(std::vector<double> radii,double originX,double originY) ;
	LayeredInclusion(std::vector<double> radii,const Point center) ;
	LayeredInclusion(double r,const Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };

	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree_3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a layered inclusion" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}

	virtual Form * getBehaviour(const Point & p) ;
// 	virtual Form * getBehaviour(const Point & p, bool final) ;
	virtual void setBehaviour(Form * b) ;
	virtual void setBehaviours(std::vector<Form *> b) ;


	
public:
	
	GEO_DERIVED_OBJECT(LayeredCircle) ;

	virtual void sample(size_t n) ;
	
} ;

class VirtualLayer :  public Circle,  public Feature
{
public:
	VirtualLayer(LayeredInclusion *father, double radius, double x, double y) ;
	VirtualLayer(LayeredInclusion *father, double radius,  Point center) ;
	virtual ~VirtualLayer() { } ;
	virtual void addSamplePoints(PointSet * po ) { };

	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree_3D * dt) ;
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a virtual layer" << std::endl ;
	}
	
	virtual Form * getBehaviour(const Point & p) ;

	virtual bool isVoid( const Point &) const {return false ;}

	virtual bool inBoundary(const Point & v) const ;
	virtual bool inBoundary(const Point *v) const ;
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;

	virtual void sample(size_t n) ;
	
	
} ;

} ;

#endif
