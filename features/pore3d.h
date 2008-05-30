//
// C++ Interface: pore3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __PORE3D_H__
#define __PORE3D_H__

#include "features.h"

namespace Mu
{

class Pore3D : public Feature, public Sphere
{
public:
	Pore3D(Feature *father, double radius, double x, double y, double z) ;
	Pore3D(Feature *father, double radius, Point center) ;
	Pore3D(double radius, double x, double y, double z) ;
	Pore3D(double radius, Point center) ;
	
	virtual void addSamplePoints(PointSet * po ) { };
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt)  ;
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void computeCenter()
	{
		return this->Sphere::computeCenter() ;
	}
	
	virtual void print() const
	{
		std::cout << "I am a pore" << std::endl ;
	}
	
	virtual bool isVoid( const Point & p) const {return squareDist(p, getCenter()) < 0.90*radius ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	/** Sample the bounding Surface
	 * 
	 * @param n number of points for the sampling.
	 */
	virtual void sample(size_t n) ;
	
	
} ;











} ;

#endif
