//
// C++ Interface: sample
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include "features.h"

namespace Mu
{

class Sample :  public Rectangle,  public Feature
{
public:
	Sample(Feature *father, double x, double y, double originX, double originY) ;
	Sample(double x, double y, double originX, double originY) ;

	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t) const { return std::vector<Geometry *>() ;}
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
	virtual void computeCenter()
	{
		return this->Rectangle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a sample" << std::endl ;
	}
	virtual bool isVoid( const Point &p) const {return false;}
	
public:
	
	GEO_DERIVED_OBJECT(Rectangle) ;
	
	virtual void sample(size_t n)
	{
		this->sampleSurface(n) ;
	}

} ;


} ;

#endif
