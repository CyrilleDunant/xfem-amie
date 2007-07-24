//
// C++ Interface: sample3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SAMPLE3D_H__
#define __SAMPLE3D_H__

#include "features.h"

namespace Mu
{

class Sample3D :  public Hexahedron,  public Feature
{
public:
	Sample3D(Feature *father, double x, double y, double z, double originX, double originY, double originZ) ;
	Sample3D(double x, double y, double z,double originX, double originY, double originZ) ;
	
	virtual bool interacts(Feature * f) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t) const { return std::vector<Geometry *>() ;}
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt)  ;
	
	virtual void computeCenter()
	{
		return this->Hexahedron::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a sample" << std::endl ;
	}
	virtual bool isVoid( const Point & p) const {return false;}
	
public:
	
	GEO_DERIVED_OBJECT(Hexahedron) ;
	
	virtual void sample(size_t n)
	{
		this->sampleSurface(n) ;
	}
	
} ;


} ;

#endif
