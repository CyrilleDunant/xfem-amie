
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __ENR_INCLUSION_H__
#define __ENR_INCLUSION_H__

#include "features.h"

namespace Mu
{

/** Enrich the mesh to simulate a circular inclusion.
 * This enrichment feature will add a soft discontinuity 
 * enrichment to the elements it crosses. If the inclusion
 * is smaller than, or fully contained by an element,
 * The enrichment will still be performed, yielding a 
 * composite field within the element.
 */
class EnrichmentInclusion :  public EnrichmentFeature,  public Circle
{
protected:
	bool updated ;
	std::vector<DelaunayTriangle *> cache ;
public:

	EnrichmentInclusion(Feature *father, double radius, double x, double y) ;
	EnrichmentInclusion(double radius, double x, double y) ;
	virtual ~EnrichmentInclusion() ;
	
	virtual bool enrichmentTarget(DelaunayTriangle * t) ;
	virtual void enrich(size_t &,  DelaunayTree * dtree) ;
	
	virtual bool interacts(Feature * f) const ;
	virtual void snap(DelaunayTree * dtree) ;
	
	virtual bool inBoundary(const Point v) const ;
	virtual bool inBoundary(const Point *v) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	std::vector<DelaunayTriangle *> getIntersectingTriangles( DelaunayTree * dt) ;
	
	virtual void setInfluenceRadius(double r) ;
	
	virtual std::vector<Point *> getSamplingPoints() const { return std::vector<Point *>(0) ; }
	
	virtual Point * pointAfter(size_t i) { return NULL ; }
	
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const
	{
		std::cout << "I am an enriched inclusion" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}

	virtual void setRadius(double newR) ;
	
// 	virtual void setSingularityHints(const Point & i, const Point & s, std::vector<Point> * hints) const ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
	
	/** Sample the surface. Does nothing.
	 * 
	 * This is necessary as we need to implement the interface. Of course, for a segmented line, it makes no sense.
	 * 
	 * @param n 
	 */
	virtual void sample(size_t n)
	{
		
	}
	
	virtual void step(double dt, std::valarray<double> *, const DelaunayTree * dtree);
	
	virtual bool moved() const ;

	void update(DelaunayTree * dtree) ;

protected:
	virtual void computeCenter() { };
	
	bool changed ;
} ;

}

#endif
