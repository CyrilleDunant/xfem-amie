
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __ENR_RING_H__
#define __ENR_RING_H__

#include "features.h"

namespace Mu
{

/** \brief Enrich the mesh to simulate a circular ring.
 *
 * This enrichment feature will add a soft discontinuity 
 * enrichment to the elements it crosses. If the inclusion
 * is smaller than, or fully contained by an element,
 * The enrichment will still be performed, yielding a 
 * composite field within the element.
 * 
 * This Inclusion will set no behaviour in the elements
 * thus derived classes should be used for actual simulation
 * see for example expansiveRing.h
 */
class EnrichmentRing :  public EnrichmentFeature,  public Circle
{
protected:
	bool updated ;
	Circle self ;
	std::vector<DelaunayTriangle *> cache ;
public:

/** \brief Construct the enrichment ring from a supporting feature, two radii and center coordinates
*
* @param father supporting Feature, has no effect
* @param radius external radius
* @param inRadius internal radius
* @param x center x
* @param y center y 
*/
	EnrichmentRing(Feature *father, double radius, double inradius, double x, double y) ;

/** \brief Construct the enrichment ring from two radii and center coordinates
*
* @param father supporting Feature, has no effect
* @param radius external radius
* @param inRadius internal radius
* @param x center x
* @param y center y 
*/
	EnrichmentRing(double radius, double inradius, double x, double y) ;
	virtual ~EnrichmentRing() ;
	
/** \brief return true if the argument intersects one of the features*/
	virtual bool enrichmentTarget(DelaunayTriangle * t) ;

/** \brief Enrich the elements cut with either of the two circles with a soft discontinuity*/
	virtual void enrich(size_t &,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
/** \brief return false*/
	virtual bool interacts(Feature * f, double d) const ;

/** \brief do nothing*/
	virtual void snap(DelaunayTree * dtree) ;
	
/** \brief return false*/
	virtual bool inBoundary(const Point v) const ;

/** \brief return false*/
	virtual bool inBoundary(const Point *v) const ;
	
/** \brief return the triangles in the tree intersecting with either of the circles*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D(FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
/** \brief return the triangles in the tree intersecting with either of the circles*/
	std::vector<DelaunayTriangle *> getBoundingElements2D( FeatureTree * dt) ;
	
/** \brief do nothing*/
	virtual void setInfluenceRadius(double r) ;
	
/** \brief return empty vector*/
	virtual std::vector<Point *> getSamplingPoints() const { return std::vector<Point *>(0) ; }
	
/** \brief do nothing, return NULL*/
	virtual Point * pointAfter(size_t i) { return NULL ; }
	
/** \brief return empty vector*/
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const
	{
		std::cout << "I am an enriched inclusion" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}

/** \brief set new external radius*/
	virtual void setRadius(double newR) ;

/** \brief set new internal radius*/
	virtual void setInRadius(double newR) ;
	
// 	virtual void setSingularityHints(const Point & i, const Point & s, std::vector<Point> * hints) const ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
	
	/** \brief Sample the surface. Does nothing.
	 * 
	 * @param n 
	 */
	virtual void sample(size_t n)
	{
		
	}
	
/** \brief do nothing */
	virtual void step(double dt, std::valarray<double> *, const Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree);

/** \brief return true if either radii changed*/
	virtual bool moved() const ;

/** \brief precompute and cache the elements intersected*/
	void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

	double getInRadius() const ;

protected:
	
	bool changed ;
} ;

}

#endif
