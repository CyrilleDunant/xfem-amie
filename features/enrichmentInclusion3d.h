
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __ENR_INCLUSION_3D_H__
#define __ENR_INCLUSION_3D_H__

#include "features.h"

namespace Mu
{

/** \brief Enrich the mesh to simulate a circular inclusion.
 *
 * This enrichment feature will add a soft discontinuity 
 * enrichment to the elements it crosses. If the inclusion
 * is smaller than, or fully contained by an element,
 * The enrichment will still be performed, yielding a 
 * composite field within the element.
 * 
 * This Inclusion will set no behaviour in the elements
 * thus derived classes should be used for actual simulation
 * see for example expansiveZone.h
 */
class EnrichmentInclusion3D :  public EnrichmentFeature,  public Sphere
{
protected:
	bool updated ;
	std::vector<DelaunayTetrahedron *> cache ;
public:

/** \brief Constructor. Construct the inclusion from a supporting feature, a radius and center coordinates */
	EnrichmentInclusion3D(Feature *father, double radius, double x, double y, double z) ;

/** \brief Constructor. Construct the inclusion from a radius and center coordinates */
	EnrichmentInclusion3D(double radius, double x, double y, double z) ;
	virtual ~EnrichmentInclusion3D() ;
	
/** \brief return true if the argument is cut by the sphere*/
	virtual bool enrichmentTarget(DelaunayTetrahedron * t) ;

/** \brief Enrich elements cut by the feature*/
	virtual void enrich(size_t& counter, Mu::Mesh< Mu::DelaunayTetrahedron, Mu::DelaunayTreeItem3D >* dtree) ;
	
/** \brief return false*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return false*/
	virtual bool inBoundary(const Point & v) const ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt) {return std::vector<DelaunayTriangle *>(0) ;} 

/** \brief return the list of elements cut by the feature*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D(FeatureTree * dt)  ;
	
	
/** \brief return list of elements cut by the feature*/
	std::vector<DelaunayTetrahedron *> getBoundingElements3D( FeatureTree *  dt) ;

/** \brief do nothing*/
	virtual void setInfluenceRadius(double r) ;
	
/** \brief return empty vector*/
	virtual std::vector<Point *> getSamplingPoints() const { return std::vector<Point *>(0) ; }
	
/** \brief do nothing, return null*/
	virtual Point * pointAfter(size_t i) { return NULL ; }
	
/** \brief return empty vector*/
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const
	{
		std::cout << "I am an enriched inclusion 3D" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}

/** \brief set the radius*/
	virtual void setRadius(double newR) ;
	
// 	virtual void setSingularityHints(const Point & i, const Point & s, std::vector<Point> * hints) const ;
	
public:
	GEO_DERIVED_OBJECT(Sphere) ;
	
	
	/** \brief Sample the surface. Does nothing.
	 * 
	 * This is necessary as we need to implement the interface. Of course, for a segmented line, it makes no sense.
	 * 
	 * @param n 
	 */
	virtual void sample(size_t n)
	{
		
	}
	
/** \brief do nothing*/
	virtual void step(double dt, std::valarray<double> *, const DelaunayTree * dtree);
	
/** \brief return true if the radius has changed*/
	virtual bool moved() const ;

/** \brief compute and cache the elements to enrich*/
	void update(Mu::Mesh< Mu::DelaunayTetrahedron, Mu::DelaunayTreeItem3D >* dtree) ;

protected:
	bool changed ;
} ;

}

#endif
