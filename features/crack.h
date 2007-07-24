//
// C++ Interface: crack
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __CRACK_H__
#define __CRACK_H__

#include "features.h"

namespace Mu
{


class Crack :  public EnrichmentFeature,  public SegmentedLine
{
	typedef enum
	{
		VOID_ENRICHMENT, 
		SPLIT_ENRICHMENT,
		SINGULAR_ENRICHMENT,
		KINKS_ENRICHMENT,
		HEAD_ENRICHMENT,
		TAIL_ENRICHMENT,
		HEAD_AND_TAIL_ENRICHMENT,
		DONE_ENRICHMENT
	} EnrichmentType ;
	
	class EnrichmentData
	{
	protected:
		std::vector<size_t> iD ;
		EnrichmentType type ;
		std::map<std::pair<std::pair<int, int>, int >, bool> state ;
	public:
		EnrichmentData(std::vector<size_t> id, EnrichmentType t) ;
		EnrichmentData(size_t id, EnrichmentType t) ;
		EnrichmentData() ;
		std::vector< std::pair<Point *, Function> > getEnrichment(const DelaunayTriangle *, const std::vector<Segment>) ;
		
		std::vector<size_t> getID() const ;
		void setID(const std::vector<size_t>) ;
		
		EnrichmentType getType() const ;
		void setType(EnrichmentType)  ;
		
		bool enriched(DelaunayTriangle *)  ;
		void setEnriched(DelaunayTriangle *,bool) ;
		
		bool operator !=(const EnrichmentData & e) const { return e.getType() != type ;}
	} ;
	
	class EnrichmentMap
	{
	protected:
		std::map<size_t, EnrichmentData> props ;
	public:
		EnrichmentMap() ;
		EnrichmentData getEnrichment( size_t) ;
		
		void update(std::vector<DelaunayTriangle *> * my_triangles, Crack * myself, size_t &start) ;
		
		bool inMap(size_t) const ;
		
	} ;
	
	EnrichmentMap map ;
	
public:

	Crack(Feature *father, const std::valarray<Point *> &points, double radius) ;
	Crack(const std::valarray<Point *> &points, double radius) ;
	virtual ~Crack() ;
	
	virtual bool enrichmentTarget(DelaunayTriangle * t) ;
	virtual void enrich(size_t &,  DelaunayTree * dtree) ;
	
	virtual bool interacts(Feature * f) const ;
	virtual void snap(DelaunayTree * dtree) ;
	
	virtual bool inBoundary(const Point v) const ;
	virtual bool inBoundary(const Point *v) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(const DelaunayTree_3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	std::vector<DelaunayTriangle *> getIntersectingTriangles( DelaunayTree * dt) ;
	
	virtual void setInfluenceRadius(double r) ;
	
	virtual std::vector<Point *> getSamplingPoints() const ;
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const
	{
		std::cout << "I am a crack" << std::endl ;
	}
	
	virtual bool isVoid( const Point &) const {return false ;}
	
// 	virtual void setSingularityHints(const Point & i, const Point & s, std::vector<Point> * hints) const ;
	
public:
	GEO_DERIVED_OBJECT(SegmentedLine) ;
	
	
	/** Sample the surface. Does nothing.
	 * 
	 * This is necessary as we need to implement the interface. Of course, for a segmented line, it makes no sense.
	 * 
	 * @param n 
	 */
	virtual void sample(size_t n)
	{
		this->sampleSurface(n) ;
	}
	
	virtual void step(double dt, std::valarray<double> *, const DelaunayTree * dtree);
	
	virtual bool moved() const ;

protected:
	virtual void computeCenter() { };
	std::pair<double, double> computeJIntegralAtHead(double dt, const DelaunayTree * dtree) ;
	std::pair<double, double> computeJIntegralAtTail(double dt, const DelaunayTree * dtree) ;
	
	bool changed ;
} ;












} ;

#endif

