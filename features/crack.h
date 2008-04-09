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

/** Branching cracks
 * Branches are all the segmented lines that form the non-overlaping structure of the crack. 
 * A branch contains no information on tips or bifurcations.
 * Tips are the singularities of the cracks. That is, all the branch ends that are not
 * bifurcations.
 * Forks are all the possible dual domain decomposition frontiers (as segmented lines) in 
 * elements containing branches. The affected element should be obtained from conflicts(center)
 *
 * The enrichment strategy is the following: first the tips are enriched. Then the forks are enriched.
 * Then, if there remain unenriched elements crossed by the branches, those are enriched.
 */
class BranchedCrack : public EnrichmentFeature,  public SegmentedLine
{
protected:
	std::vector<SegmentedLine *> branches ;
	std::vector<Point * > tips ;
	std::vector<SegmentedLine * > forks ;
	
	void enrichTips(size_t &, DelaunayTree * dt) ;
	void enrichTip(size_t &, DelaunayTree * dt, const Point * tip, double angle) ;

	void enrichForks(size_t &, DelaunayTree * dt) ;
	void enrichBranches(size_t &, DelaunayTree * dt) ;
	void enrichSegmentedLine(size_t &, DelaunayTree * dt, const SegmentedLine * line) ;
	
	typedef enum
	{
		NOT_ENRICHED,
		TIP_ENRICHED,
		BIFURCATION_ENRICHED,
		LINE_ENRICHED
	} EnrichementState;

	std::map<std::pair<DelaunayTriangle *, Point *>, EnrichementState > enrichmentMap ;
	
	double enrichementRadius ;
	bool changed ;

	std::set<Point *> donePoints ;
	
public:
	
	/** Instantiate a branching crack.
	 *
	 * The crack is instantiated as a segment, given two points.
	 * 
	 * @param father father feature. This is used to determin whether the points given are tips or lie outside of the domain.
	 * @param a first point
	 * @param b second point
	 */
	BranchedCrack(Feature *father, Point * a, Point * b) ;

	/** Instantiate a branching crack.
	 *
	 * The crack is instantiated as a segment, given two points. There is no check of wheter the points lie within a given domain
	 * therefore they are assumed to be tips.
	 *
	 * @param a first point
	 * @param b second point
	 */
	BranchedCrack(Point * a, Point * b) ;

	double getEnrichementRadius() const ;
	const std::vector<SegmentedLine *> & getBranches() const;
	const std::vector<Point * > & getTips() const;
	const std::vector<SegmentedLine * > & getForks() const;

	std::vector<SegmentedLine *> & getBranches();
	std::vector<Point * > & getTips();
	std::vector<SegmentedLine * > & getForks();
	void setEnrichementRadius(double newRadius) ;

	/** Branch the crack from the given tip, provided the two new tips.
	 * The branching operation will create a new branch, expand the original branche and add a fork.
	 *
	 * @param fromTip the original tip, which will be removed
	 * @param newTip0 first new tip, used to extend the affected branch
	 * @param newTip1 second new tip. forms the new branch with the original tip
	 */
	void branch(Point* fromTip, Point * newTip0, Point * newTip1 ) ;

	/** Merge two branchedCracks.
	 * the merging will locate the nearest (segment, tip) couple from the (original, argument) cracks and create a new intersection point.
	 * 
	 * the affected branch will have an additional point at the intersection. A fork will be added, and all the branches, forks and tips from the target will be merged into this crack. The crack given as an argument will be emptied.
	 * 
	 * @param newSet crack to merge
	 */
	void merge(BranchedCrack & newSet) ;

	/** Is this crack still existing.
	 * This is useful to know if the crack has already been merged.
	 * 
	 * @return false if there is at least a branch
	 */
	bool isEmpty() const ;

	virtual void enrich(size_t &,  DelaunayTree * dtree) ;

	virtual void computeCenter() ;
	virtual std::vector<DelaunayTriangle*> getTriangles(DelaunayTree*) ;
	virtual std::vector<DelaunayTetrahedron*> getTetrahedrons(DelaunayTree_3D*) ;
	virtual bool interacts(Feature*) const ;
	virtual Point* pointAfter(size_t) ;
	virtual std::vector<Mu::Geometry*> getRefinementZones(size_t) const ;
	virtual void print() const ;
	//	virtual void printFile(const std::string& filename) const;//SB
	virtual void sample(size_t) ;
	virtual bool isVoid(const Point&) const ;
	virtual std::vector<Point*> getSamplingPoints() const ;
	virtual bool enrichmentTarget(DelaunayTriangle*) ;
	virtual void step(double, Vector*, const DelaunayTree*) ;
	virtual void snap(DelaunayTree*) ;
	virtual bool moved() const ;

	virtual ~BranchedCrack() ;
		
	
public:
	GEO_DERIVED_OBJECT(SegmentedLine) ;
} ;

class Crack :  public EnrichmentFeature,  public SegmentedLine
{
	double stepLength ;
	double criticalJ ;
	typedef enum
	{
		VOID_ENRICHMENT, 
		SPLIT_ENRICHMENT,
		MULTI_SPLIT_ENRICHEMENT,
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
	// SB
	Crack ( Feature * father, const std::valarray<Point *> & points, double radius, double gfa, double critJ );
	Crack ( const std::valarray<Point *> & points, double radius, double gfac, double critJ );
	//end SB
	virtual ~Crack() ;
	
	virtual bool enrichmentTarget(DelaunayTriangle * t) ;
	virtual void enrich(size_t &,  DelaunayTree * dtree) ;
	
	virtual bool interacts(Feature * f) const ;
	virtual void snap(DelaunayTree * dtree) ;
	
	virtual bool inBoundary(const Point & v) const ;
	virtual bool inBoundary(const Point *v) const ;
	
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree_3D * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	void setCriticalJ(double newJ) ;

	double getCriticalJ() const ;

	std::vector<DelaunayTriangle *> getIntersectingTriangles( DelaunayTree * dt) ;
	
	virtual void setInfluenceRadius(double r) ;
	// SB	
	virtual void setParams ( double r, double gfac, double critJ);
	  // end SB
	virtual std::vector<Point *> getSamplingPoints() const ;
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual std::vector<Geometry *> getRefinementZones( size_t level) const ;
	
	virtual void print() const ;
	virtual void printFile(const std::string& filename) const;//SB
	
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

