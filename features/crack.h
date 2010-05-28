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



/** \brief Branching cracks
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
	std::vector<std::pair<Point *, double> > tips ;
	std::vector<SegmentedLine * > forks ;
	std::set<size_t> freeIds ;
	
	void enrichTips(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt) ;
	void enrichTip(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt, const std::pair<Point *, double> & tip) ;

	void enrichForks(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt) ;
	void enrichBranches(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt) ;
	void enrichSegmentedLine(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt, const SegmentedLine * line) ;
	void enrichSegmentedLine(size_t &, Mesh<DelaunayTriangle, DelaunayTreeItem> * dt, const SegmentedLine * line, const Point *) ;
	
	double influenceRadius ;
	
	typedef enum
	{
		NOT_ENRICHED,
		TIP_ENRICHED,
		BIFURCATION_ENRICHED,
		LINE_ENRICHED
	} EnrichementState;

	std::set<DelaunayTriangle *> enrichmentMap ;
	std::set<DelaunayTriangle *> forkEnrichmentMap ;
	std::set<DelaunayTriangle *> tipEnrichmentMap ;
	
	double enrichementRadius ;
	bool changed ;
	double propagationAngleFromTip(const std::pair< Mu::Point*, double >& tip, Mu::Mesh< Mu::DelaunayTriangle, Mu::DelaunayTreeItem >* dtree) ;
	std::pair<double, double> computeJIntegralAtTip ( std::pair< Mu::Point*, double >& tip, Mu::Mesh< Mu::DelaunayTriangle, Mu::DelaunayTreeItem >* dtree ) ;
	std::set<Point *> donePoints ;
	
public:
	
	/** \brief Instantiate a branching crack.
	 *
	 * The crack is instantiated as a segment, given two points.
	 * 
	 * @param father father feature. This is used to determin whether the points given are tips or lie outside of the domain.
	 * @param a first point
	 * @param b second point
	 */
	BranchedCrack(Feature *father, Point * a, Point * b) ;

	/** \brief Instantiate a branching crack.
	 *
	 * The crack is instantiated as a segment, given two points. There is no check of wheter the points lie within a given domain
	 * therefore they are assumed to be tips.
	 *
	 * @param a first point
	 * @param b second point
	 */
	BranchedCrack(Point * a, Point * b) ;

/** \brief return the radius in which elements are enriched with tip enrichment around each tip*/
	double getEnrichementRadius() const ;

/** \brief return the vector of SegmentedLine representing the geometry of the crack*/
	const std::vector<SegmentedLine *> & getBranches() const;

/** \brief return the tips and their angle*/
	const std::vector<std::pair<Point *, double> > & getTips() const;

/** \brief return the set of SegmentedLine s representing the alternative paths to describe the crack geometry.*/
	const std::vector<SegmentedLine * > & getForks() const;

/** \brief return the vector of SegmentedLine representing the geometry of the crack*/
	std::vector<SegmentedLine *> & getBranches();

/** \brief return the tips and their angle*/
	std::vector<std::pair<Point *, double> > & getTips();

/** \brief return the vector of SegmentedLine representing the geometry of the crack*/
	std::vector<SegmentedLine * > & getForks();
	void setEnrichementRadius(double newRadius) ;

	/** \brief Branch the crack from the given tip, provided the two new tips.
	 * The branching operation will create a new branch, expand the original branche and add a fork.
	 *
	 * @param fromTip the original tip, which will be removed
	 * @param newTip0 first new tip, used to extend the affected branch
	 * @param newTip1 second new tip. forms the new branch with the original tip
	 */
	void branch(Point* fromTip, Point * newTip0, Point * newTip1 ) ;
	
	/** \brief Branch the crack from the given tip, provided the two new tips.
	 * The branching operation will create a new branch, expand the original branche and add n forks.
	 *
	 * @param fromTip the original tip, which will be removed
	 * @param newTips first new tip, used to extend the affected branch
	 */
	void branch(Point* fromTip, const std::vector<Point *> & newTip) ;

	/** \brief Merge two branchedCracks.
	 * the merging will locate the nearest (segment, tip) couple from the (original, argument) cracks and create a new intersection point.
	 * 
	 * the affected branch will have an additional point at the intersection. A fork will be added, and all the branches, forks and tips from the target will be merged into this crack. The crack given as an argument will be emptied.
	 * 
	 * @param newSet crack to merge
	 */
	void merge(BranchedCrack & newSet) ;

	/** \brief Grow the crack from the given tip, provided the new tip.
	 * The growing operation will expand the original branch.
	 *
	 * @param fromTip the original tip, which will be removed
	 * @param newTip first new tip, used to extend the affected branch
	 */
	void grow( Point* fromTip, Point* newTip) ;
	
	/** \brief Is this crack still existing.
	 * This is useful to know if the crack has already been merged.
	 * 
	 * @return false if there is at least a branch
	 */
	bool isEmpty() const ;

/** \brief Enrich the elements contained in the argument which interact with the crack*/
	virtual void enrich(size_t &,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;

	virtual void computeCenter() ;

/** \brief return the list of elements interacting with the crack*/
	virtual std::vector<DelaunayTriangle*> getElements(Mesh<DelaunayTriangle, DelaunayTreeItem>*) ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron*> getElements(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>*) ;

/** \brief Return fales*/
	virtual bool interacts(Feature*, double) const ;

/** \brief insert a point in the main branch after its ith point.*/
	virtual Point* pointAfter(size_t i) ;

/** \brief Return the vector of geometries in which the mesh should be refined around this Feature*/
	virtual std::vector<Mu::Geometry*> getRefinementZones(size_t) const ;
	virtual void print() const ;
	//	virtual void printFile(const std::string& filename) const;//SB

/** \brief do nothing*/
	virtual void sample(size_t) ;

/** \brief return false*/
	virtual bool isVoid(const Point&) const ;

/** \brief return an empty vector*/
	virtual std::vector<Point*> getSamplingPoints() const ;

/** \brief return true if this element should be enriched by this crack*/
	virtual bool enrichmentTarget(DelaunayTriangle*) ;

/** \brief Grow this crack
* 
* @param dt timestep
* @param v ignore
* @param tree Element tree
*
*/
	virtual void step(double dt , Vector* v, Mesh<DelaunayTriangle, DelaunayTreeItem>* tree) ;

/** \brief Do nothing*/
	virtual void snap(Mesh<DelaunayTriangle, DelaunayTreeItem>*) ;

/** \brief if the crack grew this simulation step, return true*/
	virtual bool moved() const ;

	virtual ~BranchedCrack() ;
		
	
public:
	GEO_DERIVED_OBJECT(SegmentedLine) ;
} ;


} ;

bool operator ==(const std::pair<Mu::Point*, double> & a, const Mu::Point* b) ;


#endif

