//
// C++ Interface: delaunay
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __DELAUNAY_H_
#define  __DELAUNAY_H_

#include "../geometry/geometry_2D.h"
#include "../utilities/samplingcriterion.h"
#include "../elements/elements.h"
#include "mesh.h"
#include <vector>
#include <complex>
#include <list>
#include <iostream>

namespace Mu
{

class Star ;
class SamplingCriterion ;
class DelaunayTriangle ;
class TriElement ;
class DelaunayTree ;

/*! \brief Base class of the delaunay tree. 

It defines the structure: an item has a neighbourhood, father, sons and stepsons. It also has a creator and a killer 
*/
class DelaunayTreeItem
{
	
protected:
	bool dead ; //!< Marker. is false when the item is isVertex the current triangulation no more.

	const Point * m_k ; //!< Point killer.
	const Point * m_c ; //!< Point creator.
	
public:	

	unsigned int index ;
	Mesh<DelaunayTriangle, DelaunayTreeItem> * tree ;
	DelaunayTreeItem * father ; //!< Item destroyed by the insertion of the creator point.
	DelaunayTreeItem * stepfather ; //!< Still-alive neighbour of the father.
	
	Point * first ; //!< Defining point. Is always <em>isVertex</em> the item.
	Point * second ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * third ; //!<  Defining point. Function differs if item is a triangle or point.
	
	bool isPlane  ;//!< Marker. This allows for a bit of reflectivity, cheaper than using casts.
	bool isTriangle ;//!< Marker. This allows for a bit of reflectivity, cheaper than using casts.
	bool isDeadTriangle ;
	
	bool visited ;//!< Marker. Useful not to lose ourselves isVertex the tree.
	bool erased ;
	std::valarray<unsigned int> stepson ; ;//!< neighbours created later than ourselves
	std::valarray<unsigned int> neighbour ; //!< neighbours. three for triangles, any number for planes.
	std::valarray<unsigned int> son ;//!< items created by our destruction.
	
	DelaunayTreeItem * getNeighbour(size_t i) const ; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem * getSon(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem * getStepson(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	//! Constructor, takes the father and creator point as arguments
	/*! \a father is the father. Needed for the maintenance of the tree.
		\a c is the Creator Point. It is useful when building neighbourhood relationships. Also, it allowfor the removal of elements from the tree.
	 */
	DelaunayTreeItem( Mesh<DelaunayTriangle, DelaunayTreeItem> *tree, DelaunayTreeItem * father, const Point * c) ;
	
	virtual ~DelaunayTreeItem() ;
	
	const Point * killer() const ; //!< Accessor. Returns the killer.
	const Point * creator() const ; //!< Accessor. Returns the creator.
	
	void setCreator(const Point * p) ; //!< Accessor. sets the creator.
	
	void removeNeighbour(DelaunayTreeItem * t) ; //!< Utility removes neighbour. Is saf
	void addNeighbour(DelaunayTreeItem * t) ; //!< Utility adds neighbour. Is safe.	

	virtual void kill(const Point * p) ; //!< kill and update the neighbourhood (livings do not neighbour the deads).
	virtual void erase(const Point * p) ;//!< kill and don't update the neighbourhood (do not use).
	
	bool isAlive() const ; //!< Accessor. Are we dead ?
	
	void addStepson(DelaunayTreeItem * s) ;  //!< Utility adds stepson. Is safe.
	void removeStepson(DelaunayTreeItem * s) ;  //!< Utility removes stepson. Is safe.
	
	void addSon(DelaunayTreeItem * s) ;//!< Utility adds son. Is safe.
	void removeSon(DelaunayTreeItem * s) ;//!< Utility removes son. Is safe.
	
	void setStepfather(DelaunayTreeItem * s) ;  //!< Accessor. sets the stepfather.
	
	void clearVisited() ; //!< Accessor. We are not marked visited anymore.
	
	virtual bool isVertex(const Point *p) const = 0 ; //!< Test. Is this point \a p isVertex ?
	virtual std::pair< Point*,  Point*> nearestEdge(const Point & p) const = 0;  //!< What is the nearest edge from this point \a p.
	virtual std::pair< Point*,  Point*> commonEdge(const DelaunayTreeItem * t) const  = 0; //!< What is the common edge with this item. returns a null pair if none.
	virtual bool inCircumCircle(const Point & p) const = 0 ; //!< Test. Are we isVertex conflict with the point ?
	virtual bool isNeighbour( const DelaunayTreeItem *) const = 0 ;  //!< Test. Are we a neighbour ?
	virtual void insert(std::vector<DelaunayTreeItem *> &, Point *p,  Star *s) = 0 ; //!< Insert the point isVertex the Neighbourhood given by \a s. Returns the new elements
	virtual void conflicts(std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > &,const Point *p) ; //!< Test. Recursively give all elements isVertex conflict with \a p.
	virtual void conflicts(std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > &, const Geometry *g) ;
	virtual void conflicts(std::pair<std::vector<Point *>, std::vector<DelaunayTreeItem *> > &, const Segment *s) ;
	std::vector<DelaunayTreeItem *> flatConflicts(std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > & ret, const Geometry *g) ;
	std::vector<DelaunayTreeItem *> flatConflicts(std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> >& ret,const Point *p) ;
	virtual void print() const { } ;
	
	virtual bool in( const Point & p) const  = 0;
	
	size_t numberOfCommonVertices(const DelaunayTreeItem * s) const;
	
	virtual bool isConflicting(const Geometry * g) const ;

	bool isDuplicate(const DelaunayTreeItem * t) const ;
	
} ;


//! \brief Triangle item the tree, defined by three points. 
/*!The points are also stored as a valarray of points(inherited from \c Triangle ). Those are stored clockwise but should only be used when creating the final mesh. They insure the proper orientation of the triangles.
*/
class DelaunayTriangle : virtual public TriElement, public DelaunayTreeItem
{
	std::valarray<std::pair<Point, double> > * cachedGPs ;
	std::vector<Point * > getIntegrationHints() const ;
public:
	
	GEO_DERIVED_OBJECT(Triangle) ;
	
	DelaunayTriangle( Mesh<DelaunayTriangle, DelaunayTreeItem> *tree, DelaunayTreeItem * father,   Point *p0,   Point *p1,   Point *p2,  Point * c) ;
	DelaunayTriangle() ;
	
	
	virtual ~DelaunayTriangle() ;

	virtual bool isVertex(const Point * p) const ;
	bool isVertexByID(const Point * p) const ;
	bool hasVertexByID(const std::valarray<Point *> * p) const ;
	
	virtual std::pair< Point*,  Point*> nearestEdge(const Point &p) const ;
	std::pair< Point*,  Point*> commonEdge(const DelaunayTreeItem * t) const ;
	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the triangle's CircumCircle and false if we are on or outside. 
	 */
	virtual bool inCircumCircle(const Point & p) const ;
	virtual bool isNeighbour( const DelaunayTreeItem * t) const;
	
	void insert( std::vector<DelaunayTreeItem *> &, Point *p,   Star *s) ;
	
	void print() const;
	virtual void refresh(const TriElement *) ;
	
	virtual std::vector<std::vector<Matrix> > & getElementaryMatrix() ;
	virtual Vector getForces() ;
	virtual std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix()  ;
	virtual Vector getNonLinearForces() ;
	virtual GaussPointArray getSubTriangulatedGaussPoints() const ;

	std::valarray<unsigned int> neighbourhood ;
	DelaunayTriangle * getNeighbourhood(size_t i) const ;
	bool isInNeighbourhood(const DelaunayTriangle * t) const ;

	void addNeighbourhood(DelaunayTriangle * t) ;
	void removeNeighbourhood(DelaunayTriangle *) ;
	
	virtual bool isConflicting(const Geometry * g) const ;
	
} ;


//! \brief Demi-plane isVertex the tree, defined by three points. 
/*! The two first points form the frontier segment, whereas the last is chosen <em>outside</em> the demi-plane
*/
class DelaunayDemiPlane :  public DelaunayTreeItem
{
protected:
	Point vector ; //!< Frontier vector. Precalculated for performace reasons
	double direction ;//!< test vector. Precalculated for performace reasons
public:
	
	DelaunayDemiPlane( Mesh<DelaunayTriangle, DelaunayTreeItem> * tree, DelaunayTreeItem * father,   Point  * _begin,   Point  * _end,    Point  * p,   Point * c) ;
	
	virtual ~DelaunayDemiPlane() ;
	
	virtual std::pair< Point*,  Point*> nearestEdge(const Point & p) const ;
	std::pair< Point*,  Point*> commonEdge(const DelaunayTreeItem * t) const ;
	
	/** \brief Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the demi plane, false if we are outside or on the limit.
	 */
	virtual bool inCircumCircle(const Point & p) const ;
	virtual bool isNeighbour( const DelaunayTreeItem * t) const;
	virtual bool isVertex(const Point *p) const ;
	
	void merge(DelaunayDemiPlane * p) ;//!<  Merge two planes. If the planes are found to form a partition of the universe, they are both killed. If they are found to define the same demiplane, their families are merged, and one of them is killed.
	
	//virtual void kill(Point * p) ;
	
	virtual void insert(std::vector<DelaunayTreeItem *> &,Point *p, Star *s) ;
	
	virtual bool in( const Point & p) const
	{
		return isOnTheSameSide( p, *third, *first, *second) ;
	}

	virtual bool isConflicting(const Geometry * g) const ;
	
	void print() const;
	
} ;


//! \brief Root of the tree.
/*! Neither Plane nor triangle. The constructor should be extended to provide valid starting points for all sorts of geometries.
 */
class DelaunayRoot :  public DelaunayTreeItem
{
public:
	DelaunayRoot( Mesh<DelaunayTriangle, DelaunayTreeItem> * tree, Point * p0,  Point * p1,  Point * p2) ;
	
	virtual ~DelaunayRoot() { };
	
	virtual bool isVertex(const Point *p) const ;
	
	virtual bool inCircumCircle(const Point & p) const ;

	virtual std::pair< Point*,  Point*> nearestEdge(const Point& p) const ;
	
	virtual std::pair< Point*,  Point*> commonEdge(const DelaunayTreeItem * t) const { return std::pair< Point*,  Point*>(NULL, NULL) ; } 
	
	virtual bool isNeighbour( const DelaunayTreeItem *) const { return false ; } 
	
	virtual void insert(std::vector<DelaunayTreeItem *> &, Point *p,   Star *s) ;
	
	virtual void conflicts(std::pair<std::vector<DelaunayTreeItem *>, std::vector<DelaunayTreeItem *> > &, const Point *p )  ;
	
	virtual void conflicts(std::pair<std::vector<DelaunayTriangle *>, std::vector<DelaunayTreeItem *> > &, const Geometry *g)  ;
	
	virtual void print() const ;
	
	virtual bool in( const Point & p) const
	{
		return true ;
	}
} ;


//! \brief Neighbourhood of conflicing triangles.
/*! The star knows about all the conflicting elements when inserting a point. Thus, It has all the information available to create the correct neighbourhood relationships.
*/

class Star
{
protected:
	std::vector<const Point *> edge ;
	std::vector<DelaunayTreeItem *> treeitem ;
	
public:
	Star(std::vector<DelaunayTreeItem *> *t, const Point *p) ;
	
	size_t size() ;
	
	const Point * getEdge(size_t i) const;
	
	void updateNeighbourhood() ;
} ;

/** \brief light-weight dead triangle. Serves to optimise memory usage*/
class DelaunayDeadTriangle : public DelaunayTreeItem
{
protected:
	Point center ; 
	double radius ;
	double sqradius ;

public:
	
	DelaunayDeadTriangle(DelaunayTriangle * father) ;
	
	virtual ~DelaunayDeadTriangle() ;
	
	virtual std::pair<Point*, Point*> commonEdge(const DelaunayTreeItem * t) const ;
	virtual std::pair< Point*,  Point*> nearestEdge(const Point& p) const ;

	bool inCircumCircle(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem * t) const ;
	bool isVertex(const Point *p) const ;
	bool isVertexByID(const Point *p) const ;
	
	void insert(std::vector<DelaunayTreeItem *>& , Point *p, Star *s) { };

	bool isConflicting(const Geometry * g) const ;
	
	const Point * getCircumCenter() const ;
	double getRadius() const ;
	
	bool in( const Point & p) const;

	void print() const;
	
} ;

/** \brief Mesh container. provides log n search facilities*/
class DelaunayTree : public Mesh<DelaunayTriangle, DelaunayTreeItem>
{
friend class FeatureTree ;
friend class Geometry ;
	
protected:
	size_t global_counter ;
	bool neighbourhood ;
	std::vector<Point * > additionalPoints ;
	std::vector<DelaunayTreeItem *> tree ;
	
public:
	virtual std::vector< DelaunayTriangle* > getElements() {return getTriangles() ;};
	virtual std::vector< DelaunayTriangle* > getConflictingElements(const Mu::Point* p)
	{
		std::vector< DelaunayTreeItem* > targets = conflicts(p) ;
		std::vector<DelaunayTriangle*> ret ;
		for(size_t i = 0 ; i < targets.size() ; i++)
		{
			if(targets[i]->isAlive() && targets[i]->isTriangle && !targets[i]->isDeadTriangle )
				ret.push_back(static_cast< DelaunayTriangle *> (targets[i])) ;
		}
		return ret ;
	};
	virtual std::vector< DelaunayTriangle* > getConflictingElements(const Geometry* g)
	{
		return conflicts(g) ;
	}
	virtual std::vector<DelaunayTreeItem *> & getTree() { return tree ;}
	virtual const std::vector<DelaunayTreeItem *> & getTree() const { return tree ;}
public:

	size_t & getLastNodeId() {return global_counter ;};
	const size_t & getLastNodeId() const {return global_counter ;}  ;
	std::vector<DelaunayDemiPlane *> plane ;
	
	DelaunayTree( Point * p0,  Point *p1,  Point *p2) ;
	
	virtual ~DelaunayTree() ;
	
	void insert( Point *p) ;
	void insert( Segment *s) ;
	
	/** \brief Conditionnal insertion of points.
	 * 
	 * @param p Point to insert.
	 * @param v vector of samplingCriterions.
	 * @param minScore minimum score to insert the point. The score is calculated as : 
	 *	\f$ \frac{N_{tot}}{N_{pass}}\f$
	 */
	void insertIf( Point *p, std::vector<SamplingCriterion *> v, double minScore ) ;
	
	/** \brief Get the list of triangles in comflict with a point. This method is O(log(n)), except when construction of the mesh has led to circular stepparenthoods. In that case, it is O(log(n)) in general, and O(n) 5% of the times.
	 * 
	 * @param p Point to check.
	 * @return the list of triangles in conflict with p. A triangle is in conflict if the point is tricly in the circumcircle.
	 */
	std::vector<DelaunayTreeItem *> conflicts( const Point *p) const ;
	
	/** \brief Get the list of triangles in conflict with a given geometry.
	 * 
	 * @param g test geometry.
	 * @return all the triangles for which at least one node is in the geometry.
	 */
	std::vector<DelaunayTriangle *> conflicts( const Geometry *g) const ;
	
	
	/** \brief Get the list of Points in conflict with a given segment.
	 * 
	 * @param s test segment.
	 * @return all the intersection points of triangles and the segment.
	 */
	std::vector<Point *> conflicts( const Segment *s) const ;
	
	
	/** \brief  Find the boundaries of the triangulation
	 * 
	 * @return the planes bordering the triangulated domain.
	 */
	std::vector<DelaunayDemiPlane *> * getConvexHull() ;
	
	/** \brief Return the result of the trianglulation.
	 * 
	 * @return all the living triangles resulting from the triangulation.
	 */
	std::vector<DelaunayTriangle *> getTriangles(bool buildNeighbourhood = true) ;
	
	void addSharedNodes(size_t nodes_per_side, size_t time_planes = 1, double timestep = 2, const TriElement * father = NULL) ; 
	virtual void setElementOrder(Order elemOrder) ;
	void refresh(TriElement *father) ;
	
	size_t numPoints() const;
	
	void print() const;
} ;

/** \brief functor for usage with STL containers. Compares two delaunayTreeItems*/
struct EqItems
{
	double tol ;
	DelaunayTreeItem * it ;
	EqItems(DelaunayTreeItem * i, double t) : tol(t), it(i) { } ;
	bool operator()(const DelaunayTreeItem * a) const
	{
		return squareDist2D(a->first, it->first) + squareDist2D(a->second, it->second) + squareDist2D(a->third, it->third) < tol ;
	}
} ;

std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > quad(const DelaunayTriangle * t) ;

} ;

//! Make \a t0 and \a t1 Neighbours. Safe.
void makeNeighbours( Mu::DelaunayTreeItem *t0, Mu::DelaunayTreeItem *t1 ) ;

#endif  //__DELAUNAY_H_
