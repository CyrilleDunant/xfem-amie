
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


#ifndef __DELAUNAY_3D_H_
#define  __DELAUNAY_3D_H_
#include "geometry/geometry_base.h"
#include "geometry/geometry_3D.h"
#include "elements/elements.h"
#include <vector>
#include <complex>
#include <list>
#include <set>
#include <iostream>

namespace Mu
{

class Star_3D ;
class DelaunayTetrahedron ;
class DelaunayTree_3D ;

/*! Base class of the delaunay tree. 

It defines the structure: an item has a neighbourhood, father, sons and stepsons. It also has a creator and a killer 
*/
class DelaunayTreeItem_3D
{
	
protected:
	bool dead ; //!< Marker. is false when the item is isVertex the current triangulation no more.
	
	
	const Point * m_k ; //!< Point killer.
	const Point * m_c ; //!< Point creator.
	
public:
	unsigned int index ;
	bool erased ;
	DelaunayTree_3D * tree ;
	DelaunayTreeItem_3D * father ; //!< Item destroyed by the insertion of the creator point.
	DelaunayTreeItem_3D * stepfather ; //!< Still-alive neighbour of the father.
	
	Point * first ; //!< Defining point. Is always <em>isVertex</em> the item.
	Point * second ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * third ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * fourth ; //!<  Defining point. Function differs if item is a triangle or point.
	
	bool isSpace ;  //newly added
	bool isTetrahedron ; // newly added
	bool isDeadTetrahedron ;
	bool visited ;//!< Marker. Useful not to lose ourselves isVertex the tree.
	
	std::valarray<unsigned int> stepson ; ;//!< neighbours created later than ourselves
	std::valarray<unsigned int> neighbour ; //!< neighbours. three for triangles, any number for planes.
	std::valarray<unsigned int> son ;//!< items created by our destruction.
	
	//! Constructor, takes the father and creator point as arguments
	/*! \a father is the father. Needed for the maintenance of the tree.
		\a c is the Creator Point. It is useful when building neighbourhood relationships. Also, it allowfor the removal of elements from the tree.
	 */
	DelaunayTreeItem_3D( DelaunayTree_3D *tree, DelaunayTreeItem_3D * father, const Point * c) ;
	
	virtual ~DelaunayTreeItem_3D() ;
	
	const Point * killer() const ; //!< Accessor. Returns the killer.
	const Point * creator() const ; //!< Accessor. Returns the creator.
	
	void setCreator(const Point * p) ; //!< Accessor. sets the creator.
	
	void removeNeighbour(DelaunayTreeItem_3D * t) ; //!< Utility removes neighbour. Is safe.
	
	void addNeighbour(DelaunayTreeItem_3D * t) ; //!< Utility adds neighbour. Is safe.
	
	DelaunayTreeItem_3D * getNeighbour(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem_3D * getSon(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem_3D * getStepson(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	
	virtual void kill(const Point * p) ; //!< kill and update the neighbourhood (livings do not neighbour the deads).
	virtual void erase(const Point * p) ;//!< kill and don't update the neighbourhood (do not use).
	
	bool isAlive() const ; //!< Accessor. Are we dead ?
	
	void addStepson(DelaunayTreeItem_3D * s) ;  //!< Utility adds stepson. Is safe.
	void removeStepson(DelaunayTreeItem_3D * s) ;  //!< Utility removes stepson. Is safe.
	
	void addSon(DelaunayTreeItem_3D * s) ;//!< Utility adds son. Is safe.
	void removeSon(DelaunayTreeItem_3D * s) ;//!< Utility removes son. Is safe.
		
	void setStepfather(DelaunayTreeItem_3D * s) ;  //!< Accessor. sets the stepfather.
	void setFather(DelaunayTreeItem_3D * s) ;
	
	void clearVisited() ; //!< Accessor. We are not marked visited anymore.
	
	virtual bool isVertex(const Point *p) const = 0 ; //!< Test. Is this point \a p isVertex ?

	virtual std::vector< Point*> commonSurface(const DelaunayTreeItem_3D *t) const  = 0; //!< What is the common edge with this item. returns a null pair if none
	
	virtual bool inCircumSphere(const Point &p) const = 0 ; //!< Test. Are we isVertex conflict with the point ?
	//newly aded
	virtual bool isNeighbour( const DelaunayTreeItem_3D *) const = 0 ;  //!< Test. Are we a neighbour ?
	virtual void insert(std::vector<DelaunayTreeItem_3D *> &, Point *p,  Star_3D *s) = 0 ; //!< Insert the point isVertex the Neighbourhood given by \a s. Returns the new elements
	virtual  void conflicts(std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > &,const Point *p) ; //!< Test. Recursively give all elements isVertex conflict with \a p.
	virtual  void conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > &, const Geometry *g) ;
	virtual void print() const = 0 ;
	
	virtual bool in( const Point & p) const  = 0;
	
	size_t numberOfCommonVertices(const DelaunayTreeItem_3D * s) const;
	
	bool isDuplicate(const DelaunayTreeItem_3D * t) const ;
	
} ;


//! Tetrahedron item the tree, defined by four points. 
/*!The points are also stored as a valarray of points(inherited from \c Triangle ). Those are stored clockwise but should only be used when creating the final mesh. They insure the proper orientation of the triangles.
*/
class DelaunayTetrahedron : virtual public TetrahedralElement, public DelaunayTreeItem_3D
{
public:
	
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	DelaunayTetrahedron(DelaunayTree_3D *tree,  DelaunayTreeItem_3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,  Point * c) ;
	DelaunayTetrahedron(DelaunayTree_3D *tree,  DelaunayTreeItem_3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,   
	                     Point *p4,   Point *p5,   Point *p6, Point *p7,
	                     Point * c) ;
	DelaunayTetrahedron() ;
	
	
	virtual ~DelaunayTetrahedron() ;

	bool isVertex(const Point * p) const ;
	bool isVertexByID(const Point * p) const ;
	bool hasVertexByID(const std::valarray<Point *> * p) const ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem_3D * t) const ;

	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the triangle's CircumCircle and false if we are on or outside. 
	 */
	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem_3D * t) const ;
	void kill(const Point * p) ;
	void insert(std::vector<DelaunayTreeItem_3D *> &, Point *p,   Star_3D *s) ;
	
	void print() const;
	
	//virtual std::vector<DelaunayTetrahedron *> withCommonPoint(Point * p) ;
	//	virtual std::vector<DelaunayTetrahedron *> withCommonPoints() ;
	void refresh(const TetrahedralElement *) ;
	
	std::vector<std::vector<Matrix> > getElementaryMatrix() const ;
	std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix()  ;
	Vector getForces() const ;
	Vector getNonLinearForces() const ;
	
	GaussPointArray getSubTriangulatedGaussPoints() const ;
	
	std::valarray<unsigned int> neighbourhood ;
	DelaunayTetrahedron * getNeighbourhood(size_t i) const ;
	void addNeighbourhood(DelaunayTetrahedron * t) ;
	void removeNeighbourhood(DelaunayTetrahedron *) ;
	bool isInNeighbourhood(const DelaunayTetrahedron * t) const ;
		
	
} ;
//! Demi-plane isVertex the tree, defined by three points. 
/*! The two first points form the frontier segment, whereas the last is chosen <em>outside</em> the demi-plane
*/

//newly added class
class DelaunayDemiSpace : public DelaunayTreeItem_3D
{
protected:
	Point vector1 ; //!< Frontier vector. Precalculated for performance reasons
	Point vector2 ; //!< Frontier vector. Precalculated for performance reason
	double direction ;//!< test vector. Precalculated for performance reasons
public:
	Point pseudonormal ;
	
public:
	
	DelaunayDemiSpace(DelaunayTree_3D *tree,  DelaunayTreeItem_3D * father,   Point  * _one,   Point  * _two, Point  * _three,   Point  * p,   Point * c) ;
	
	virtual ~DelaunayDemiSpace() ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem_3D * t) const ;	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the demi plane, false if we are outside or on the limit.
	 */
	
	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem_3D * t) const ;
	bool isVertex(const Point *p) const ;
	
	void merge(DelaunayDemiSpace * p) ;//!< ????
	
	//virtual void kill(Point * p) ;
	
	void insert(std::vector<DelaunayTreeItem_3D *>& , Point *p, Star_3D *s) ;
	
	bool in( const Point & p) const
	{
		return isOnTheSameSide( &p, third, first, second, fourth) ;
	}

	void print() const;
	
} ;


class DelaunayDeadTetrahedron : public DelaunayTreeItem_3D
{
protected:
	Point center ; 
	double radius ;
	double sqradius ;

public:
	
	DelaunayDeadTetrahedron(DelaunayTetrahedron * father) ;
	
	virtual ~DelaunayDeadTetrahedron() ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem_3D * t) const ;	

	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem_3D * t) const ;
	bool isVertex(const Point *p) const ;
	bool isVertexByID(const Point *p) const ;
	
	void insert(std::vector<DelaunayTreeItem_3D *>& , Point *p, Star_3D *s) { };
	
	const Point * getCircumCenter() const ;
	double getRadius() const ;
	
	bool in( const Point & p) const;

	void print() const;
	
} ;

//! Root of the tree.
/*! Neither Plane nor triangle. The constructor should be extended to provide valid starting points for all sorts of geometries.
 */
class DelaunayRoot_3D : public DelaunayTreeItem_3D
{
public:
	DelaunayRoot_3D(DelaunayTree_3D *tree, Point * p0,  Point * p1,  Point * p2, Point * p3) ;
	
	bool isVertex(const Point *p) const ;
	
	bool inCircumSphere(const Point & p) const ;
			
	std::vector< Point*> commonSurface(const DelaunayTreeItem_3D *t) const ; //!< What is the common edge with this item. returns a null pair if none
	
	bool isNeighbour( const DelaunayTreeItem_3D *) const { return false ; } 
	
	void insert(std::vector<DelaunayTreeItem_3D *> &, Point *p,   Star_3D *s) ;
	
	void conflicts(std::pair<std::vector<DelaunayTreeItem_3D *>, std::vector<DelaunayTreeItem_3D *> > &, const Point *p )  ;
	
	void conflicts( std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem_3D *> > &, const Geometry *g)  ;
	
	void print() const ;
	
	bool in( const Point & p) const
	{
		return true ;
	}
} ;


//! Neighbourhood of conflicing triangles.

/*! The star knows about all the conflicting elements when inserting a point. Thus, It has all the information available to create the correct neighbourhood relationships.
*/

class Star_3D
{
protected:
	std::vector<const Point *> edge ;
	std::vector<DelaunayTreeItem_3D *> treeitem ;
	
public:
	std::vector<DelaunayTreeItem_3D *> cleanup ;
	~Star_3D() ;
	Star_3D(std::vector<DelaunayTreeItem_3D *> *t, const Point *p) ;
	
	size_t size() ;
	
	const Point * getEdge(size_t i) const;
	
	void updateNeighbourhood() ;
	
	const Point *creator ;
} ;
//changed

class DelaunayTree_3D 
{
friend class FeatureTree ;
friend class Geometry ;
	
protected:
	size_t global_counter ;
	bool neighbourhood ;
	
public:

	std::vector<DelaunayTreeItem_3D *> tree ;
	std::vector<DelaunayDemiSpace *> space ;
	
	DelaunayTree_3D ( Point * p0,  Point *p1,  Point *p2, Point *p3) ;
	
	virtual ~DelaunayTree_3D() ;
	
	void insert( Point *p) ;
	void insert( Segment *s) ;
	
	/** Conditionnal insertion of points.
	 * 
	 * @param p Point to insert.
	 * @param v vector of samplingCriterions.
	 * @param minScore minimum score to insert the point. The score is calculated as : 
	 *	\f$ \frac{N_{tot}}{N_{pass}}\f$
	 */
// 	void insertIf( Point *p, std::vector<SamplingCriterion_3D *> v, double minScore ) ;
	
	/** Get the list of triangles in comflict with a point. This method is O(log(n)), except when construction of the mesh has led to circular stepparenthoods. In that case, it is O(log(n)) in general, and O(n) 5% of the times.
	 * 
	 * @param p Point to check.
	 * @return the list of triangles in conflict with p. A triangle is in conflict if the point is tricly in the circumcircle.
	 */
	std::vector<DelaunayTreeItem_3D *> conflicts( const Point *p) const ;
	
	/** Get the list of triangles in conflict with a given geometry.
	 * 
	 * @param g test geometry.
	 * @return all the triangles for which at least one node is in the geometry.
	 */
	std::vector<DelaunayTetrahedron *> conflicts( const Geometry *g) const ;

	
	/** Find the boundaries of the triangulation
	 * 
	 * @return the planes bordering the triangulated domain.
	 */
	std::vector<DelaunayDemiSpace *> getConvexHull() ;
	
	/** Return the result of the trianglulation.
	 * 
	 * @return all the living triangles resulting from the triangulation.
	 */
	std::vector<DelaunayTetrahedron *> getTetrahedrons(bool buildNeighbourhood = true) ;
	
	void addSharedNodes(size_t nodes_per_side, const TetrahedralElement * father = NULL) ; 
	void addSharedNodes(size_t nodes_per_side, size_t time_planes = 2, double timestep = 2, const TetrahedralElement * father = NULL) ;
	
	void refresh(const TetrahedralElement *father) ;
	
	size_t numPoints() const;
	
	void print() const;
} ;

} ;

//! Make \a t0 and \a t1 Neighbours. Safe.
void makeNeighbours( Mu::DelaunayTreeItem_3D *t0, Mu::DelaunayTreeItem_3D *t1 ) ;

#endif  //__DELAUNAY_3D_H_
