
//
// C++ Interface: delaunay
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef __DELAUNAY_3D_H_
#define  __DELAUNAY_3D_H_
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_3D.h"
#include "../elements/elements.h"
#include <vector>
#include <complex>
#include <list>
#include <set>
#include <iostream>
#include <bitset>
#include "mesh.h"

namespace Mu
{

class Star3D ;
class DelaunayTetrahedron ;
class DelaunayTree3D ;

/*! \brief Base class of the delaunay tree. 

It defines the structure: an item has a neighbourhood, father, sons and stepsons. It also has a creator and a killer 
*/
class DelaunayTreeItem3D
{
	
protected:
	bool m_dead ;
	bool m_isSpace ;
	bool m_isTetrahedron ;
	bool m_isDeadTetrahedron ;
	bool m_visited ;
	bool m_erased ;
// 	std::bitset<6> state ;
// 	std::bitset<6>::reference dead() ;
	bool dead() const ;
	bool & dead() ;
public:
	const Point * killer ; //!< Point killer.
	const Point * creator ; //!< Point creator.
	int index ;
	
	Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * tree ;
	int father ; //!< Item destroyed by the insertion of the creator point.
	int stepfather ; //!< Still-alive neighbour of the father.
	
	Point * first ; //!< Defining point. Is always <em>isVertex</em> the item.
	Point * second ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * third ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * fourth ; //!<  Defining point. Function differs if item is a triangle or point.
	
	
// 	std::bitset<6>::reference isSpace() ;
// 	std::bitset<6>::reference isTetrahedron() ;
// 	std::bitset<6>::reference isDeadTetrahedron() ;
//	std::bitset<6>::reference visited() ;//!< Marker. Useful not to lose ourselves isVertex the tree.
// 	std::bitset<6>::reference erased() ;
	bool & isSpace() ;
	bool & isTetrahedron() ;
	bool & isDeadTetrahedron() ;
	bool & visited() ;//!< Marker. Useful not to lose ourselves isVertex the tree.
	bool & erased()  ;
	bool isSpace() const;
	bool isTetrahedron() const;
	bool isDeadTetrahedron() const;
	bool visited() const;//!< Marker. Useful not to lose ourselves isVertex the tree.
	bool erased() const ;
	bool isAlive() const ; //!< Accessor. Are we dead ?

	std::valarray<unsigned int> stepson ; ;//!< neighbours created later than ourselves
	std::valarray<unsigned int> neighbour ; //!< neighbours. three for triangles, any number for planes.
	std::valarray<unsigned int> son ;//!< items created by our destruction.
	
	//! \brief Constructor, takes the father and creator point as arguments
	/*! \a father is the father. Needed for the maintenance of the tree.
		\a c is the Creator Point. It is useful when building neighbourhood relationships. Also, it allowfor the removal of elements from the tree.
	 */
	DelaunayTreeItem3D( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *tree, DelaunayTreeItem3D * father, const Point * c) ;
	
	virtual ~DelaunayTreeItem3D() ;
		
	void removeNeighbour(DelaunayTreeItem3D * t) ; //!< Utility removes neighbour. Is safe.
	void addNeighbour(DelaunayTreeItem3D * t) ; //!< Utility adds neighbour. Is safe.
	
	DelaunayTreeItem3D * getNeighbour(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem3D * getSon(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	DelaunayTreeItem3D * getStepson(size_t i) const; //!< Accessor. returns the i<sup>th</sup> Neighbour. Safe
	
	virtual void kill(const Point * p) ; //!< kill and update the neighbourhood (livings do not neighbour the deads).
	virtual void erase(const Point * p) ;//!< kill and don't update the neighbourhood (do not use).
	

	
	void addStepson(DelaunayTreeItem3D * s) ;  //!< Utility adds stepson. Is safe.
	void removeStepson(DelaunayTreeItem3D * s) ;  //!< Utility removes stepson. Is safe.
	
	void addSon(DelaunayTreeItem3D * s) ;//!< Utility adds son. Is safe.
	void removeSon(DelaunayTreeItem3D * s) ;//!< Utility removes son. Is safe.
		
	void setStepfather(DelaunayTreeItem3D * s) ;  //!< Accessor. sets the stepfather.
	void setFather(DelaunayTreeItem3D * s) ;
	
	void clearVisited() ; //!< Accessor. We are not marked visited anymore.
	
	virtual bool isVertex(const Point *p) const = 0 ; //!< Test. Is this point \a p isVertex ?

	virtual std::vector< Point*> commonSurface(DelaunayTreeItem3D *t)  = 0; //!< What is the common edge with this item. returns a null pair if none
	
	virtual bool inCircumSphere(const Point &p) const = 0 ; //!< Test. Are we isVertex conflict with the point ?
	virtual bool onCircumSphere(const Point &p) const = 0 ;
	virtual bool isNeighbour(const  DelaunayTreeItem3D *) const = 0 ;  //!< Test. Are we a neighbour ?
	virtual void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,  Star3D *s) = 0 ; //!< Insert the point isVertex the Neighbourhood given by \a s. Returns the new elements
	virtual  void conflicts(std::valarray<bool> & visitedItems,  std::vector<DelaunayTreeItem3D *>  & ret ,const Point *p) ; //!< Test. Recursively give all elements isVertex conflict with \a p.
	virtual  void conflicts(std::valarray<bool> & visitedItems, std::vector<DelaunayTreeItem3D *> & ret, const Geometry *g) ;
	virtual void flatConflicts(std::valarray<bool> & visitedItems,  std::vector<DelaunayTreeItem3D *> &toTest, std::vector<DelaunayTreeItem3D *> & ret, const Geometry *g) ;
	virtual void flatConflicts(std::valarray<bool> & visitedItems,  std::vector<DelaunayTreeItem3D *> &toTest , std::vector<DelaunayTreeItem3D *> & ret,const Point *p, int threadid = -1) ;

	virtual void print() const = 0 ;

	virtual bool normalisedIn( const Point & p) const  = 0;
	
	size_t numberOfCommonVertices(const DelaunayTreeItem3D * s) const;
	
	bool isDuplicate( const DelaunayTreeItem3D * t) const ;
	
} ;


//! \brief Tetrahedron item the tree, defined by four points. 
/*!The points are also stored as a valarray of points(inherited from \c Triangle ). Those are stored clockwise but should only be used when creating the final mesh. They insure the proper orientation of the triangles.
*/
class DelaunayTetrahedron : virtual public TetrahedralElement, public DelaunayTreeItem3D
{
public:
	
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	DelaunayTetrahedron(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *tree,  DelaunayTreeItem3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,  Point * c) ;
	DelaunayTetrahedron(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *tree,  DelaunayTreeItem3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,   
	                     Point *p4,   Point *p5,   Point *p6, Point *p7,
	                     Point * c) ;
	DelaunayTetrahedron() ;
	
	
	virtual ~DelaunayTetrahedron() ;

	virtual bool isVertex(const Point * p) const ;
	bool isVertexByID(const Point * p) const ;
	bool hasVertexByID(const std::valarray<Point *> * p) const ;
	
	virtual std::vector< Point*> commonSurface(DelaunayTreeItem3D * t) ;

	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the triangle's CircumCircle and false if we are on or outside. 
	 */
	virtual bool inCircumSphere(const Point & p) const ;
	virtual bool onCircumSphere(const Point & p) const ;
	virtual bool isNeighbour( const DelaunayTreeItem3D * t) const;
	virtual void kill(const Point * p) ;
	virtual void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,   Star3D *s) ;
	std::vector<Point *> getIntegrationHints() const ;

	virtual void print() const;
	virtual bool normalisedIn(const Point & p) const ;

	void refresh(const TetrahedralElement *) ;
	
	std::valarray<std::valarray<Matrix> > & getElementaryMatrix() ;
	std::valarray<std::valarray<Matrix> > & getViscousElementaryMatrix() ;
	void clearElementaryMatrix() ;
	std::valarray<std::valarray<Matrix> > getNonLinearElementaryMatrix()  ;
// 	Vector getForces() ;
	Vector getNonLinearForces()  ;
	
	const GaussPointArray & getSubTriangulatedGaussPoints()  ;
	
	std::valarray<unsigned int> neighbourhood ;
	DelaunayTetrahedron * getNeighbourhood(size_t i) const ;
	void addNeighbourhood(DelaunayTetrahedron * t) ;
	void removeNeighbourhood(DelaunayTetrahedron *) ;
	
	virtual Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh() const {return nullptr ; } ;
	virtual Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh() const {return tree ; } ;
	
} ;
//! \brief Demi-space in the tree, defined by three points and a control point out of the demi-space. 
/*! The three first points form the frontier plane, whereas the last is chosen <em>outside</em> the demi-plane
*/

class DelaunayDemiSpace : public DelaunayTreeItem3D
{

public:
	Point pseudonormal ;
	
public:
	
	DelaunayDemiSpace(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *tree,  DelaunayTreeItem3D * father,   Point  * _one,   Point  * _two, Point  * _three,   Point  * p,   Point * c) ;
	
	virtual ~DelaunayDemiSpace() ;
	
	virtual std::vector< Point*> commonSurface(DelaunayTreeItem3D * t)  ;	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the demi plane, false if we are outside or on the limit.
	 */
	
	virtual bool inCircumSphere(const Point & p) const ;
	virtual bool onCircumSphere(const Point & p) const ;
	virtual bool isNeighbour( const DelaunayTreeItem3D * t) const ;
	virtual bool isVertex(const Point *p) const ;
	
	void merge(DelaunayDemiSpace * p) ;//!< ????
	
	//virtual void kill(Point * p) ;
	
	void insert(std::vector<DelaunayTreeItem3D *>& , Point *p, Star3D *s) ;
	
	virtual bool normalisedIn( const Point & p) const
	{
		return isOnTheSameSide( &p, third, first, second, fourth, tree->getInternalScale()) ;
	}


	virtual void print() const;

} ;

/** \brief light-weight class to replace dead Tetrahedrons*/
class DelaunayDeadTetrahedron : public DelaunayTreeItem3D
{
public:
	
	Point center ;
	double radius ;

	DelaunayDeadTetrahedron(DelaunayTetrahedron * father) ;
	
	virtual ~DelaunayDeadTetrahedron() ;
	
	virtual std::vector< Point*> commonSurface(DelaunayTreeItem3D * t) ;	

	virtual bool inCircumSphere(const Point & p) const ;
	virtual bool onCircumSphere(const Point & p) const ;
	virtual bool isNeighbour(const DelaunayTreeItem3D * t) const;
	virtual bool isVertex(const Point *p) const ;
	bool isVertexByID(const Point *p) ;
	
	virtual void insert(std::vector<DelaunayTreeItem3D *>& , Point *p, Star3D *s) { };
	const Point * getCircumCenter() const ;
	double getRadius() const ;
	
	virtual bool normalisedIn( const Point & p) const;

	void print() const;

} ;

//! \brief Root of the tree.
/*! Neither Plane nor triangle. The constructor should be extended to provide valid starting points for all sorts of geometries.
 */
class DelaunayRoot3D : public DelaunayTreeItem3D
{
public:
	DelaunayRoot3D(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *tree, Point * p0,  Point * p1,  Point * p2, Point * p3) ;
	
	virtual bool isVertex(const Point *p) const ;
	
	virtual bool inCircumSphere(const Point & p) const ;
	virtual bool onCircumSphere(const Point & p) const ;
			
	virtual std::vector< Point*> commonSurface(DelaunayTreeItem3D *t) ; //!< What is the common edge with this item. returns a null pair if none
	
	virtual bool isNeighbour(const  DelaunayTreeItem3D *) const { return false ; } 
	
	virtual void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,   Star3D *s) ;
	
	virtual void conflicts(std::valarray<bool> & visitedItems, std::vector<DelaunayTreeItem3D *> &, const Point *p )  ;
	
	virtual void conflicts( std::valarray<bool> & visitedItems, std::vector<DelaunayTreeItem3D *> &, const Geometry *g)  ;

	virtual void print() const ;

	virtual bool normalisedIn( const Point & p) const
	{
		return true ;
	}
} ;


//! \brief Neighbourhood of conflicing triangles.

/*! The star knows about all the conflicting elements when inserting a point. Thus, It has all the information available to create the correct neighbourhood relationships.
*/
class Star3D
{
protected:
	std::vector<const Point *> edge ;
	std::vector<DelaunayTreeItem3D *> treeitem ;
	
public:
	std::vector<DelaunayTreeItem3D *> cleanup ;
	~Star3D() ;
	Star3D(std::vector<DelaunayTreeItem3D *> *t, const Point *p) ;
	
	size_t size() ;
	
	const Point * getEdge(size_t i) const;
	
	void updateNeighbourhood() ;
	
	const Point *creator ;
} ;
//changed

/** \brief Mesh container. provides log n search facilities*/
class DelaunayTree3D :public Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>
{
friend class FeatureTree ;
friend class Geometry ;
	
std::valarray<bool> visitedItems ;

protected:
	double internalScale ;
	size_t global_counter ;
	bool neighbourhood ;
	bool fixedScale ;
	std::vector<Point *> additionalPoints ;
	
public:
	virtual double getInternalScale() const {return internalScale ;} ;
	virtual std::vector<Point * > & getAdditionalPoints() {return additionalPoints ; };
	virtual const std::vector<Point * > & getAdditionalPoints() const {return additionalPoints ;};

	virtual std::vector< DelaunayTetrahedron* > getElements() {return getTetrahedrons() ;};
	virtual std::vector< DelaunayTetrahedron* > getConflictingElements(const Mu::Point* p) 
	{
		std::vector< DelaunayTreeItem3D* > targets = conflicts(p) ;
		std::vector<DelaunayTetrahedron*> ret ;
		for(size_t i = 0 ; i < targets.size() ; i++)
		{
			if(targets[i]->isAlive() && targets[i]->isTetrahedron() && !targets[i]->isDeadTetrahedron() )
				ret.push_back(static_cast< DelaunayTetrahedron *> (targets[i])) ;
		}
		return ret ;
	};
	virtual std::vector< DelaunayTetrahedron* > getConflictingElements(const Geometry* g) 
	{
		return conflicts(g) ;
	}
	
	std::vector<DelaunayTreeItem3D *> tree ;
public:
	const size_t & getLastNodeId() const {return global_counter ;}  ;
	std::vector<DelaunayDemiSpace *> space ;
	
	DelaunayTree3D ( Point * p0,  Point *p1,  Point *p2, Point *p3) ;
	
	virtual ~DelaunayTree3D() ;
	
	virtual size_t addToTree(DelaunayTreeItem3D * it)
	{
		tree.push_back(it);
		return tree.size()-1;
	}
	
	virtual DelaunayTreeItem3D * getInTree(int index) 
	{
		return tree[std::abs(index)] ;
	}
	
// 	virtual std::vector<DelaunayTreeItem3D *> & getTree() { return tree ;}
// 	virtual const std::vector<DelaunayTreeItem3D *> & getTree() const { return tree ;}
	void insert( Point *p) ;
	void insert( Segment *s) ;
	
	/** \brief Get the list of triangles in comflict with a point. This method is O(log(n)), except when construction of the mesh has led to circular stepparenthoods. In that case, it is O(log(n)) in general, and O(n) 5% of the times.
	 * 
	 * @param p Point to check.
	 * @return the list of triangles in conflict with p. A triangle is in conflict if the point is tricly in the circumcircle.
	 */
	std::vector<DelaunayTreeItem3D *> conflicts( const Point *p) ;

	virtual void extrude(double dt) ;
	
	virtual void extrude(const Vector & dt) ;
	
	/** \brief Get the list of triangles in conflict with a given geometry.
	 * 
	 * @param g test geometry.
	 * @return all the triangles for which at least one node is in the geometry.
	 */
	std::vector<DelaunayTetrahedron *> conflicts( const Geometry *g) ;
	
	/** \brief Find the boundaries of the triangulation
	 * 
	 * @return the planes bordering the triangulated domain.
	 */
	std::vector<DelaunayDemiSpace *> getConvexHull() ;
	
	/** \brief Return the result of the trianglulation.
	 * 
	 * @return all the living triangles resulting from the triangulation.
	 */
	std::vector<DelaunayTetrahedron *> getTetrahedrons(bool buildNeighbourhood = true) ;
	
	void addSharedNodes(size_t nodes_per_side, const TetrahedralElement * father = nullptr) ; 
	void addSharedNodes(size_t nodes_per_side, size_t time_planes = 2, double timestep = 2, const TetrahedralElement * father = nullptr) ;
	void addSharedNodes(DelaunayTree3D * dt) ;
	virtual void setElementOrder(Order elemOrder, double dt = 0) ;
	
	void addElements(std::vector<DelaunayTreeItem3D *> & cons, Point * p) ;
	void refresh(const TetrahedralElement *father) ;
	
	size_t numPoints() const;

	void purge() ;

	void print() const;

} ;

std::pair<std::vector<DelaunayTetrahedron *>, std::vector<Point *> > quad(const DelaunayTetrahedron * t) ;

} ;



//! Make \a t0 and \a t1 Neighbours. Safe.
void makeNeighbours( Mu::DelaunayTreeItem3D *t0, Mu::DelaunayTreeItem3D *t1 ) ;

#endif  //__DELAUNAY_3D_H_
