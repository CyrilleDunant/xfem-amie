
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
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_3D.h"
#include "../elements/elements.h"
#include <vector>
#include <complex>
#include <list>
#include <set>
#include <iostream>
#include <bitset>

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
	std::bitset<6> state ;
	std::bitset<6>::reference dead() ;
	bool dead() const ;
public:
	const Point * killer ; //!< Point killer.
	const Point * creator ; //!< Point creator.
	int index ;
	
	DelaunayTree3D * tree ;
	int father ; //!< Item destroyed by the insertion of the creator point.
	int stepfather ; //!< Still-alive neighbour of the father.
	
	Point * first ; //!< Defining point. Is always <em>isVertex</em> the item.
	Point * second ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * third ; //!<  Defining point. Is always <em>isVertex</em> the item.
	Point * fourth ; //!<  Defining point. Function differs if item is a triangle or point.
	
	
	std::bitset<6>::reference isSpace() ;
	std::bitset<6>::reference isTetrahedron() ;
	std::bitset<6>::reference isDeadTetrahedron() ;
	std::bitset<6>::reference visited() ;//!< Marker. Useful not to lose ourselves isVertex the tree.
	std::bitset<6>::reference erased() ;
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
	DelaunayTreeItem3D( DelaunayTree3D *tree, DelaunayTreeItem3D * father, const Point * c) ;
	
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

	virtual std::vector< Point*> commonSurface(const DelaunayTreeItem3D *t) const  = 0; //!< What is the common edge with this item. returns a null pair if none
	
	virtual bool inCircumSphere(const Point &p) const = 0 ; //!< Test. Are we isVertex conflict with the point ?
	//newly aded
	virtual bool isNeighbour( const DelaunayTreeItem3D *) const = 0 ;  //!< Test. Are we a neighbour ?
	virtual void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,  Star3D *s) = 0 ; //!< Insert the point isVertex the Neighbourhood given by \a s. Returns the new elements
	virtual  void conflicts(std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > &,const Point *p) ; //!< Test. Recursively give all elements isVertex conflict with \a p.
	virtual  void conflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > &, const Geometry *g) ;
	std::vector<DelaunayTreeItem3D *> flatConflicts(std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > & ret, const Geometry *g) ;
	std::vector<DelaunayTreeItem3D *> flatConflicts(std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> >& ret,const Point *p) ;

	virtual void print() const = 0 ;

	virtual bool in( const Point & p) const  = 0;
	
	size_t numberOfCommonVertices(const DelaunayTreeItem3D * s) const;
	
	bool isDuplicate(const DelaunayTreeItem3D * t) const ;
	
} ;


//! \brief Tetrahedron item the tree, defined by four points. 
/*!The points are also stored as a valarray of points(inherited from \c Triangle ). Those are stored clockwise but should only be used when creating the final mesh. They insure the proper orientation of the triangles.
*/
class DelaunayTetrahedron : virtual public TetrahedralElement, public DelaunayTreeItem3D
{
public:
	
	GEO_DERIVED_OBJECT(Tetrahedron) ;
	
	DelaunayTetrahedron(DelaunayTree3D *tree,  DelaunayTreeItem3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,  Point * c) ;
	DelaunayTetrahedron(DelaunayTree3D *tree,  DelaunayTreeItem3D * father,   Point *p0,   Point *p1,   Point *p2, Point *p3,   
	                     Point *p4,   Point *p5,   Point *p6, Point *p7,
	                     Point * c) ;
	DelaunayTetrahedron() ;
	
	
	virtual ~DelaunayTetrahedron() ;

	bool isVertex(const Point * p) const ;
	bool isVertexByID(const Point * p) const ;
	bool hasVertexByID(const std::valarray<Point *> * p) const ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem3D * t) const ;

	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the triangle's CircumCircle and false if we are on or outside. 
	 */
	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem3D * t) const ;
	void kill(const Point * p) ;
	void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,   Star3D *s) ;
	std::vector<Point *> getIntegrationHints() const ;

	void print() const;

	void refresh(const TetrahedralElement *) ;
	
	std::vector<std::vector<Matrix> > & getElementaryMatrix() ;
	void clearElementaryMatrix() ;
	std::vector<std::vector<Matrix> > getNonLinearElementaryMatrix()  ;
	Vector getForces() const ;
	Vector getNonLinearForces() const ;
	
	GaussPointArray getSubTriangulatedGaussPoints() const ;
	
	std::valarray<unsigned int> neighbourhood ;
	DelaunayTetrahedron * getNeighbourhood(size_t i) const ;
	void addNeighbourhood(DelaunayTetrahedron * t) ;
	void removeNeighbourhood(DelaunayTetrahedron *) ;
	
} ;
//! \brief Demi-space in the tree, defined by three points and a control point out of the demi-space. 
/*! The three first points form the frontier plane, whereas the last is chosen <em>outside</em> the demi-plane
*/

class DelaunayDemiSpace : public DelaunayTreeItem3D
{

public:
	Point pseudonormal ;
	
public:
	
	DelaunayDemiSpace(DelaunayTree3D *tree,  DelaunayTreeItem3D * father,   Point  * _one,   Point  * _two, Point  * _three,   Point  * p,   Point * c) ;
	
	virtual ~DelaunayDemiSpace() ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem3D * t) const ;	
	/** Check for point location.
	 * 
	 * @param p Point to check.
	 * @return true if we are in the demi plane, false if we are outside or on the limit.
	 */
	
	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem3D * t) const ;
	bool isVertex(const Point *p) const ;
	
	void merge(DelaunayDemiSpace * p) ;//!< ????
	
	//virtual void kill(Point * p) ;
	
	void insert(std::vector<DelaunayTreeItem3D *>& , Point *p, Star3D *s) ;
	
	bool in( const Point & p) const
	{
		return isOnTheSameSide( &p, third, first, second, fourth) ;
	}


	void print() const;

} ;

/** \brief light-weight class to replace dead Tetrahedrons*/
class DelaunayDeadTetrahedron : public DelaunayTreeItem3D
{
public:
	
	double x ;
	double y ;
	double z ;
	double radius ;

	DelaunayDeadTetrahedron(DelaunayTetrahedron * father) ;
	
	virtual ~DelaunayDeadTetrahedron() ;
	
	std::vector< Point*> commonSurface(const DelaunayTreeItem3D * t) const ;	

	bool inCircumSphere(const Point & p) const ;
	bool isNeighbour( const DelaunayTreeItem3D * t) const ;
	bool isVertex(const Point *p) const ;
	bool isVertexByID(const Point *p) const ;
	
	void insert(std::vector<DelaunayTreeItem3D *>& , Point *p, Star3D *s) { };
	
	bool in( const Point & p) const;

	void print() const;

} ;

//! \brief Root of the tree.
/*! Neither Plane nor triangle. The constructor should be extended to provide valid starting points for all sorts of geometries.
 */
class DelaunayRoot3D : public DelaunayTreeItem3D
{
public:
	DelaunayRoot3D(DelaunayTree3D *tree, Point * p0,  Point * p1,  Point * p2, Point * p3) ;
	
	bool isVertex(const Point *p) const ;
	
	bool inCircumSphere(const Point & p) const ;
			
	std::vector< Point*> commonSurface(const DelaunayTreeItem3D *t) const ; //!< What is the common edge with this item. returns a null pair if none
	
	bool isNeighbour( const DelaunayTreeItem3D *) const { return false ; } 
	
	void insert(std::vector<DelaunayTreeItem3D *> &, Point *p,   Star3D *s) ;
	
	void conflicts(std::pair<std::vector<DelaunayTreeItem3D *>, std::vector<DelaunayTreeItem3D *> > &, const Point *p )  ;
	
	void conflicts( std::pair<std::vector<DelaunayTetrahedron *>, std::vector<DelaunayTreeItem3D *> > &, const Geometry *g)  ;

	void print() const ;

	bool in( const Point & p) const
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
class DelaunayTree3D 
{
friend class FeatureTree ;
friend class Geometry ;
	
protected:
	size_t global_counter ;
	bool neighbourhood ;
	std::vector<Point *> additionalPoints ;
public:

	std::vector<DelaunayTreeItem3D *> tree ;
	std::vector<DelaunayDemiSpace *> space ;
	
	DelaunayTree3D ( Point * p0,  Point *p1,  Point *p2, Point *p3) ;
	
	virtual ~DelaunayTree3D() ;
	
	void insert( Point *p) ;
	void insert( Segment *s) ;
	
	/** \brief Conditionnal insertion of points.
	 * 
	 * @param p Point to insert.
	 * @param v vector of samplingCriterions.
	 * @param minScore minimum score to insert the point. The score is calculated as : 
	 *	\f$ \frac{N_{tot}}{N_{pass}}\f$
	 */
// 	void insertIf( Point *p, std::vector<SamplingCriterion_3D *> v, double minScore ) ;
	
	/** \brief Get the list of triangles in comflict with a point. This method is O(log(n)), except when construction of the mesh has led to circular stepparenthoods. In that case, it is O(log(n)) in general, and O(n) 5% of the times.
	 * 
	 * @param p Point to check.
	 * @return the list of triangles in conflict with p. A triangle is in conflict if the point is tricly in the circumcircle.
	 */
	std::vector<DelaunayTreeItem3D *> conflicts( const Point *p) const ;
	
	/** \brief Get the list of triangles in conflict with a given geometry.
	 * 
	 * @param g test geometry.
	 * @return all the triangles for which at least one node is in the geometry.
	 */
	std::vector<DelaunayTetrahedron *> conflicts( const Geometry *g) const ;
	
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
	
	void addSharedNodes(size_t nodes_per_side, const TetrahedralElement * father = NULL) ; 
	void addSharedNodes(size_t nodes_per_side, size_t time_planes = 2, double timestep = 2, const TetrahedralElement * father = NULL) ;
	
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
