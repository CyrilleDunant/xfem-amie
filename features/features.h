// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __FEATURES_H_
#define __FEATURES_H_

#include "../geometry/geometry_base.h"
#include "../elements/elements.h" 
#include "../delaunay.h"
#include "../delaunay_3d.h"
#include "../samplingcriterion.h"
#include "../solvers/assembly.h"

#include <valarray>
#include <deque>
#include <iostream>
#include <algorithm>

static const double DEFAULT_BOUNDARY = 1 ;
static const size_t DEFAULT_POINTNUMBER = 1000 ;

namespace Mu
{
/** Feature. 
 * a feature is the essential unit of description of a given test setup. A feature can be the sample itself, 
 * a crack, an inclusion, or a hole. Features are defined both by a constitutive law of comportment, 
 * given by the Cauchy-Green stress tensor, and by a geometrical definition.
 * 
 * Features are responsible for providing geometry-geometry interactions between them.
 */
class Feature : public Geometry
{
protected:
	/** Boundary for the sampling of the feature.
	 * This boundary is necessary, as the triangulation wich would allow for finer control does 
	 * not exist at the sampling stage. As most features are convex, if a distance is kept between 
	 * the points of the parent and the points of the children, the boundaries of the actual
	 * geometry defining the sample will automatically be respected by Delaunay triangulation:
	 * Delaunay keeps convex sets.
	 * 
	 */
	Geometry * boundary ;
	Geometry * boundary2 ;
	
	/** Influence radius. deprecated.
	 */
	double infRad ;
	
	/** Children of the feature.
	 * 
	 * Features are arranged in a tree of features, a parent containing all its children, from a geometrical 
	 * point of view. In fact, children can go beyond their parents boundaries, and this may be used to 
	 * produce various irregular shapes. This generally works well, as long as one doesn't wish to use auto-refinement.
	 * 
	 */
	std::vector<Feature *> m_c ;
	
	/** Father of the feature. 
	 */
	Feature * m_f ;
	
	/** Constitutive tensor describing the elastic material of the feature.
	 * 
	 * This can (and should) be extended to allow for various non-linear laws of behaviour. This could be done by 
	 * including physics.h and having a pointer to the law of behaviour and not simply the C-G tensor.
	 */
	Form * behaviour ;
	
	double now ;

public:
	
	bool isEnrichmentFeature ;
	bool isCompositeFeature ;
	bool isVirtualFeature ;
	
	/** Feature constructor.
	 * 
	 * @param father sets the father.
	 */
	Feature(Feature *father) ;

	Feature() ;
	
	/** Feature constructor.
	 * 
	 * @param father sets the father.
	 * @param b      sets the boundary.
	 */
	Feature(Feature *father, Geometry * b) ;
	
	virtual ~Feature() ;
	
	virtual void setBoundary(Geometry * g) ;
	virtual const Geometry * getBoundary() const ;
	virtual Geometry * getBoundary() ;
	/**  Is Point in boundary?
	 * 
	 * @param v  point to check. 
	 * @return   true if in boundary, or in boundary of one of the descendants.
	 */
	virtual bool inBoundary(const Point & v) const ;
	
	 /**  Is Point in boundary?
	 * 
	 * @param v  point to check. 
	 * @return   true if in boundary, or in boundary of one of the descendants.
	 */
	virtual bool inBoundary(const Point *v) const ;
	
	virtual bool inBoundaryLayer(const Point *v) const ;
	 
	/** Sets Influence radius. deprecated.
	 */
	virtual void setInfluenceRadius(double r) ;
	 
	
	/** Set the Cauchy-Green Strain Tensor.
	 * 
	 * @param m Matrix to point to. Be careful not to delete the matrix while it is in use by any feature!
	 */
	virtual void setBehaviour(Form * b) ;
	
	
	/** Get the Cauchy-Green Strain Tensor.
	 * 
	 * @return the Cauchy-Green Strain tensor.
	 */
	virtual Form * getBehaviour( const Point & p ) ;
	
	/** Add a children to the feature.
	 * 
	 * @param f Children to add.
	 */
	virtual void addChild(Feature *f) ;
	virtual void removeChild(Feature *f) ;
	
	
	/** Return i<sup>th</sup> child.
	 * 
	 * @param i index of the child to return.
	 * @return The i<sup>th</sup> child.
	 */
	virtual Feature * getChild(size_t i) const ;
	
	
	/**Return the father.
	 * 
	 * @return thefather Feature. It can be NULL !
	 */
	virtual Feature * getFather() const;
	
	
	/** Return the children
	 * 
	 * @return a pointer to the m_c member. This is a vector containing the pointers to the children.
	 */
	virtual const std::vector<Feature *> & getChildren() const;
	virtual std::vector<Feature *> & getChildren();
	virtual std::vector<Feature *> getDescendants() const ;
	
	/** Reparent the feature.
	 * 
	 * @param f the new father. The feature is automatically added to the children of the new parent.
	 */
	virtual void setFather(Feature *f) ;
	
	/** Return all the triangle stricly <em>in</em> the feature.
	 * 
	 * The DealunayTree::conflict() function returns all the triangles having <em>at least</em> on point in the boundary of the feature. This function returns only those triangles whose <em>center</em> lies <em>in</em> the feature.
	 * 
	 * @param dt the DelaunayTree from which to get the triangles.
	 * @return the vector of triangles satisfying the condition center \f$ \in \f$ Feature. 
	 */
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  = 0;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(DelaunayTree_3D * dt)  = 0;
	virtual std::vector<DelaunayTriangle *> getBoundingTriangles( DelaunayTree * dt) ;
	virtual std::vector<DelaunayTetrahedron *> getBoundingTetrahedrons(DelaunayTree_3D * dt) ;
	
	/** Check for interaction.
	 * 
	 * Two features are said to be interactiong if ther boundaries are intersecting.
	 * 
	 * @param f Feature for which to check for the interaction.
	 * @return true if interacting.
	 */
	virtual bool interacts(Feature * f) const = 0;
	
	
	/** Insert a point on the bounding surface of the feature.
	 * 
	 * @param i index of the Point after which to do the insertion.
	 * @return a pointer to the Point just inserted.
	 */
	virtual Point * pointAfter(size_t i) = 0 ;
	
	
	/** Define a vector of geometries to use for sucessive refinement.
	 * 
	 * @return the vector of geometries.
	 */
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const = 0 ;
	
	/** Projects a point on the boundary of the feature.
	 * 
	 * the point is not copied, its coordinates are simply changed.
	 * 
	 * @param  p Point to project.
	 */
	virtual void project(Point * p) const = 0 ;
	
// 	virtual std::vector<Point> intersection(const Geometry * f) const = 0;
// 	virtual bool intersects(const Feature * f) const = 0;
	
	
	virtual void print() const = 0 ;
	
	virtual const std::valarray<Point *> & getBoundingPoints() const = 0 ;
	virtual std::valarray<Point *> & getBoundingPoints() = 0 ;
	virtual std::vector<Point *> doubleSurfaceSampling() ;
	virtual const Point & getBoundingPoint(size_t) const = 0 ;
	virtual Point & getBoundingPoint(size_t)  = 0 ;
	virtual std::valarray<Point *> & getInPoints() = 0 ;
	virtual const std::valarray<Point *> & getInPoints() const = 0 ;
	virtual const Point & getInPoint(size_t) const = 0 ;
	virtual Point & getInPoint(size_t) = 0 ;
	virtual bool in( const Point & ) const = 0 ;
	virtual const Point & getCenter() const = 0 ;
	virtual double getRadius() const = 0 ;
	virtual double area() const = 0 ;
	virtual void sample(size_t n) = 0 ;
	
	virtual bool isVoid( const Point &) const = 0 ;
	
} ;

class VirtualFeature : virtual public Feature
{
public:

	VirtualFeature(Feature *father) : Feature(father) { isVirtualFeature = true ;}
	
	VirtualFeature(Feature *father, Geometry * b) : Feature(father, b)  { isVirtualFeature = true ;}
	
	virtual ~VirtualFeature() { };

	virtual void print() const = 0 ;
	virtual Feature * getSource() const = 0;
} ;

class CompositeFeature : virtual public Feature
{
protected:
	std::vector<VirtualFeature *> components ;
public:

	CompositeFeature(Feature *father) : Feature(father) { isCompositeFeature = true ;}
	
	CompositeFeature(Feature *father, Geometry * b) : Feature(father, b)  { isCompositeFeature = true ;}
	
	virtual ~CompositeFeature() ;

	std::vector<VirtualFeature *> & getComponents() ;
	const std::vector<VirtualFeature *> & getComponents() const ;

	virtual void print() const = 0 ;
} ;



class EnrichmentFeature : virtual public Feature
{
public:
	
	/** Feature constructor.
	 * 
	 * @param father sets the father.
	 */
	EnrichmentFeature(Feature *father) : Feature(father) { this->isEnrichmentFeature = true ;}
	
	/** Feature constructor.
	 * 
	 * @param father sets the father.
	 * @param b      sets the boundary.
	 */
	EnrichmentFeature(Feature *father, Geometry * b) : Feature(father, b) { this->isEnrichmentFeature = true ;}
	
	virtual ~EnrichmentFeature() { };
	
	virtual std::vector<Point *> getSamplingPoints() const = 0 ;
	
	virtual bool enrichmentTarget(DelaunayTriangle * t) = 0;
	virtual void enrich(size_t & , DelaunayTree * dtree) = 0 ;
	
	virtual void step(double dt, Vector *, const DelaunayTree * dtree) = 0;
	virtual void snap(DelaunayTree * dtree) = 0 ;
	
	virtual bool moved() const = 0;
	
} ;

class Pixel
{
protected:
	std::vector<Feature *> features ;
	Point tl ;
	Point tr ;
	Point bl ;
	Point br ;
	bool filled ;
	std::valarray<Pixel> pixels ;
	void refine() ;
	int computeFillFactor() const;
// 	const short level ;
// 	const short levels ;
public:
	Pixel();
	~Pixel();
	Pixel(double x, double y, double s) ;

	const std::vector<Feature *> & getFeatures() const;
	
	std::vector<Feature *> & getFeatures();
	
	bool in(const Point & p) const;

	bool coOccur(const Geometry * inc) const;
	bool coOccur(const Point & p) const;
	void coOccuringFeatures(std::vector<Feature *>&, const Geometry * inc) const ;
	void coOccuringFeatures(std::vector<Feature *>&, const Point & p) const ;
	void remove(Feature * inc);
	
	bool add(Feature * inc);
	void forceAdd(Feature * inc) ;

	void print() const ;

} ;

class Voxel
{
protected:
	std::vector<Feature *> features ;
	Point tlf ;
	Point trf ;
	Point blf ;
	Point brf ;
	Point tlb ;
	Point trb;
	Point blb ;
	Point brb ;
	bool filled ;

	std::vector<Voxel *> pixels ;
	void refine() ;
	int computeFillFactor() const;
// 	const short level ;
// 	const short levels ;
public:
	Voxel();
	
	Voxel(double x, double y, double z ,double s) ;

	~Voxel();

	const std::vector<Feature *> & getFeatures() const;
	
	std::vector<Feature *> & getFeatures();
	
	bool in(const Point & p) const;

	bool coOccur(const Geometry * inc) const;
	bool coOccur(const Point & p) const;
	void coOccuringFeatures(std::vector<Feature *>&, const Geometry * inc) const ;
	void coOccuringFeatures(std::vector<Feature *>&, const Point & p) const ;
	void remove(Feature * inc);
	
	bool add(Feature * inc);
	void forceAdd(Feature * inc) ;

	void print() const ;
	
	bool isFilled() const ;
	
	Point center() const ;

} ;

class Grid
{
protected:
	std::valarray< std::valarray<Pixel *> > pixels;
	double x ;
	double y ;
	size_t lengthX ;
	size_t lengthY ;
	
	double psize ;
public:
		
	Grid(double sizeX, double sizeY, int div );
	
	~Grid() ;
	bool add(Feature * inc);
	void forceAdd(Feature * inc) ;
	std::vector<Feature *> coOccur(const Geometry * geo) const ;
	std::vector<Feature *> coOccur(const Point & p) const ;
} ;

class Grid3D
{
protected:
	std::valarray<std::valarray< std::valarray<Voxel *> > > pixels;
	std::vector<Voxel *> unfilledpixel ;
	double x ;
	double y ;
	double z ;
	size_t lengthX ;
	size_t lengthY ;
	size_t lengthZ ;
	
	double psize ;
	int dirtyCounter ;
public:
		
	Grid3D(double sizeX, double sizeY, double sizeZ, int div );
	Point randomFreeCenter() const ;
	~Grid3D() ;
	bool add(Feature * inc);
	void forceAdd(Feature * inc) ;
	std::vector<Feature *> coOccur(const Geometry * geo) const ;
	std::vector<Feature *> coOccur(const Point & p) const;
	double fraction() const ;
} ;

/** Container for the features defining the setup.
 * 
 * The feature tree is responsible for all global operations: meshing, matrix assembly, 
 * applying boundary conditions.
 * 
 * Postprocessing the results is also done here.
 * 
 */
class FeatureTree 
{
	
protected:
	/** Contains all the features. */
	std::vector<Feature *> tree ;

	/**For fast Access*/
	Grid * grid ;
	Grid3D * grid3d ;
	
	/** Contains the mesh in the form of a delaunay tree. 
	 * The mesh is generated with linear triangles, and when it is final, midpoints are added and 
	 * projected. No operations should add midpoints before meshing is complete.
	 */
	DelaunayTree * dtree ;
	DelaunayTree_3D * dtree3D ;
	
	TetrahedralElement *father3D  ;
	TriElement *father2D  ;

	bool meshChange ;
	bool solverConvergence ;
	bool enrichmentChange ;
	
	/** List of points used for the mesh.
	 * Each point is associated with the feature from whose discretiation it was generated.
	 */
	std::deque<std::pair<Point *, Feature *> > meshPoints;
	
	/** List of the elements.
	 */
	std::vector<DelaunayTetrahedron * > elements3D;
	
	/** Assembly used for the generation of the stiffness matrix and the solving of the problem.
	 */
	Assembly * K ;
	
	//Assembly2D
	/** Order to which Elements should be brought once all operations are accomplished.
	 */
	Order elemOrder ;
	
	/** Insert triangle midpoints. */
	void makeToOrder() ;
	
	/** Project all points on their respective boundaries.*/
	void stitch() ;
	
	void renumber() ;
public:
		/** Generate the triangulation.
	 * Once the sampling is done, the sampling points are fed into a Delaunay Triangulation algorithm, 
	 * which generates the triangles. The mesh is still composed of 3 points triangles at this point.
	 * 
	 * @param correctionSteps additional steps where points are inserted in incorrect tetrahedrons.
	 */
	void generateElements( size_t correctionSteps = 0) ;
	
	/** Perform the assembly.
	 * 
	 * The mesh is completed with eventual intermediate points (if higher order elements were asked for), and elements 
	 * are assembled in K.
	 * 
	 */
	void assemble() ;


	void print() const;
	void printForFeature(const Feature *f) const;
	
	double now ;
	
	size_t numdofs ;
	
	bool stitched ;
	bool renumbered ;
	bool needAssembly ;
	bool initialized ;
	
	Point * checkElement( const DelaunayTetrahedron * t ) const;
	Point * checkElement( const DelaunayTriangle * t ) const ;
	Feature * getFeatForTetra( const DelaunayTetrahedron * t ) const;

	bool solverConverged() const ;
	bool meshChanged() const ;
	bool enrichmentChanged() const ;
	
public:
	
	/** Initialise the feature tree with the parent feature.
	 * 
	 * @param first Parent of all features. This shoud typically be the sample itself.
	 * @return 
	 */
	FeatureTree(Feature *first) ;
	virtual ~FeatureTree() ;
	
	/** Add feature as the daghter of another.
	 * 
	 * A feature being the daughter of another typically implies that it is fully contained therein. 
	 * This is however not necessarilly the case.
	 * 
	 * @param father Parent feature.
	 * @param t daughter feature.
	 */
	void addFeature(Feature *father, Feature * t) ;

	void twineFeature(CompositeFeature * father, CompositeFeature * f) ;
	Vector getDisplacements() const ;
	
	/** Generate the sample points for all the features. The features are passed a sampling 
	 * argument proportionnal to their area compared with the area of the root feature. 
	 * If the number is lower than 10, than the argument passed is 10.
	 * 
	 * @param npoints number of sampling points on the boundary of the sample. This number must be 
	 * divisable by four if the sample is a rectangle.
	 */
	void sample(size_t npoints) ;
	
	/** Attempt to enhance the mesh, based on a sampling citerion.
	 * 
	 * Criterions are typically made to check the adequacy of the geometry of the triangles generated 
	 * during the triangulation phase. They also provide a hint as to the placement of new points which could 
	 * enhance the mesh quality.
	 * 
	 * @param nit Number of refinement iterations to perform.
	 * @param cri Citerion to use.
	 */
	void refine(size_t nit, SamplingCriterion *cri) ;
	
	
	/** Refine the mesh around the features.
	 * 
	 * Features provide a set of geometries which are targets for successive refinement. Refinement is 
	 * done by inserting a point in the center of each triangle in the zone.
	 * 
	 */
	void refine(size_t level) ;
	
	
	/** Set the constitutive tensor of the given feature.
	 * 
	 * <b>Beware!</b> this function <i>deletes</i> the previous tensor. Make sure is is not in use by another feature. 
	 * 
	 * @param m pointer to the tensor. 
	 * @param f Feature to affect.
	 */
	void setStrainTensor(Matrix * m, Feature * f) ;
	
	/** set the target order of Elements
	 * 
	 * @param ord order of elements to use.
	 * 
	 * \todo Make it work for values other than 1 and 2.
	 */
	void setOrder(Order ord) ;
	
	/** Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * strain values at the mesh points.
	 * 
	 * \todo make it cleaner and element-order independant.
	 * 
	 * \todo chose an apropriate format for storing the result. This is necessary for visco-elastic behaviour.
	 * 
	 * @param disp displacements
	 * @return strain.
	 */
	Vector strainFromDisplacements() const ;
	
	/** Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * stress values at the mesh points.
	 * 
	 * \todo chose an apropriate format for storing the result. This is necessary for visco-elastic behaviour.
	 * 
	 * @param disp displacements
	 * @return stress.
	 */
	Vector stressFromDisplacements() const ;
	
	std::pair<Vector , Vector > getStressAndStrain() ;
	
	size_t numPoints() const ;
	
	bool step(double dt) ;
	void stepBack() ;
	void elasticStep() ;

	std::deque<std::pair<Point *, Feature *> >::iterator  begin() ;
	std::deque<std::pair<Point *, Feature *> >::iterator end() ;
	
	Assembly * getAssembly() ;
	std::vector<DelaunayTriangle *> getTriangles();
	std::vector<DelaunayTetrahedron *> getTetrahedrons() ;
		
	std::vector<DelaunayTriangle *> getBoundingTriangles(Feature * f = NULL) ;	
	
	Form * getElementBehaviour(const DelaunayTriangle *) const ;
	Form * getElementBehaviour(const DelaunayTetrahedron *) const ;
	
	void insert(Point * p ) ;
	
	DelaunayTree * getDelaunayTree() ;
	DelaunayTree_3D * getDelaunayTree3D() ;
	
	bool inRoot(const Point &p) const ;
	
	Feature * getFeature(size_t i)
	{
		return tree[i] ;
	}
	
	void initializeElements() ;
	
	double getMaximumDisplacement() const ;
	double getMinimumDisplacement() const ;
	
	bool is3D() const ;
	bool is2D() const ;

	std::vector<DelaunayTriangle> getSnapshot2D() const ;
	
} ;



struct PairPointFeatureEqual
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first == *p2.first ;
	}
} ;

struct PairPointFeatureLess_Than
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first < *p2.first ;
	}
} ;



} ;

#endif // __FEATURES_H_
