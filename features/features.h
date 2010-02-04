// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __FEATURES_H_
#define __FEATURES_H_

#include "../geometry/geometry_base.h"
#include "../elements/elements.h" 
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "../utilities/samplingcriterion.h"
#include "../solvers/assembly.h"

#include <valarray>
#include <deque>
#include <iostream>
#include <algorithm>

static const double DEFAULT_BOUNDARY = 1 ;
static const size_t DEFAULT_POINTNUMBER = 1000 ;

namespace Mu
{
/** \brief Feature.
 * 
 * a feature is the essential unit of description of a given test setup. A feature can be the sample itself, 
 * a crack, an inclusion, or a hole. Features are defined both by a constitutive law of comportment, 
 * given by the Cauchy-Green stress tensor, and by a geometrical definition.
 * 
 * Features are responsible for providing geometry-geometry interactions between them.
 */
class Feature : public Geometry
{
protected:
	/** \brief Boundary for the sampling of the feature.
	 *
	 * This boundary is necessary, as the triangulation wich would allow for finer control does 
	 * not exist at the sampling stage. As most features are convex, if a distance is kept between 
	 * the points of the parent and the points of the children, the boundaries of the actual
	 * geometry defining the sample will automatically be respected by Delaunay triangulation:
	 * Delaunay keeps convex sets.
	 * 
	 */
	Geometry * boundary ;
	Geometry * boundary2 ;
	
	/** \brief Influence radius. deprecated.
	 */
	double infRad ;
	
	/** \brief Children of the feature.
	 * 
	 * Features are arranged in a tree of features, a parent containing all its children, from a geometrical 
	 * point of view. In fact, children can go beyond their parents boundaries, and this may be used to 
	 * produce various irregular shapes. This generally works well, as long as one doesn't wish to use auto-refinement.
	 * 
	 */
	std::vector<Feature *> m_c ;
	
	/** \brief Father of the feature. 
	 */
	Feature * m_f ;
	
	/** \brief Constitutive Law describing the behaviour of the feature.
	 * 
	 */
	Form * behaviour ;
	
	double now ;

public:
	
	bool isEnrichmentFeature ;
	bool isCompositeFeature ;
	bool isVirtualFeature ;
	
	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 */
	Feature(Feature *father) ;

	Feature() ;
	
	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 * @param b      sets the boundary.
	 */
	Feature(Feature *father, Geometry * b) ;
	
	virtual ~Feature() ;
	
	virtual void setBoundary(Geometry * g) ;
	virtual const Geometry * getBoundary() const ;
	virtual Geometry * getBoundary() ;

	/** \brief Is Point in boundary?
	 * 
	 * @param v  point to check. 
	 * @return   true if in boundary, or in boundary of one of the descendants.
	 */
	virtual bool inBoundary(const Point & v) const ;
	
	 /** \brief Is Point in boundary?
	 * 
	 * @param v  point to check. 
	 * @return   true if in boundary, or in boundary of one of the descendants.
	 */
	virtual bool inBoundary(const Point *v) const ;
	
	/** \brief If the feature contains an internal frontier, return true if the argument is in the boundary, but out of the interanl frontier*/
	virtual bool inBoundaryLayer(const Point *v) const ;
	 
	/** \brief Sets Influence radius.
	 */
	virtual void setInfluenceRadius(double r) ;
	 
	
	/** \brief Set the Behaviour
	 * 
	 * @param b Behaviour of the feature. getCopy() from this behaviour will be called to generate the behaviour of elements depending from this feature.
	 */
	virtual void setBehaviour(Form * b) ;
	
	/** \brief Get the Cauchy-Green Strain Tensor.
	 * 
	 * @return the Cauchy-Green Strain tensor.
	 */
	virtual Form * getBehaviour( const Point & p ) ;
	
	/** \brief Add a child to the feature.
	 * 
	 * @param f Children to add.
	 */
	virtual void addChild(Feature *f) ;

	/** \brief Remove a child of the feature.
	 * 
	 * @param f Children to add.
	 */
	virtual void removeChild(Feature *f) ;
	
	
	/** \brief Return i<sup>th</sup> child.
	 * 
	 * @param i index of the child to return.
	 * @return The i<sup>th</sup> child.
	 */
	virtual Feature * getChild(size_t i) const ;
	
	
	/** \brief Return the father.
	 * 
	 * @return thefather Feature. It can be NULL !
	 */
	virtual Feature * getFather() const;
	
	
	/** \brief Return the children
	 * 
	 * @return a pointer to the m_c member. This is a vector containing the pointers to the children.
	 */
	virtual const std::vector<Feature *> & getChildren() const;

	/** \brief Return the children
	 * 
	 * @return a pointer to the m_c member. This is a vector containing the pointers to the children.
	 */
	virtual std::vector<Feature *> & getChildren();

	/** \brief Return the children and their children recursively
	 * 
	 * @return all descendants
	 */
	virtual std::vector<Feature *> getDescendants() const ;
	
	/** \brief Reparent the feature.
	 * 
	 * @param f the new father. The feature is automatically added to the children of the new parent.
	 */
	virtual void setFather(Feature *f) ;
	
	/** \brief Return all the triangle <em>in</em> the feature.
	 * 
	 * The DealunayTree::conflict() function returns all the triangles having <em>at least</em> on point in the boundary of the feature. This function returns only those triangles whose <em>center</em> lies <em>in</em> the feature.
	 * 
	 * @param dt the DelaunayTree from which to get the triangles.
	 * @return the vector of triangles satisfying the condition center \f$ \in \f$ Feature. 
	 */
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  = 0;
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons(DelaunayTree3D * dt)  = 0;

/** \brief return triangles intersecting the feature*/
	virtual std::vector<DelaunayTriangle *> getBoundingTriangles( DelaunayTree * dt) ;

/** \brief return tetrahedrons intersecting the feature*/
	virtual std::vector<DelaunayTetrahedron *> getBoundingTetrahedrons(DelaunayTree3D * dt) ;
	
	/** \brief Check for interaction.
	 * 
	 * Two features are said to be interactiong if ther boundaries are intersecting.
	 * 
	 * @param f Feature for which to check for the interaction.
	 * @return true if interacting.
	 */
	virtual bool interacts(Feature * f) const = 0;
	
	
	/** \brief Insert a point on the bounding surface of the feature.
	 * 
	 * @param i index of the Point after which to do the insertion.
	 * @return a pointer to the Point just inserted.
	 */
	virtual Point * pointAfter(size_t i) = 0 ;
	
	
	/** \brief Define a vector of geometries to use for sucessive refinement.
	 * 
	 * @return the vector of geometries.
	 */
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const = 0 ;
	
	/** \brief Projects a point on the boundary of the feature.
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

/** \brief Class which induces behaviou in elements but does not affect the constitution of the mesh*/
class VirtualFeature : virtual public Feature
{
public:

	VirtualFeature(Feature *father) : Feature(father) { isVirtualFeature = true ;}
	
	VirtualFeature(Feature *father, Geometry * b) : Feature(father, b)  { isVirtualFeature = true ;}
	
	virtual ~VirtualFeature() { };

	virtual void print() const = 0 ;
	virtual Feature * getSource() = 0;
} ;

/** \brief Feature composed of sub-features*/
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

/** \brief Feature which has no effect on the mesh, but enriches elements interacting with it*/
class EnrichmentFeature : virtual public Feature
{
public:
	
	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 */
	EnrichmentFeature(Feature *father) : Feature(father) { this->isEnrichmentFeature = true ;}
	
	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 * @param b      sets the boundary.
	 */
	EnrichmentFeature(Feature *father, Geometry * b) : Feature(father, b) { this->isEnrichmentFeature = true ;}
	
	virtual ~EnrichmentFeature() { };
	
	/** \brief Typically not used*/
	virtual std::vector<Point *> getSamplingPoints() const = 0;
	
	/** \brief return true if the argument should be enriched*/
	virtual bool enrichmentTarget(DelaunayTriangle * t) { return false ; };

	/** \brief return true if the argument should be enriched*/
	virtual bool enrichmentTarget(DelaunayTetrahedron * t) { return false ;};

	/** \brief enrich the mesh*/
	virtual void enrich(size_t & , DelaunayTree * dtree) { } ;

	/** \brief enrich the mesh*/
	virtual void enrich(size_t & , DelaunayTree3D * dtree) { } ;
	
	/** \brief update enrichment geometry*/
	virtual void step(double dt, Vector *, const DelaunayTree * dtree) { };

	/** \brief update enrichment geometry*/
	virtual void step(double dt, Vector *, const DelaunayTree3D * dtree) { };
	
	virtual void snap(DelaunayTree * dtree) { } ;
	virtual void snap(DelaunayTree3D * dtree) { } ;
	
	/** \brief return true if enrichment geometry has changed*/
	virtual bool moved() const = 0;
	
} ;

/** \brief Abstract boundary condition object for usage in multigrid solver. Work in Progress*/
class BoundaryCondition
{	
protected:
	LagrangeMultiplierType condition;
	std::vector<double> data ;

public:
	BoundaryCondition(LagrangeMultiplierType t, const std::vector<double> & d) ;
	virtual void apply(Assembly * a, DelaunayTree * t) const = 0 ;
	virtual void apply(Assembly * a, DelaunayTree3D * t) const = 0 ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class ProjectionDefinedBoundaryCondition : public BoundaryCondition
{
private:
	Point direction ;
	Point from ;

public:
	ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t, const std::vector<double> & d,const Point & direction, const Point & from) ;
	virtual void apply(Assembly * a, DelaunayTree * t) const ;
	virtual void apply(Assembly * a, DelaunayTree3D * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class GeometryDefinedBoundaryCondition : public BoundaryCondition
{
private:
	Geometry * domain ;

public:
	GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, const std::vector<double> & d, Geometry * source) ;
	virtual void apply(Assembly * a, DelaunayTree * t) const ;
	virtual void apply(Assembly * a, DelaunayTree3D * t)  const ;
} ;

/** \brief Boundary condition object for usage in multigrid solver. Work in Progress*/
class GeometryProjectedBoundaryCondition : public BoundaryCondition
{
private:
	Geometry * domain ;
	Point from ;
	Point direction ;
public:
	GeometryProjectedBoundaryCondition(LagrangeMultiplierType t, const std::vector<double> & d, Geometry * source, const Point & from,  const Point & direction ) ;
	virtual void apply(Assembly * a, DelaunayTree * t) const ;
	virtual void apply(Assembly * a, DelaunayTree3D * t)  const ;
} ;

/** \brief Member of an acess grid. Contains references to finer grid level and/or Feature s*/
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
	/** \brief default constructor*/
	Pixel();
	~Pixel();
	/** \brief Construct a pixel centred at x,y of size s
	 *
	 * @param x center x
	 * @param y center y
	 * @param s size
	*/
	Pixel(double x, double y, double s) ;

	/** \brief return features stored in this pixel*/
	const std::vector<Feature *> & getFeatures() const;
	
	/** \brief return features stored in this pixel*/
	std::vector<Feature *> & getFeatures();
	
	/** \brief return true if the argument lies in the pixel*/
	bool in(const Point & p) const;

	/** \brief return true if the Geometry overlaps the pixel*/
	bool coOccur(const Geometry * inc) const;
	
	/** \brief return true if the argument lies in the pixel*/
	bool coOccur(const Point & p) const;

	/** \brief Given a Geometry, return a list of stored Feature pointers overlapping the Geometry
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Geometry to check for overlaps
*/
	void coOccuringFeatures(std::vector<Feature *>& ret, const Geometry * inc) const ;

	/** \brief Given a Point, return a list of stored Feature pointers in which the point lies
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Point to check for overlaps
*/
	void coOccuringFeatures(std::vector<Feature *>& ret , const Point & p) const ;

/** \brief remove the argument from the Feature list*/
	void remove(Feature * inc);
	
/** \brief add the argument to the Feature list if it does not overlap with another already present Feature*/
	bool add(Feature * inc);

/** \brief add the argument to the Feature list unconditionnally*/
	void forceAdd(Feature * inc) ;

	void print() const ;

} ;


/** \brief Member of an acess grid. Contains references to finer grid level and/or Feature s*/
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
	/** \brief default constructor*/
	Voxel();

		/** \brief Construct a voxel centred at x,y,z of size s
	 *
	 * @param x center x
	 * @param y center y
	 * @param z center z
	 * @param s size
	*/
	Voxel(double x, double y, double z ,double s) ;

	~Voxel();

	/** \brief return features stored in this voxel*/
	const std::vector<Feature *> & getFeatures() const;
	
	/** \brief return features stored in this voxel*/
	std::vector<Feature *> & getFeatures();
	
	/** \brief return true if the argument lies in the voxel*/
	bool in(const Point & p) const;

	/** \brief return true if the Geometry overlaps the voxel*/
	bool coOccur(const Geometry * inc) const;

	/** \brief return true if the argument lies in the voxel*/
	bool coOccur(const Point & p) const;

	/** \brief Given a Geometry, return a list of stored Feature pointers overlapping the Geometry
*
* The search goes recursively through all the sub-voxels
* @param  ret Vector in which the result should be stored
* @param inc Geometry to check for overlaps
*/
	void coOccuringFeatures(std::vector<Feature *>&, const Geometry * inc) const ;

	/** \brief Given a Point, return a list of stored Feature pointers in which the point lies
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Point to check for overlaps
*/
	void coOccuringFeatures(std::vector<Feature *>&, const Point & p) const ;

/** \brief remove the argument from the Feature list*/
	void remove(Feature * inc);
	
/** \brief add the argument to the Feature list if it does not overlap with another already present Feature*/
	bool add(Feature * inc);

/** \brief add the argument to the Feature list unconditionnally*/
	void forceAdd(Feature * inc) ;

	void print() const ;
	
	bool isFilled() const ;
	
	Point center() const ;

} ;

/** \brief access grid for Features*/
class Grid
{
protected:
	std::valarray< std::valarray<Pixel *> > pixels;
	double x ;
	double y ;
	Point c ;
	size_t lengthX ;
	size_t lengthY ;
	
	double psize ;
public:
		
/** \brief Copnstruct a grid from a size, an initial number of divisions and a center
*
* @param sizeX Length of the access grid
* @param sizeY width of the access grid
* @param div maximum number of spacial divisions to use
* @param center center of the grid
 */
	Grid(double sizeX, double sizeY, int div, const Point & center );
	
	~Grid() ;

/** \brief Add a Feature if it does not overlap with another allready present Feature
*
* @param inc Feature to add
* @return true if insertion was successful
*/
	bool add(Feature * inc);

/** \brief Add a Feature unconditionnally*/
	void forceAdd(Feature * inc) ;

/** \brief Return the lis of Features overlapping the argument*/
	std::vector<Feature *> coOccur(const Geometry * geo) const ;

/** \brief Return the list of Features containing the argument*/
	std::vector<Feature *> coOccur(const Point & p) const ;

/** \brief Get a new grid, given a number of divisions*/
	Grid getGrid(int div) const;
} ;

/** \brief access grid for Features*/
class Grid3D
{
protected:
	std::valarray<std::valarray< std::valarray<Voxel *> > > pixels;
	std::vector<Voxel *> unfilledpixel ;
	double x ;
	double y ;
	double z ;
	Point c ;
	size_t lengthX ;
	size_t lengthY ;
	size_t lengthZ ;
	
	double psize ;
	int dirtyCounter ;
public:
		
/** \brief Copnstruct a grid from a size, an initial number of divisions and a center
*
* @param sizeX Length of the access grid
* @param sizeY width of the access grid
* @param sizeZ breadth of the access grid
* @param div maximum number of spacial divisions to use
* @param center center of the grid
 */
	Grid3D(double sizeX, double sizeY, double sizeZ, int div, const Point & center );

/** \brief Return a Point contained in no allready placed Feature in the grid*/
	Point randomFreeCenter() const ;
	~Grid3D() ;

/** \brief Add a Feature if it does not overlap with another allready present Feature
*
* @param inc Feature to add
* @return true if insertion was successful
*/
	bool add(Feature * inc);

/** \brief Add a Feature unconditionnally*/
	void forceAdd(Feature * inc) ;

/** \brief Return the lis of Features overlapping the argument*/
	std::vector<Feature *> coOccur(const Geometry * geo) const ;

/** \brief Return the list of Features containing the argument*/
	std::vector<Feature *> coOccur(const Point & p) const;
	double fraction() const ;

/** \brief Get a new grid, given a number of divisions*/
	Grid3D getGrid(int div) const ;
} ;

/** \brief Container for the features defining the setup.
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
	/** \brief Contains all the features. */
	std::vector<Feature *> tree ;

	/** \brief For fast Access*/
	Grid * grid ;
	Grid3D * grid3d ;
	
	/** \brief  Contains the mesh in the form of a delaunay tree. 
	 * The mesh is generated with linear triangles, and when it is final, midpoints are added and 
	 * projected. No operations should add midpoints before meshing is complete.
	 */
	DelaunayTree * dtree ;
	DelaunayTree3D * dtree3D ;
	
	TetrahedralElement *father3D  ;
	TriElement *father2D  ;

	bool meshChange ;
	bool solverConvergence ;
	bool enrichmentChange ;
	bool setBehaviours ;
	
	/** \brief  List of points used for the mesh.
	 * Each point is associated with the feature from whose discretiation it was generated.
	 */
	std::deque<std::pair<Point *, Feature *> > meshPoints;
	std::vector<Point *> additionalPoints ;
	
	/** \brief  List of the elements.
	 */
	std::vector<DelaunayTetrahedron * > elements3D;
	
	/** \brief  Assembly used for the generation of the stiffness matrix and the solving of the problem.
	 */
	Assembly * K ;
	
	//Assembly2D
	/** \brief  Order to which Elements should be brought once all operations are accomplished.
	 */
	Order elemOrder ;
	
	/** \brief  Insert triangle midpoints. */
	void makeToOrder() ;
	
	/** \brief  Project all points on their respective boundaries.*/
	void stitch() ;
	
	void renumber() ;

	void setElementBehaviours() ;
public:
		/** \brief  Generate the triangulation.
	 * Once the sampling is done, the sampling points are fed into a Delaunay Triangulation algorithm, 
	 * which generates the triangles. The mesh is still composed of 3 points triangles at this point.
	 * 
	 * @param correctionSteps additional steps where points are inserted in incorrect tetrahedrons.
	 */
	void generateElements( size_t correctionSteps = 0, bool computeIntersections = true) ;
	
	/** \brief  Perform the assembly.
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

	double crackedVolume ;
	double damagedVolume ;
	
public:
	
	/** \brief  Initialise the feature tree with the parent feature.
	 * 
	 * @param first Parent of all features. This shoud typically be the sample itself.
	 * @return 
	 */
	FeatureTree(Feature *first) ;
	virtual ~FeatureTree() ;
	
	/** \brief  Add feature as the daughter of another.
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
	
	/** \brief  Generate the sample points for all the features. The features are passed a sampling 
	 * argument proportionnal to their area compared with the area of the root feature. 
	 * If the number is lower than 10, than the argument passed is 10.
	 * 
	 * @param npoints number of sampling points on the boundary of the sample. This number must be 
	 * divisable by four if the sample is a rectangle.
	 */
	void sample(size_t npoints) ;
	
	/** \brief  Attempt to enhance the mesh, based on a sampling citerion.
	 * 
	 * Criterions are typically made to check the adequacy of the geometry of the triangles generated 
	 * during the triangulation phase. They also provide a hint as to the placement of new points which could 
	 * enhance the mesh quality.
	 * 
	 * @param nit Number of refinement iterations to perform.
	 * @param cri Citerion to use.
	 */
	void refine(size_t nit, SamplingCriterion *cri) ;
	
	
	/** \brief  Refine the mesh around the features.
	 * 
	 * Features provide a set of geometries which are targets for successive refinement. Refinement is 
	 * done by inserting a point in the center of each triangle in the zone.
	 * 
	 */
	void refine(size_t level) ;
	
	
	/** \brief  Set the constitutive law of the given feature.
	 * 
	 * <b>Beware!</b> this function <i>deletes</i> the previous tensor. Make sure is is not in use by another feature. 
	 * 
	 * @param m pointer to the tensor. 
	 * @param f Feature to affect.
	 */
	void setStrainTensor(Matrix * m, Feature * f) ;
	
	/** \brief  set the target order of Elements
	 * 
	 * @param ord order of elements to use.
	 * 
	 * \todo Make it work for values other than 1 and 2.
	 */
	void setOrder(Order ord) ;
	
	/**  \brief  Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * strain values at the mesh points.
	 * 
	 * @return strain.
	 */
	Vector strainFromDisplacements() const ;
	
	/**  \brief  Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * stress values at the mesh points.
	 * @return stress.
	 */
	Vector stressFromDisplacements() const ;
	
/** \brief Return the stress and strain of a vector of Tetrahedrons*/
	std::pair<Vector , Vector > getStressAndStrain(const std::vector<DelaunayTetrahedron *> &) ;

/** \brief Return the stress and strain of the elements of the current mesh*/
	std::pair<Vector , Vector > getStressAndStrain() ;
	
	size_t numPoints() const ;
	
/** \brief Step in time
* @param dt timestep
*/
	bool step(double dt) ;

/** \brief annul the last timestep*/
	void stepBack() ;

/** \brief Perform a time step, but do not update the features*/
	void elasticStep() ;

	std::deque<std::pair<Point *, Feature *> >::iterator begin() ;
	std::deque<std::pair<Point *, Feature *> >::iterator end() ;
	
/** \brief Return the Assembly*/
	Assembly * getAssembly() ;

/** \brief return the triangles of the mesh*/
	std::vector<DelaunayTriangle *> getTriangles();

/** \brief return the tetrahedrons of the mesh*/
	std::vector<DelaunayTetrahedron *> getTetrahedrons() ;
		
/** \brief return the triangles lying next to a mesh border*/
	std::vector<DelaunayTriangle *> getBoundingTriangles(Feature * f = NULL) ;	
	
/** \brief return the Behaviour of the argument, deduced from the Feature s*/
	Form * getElementBehaviour(const DelaunayTriangle *) const ;

/** \brief return the Behaviour of the argument, deduced from the Feature s*/
	Form * getElementBehaviour(const DelaunayTetrahedron *) const ;
	
/** \brief insert a point in the mesh*/
	void insert(Point * p ) ;
	
/** \brief return the 2D mesh*/
	DelaunayTree * getDelaunayTree() ;

/** \brief return the 3D mesh*/
	DelaunayTree3D * getDelaunayTree3D() ;
	
/** \brief Return true id the argument lies in the root feature*/
	bool inRoot(const Point &p) const ;
	
	Feature * getFeature(size_t i)
	{
		return tree[i] ;
	}
	
/** \brief initialise the element states*/
	void initializeElements() ;
	
	double getMaximumDisplacement() const ;
	double getMinimumDisplacement() const ;
	
/** \brief return true if the currently defined problem is 3D*/
	bool is3D() const ;

/** \brief return true if the currently defined problem is 3D*/
	bool is2D() const ;

	std::vector<DelaunayTriangle> getSnapshot2D() const ;

} ;


/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureEqual
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first == *p2.first ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first < *p2.first ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_x
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->x < p2.first->x ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_y
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->y < p2.first->y ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_z
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->z < p2.first->z ;
	}
} ;

} ;

#endif // __FEATURES_H_
