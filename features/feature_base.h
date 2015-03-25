// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef FEATURE_BASE_H
#define FEATURE_BASE_H
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_2D.h"
#include "../elements/integrable_entity.h"
#include "../mesher/mesh.h"
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "../polynomial/vm_function_base.h"


namespace Amie
{
	static const double DEFAULT_BOUNDARY = 1 ;
	static const size_t DEFAULT_POINTNUMBER = 1000 ;
	class DelaunayTriangle ;
	class DelaunayTetrahedron ;
	class FeatureTree ;
	
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
	Feature * behaviourSource ;
	
	/** \brief Constitutive Law describing the behaviour of the feature.
	 * 
	 */
	Form * behaviour ;
	
	double now ;
	
	int layer ;
	double fraction ;

public:
	
	bool isEnrichmentFeature ;
	bool isCompositeFeature ;
	bool isVirtualFeature ;
	bool isUpdated ;
	
	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 */
	Feature(Feature * const father = nullptr) ;

	/** \brief Feature constructor.
	 * 
	 * @param father sets the father.
	 * @param b      sets the boundary.
	 */
	Feature(Feature * const father, Geometry * const b) ;
	
	virtual ~Feature() ;
	
	Feature * getBehaviourSource() ;
	const Feature * getBehaviourSource() const ;
	void setBehaviourSource( Feature * const f) ;
	 
	bool inBoundary(const Point &p, double d) const ;
	bool onBoundary(const Point &p, double d) const ;
	
	int getLayer() const {return layer ;}
	void setLayer(int l) {layer = l;}
	double getFraction() const {return fraction ;}
	void setFraction(double f) {fraction = f;}
	
	/** \brief Set the Behaviour
	 * 
	 * @param b Behaviour of the feature. getCopy() from this behaviour will be called to generate the behaviour of elements depending from this feature.
	 */
	virtual void setBehaviour(Form * const b) ;
	
	/** \brief Get the Behaviour
	 * 
	 * @return the Behaviour
	 */
	virtual Form * getBehaviour( const Point & p ) const ;
	
	virtual Form * getBehaviour() const ;
	
	/** \brief Add a child to the feature.
	 * 
	 * @param f Children to add.
	 */
	virtual void addChild(Feature * const f) final;

	/** \brief Remove a child of the feature.
	 * 
	 * @param f Children to add.
	 */
	virtual void removeChild(Feature * const f) ;
	
	
	/** \brief Return i<sup>th</sup> child.
	 * 
	 * @param i index of the child to return.
	 * @return The i<sup>th</sup> child.
	 */
	virtual Feature * getChild(size_t i) const ;
	
	
	/** \brief Return the father.
	 * 
	 * @return thefather Feature. It can be nullptr !
	 */
	virtual const Feature * getFather() const;
	virtual  Feature * getFather() final;
	
	
	/** \brief Return the children
	 * 
	 * @return a pointer to the m_c member. This is a vector containing the pointers to the children.
	 */
	virtual const std::vector< Feature *> & getChildren() const;

	/** \brief Return the children
	 * 
	 * @return a pointer to the m_c member. This is a vector containing the pointers to the children.
	 */
	virtual std::vector< Feature *> & getChildren() final;

	virtual void addMeshPointsInFather() { } ;

	/** \brief Return the children and their children recursively
	 * 
	 * @return all descendants
	 */
	virtual std::vector<Feature *> getDescendants() const ;
	
	/** \brief Reparent the feature.
	 * 
	 * @param f the new father. The feature is automatically added to the children of the new parent.
	 */
	virtual void setFather(Feature * const f) ;
	
	/** \brief Return all the triangle <em>in</em> the feature.
	 * 
	 * The DealunayTree::conflict() function returns all the triangles having <em>at least</em> on point in the boundary of the feature. This function returns only those triangles whose <em>center</em> lies <em>in</em> the feature.
	 * 
	 * @param dt the DelaunayTree from which to get the triangles.
	 * @return the vector of triangles satisfying the condition center \f$ \in \f$ Feature. 
	 */
	virtual std::vector<DelaunayTriangle *> getElements2D(  FeatureTree* dt)  = 0;
	virtual std::vector<DelaunayTetrahedron *> getElements3D(  FeatureTree* dt)  = 0;

/** \brief return triangles intersecting the feature*/
	virtual std::vector<DelaunayTriangle *> getBoundingElements2D( FeatureTree* dt) const ;

/** \brief return tetrahedrons intersecting the feature*/
	virtual std::vector<DelaunayTetrahedron *> getBoundingElements3D(FeatureTree* dt) const ;
	
	/** \brief Check for interaction.
	 * 
	 * Two features are said to be interactiong if ther boundaries are intersecting.
	 * 
	 * @param f Feature for which to check for the interaction.
	 * @return true if interacting.
	 */
	virtual bool interacts(Feature * f, double d) const = 0;

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
	
	
	virtual void print() const { std::cout << "i'm a feature" << std::endl ; } ;
	
	virtual const Amie::PointArray & getBoundingPoints() const = 0 ;
	virtual Amie::PointArray & getBoundingPoints() = 0 ;
	virtual std::vector<Point *> doubleSurfaceSampling() ;
	virtual const Point & getBoundingPoint(size_t) const = 0 ;
	virtual Point & getBoundingPoint(size_t)  = 0 ;
	virtual Amie::PointArray & getInPoints() = 0 ;
	virtual const Amie::PointArray & getInPoints() const = 0 ;
	virtual const Point & getInPoint(size_t) const = 0 ;
	virtual Point & getInPoint(size_t) = 0 ;
	virtual bool in( const Point & ) const = 0 ;
	virtual const Point & getCenter() const = 0 ;
	virtual double getRadius() const = 0 ;
	virtual double area() const = 0 ;
	virtual void sample(size_t n) = 0 ;
	
	virtual bool isVoid( const Point &) const = 0 ;
	
} ;

/** \brief Class which induces behaviour in elements but does not affect the constitution of the mesh*/
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
	virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) { } ;

	/** \brief enrich the mesh*/
	virtual void enrich(size_t & , Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) { } ;
	
	/** \brief update enrichment geometry*/
	virtual void step(double dt, Vector *, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) { };

	/** \brief update enrichment geometry*/
	virtual void step(double dt, Vector *, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) { };
	
	virtual void update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) { } ;
	virtual void update(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) { } ;
	
	/** \brief return true if enrichment geometry has changed*/
	virtual bool moved() const = 0;
        
        virtual double grade() { return -1 ;} ;
	
} ;


	
} 

#endif //FEATURE_BASE_H
