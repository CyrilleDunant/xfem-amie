//
// C++ Interface: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LAYERED_INCLUSION_H__
#define __LAYERED_INCLUSION_H__

#include "features.h"

namespace Mu
{

/** \brief Layered Inclusion */
class LayeredInclusion :public LayeredCircle, virtual public CompositeFeature
{
protected:
	std::vector<Form *> layeredBehaviour ;
public:
	/** \brief Construct an Inclusion.
	 * 
	 * @param father father of the current Feature 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	LayeredInclusion(Feature *father, std::vector<double> radii,double originX,double originY) ;

	/** \brief Construct an Inclusion.
	 * 
	 * @param father father of the current Feature 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param center center.
	 */
	LayeredInclusion(Feature *father,std::vector<double> radii, const Point center) ;

	/** \brief Construct an Inclusion.
	 * 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param originX global x coordinate of the center.
	 * @param originY global y coordinate of the center.
	 */
	LayeredInclusion(std::vector<double> radii,double originX,double originY) ;

	/** \brief Construct an Inclusion.
	 * 
	 * @param radii radii of the circles to construct. Will be reordered.
	 * @param center center.
	 */
	LayeredInclusion(std::vector<double> radii,const Point center) ;

	/** \brief Construct an Inclusion.
	 * 
	 * @param r radius of the circle to construct.
	 * @param center center.
	 */
	LayeredInclusion(double r,const Point center) ;
	
/** \brief do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in the mesh with at least one vertex in the outermose circle*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle, DelaunayTreeItem> * dt);

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
/** \brief Do nothing*/
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const ;
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}

/** \brief return the behaviour of the layer in which the argument lies*/
	virtual Form * getBehaviour(const Point & p) ;

/** \brief set all layers behaviour to the given Behaviour*/
	virtual void setBehaviour(Form * b) ;

/** \brief set the layer behaviour*/
	virtual void setBehaviours(std::vector<Form *> b) ;


	
public:
	
	GEO_DERIVED_OBJECT(LayeredCircle) ;

	virtual void sample(size_t n) ;
	
} ;

/** \brief Helper class for the LayeredInclusion*/
class VirtualLayer :public Circle, virtual public VirtualFeature
{
	LayeredInclusion * source ;
public:
	VirtualLayer(LayeredInclusion *father, double radius, double x, double y) ;
	VirtualLayer(LayeredInclusion *father, double radius,Point center) ;
	virtual ~VirtualLayer() { } ;
	virtual void addSamplePoints(PointSet * po ) { };

	virtual bool interacts(Feature * f, double d) const ;
	
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle, DelaunayTreeItem> * dt);
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dt) ;
	
	virtual void computeCenter()
	{
		return this->Circle::computeCenter() ;
	}
	
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const ;
	
	virtual Form * getBehaviour(const Point & p) ;

	virtual bool isVoid( const Point &) const {return false ;}

	virtual Feature * getSource() ;
	
	
public:
	
	virtual std::vector<Point> getSamplingBoundingPoints(size_t num_points) const
	{
		return this->Circle::getSamplingBoundingPoints(num_points) ;
	}
	virtual const PointArray & getBoundingPoints() const 
	{
		return this->source->getBoundingPoints() ;
	}
	virtual PointArray & getBoundingPoints()
	{
		return this->source->getBoundingPoints() ;
	}
	virtual GeometryType getGeometryType() const
	{
		return this->Circle::gType ;
	}
	virtual const Point & getBoundingPoint(size_t i) const 
	{
		return this->source->getBoundingPoint(i) ;
	}
	virtual Point & getBoundingPoint(size_t i)
	{
		return this->source->getBoundingPoint(i) ;
	}
	virtual std::vector<Point> getBoundingBox() const
	{
		return this->Circle::getBoundingBox() ;
	}
	virtual void setBoundingPoint(size_t i, Point * p)
	{
		this->Circle::setBoundingPoint(i,p) ; 
	}
	virtual void setBoundingPoints(const PointArray & nb)
	{
		this->Circle::setBoundingPoints(nb) ; 
	}
	virtual void setInPoints(const PointArray & nb)
	{
		this->Circle::setInPoints(nb) ;
	}
	virtual void project(Point * p) const
	{
		this->Circle::project(p) ;
	}
	virtual double getRadius() const
	{
		return this->Circle::getRadius() ; 
	}
	virtual void sampleBoundingSurface(size_t n)
	{
		this->Circle::sampleBoundingSurface(n) ;
	}
	virtual void sampleSurface(size_t n) 
	{
		this->Circle::sampleSurface(n) ;
	}
	virtual SpaceDimensionality spaceDimensions() const 
	{
		return this->Circle::spaceDimensions() ;
	}
	virtual bool intersects(const Geometry * g) const
	{
		return this->Circle::intersects(g) ;
	}
	virtual std::vector<Point> intersection(const Geometry * g) const
	{
		return this->Circle::intersection(g) ;
	}
	virtual bool in(const Point & p) const
	{
		return this->Circle::in(p) ;
	}
	const Point & getPoint(size_t i) const
	{
		return this->Circle::getPoint(i) ; 
	}
	Point & getPoint(size_t i) 
	{
		return this->Circle::getPoint(i) ; 
	}
	virtual const Point &getInPoint(size_t i) const
	{
		return this->source->getInPoint(i) ;
	}
	virtual Point &getInPoint(size_t i)
	{
		return this->source->getInPoint(i) ;
	}
	virtual const std::valarray<Mu::Point*> & getInPoints() const
	{
		return this->source->getInPoints() ;
	}
	virtual std::valarray<Mu::Point*> & getInPoints() 
	{
		return this->source->getInPoints() ;
	}
	virtual size_t size() const
	{
		return this->source->size() ;
	}
	virtual size_t sides() const
	{
		return this->Circle::sides() ;
	}
	virtual Point & getCenter()
	{
		return this->Circle::getCenter() ;
	}
	virtual const Point & getCenter() const
	{
		return this->Circle::getCenter() ;
	}
	virtual double area() const
	{
		return this->Circle::area() ;
	}
	virtual double volume() const
	{
		return this->Circle::volume() ;
	}
	virtual void sample(size_t n) ;
	
	
} ;

} ;

#endif
