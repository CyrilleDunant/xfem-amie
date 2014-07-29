//
// C++ Interface: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __LAYERED_INCLUSION_H__
#define __LAYERED_INCLUSION_H__

#include "features.h"

namespace Amie
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

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in the mesh with at least one vertex in the outermose circle*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt);

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	
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
	
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt);
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) ;
	
	virtual void print() const ;
	
	virtual Form * getBehaviour(const Point & p) ;

	virtual bool isVoid( const Point &) const {return false ;}

	virtual Feature * getSource() ;
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;

	virtual void sample(size_t n) ;
	
	
} ;

} ;

#endif
