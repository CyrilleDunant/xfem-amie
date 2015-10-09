//
// C++ Interface: inclusion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
// Author: Cyrille Dunant <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCLUSION_H__
#define __INCLUSION_H__

#include "features.h"
#include "../geometry/level_set.h" 

namespace Amie
{

/** \brief Circular inclusion*/
class Inclusion :  virtual public Circle,   public Feature
{
public:
/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Inclusion(Feature *father, double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	Inclusion(Feature *father, double radius,  Point center) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
*/
	Inclusion(double radius, double x, double y) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	Inclusion(double radius, Point center) ;

	Inclusion(Circle c) ;

/** \brief do nothing */
	virtual void addSamplePoints(PointSet * po ) { };

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
/** \brief return false */
	virtual bool isVoid( const Point &) const {return false ;}
	
	virtual void setRadius(double newR) ;
	
	
public:
	
	GEO_DERIVED_OBJECT(Circle) ;

	virtual void sample(size_t n) ;
	
} ;

/** \brief Triangular inclusion*/
class TriangularInclusion :  public Triangle,  public Feature
{
public:

/** \brief construct a triangular inclusion from three points
*
* @param father father feature
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularInclusion(Feature *father, const Point & a, const Point & b, const Point & c) ;

/** \brief construct a triangular inclusion from three points
*
* @param a first vertex
* @param b second vertex
* @param c third vertex
*/
	TriangularInclusion(const Point & a, const Point & b, const Point & c) ;


/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
		
	virtual void print() const
	{
		std::cout << "I am a triangular inclusion" << std::endl ;
	}
	
/** \brief return false */
	virtual bool isVoid(const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Triangle) ;

	virtual void sample(size_t n) ;
	
} ;

/** \brief Ellipsoidal inclusion*/
class EllipsoidalInclusion :  virtual public Ellipse,   public Feature
{
public:
/** \brief construct an ellipsoidal inclusion (see Ellipse class for Ellipse constructors)
*
* @param father father feature
*/
	EllipsoidalInclusion(Feature *father, Point center, Point a, Point b) ;

/** \brief construct an ellipsoidal inclusion (see Ellipse class for Ellipse constructors)
*
* @param father father feature
*/
	EllipsoidalInclusion(Feature *father, Point center, Point a, double b) ;

/** \brief construct an ellipsoidal inclusion (see Ellipse class for Ellipse constructors) */
	EllipsoidalInclusion(Point center, Point a, Point b) ;

/** \brief construct an ellipsoidal inclusion (see Ellipse class for Ellipse constructors) */
	EllipsoidalInclusion(Point center, Point a, double b) ;

/** \brief Do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;

/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
		
	virtual void print() const
	{
		std::cout << "I am an ellipsoidal inclusion" << std::endl ;
	}
	
/** \brief return false */
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Ellipse) ;

	virtual void sample(size_t n) ;
	
} ;

/** \brief Triangular inclusion*/
class RectangularInclusion :  public OrientedRectangle,  public Feature
{
public:
	
	/** \brief construct a rectangular inclusion from three points
	 * 
	 * @param father father feature
	 * @param a first vertex
	 * @param b second vertex
	 * @param c third vertex
	 * @param d fourth vertex
	 */
	RectangularInclusion(Feature *father, const Point & a, const Point & b, const Point & c, const Point & d) ;
	
	/** \brief construct a rectangular inclusion from three points
	 * 
	 * @param a first vertex
	 * @param b second vertex
	 * @param c third vertex
	 * @param d fourth vertex
	 */
	RectangularInclusion(const Point & a, const Point & b, const Point & c, const Point & d) ;
	
	
	/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
	/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
	/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;
	
	/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 
	
	virtual void print() const
	{
		std::cout << "I am a rectangular inclusion" << std::endl ;
	}
	
	/** \brief return false */
	virtual bool isVoid(const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(OrientedRectangle) ;
	
	virtual void sample(size_t n) ;
	
} ;


/** \brief Interface Transition Zone*/
class ITZFeature : public LevelSet, public VirtualFeature
{
protected:
	double length ;
	double min ;
	double max ;
	std::vector<Feature *> sources ;
	virtual void computeCenter() {} ;

public:
	ITZFeature(Feature *father, Feature * g, double l, double a = 0., double b = 1.) ;
	ITZFeature(Feature *father, std::vector<Inclusion *> & g, double l, double a = 0., double b = 1.) ;

	virtual void print() {} ;
	virtual Form * getBehaviour( const Point & p ) const ;
	virtual Feature * getSource() {return sources[0] ; } ;
	virtual const Feature * getSource() const {return sources[0] ; } ;
	virtual bool in( const Point &) const ;

	double getLength() const {return length ; } ;
	void setLength(double l) {length = l ;} ;

//	virtual std::vector<DelaunayTriangle*> getTriangles(Mu::DelaunayTree*) { return std::vector<Mu::DelaunayTriangle*>() ;}
	
//	virtual std::vector<DelaunayTetrahedron*> getTetrahedrons(Mu::DelaunayTree3D*) {return std::vector<Mu::DelaunayTetrahedron*>() ;}
 	virtual bool interacts(Amie::Feature*, double) const {return false ;}
	virtual Point* pointAfter(size_t) {return nullptr ;}
	virtual std::vector<Amie::Geometry*> getRefinementZones(size_t) const {return std::vector<Amie::Geometry*>() ;}
	virtual void print() const {std::cout << "ITZ !" << std::endl;}
	virtual void sample(size_t) {} ;
	virtual bool isVoid(const Amie::Point&) const {return false ;}

	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  { return std::vector<Amie::DelaunayTriangle*>() ;} ;
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt) {return std::vector<DelaunayTetrahedron *>(0) ;} 

public:

   LEVEL_SET_DERIVED_OBJECT() ;


} ;



} 

#endif
