//
// C++ Interface: inclusion3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __INCLUSION3D_H__
#define __INCLUSION3D_H__

#include "features.h"
#include "../utilities/xml.h"

namespace Mu
{

/** \brief 3D spherical inclusion*/
class Inclusion3D :  public Sphere,  public Feature
{
	friend class Sphere ;
public:

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	Inclusion3D(Feature *father, double radius, double x, double y, double z) ;

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	Inclusion3D(Feature *father, double radius,  Point center) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	Inclusion3D(double radius, double x, double y, double z) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	Inclusion3D(double radius, Point center) ;
	
	virtual XMLTree * toXML() ;

/** \brief do nothing */
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle, DelaunayTreeItem> * dt)  ;
	
/** \brief return all tets in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dt)  ;
	
/** \brief do nothing */
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	virtual void sample(size_t n) ;
	
} ;


/** \brief Regular octahedron inclusion*/
class OctahedralInclusion :  public RegularOctahedron,  public Feature
{
	friend class RegularOctahedron ;
public:

/** \brief Construct inclusion from length and center position
*
* @param father Father feature 
* @param length of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	OctahedralInclusion(Feature *father, double radius, double x, double y, double z) ;

/** \brief Construct inclusion from length and center position
*
* @param father Father feature 
* @param length of the inclusion
* @param center center
*/
	OctahedralInclusion(Feature *father, double radius,  Point center) ;

/** \brief Construct inclusion from length and center position
*
* @param length of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	OctahedralInclusion(double radius, double x, double y, double z) ;

/** \brief Construct inclusion from length size and center position
*
* @param length of the inclusion
* @param center center
*/
	OctahedralInclusion(double length, Point center) ;
	
/** \brief do nothing */
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
/** \brief return all tets in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt)  ;
		
/** \brief do nothing */
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(RegularOctahedron) ;
	
	virtual void sample(size_t n) ;
	
} ;

/** \brief Inclusion of spherical shape used only to attribute behaviour. This feature has no effect on the mesh*/
class VirtualInclusion3D :  public Sphere,  public VirtualFeature
{
	friend class Sphere ;
public:

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	VirtualInclusion3D(Feature *father, double radius, double x, double y, double z) ;

/** \brief Construct inclusion from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	VirtualInclusion3D(Feature *father, double radius,  Point center) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	VirtualInclusion3D(double radius, double x, double y, double z) ;

/** \brief Construct inclusion from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	VirtualInclusion3D(double radius, Point center) ;
	
/** \brief do nothing*/
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return false*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return empty vector*/
	virtual std::vector<Geometry *> getRefinementZones(size_t ) const ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTriangle *> getTriangles( DelaunayTree * dt)  ;
	
/** \brief return all tets in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTetrahedron *> getTetrahedrons( DelaunayTree3D * dt)  ;
		
/** \brief return tthis*/
	virtual Feature * getSource() ;

/** \brief return false */
	virtual bool inBoundary(const Point *) const ;
	
/** \brief return false*/
	virtual bool inBoundary(const Point) const ;
	
/** \brief do nothing, return NULL*/
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am an inclusion" << std::endl ;
	}
	
/** \brief return false*/
	virtual bool isVoid( const Point &) const {return false ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	virtual void sample(size_t n) ;
	
} ;

} ;

#endif
