//
// C++ Interface: pore3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __PORE3D_H__
#define __PORE3D_H__

#include "features.h"

namespace Mu
{

/** \brief 3D spherical pore. Has a VoidBehaviour by default*/
class Pore3D : public Feature, public Sphere
{
public:

/** \brief Construct pore from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	Pore3D(Feature *father, double radius, double x, double y, double z) ;

/** \brief Construct pore from radius and center position
*
* @param father Father feature 
* @param radius of the inclusion
* @param center center
*/
	Pore3D(Feature *father, double radius, Point center) ;

/** \brief Construct pore from radius and center position
*
* @param radius of the inclusion
* @param x center x
* @param y center y
* @param z center z
*/
	Pore3D(double radius, double x, double y, double z) ;

/** \brief Construct pore from radius and center position
*
* @param radius of the inclusion
* @param center center
*/
	Pore3D(double radius, Point center) ;
	
/** \brief do nothing */
	virtual void addSamplePoints(PointSet * po ) { };
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief get list of refinement zones*/
	virtual std::vector<Geometry *> getRefinementZones(size_t) const ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle, DelaunayTreeItem> * dt)  ;
	
/** \brief return all tets in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dt)  ;
	
/** \brief do nothing*/
	virtual Point * pointAfter(size_t i) ;
	
	virtual void computeCenter()
	{
		return this->Sphere::computeCenter() ;
	}
	
	virtual void print() const
	{
		std::cout << "I am a pore" << std::endl ;
	}
	
	virtual bool isVoid( const Point & p) const {return squareDist(p, getCenter()) < 0.90*radius ;}
	
	
public:
	
	GEO_DERIVED_OBJECT(Sphere) ;
	
	/** Sample the bounding Surface
	 * 
	 * @param n number of points for the sampling.
	 */
	virtual void sample(size_t n) ;
	
	
} ;











} ;

#endif
