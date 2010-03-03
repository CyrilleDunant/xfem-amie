//
// C++ Interface: sample
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2006-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include "features.h"

namespace Mu
{

/** \brief rectangular Feature with sides aligned to the axes*/
class Sample :  public Rectangle,  public Feature
{
public:

/** \brief Construct Sample from dimensions and center
*
* @param father Father feature 
* @param x size x
* @param y size y
* @param originX center x
* @param originY center y
*/
	Sample(Feature *father, double x, double y, double originX, double originY) ;

/** \brief Construct Sample from dimensions and center
*
* @param x size x
* @param y size y
* @param originX center x
* @param originY center y
*/
	Sample(double x, double y, double originX, double originY) ;

/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return empty list*/
	virtual std::vector<Geometry *> getRefinementZones(size_t) const { return std::vector<Geometry *>() ;}
	
/** \brief return all triangles in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTriangle *> getElements( Mesh<DelaunayTriangle> * dt)  ;
	
/** \brief return empty vector*/
	virtual std::vector<DelaunayTetrahedron *> getElements( Mesh<DelaunayTetrahedron> * dt) { return std::vector<DelaunayTetrahedron *>(0) ; }
	
	virtual void computeCenter()
	{
		return this->Rectangle::computeCenter() ;
	}
	
/** \brief do nothing */
	virtual Point * pointAfter(size_t i) ;
	
	virtual void print() const
	{
		std::cout << "I am a sample" << std::endl ;
	}

/** \brief return false */
	virtual bool isVoid( const Point &p) const {return false;}
	
public:
	
	GEO_DERIVED_OBJECT(Rectangle) ;
	
	virtual void sample(size_t n)
	{
		this->sampleSurface(n) ;
	}

} ;


} ;

#endif
