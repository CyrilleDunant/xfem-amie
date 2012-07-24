//
// C++ Interface: sample3d
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2006-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SAMPLE3D_H__
#define __SAMPLE3D_H__

#include "features.h"

namespace Mu
{

class Sample3D :  public Hexahedron,  public Feature
{
public:
	/** \brief Sample3D constructor.
	 * 
	 * @param father father feature
	 * @param x lenght
	 * @param y width
	 * @param z depth
	 * @param originX center x
	 * @param originY center y
	 * @param originZ center z
	 */
	Sample3D(Feature *father, double x, double y, double z, double originX, double originY, double originZ) ;

	/** \brief Sample3D constructor.
	 *
	 * The father feature is set to nullptr.
	 * 
	 * @param x lenght
	 * @param y width
	 * @param z depth
	 * @param originX center x
	 * @param originY center y
	 * @param originZ center z
	 */
	Sample3D(double x, double y, double z,double originX, double originY, double originZ) ;
	
/** \brief return true if the boundary overlaps that of the argument*/
	virtual bool interacts(Feature * f, double d) const ;
	
/** \brief return empty list*/
	virtual std::vector<Geometry *> getRefinementZones(size_t) const { return std::vector<Geometry *>() ;}
	
/** \brief return empty list*/
	virtual std::vector<DelaunayTriangle *> getElements2D( FeatureTree * dt)  ;
	
/** \brief return all tets in mesh with at least a vertex in this Feature*/
	virtual std::vector<DelaunayTetrahedron *> getElements3D( FeatureTree * dt)  ;
	
	virtual void print() const
	{
		std::cout << "I am a sample" << std::endl ;
	}

/** \brief return false */
	virtual bool isVoid( const Point & p) const {return false;}
	
public:
	
	GEO_DERIVED_OBJECT(Hexahedron) ;
	
	virtual void sample(size_t n)
	{
		this->sampleSurface(n) ;
	}
	
} ;


} ;

#endif
