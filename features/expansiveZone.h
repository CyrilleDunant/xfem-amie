
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __EXPANSIVE_ZONE_H__
#define __EXPANSIVE_ZONE_H__

#include "enrichmentInclusion.h"

namespace Mu
{

/** Expansive zone. 
 * This enrichement feature will introduce a 
 * soft discontinuity in the mesh, as well as 
 * attribute an imposed-strain elastic behaviour
 * to the material contained within
*/
class ExpansiveZone :  public EnrichmentInclusion
{

	Vector imposedDef ;
	Matrix cgTensor ;
public:

	ExpansiveZone(Feature *father, double radius, double x, double y, Matrix cgTensor, Vector deformation) ;
	virtual ~ExpansiveZone() ;
	
	virtual void enrich(size_t &,  DelaunayTree * dtree) ;
	
	virtual std::vector<Mu::DelaunayTetrahedron*> getTetrahedrons(Mu::DelaunayTree_3D*) { return std::vector<Mu::DelaunayTetrahedron*>() ;}
	
	virtual void print() const
	{
		std::cout << "I am an expansive zone" << std::endl ;
	}
	
	void reset() ;
	
	const Circle * getGeometry() const ;

	Circle * getGeometry() ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
} ;

}

#endif
