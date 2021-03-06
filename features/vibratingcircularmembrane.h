
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __VIBRATING_MEMBRANE_H__
#define __VIBRATING_MEMBRANE_H__

#include "enrichmentInclusion.h"

namespace Amie
{

/** \brief Vibrating zone. Incomplete
 * This enrichement feature will introduce a 
 * soft discontinuity in the mesh, as well as 
 * attribute a vibrating elastic behaviour
 * to the material contained within
*/
class VibratingMembrane :  public EnrichmentInclusion
{

	Matrix cgTensor ;
public:

	VibratingMembrane(Feature *father, double radius, double x, double y, Matrix cgTensor) ;
	virtual ~VibratingMembrane() ;
	
	virtual void enrich(size_t &,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
	virtual std::vector<Amie::DelaunayTetrahedron*> getElements3D(FeatureTree *) { return std::vector<Amie::DelaunayTetrahedron*>() ;}
	virtual std::vector<Amie::DelaunayTriangle*> getElements2D(FeatureTree *) { return std::vector<Amie::DelaunayTriangle*>() ;}
	
	virtual void print() const
	{
		std::cout << "I am an expansive zone" << std::endl ;
	}
	
	void reset() ;
	
public:
	GEO_DERIVED_OBJECT(Circle) ;
	
} ;

}

#endif
