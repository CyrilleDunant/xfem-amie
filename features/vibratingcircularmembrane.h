
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __VIBRATING_MEMBRANE_H__
#define __VIBRATING_MEMBRANE_H__

#include "enrichmentInclusion.h"

namespace Mu
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
	
	virtual std::vector<Mu::DelaunayTetrahedron*> getElements3D(FeatureTree *) { return std::vector<Mu::DelaunayTetrahedron*>() ;}
	virtual std::vector<Mu::DelaunayTriangle*> getElements2D(FeatureTree *) { return std::vector<Mu::DelaunayTriangle*>() ;}
	
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
