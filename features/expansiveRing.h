
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __EXPANSIVE_RING_H__
#define __EXPANSIVE_RING_H__

#include "enrichmentRing.h"

namespace Mu
{

/** \brief Expansive Ring. 
 * This enrichement feature will introduce a 
 * soft discontinuity in the mesh, as well as 
 * attribute an imposed-strain elastic behaviour
 * to the material contained within
*/
class ExpansiveRing :  public EnrichmentRing
{

	Vector imposedDef ;
	Matrix cgTensor ;
	std::set<DelaunayTriangle *> interfaced ;
	std::set<DelaunayTriangle *> expansive ;
	std::set<DelaunayTriangle *> biInterfaced ;
public:

/** \brief Constructor. construct the ring
*
* @param father supporting feature. This features behaviour is used for the bimaterial interfaces
* @param radius External radius
* @param inRadius Internal radius
* @param x center x
* @param y center y
* @param cgTensor Stifness of the expensive material of the ring, This is the complete experssion of the CG stress tensor
* @param deformation Vector of the imposed strain
*/
	ExpansiveRing(Feature *father, double radius, double inradius, double x, double y, const Matrix & cgTensor, Vector deformation) ;
	virtual ~ExpansiveRing() ;
	
/** \brief enrich elements and change their Behaviour if required*/
	virtual void enrich(size_t &,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) ;
	
/** \brief return empty vector*/
	virtual std::vector<Mu::DelaunayTetrahedron*> getElements3D(FeatureTree*) { return std::vector<Mu::DelaunayTetrahedron*>() ;}
	virtual std::vector<Mu::DelaunayTriangle*> getElements2D(FeatureTree*) { return std::vector<Mu::DelaunayTriangle*>() ;}
	
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
