
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

/** \brief Expansive zone. 
 *
 * This enrichement feature will introduce a 
 * soft discontinuity in the mesh, as well as 
 * attribute an imposed-strain elastic behaviour
 * to the material contained within
*/
class ExpansiveZone :  public EnrichmentInclusion
{
	std::set<DelaunayTriangle *> bimateralInterfaced ;
	std::set<DelaunayTriangle *> expansive ;
	Vector imposedDef ;
	Matrix cgTensor ;
public:

/** \brief Constructor. construct the zone
*
* @param father supporting feature. This features behaviour is used for the bimaterial interfaces
* @param radius External radius
* @param x center x
* @param y center y
* @param cgTensor Stifness of the expensive material of the ring, This is the complete experssion of the CG stress tensor
* @param deformation Vector of the imposed strain
*/
	ExpansiveZone(Feature *father, double radius, double x, double y, const Matrix & cgTensor, Vector deformation) ;
	virtual ~ExpansiveZone() ;
	
/** \brief enrich elements and change their Behaviour if required*/
	 virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree);
	
/** \brief return empty vector*/
	virtual std::vector<Mu::DelaunayTetrahedron*> getElements(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>*) { return std::vector<Mu::DelaunayTetrahedron*>() ;}
	virtual std::vector<Mu::DelaunayTriangle*> getElements(Mesh<DelaunayTriangle, DelaunayTreeItem>*) { return std::vector<Mu::DelaunayTriangle*>() ;}
	
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
