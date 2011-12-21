
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __EXPANSIVE_ZONE_3D_H__
#define __EXPANSIVE_ZONE_3D_H__

#include "enrichmentInclusion3d.h"
#include "../physics/stiffness_with_imposed_deformation.h"

namespace Mu
{

/** \brief Expansive zone. 
 *
 * This enrichement feature will introduce a 
 * soft discontinuity in the mesh, as well as 
 * attribute an imposed-strain elastic behaviour
 * to the material contained within
*/
class ExpansiveZone3D :  public EnrichmentInclusion3D
{
	std::set<DelaunayTetrahedron *> bimateralInterfaced ;
	std::set<DelaunayTetrahedron *> expansive ;
	Vector imposedDef ;
	Matrix cgTensor ;
	Form * original ;
	
public:

/** \brief Constructor. construct the zone
*
* @param father supporting feature. This features behaviour is used for the bimaterial interfaces
* @param radius External radius
* @param x center x
* @param y center y
* @param z center z
* @param cgTensor Stifness of the expensive material of the ring, This is the complete experssion of the CG stress tensor
* @param deformation Vector of the imposed strain
*/
	ExpansiveZone3D(Feature *father, double radius, double x, double y, double z, const Matrix & cgTensor, Vector deformation) ;
	ExpansiveZone3D(Feature *father, double radius, double x, double y, double z, StiffnessWithImposedDeformation * exp) ;
	virtual ~ExpansiveZone3D() ;
	
/** \brief enrich elements and change their Behaviour if required*/
	virtual void enrich(size_t &,  Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) ;
	
	virtual void print() const
	{
		std::cout << "I am an expansive zone" << std::endl ;
	}
	
	void reset() ;

	virtual std::vector<Mu::DelaunayTetrahedron*> getElements3D(FeatureTree*) { return std::vector<Mu::DelaunayTetrahedron*>() ;}
	virtual std::vector<Mu::DelaunayTriangle*> getElements2D(FeatureTree*) { return std::vector<Mu::DelaunayTriangle*>() ;}
	
public:
	GEO_DERIVED_OBJECT(Sphere) ;
	
} ;

}

#endif
