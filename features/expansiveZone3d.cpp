// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveZone3d.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"

using namespace Mu ;

ExpansiveZone3D::ExpansiveZone3D(Feature *father, double radius, double x, double y, double z, const Matrix & tensor, Vector def) : EnrichmentInclusion3D(father, radius, x, y, z),  imposedDef(def),cgTensor(tensor)
{
	
}

ExpansiveZone3D::~ExpansiveZone3D() {}
	
void ExpansiveZone3D::reset() 
{
	cache.clear() ;
	updated = true ;
}

void ExpansiveZone3D::enrich(size_t&lastId, Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* dtree)
{
	EnrichmentInclusion3D::enrich(lastId, dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTetrahedron *> & disc = EnrichmentInclusion3D::cache ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTetrahedron *> ring ;
	std::vector<DelaunayTetrahedron *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(EnrichmentInclusion3D::enrichmentTarget(disc[i]))
			ring.push_back(disc[i]) ;
		else if(this->in(disc[i]->getCenter()))
			inDisc.push_back(disc[i]) ;
	}
	std::set<DelaunayTetrahedron *> newInterface ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(bimateralInterfaced.find(ring[i]) == bimateralInterfaced.end())
		{
			ring[i]->setBehaviour(new BimaterialInterface(getPrimitive(),
					      new StiffnessWithImposedDeformation(cgTensor, imposedDef),
					      ring[i]->getBehaviour()->getCopy()
					      )) ;
			ring[i]->getBehaviour()->transform(ring[i]->getXTransform(), ring[i]->getYTransform(), ring[i]->getZTransform()) ;
		}
		newInterface.insert(ring[i]) ;
	}
	
	std::set<DelaunayTetrahedron *> newExpansive ;
	for(size_t i = 0 ; i < inDisc.size() ; i++)
	{
		if(expansive.find(inDisc[i]) == expansive.end())
			inDisc[i]->setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
		
		newExpansive.insert(inDisc[i]) ;
	}
	expansive = newExpansive ;
	
	if(disc.size() == 1)
	{
		if(bimateralInterfaced.find(disc[0]) == bimateralInterfaced.end())
		{
			disc[0]->setBehaviour(new BimaterialInterface(getPrimitive(),
					      new StiffnessWithImposedDeformation(cgTensor, imposedDef),
					      disc[0]->getBehaviour()->getCopy()
					      )) ;
			disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform(), disc[0]->getZTransform()) ;
		}
		newInterface.insert(disc[0]) ;
	}
	
	bimateralInterfaced = newInterface ;
	
}
	

