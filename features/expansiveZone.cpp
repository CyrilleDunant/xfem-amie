// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"

using namespace Mu ;

ExpansiveZone::ExpansiveZone(Feature *father, double radius, double x, double y, const Matrix & tensor, Vector def) : EnrichmentInclusion(father, radius, x, y),  imposedDef(def),cgTensor(tensor)
{
	
}

ExpansiveZone::~ExpansiveZone() {}
	
const Circle * ExpansiveZone::getGeometry() const 
{
	return static_cast<const Circle *>(this) ;
}

Circle * ExpansiveZone::getGeometry() 
{
	return static_cast<Circle *>(this) ;
}

void ExpansiveZone::reset() 
{
	cache.clear() ;
	updated = true ;
}

void ExpansiveZone::enrich(size_t & counter , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	EnrichmentInclusion::enrich(counter, dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = cache ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Triangle *>(disc[i])) )
			ring.push_back(disc[i]) ;
		else if(this->in(disc[i]->getCenter()))
			inDisc.push_back(disc[i]) ;
	}
	
	std::set<DelaunayTriangle *> newInterface ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(bimateralInterfaced.find(ring[i]) == bimateralInterfaced.end())
		{
			ring[i]->setBehaviour(new BimaterialInterface(static_cast<Circle *>(this),
														new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														ring[i]->getBehaviour()->getCopy()
														)) ;
			ring[i]->getBehaviour()->transform(ring[i]->getXTransform(), ring[i]->getYTransform()) ;
		}
		newInterface.insert(ring[i]) ;
	}
	
	std::set<DelaunayTriangle *> newExpansive ;
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
			disc[0]->setBehaviour(new BimaterialInterface(static_cast<Circle *>(this),
														new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														disc[0]->getBehaviour()->getCopy()
														)) ;
			disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform()) ;
		}
		newInterface.insert(disc[0]) ;
	}
	
	bimateralInterfaced = newInterface ;
	
}
	

