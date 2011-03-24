// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"

using namespace Mu ;

ExpansiveZone::ExpansiveZone(Feature *father, double radius, double x, double y, const Matrix & tensor, Vector def) : EnrichmentInclusion(father, radius, x, y),  imposedDef(def),cgTensor(tensor)
{
	setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;	
}

ExpansiveZone::~ExpansiveZone() {}
	

void ExpansiveZone::reset() 
{
	cache.clear() ;
	updated = true ;
}

void ExpansiveZone::enrich(size_t & lastId , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	this->setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
	EnrichmentInclusion::enrich(lastId, dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> & disc = EnrichmentInclusion::cache ;/*dtree->getConflictingElements(getPrimitive()) ;
        if(disc.size() < 2)
		return ;*/
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		Segment s0(*disc[i]->first, *disc[i]->second) ;
		Segment s1(*disc[i]->second, *disc[i]->third) ;
		Segment s2(*disc[i]->third, *disc[i]->first) ;
		if(!(s0.intersection(getPrimitive()).empty() && s1.intersection(getPrimitive()).empty() && s2.intersection(getPrimitive()).empty())&& disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			ring.push_back(disc[i]) ;
		else if(in(disc[i]->getCenter()))
			inDisc.push_back(disc[i]) ;
	}
	
	std::set<DelaunayTriangle *> newInterface ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(bimateralInterfaced.find(ring[i]) == bimateralInterfaced.end())
		{
			ring[i]->setBehaviour(new BimaterialInterface(getPrimitive(),
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
			disc[0]->setBehaviour(new BimaterialInterface(getPrimitive(),
														new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														disc[0]->getBehaviour()->getCopy()
														)) ;
			disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform()) ;
		}
		newInterface.insert(disc[0]) ;
	}
	bimateralInterfaced = newInterface ;
}
	

void ExpansiveZone::setExpansion(Vector a)
{
	imposedDef = a ;
	updated = true ;

	for( auto i = bimateralInterfaced.begin() ; i != bimateralInterfaced.end() ; i++)
	{
		BimaterialInterface* interface = static_cast<BimaterialInterface *>((*i)->getBehaviour()) ;
		static_cast<StiffnessWithImposedDeformation*>(interface->inBehaviour)->imposed = a ;
	}
	
	for( auto i = expansive.begin() ; i != expansive.end() ; i++)
	{
		static_cast<StiffnessWithImposedDeformation*>((*i)->getBehaviour())->imposed = a ;
	}
}


MaterialInclusion::MaterialInclusion(Feature *father, double radius, double x, double y, LinearForm * inclusionBehaviour) : EnrichmentInclusion(father, radius, x, y),  inclusionBehaviour(inclusionBehaviour)
{
	
}

MaterialInclusion::~MaterialInclusion() {}

void MaterialInclusion::reset() 
{
	cache.clear() ;
	updated = true ;
}

void MaterialInclusion::enrich(size_t & counter , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	EnrichmentInclusion::enrich(counter, dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = cache ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Triangle *>(disc[i])) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			ring.push_back(disc[i]) ;
		else if(this->in(disc[i]->getCenter())&& disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			inDisc.push_back(disc[i]) ;
	}
	
	std::set<DelaunayTriangle *> newInterface ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(bimateralInterfaced.find(ring[i]) == bimateralInterfaced.end())
		{
			ring[i]->setBehaviour(new BimaterialInterface(getPrimitive(),
														inclusionBehaviour->getCopy(),
														ring[i]->getBehaviour()->getCopy()
														)) ;
			ring[i]->getBehaviour()->transform(ring[i]->getXTransform(), ring[i]->getYTransform()) ;
		}
		newInterface.insert(ring[i]) ;
	}
	
	std::set<DelaunayTriangle *> newExpansive ;
	for(size_t i = 0 ; i < inDisc.size() ; i++)
	{
		if(internal.find(inDisc[i]) == internal.end())
			inDisc[i]->setBehaviour(inclusionBehaviour->getCopy()) ;
		
		newExpansive.insert(inDisc[i]) ;
	}
	internal = newExpansive ;
	
	if(disc.size() == 1)
	{
		if(bimateralInterfaced.find(disc[0]) == bimateralInterfaced.end())
		{
			disc[0]->setBehaviour(new BimaterialInterface(getPrimitive(),
														inclusionBehaviour->getCopy(),
														disc[0]->getBehaviour()->getCopy()
														)) ;
			disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform()) ;
		}
		newInterface.insert(disc[0]) ;
	}
	
	bimateralInterfaced = newInterface ;
	
}
	


