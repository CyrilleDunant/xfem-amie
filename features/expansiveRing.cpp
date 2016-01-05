// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveRing.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/triple_behaviour.h"
#include "../physics/dual_behaviour.h"

namespace Amie {

ExpansiveRing::ExpansiveRing(Feature *father, double radius,double inradius, double x, double y, const Matrix & tensor, Vector def) : EnrichmentRing(father, radius, inradius, x, y),  imposedDef(def),cgTensor(tensor)
{
	
}

ExpansiveRing::~ExpansiveRing() {}
	

void ExpansiveRing::reset() 
{
	cache.clear() ;
	updated = true ;
}

void ExpansiveRing::enrich(size_t & lastId,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{

	EnrichmentRing::enrich(lastId, dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = cache ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if((intersection(static_cast<Triangle *>(disc[i])).size() == 2 
			|| self.intersection(static_cast<Triangle *>(disc[i])).size() == 2) 
			&& disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			ring.push_back(disc[i]) ;
		else if(this->in(disc[i]->getCenter()) && !self.in(disc[i]->getCenter()) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			inDisc.push_back(disc[i]) ;
	}
	
	std::set<DelaunayTriangle *> newInterfaced ;
	std::set<DelaunayTriangle *> newExpansive ;
	std::set<DelaunayTriangle *> newBiInterfaced ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(interfaced.find(ring[i]) == interfaced.end())
		{
			ring[i]->setBehaviour(dtree, new TrimaterialInterface(&self, getPrimitive(),
														ring[i]->getBehaviour()->getCopy(),new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														ring[i]->getBehaviour()->getCopy()
														)) ;
			ring[i]->getBehaviour()->transform(ring[i]) ;
		}
		
		newInterfaced.insert(ring[i]) ;
	}
	
	
	for(size_t i = 0 ; i < inDisc.size() ; i++)
	{
		if(expansive.find(ring[i]) == expansive.end())
		{
			inDisc[i]->setBehaviour(dtree, new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
		}
		
		newExpansive.insert(inDisc[i]) ;
	}
	
	if(disc.size() == 1)
	{
		if(biInterfaced.find(disc[0]) == biInterfaced.end())
		{
			disc[0]->setBehaviour(dtree, new BimaterialInterface(getPrimitive(),
														new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														disc[0]->getBehaviour()->getCopy()
														)) ;
			disc[0]->getBehaviour()->transform(disc[0]) ;
		}
		
		newBiInterfaced.insert(disc[0]) ;
		
	}
	
	biInterfaced = newBiInterfaced ;
	interfaced = newInterfaced ;
	expansive = newExpansive ;
}

}

