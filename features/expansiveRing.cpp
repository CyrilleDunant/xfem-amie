// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveRing.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/triple_behaviour.h"
#include "../physics/dual_behaviour.h"

using namespace Mu ;

ExpansiveRing::ExpansiveRing(Feature *father, double radius,double inradius, double x, double y, const Matrix & tensor, Vector def) : EnrichmentRing(father, radius, inradius, x, y),  imposedDef(def),cgTensor(tensor)
{
	
}

ExpansiveRing::~ExpansiveRing() {}
	
const Circle * ExpansiveRing::getGeometry() const 
{
	return static_cast<const Circle *>(this) ;
}

Circle * ExpansiveRing::getGeometry() 
{
	return static_cast<Circle *>(this) ;
}

void ExpansiveRing::reset() 
{
	cache.clear() ;
	updated = true ;
}

void ExpansiveRing::enrich(size_t & ,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{

	EnrichmentRing::enrich(dtree->getLastNodeId(), dtree) ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = cache ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersection(static_cast<Triangle *>(disc[i])).size() == 2 || self.intersection(static_cast<Triangle *>(disc[i])).size() == 2)
			ring.push_back(disc[i]) ;
		else if(this->in(disc[i]->getCenter()) && !self.in(disc[i]->getCenter()) )
			inDisc.push_back(disc[i]) ;
	}
	
	std::set<DelaunayTriangle *> newInterfaced ;
	std::set<DelaunayTriangle *> newExpansive ;
	std::set<DelaunayTriangle *> newBiInterfaced ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		if(interfaced.find(ring[i]) == interfaced.end())
		{
			ring[i]->setBehaviour(new TrimaterialInterface(&self, static_cast<Circle *>(this),
														ring[i]->getBehaviour()->getCopy(),new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														ring[i]->getBehaviour()->getCopy()
														)) ;
			ring[i]->getBehaviour()->transform(ring[i]->getXTransform(), ring[i]->getYTransform()) ;
		}
		
		newInterfaced.insert(ring[i]) ;
	}
	
	
	for(size_t i = 0 ; i < inDisc.size() ; i++)
	{
		if(expansive.find(ring[i]) == expansive.end())
		{
			inDisc[i]->setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
		}
		
		newExpansive.insert(inDisc[i]) ;
	}
	
	if(disc.size() == 1)
	{
		if(biInterfaced.find(disc[0]) == biInterfaced.end())
		{
			disc[0]->setBehaviour(new BimaterialInterface(static_cast<Circle *>(this),
														new StiffnessWithImposedDeformation(cgTensor, imposedDef),
														disc[0]->getBehaviour()->getCopy()
														)) ;
			disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform()) ;
		}
		
		newBiInterfaced.insert(disc[0]) ;
		
	}
	
	biInterfaced = newBiInterfaced ;
	interfaced = newInterfaced ;
	expansive = newExpansive ;
}
	

