// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"

using namespace Mu ;

ExpansiveZone::ExpansiveZone(Feature *father, double radius, double x, double y, Matrix tensor, Vector def) : EnrichmentInclusion(father, radius, x, y),  imposedDef(def),cgTensor(tensor)
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
}

void ExpansiveZone::enrich(size_t & counter,  DelaunayTree * dtree)
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
	
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		ring[i]->setBehaviour(new BimaterialInterface(static_cast<Circle *>(this),
		                                              new StiffnessWithImposedDeformation(cgTensor, imposedDef),
		                                              ring[i]->getBehaviour()->getCopy()
		                                             )) ;
		ring[i]->getBehaviour()->transform(ring[i]->getXTransform(), ring[i]->getYTransform()) ;
	}
	
	for(size_t i = 0 ; i < inDisc.size() ; i++)
	{
		inDisc[i]->setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
	}
	
	if(disc.size() == 1)
	{
		disc[0]->setBehaviour(new BimaterialInterface(static_cast<Circle *>(this),
		                                              new StiffnessWithImposedDeformation(cgTensor, imposedDef),
		                                              disc[0]->getBehaviour()->getCopy()
		                                             )) ;
		disc[0]->getBehaviour()->transform(disc[0]->getXTransform(), disc[0]->getYTransform()) ;
	}
	
}
	

