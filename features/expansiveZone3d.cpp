// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "expansiveZone3d.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"


using namespace Mu ;

ExpansiveZone3D::ExpansiveZone3D(Feature *father, double radius, double x, double y, double z, const Matrix & tensor, Vector def) : EnrichmentInclusion3D(father, radius, x, y, z),  imposedDef(def),cgTensor(tensor)
{
	this->setBehaviour(new StiffnessWithImposedDeformation(tensor, def)) ;
}

ExpansiveZone3D::ExpansiveZone3D(Feature *father, double radius, double x, double y, double z, StiffnessWithImposedDeformation * exp) : EnrichmentInclusion3D(father, radius, x, y, z)
{
	imposedDef.resize(6) ;
	imposedDef = exp->imposed ;
	cgTensor = exp->getTensor(Point(0,0,0,0)) ;
	this->setBehaviour(new StiffnessWithImposedDeformation(cgTensor, imposedDef)) ;
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

	for( size_t i = 0 ; i < ring.size() ; i++ )
	{
		if( bimateralInterfaced.find( ring[i] ) == bimateralInterfaced.end() )
		{
			BimaterialInterface *bi = nullptr ;

			if( dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() ) )
			{
//				std::cout << "get original" << std::endl ;
				bi = new BimaterialInterface( getPrimitive(),
				                              new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
				                              dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() )->original->getCopy()
				                            ) ;
//				std::cout << dynamic_cast<NonLocalMohrCoulomb *>( dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() )->original->getFractureCriterion() )->upVal << std::endl ;
			}
			else
			{
				bi = new BimaterialInterface( getPrimitive(),
			        new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
			        ring[i]->getBehaviour()->getCopy() ) ;
			}

			delete ring[i]->getBehaviour() ;
			ring[i]->setBehaviour( bi ) ;
			bi->transform( ring[i] ) ;
			bi->setSource( getPrimitive() );
		}

		newInterface.insert( ring[i] ) ;
	}

	std::set<DelaunayTetrahedron *> newExpansive ;

	for( size_t i = 0 ; i < inDisc.size() ; i++ )
	{
		if( expansive.find( inDisc[i] ) == expansive.end() )
		{
			StiffnessWithImposedDeformation * bi = new StiffnessWithImposedDeformation( cgTensor, imposedDef ) ;
			delete inDisc[i]->getBehaviour() ;
			inDisc[i]->setBehaviour( bi ) ;
			inDisc[i]->getBehaviour()->setSource( getPrimitive() );
		}

		newExpansive.insert( inDisc[i] ) ;
	}

	expansive = newExpansive ;

	if( disc.size() == 1 )
	{
// 		if( bimateralInterfaced.find( disc[0] ) == bimateralInterfaced.end() )
// 		{
		std::vector<Feature *> brother ;
		if(getFather())
			brother = this->getFather()->getChildren() ;
		std::vector<Feature *> feat ;
		for(size_t i = 0 ; i < brother.size() ; i++)
		{
			if(disc[0]->in(brother[i]->getCenter()))
				feat.push_back(brother[i]) ;
		}
		
	
		if( dynamic_cast<HomogeneisedBehaviour *>( disc[0]->getBehaviour() ) )
		{
//			std::cout << "get original" << std::endl ;
			dynamic_cast<HomogeneisedBehaviour *>( disc[0]->getBehaviour() )->updateEquivalentBehaviour(feat, disc[0]) ;
		}
		else
		{
//			std::cout << "set as homogenized" << std::endl ;
			disc[0]->setBehaviour(new HomogeneisedBehaviour( feat, disc[0] )) ;
		}
// 			else
// 			{
// 				BimaterialInterface * bi = new BimaterialInterface( getPrimitive(),
// 				                       new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
// 				                       disc[0]->getBehaviour()->getCopy()
// 				                                              ) ;
// 				delete disc[0]->getBehaviour() ;
// 				disc[0]->setBehaviour( bi ) ;
// 			}
// 
// 			disc[0]->getBehaviour()->transform( disc[0]->getXTransform(), disc[0]->getYTransform(), disc[0]->getZTransform() ) ;
// 			disc[0]->getBehaviour()->setSource( getPrimitive() );
// 		}

		newInterface.insert( disc[0] ) ;
	}

	bimateralInterfaced = newInterface ;
	
}
	

