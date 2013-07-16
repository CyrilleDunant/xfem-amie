// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "growingExpansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/damagemodels/fractiondamage.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness.h"
#include "../physics/maxwell.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/materials/aggregate_behaviour.h"

using namespace Mu ;

GrowingExpansiveZone::GrowingExpansiveZone(Feature* father, Function & g, double x, double y, ViscoelasticityAndImposedDeformation* i) : TimeDependentEnrichmentInclusion(father,g,x,y), imp(i)
{
	changed = true ;
}

GrowingExpansiveZone::~GrowingExpansiveZone() {
//  delete imp ;
} ;

void GrowingExpansiveZone::enrich(size_t & counter, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	TimeDependentEnrichmentInclusion::enrich(counter,dtree) ;
	
	std::vector<DelaunayTriangle *> & disc = TimeDependentEnrichmentInclusion::cache ;
	
// 	if(disc.size() < 6)
// 		return ;
	
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inDisc ;

	for( size_t i = 0 ; i < disc.size() ; i++ )
	{
		if( getPrimitive()->intersects( dynamic_cast<Triangle *>(disc[i]) ) )
			ring.push_back( disc[i] ) ;
		else /*if( in( disc[i]->getCenter() ) )*/
		{
			Point p = disc[i]->getCenter() ;
			p.t = VirtualMachine().eval( disc[i]->getTTransform(), 0,0,0,0) ;
			if(getPrimitive()->in(p))
				inDisc.push_back( disc[i] ) ;
		}
	}

	std::set<DelaunayTriangle *> newInterface ;

	for( size_t i = 0 ; i < ring.size() ; i++ )
	{
		if( bimateralInterfaced.find( ring[i] ) == bimateralInterfaced.end() )
		{
			BimaterialInterface *bi = new BimaterialInterface( getPrimitive(),
				imp->getCopy(),
				ring[i]->getBehaviour()->getCopy() );

			Geometry * src =  ring[i]->getBehaviour()->getSource() ;
//			delete ring[i]->getBehaviour() ;
			ring[i]->setBehaviour( bi ) ;
			bi->transform( ring[i]) ;
			bi->setSource( src );
		}
		else
			dynamic_cast<BimaterialInterface *>(ring[i]->getBehaviour())->transform( ring[i] ) ;
		  

		newInterface.insert( ring[i] ) ;
	}
	
// 	std::cout << "\n" << ring.size() << std::endl ;
// 	std::cout << inDisc.size() <<"\n" <<  std::endl ;
// 
	std::set<DelaunayTriangle *> newExpansive ;

	for( size_t i = 0 ; i < inDisc.size() ; i++ )
	{
		if( expansive.find( inDisc[i] ) == expansive.end() )
		{
//			delete inDisc[i]->getBehaviour() ;
			inDisc[i]->setBehaviour( imp->getCopy() ) ;
			inDisc[i]->getBehaviour()->setSource( getPrimitive() );
		}

		newExpansive.insert( inDisc[i] ) ;
	}

	expansive = newExpansive ;
	bimateralInterfaced = newInterface ;
	
	bool noPrevBC = dofIdPrev.empty() ;
	
	std::vector<DelaunayTriangle *> enriched(enrichedElem.begin(), enrichedElem.end()); 
//	std::cout << enriched.size() << "\t" << counter << std::endl ;
	std::valarray<bool> done(counter) ;
	done = false ;

	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		disc[i]->clearBoundaryConditions() ;
	}
	
	for(size_t i = 0 ; i < enriched.size() ; i++)
	{
		enriched[i]->clearBoundaryConditions() ;
		if( enriched[i]->getEnrichmentFunctions().size() == 0)
		{
			continue ;
		}
// 		if(enriched[i]->hasBoundaryConditions())
// 			continue ;
		
		
		GaussPointArray gp = enriched[i]->getGaussPoints() ;
		std::valarray<Matrix> Jinv( gp.gaussPoints.size() ) ;
		
		for( size_t j = 0 ; j < gp.gaussPoints.size() ;  j++ )
		{
			enriched[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}
		
		size_t enrichedDofPerTimePlanes = enriched[i]->getEnrichmentFunctions().size()/enriched[i]->timePlanes() ;
		
		
		for(size_t j = 0 ; j < enriched[i]->getEnrichmentFunctions().size() - enrichedDofPerTimePlanes ; j++)
		{
			if(done[ enriched[i]->getEnrichmentFunction(j).getDofID() ])
				continue ;
			done[ enriched[i]->getEnrichmentFunction(j).getDofID() ] = true ;
		  
			if(noPrevBC)
			{
				for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
				{
enriched[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, enriched[i], gp, Jinv, enriched[i]->getEnrichmentFunction(j).getDofID(), 0., n) ) ;
				}
			}
			else
			{
				size_t skip = 0 ;
				size_t ndof = imp->getNumberOfDegreesOfFreedom() ;
				Vector prev = enriched[i]->getState().getEnrichedDisplacements() ;
				size_t enrichedDofPerTimePlanesPrev = prev.size()/(ndof * enriched[i]->timePlanes() ) ;
				if( dofIdPrev.find( enriched[i]->getEnrichmentFunction(j).getPoint() ) == dofIdPrev.end() )
				{
					skip++ ;
					for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					{
enriched[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, enriched[i], gp, Jinv, enriched[i]->getEnrichmentFunction(j).getDofID(), 0., n) ) ;
					}
				}
				else
				{
					for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					{
enriched[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, enriched[i], gp, Jinv, enriched[i]->getEnrichmentFunction(j).getDofID(), prev[ (j-skip+enrichedDofPerTimePlanesPrev)*ndof + n ] , n) ) ;
					}
				}
				  
				  
//disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, dofIdCurrent[disc[i]->getEnrichmentFunction(j).getPoint()], 0., n) ) ;
			}
		}
//		std::cout << std::endl ;
	}
	
	dofIdPrev = dofIdCurrent ;

}


	


