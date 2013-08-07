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


	std::vector<std::vector<size_t > > nodesIterator ;
	size_t nodesPerPlane = disc[0]->getBoundingPoints().size() / disc[0]->timePlanes() ;
	size_t nodesPerEdge = nodesPerPlane/3 ;
	for(size_t i = 0 ; i < disc[0]->timePlanes() ; i++)
	{
		std::vector<size_t> tmp ;
		tmp.push_back(nodesPerPlane*i + nodesPerEdge*0);
		tmp.push_back(nodesPerPlane*i + nodesPerEdge*1);
		tmp.push_back(nodesPerPlane*i + nodesPerEdge*2);
		nodesIterator.push_back(tmp) ;
	}
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		bool added = false ;
		double dt = disc[i]->getState().getNodalDeltaTime() ;
			
		bool bin = false ;
		bool bout = false ;
		for(size_t j = 0 ; j < nodesIterator.size() && !added ; j++)
		{
			Point A = disc[i]->getBoundingPoint( nodesIterator[j][0] ) ;
			if(in(A))
				bin= true ;
			else
				bout = true ;
			Point B = disc[i]->getBoundingPoint( nodesIterator[j][1] ) ;
			if(in(B))
				bin= true ;
			else
				bout = true ;
			Point C = disc[i]->getBoundingPoint( nodesIterator[j][2] ) ;
			if(in(C))
				bin= true ;
			else
				bout = true ;

  			if(bin && bout)
  			{
				added = true ;
				ring.push_back(disc[i]) ;
 			}

		}
		if(!added)
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
//	return ;
	expansive = newExpansive ;
	bimateralInterfaced = newInterface ;

	
	bool noPrevBC = dofIdPrev.empty() ;
	
	std::vector<DelaunayTriangle *> enriched(enrichedElem.begin(), enrichedElem.end()); 
//	std::cout << enriched.size() << "\t" << counter << std::endl ;

	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		disc[i]->clearBoundaryConditions() ;
		if( disc[i]->getEnrichmentFunctions().size() == 0)
		{
			continue ;
		}
// 		if(enriched[i]->hasBoundaryConditions())
// 			continue ;
		
		
		GaussPointArray gp = disc[i]->getGaussPoints() ;
		std::valarray<Matrix> Jinv( gp.gaussPoints.size() ) ;
		
		for( size_t j = 0 ; j < gp.gaussPoints.size() ;  j++ )
		{
			disc[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}
		
		size_t enrichedDofPerTimePlanes = disc[i]->getEnrichmentFunctions().size()/disc[i]->timePlanes() ;
		
		
		for(size_t j = 0 ; j < disc[i]->getEnrichmentFunctions().size() - enrichedDofPerTimePlanes ; j++)
		{

			if(noPrevBC)
			{
				if( pointsAndValues.find( disc[i]->getEnrichmentFunction(j).getPoint()) != pointsAndValues.end())
				{
					for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					{
						disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, disc[i]->getEnrichmentFunction(j).getDofID(), 0., n) ) ;
					}
				}
			}
			else
			{
				size_t skip = 0 ;
				size_t ndof = imp->getNumberOfDegreesOfFreedom() ;
				size_t enrichedDofPerTimePlanesPrev = disc[i]->getState().getEnrichedDisplacements().size()/(ndof * disc[i]->timePlanes() ) ;
				if( pointsAndValues.find( disc[i]->getEnrichmentFunction(j).getPoint() ) == pointsAndValues.end() )
				{
					skip++ ;
					for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					{
						disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, disc[i]->getEnrichmentFunction(j).getDofID(), 0., n) ) ;
					}
				}
				else
				{
					for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					{
						disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, disc[i]->getEnrichmentFunction(j).getDofID(), pointsAndValues[disc[i]->getEnrichmentFunction(j).getPoint()][n], n) ) ;
					}
				}
				  
				  
//disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, dofIdCurrent[disc[i]->getEnrichmentFunction(j).getPoint()], 0., n) ) ;
			}
		}
//		std::cout << std::endl ;
	}
	
	dofIdPrev = dofIdCurrent ;

}

void GrowingExpansiveZone::step(double dt, std::valarray< double >* disp, const Mu::Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
{
	TimeDependentEnrichmentInclusion::step(dt, disp, dtree) ;

	pointsAndValues.clear() ;
	for(size_t i = 0 ; i < TimeDependentEnrichmentInclusion::cache.size() ; i++)
	{
		std::cout <<  TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunctions().size() << std::endl ;
		for(  size_t j = 0 ; j < TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunctions().size() ;j++)
		{
			std::cout << i << ";" << j << std::endl ;
			if(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentSource(j) == getPrimitive())
			{
				Vector disp(imp->getNumberOfDegreesOfFreedom()) ;
				for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
					disp[n] = TimeDependentEnrichmentInclusion::cache[i]->getState().getEnrichedDisplacements()[j*imp->getNumberOfDegreesOfFreedom()+n] ;

				for(size_t k = 0 ; k < TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoints().size() ; k++)
				{
					std::cout << i << ";" << j << ";" << k << std::endl ;
					if(squareDist2D(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint(), &TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoint(k)) < POINT_TOLERANCE_2D 
						  && std::abs(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint()->t - TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoint(k).t) > POINT_TOLERANCE_2D )
					{
						pointsAndValues[TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint()] = disp ;
						break ;
					}
				}
			}
		}
	}
}
	


