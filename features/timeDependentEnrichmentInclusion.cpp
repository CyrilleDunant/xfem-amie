// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "timeDependentEnrichmentInclusion.h"
#include "inclusion.h"
#include "feature_base.h"
#include "sample.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/stiffness.h"
#include "../physics/maxwell.h"
#include "../physics/void_form.h"

using namespace Mu ;

TimeDependentEnrichmentInclusion::TimeDependentEnrichmentInclusion(Feature * father, Function & r, double x, double y) : EnrichmentInclusion( father, VirtualMachine().eval(r, Point(0,0,0,0)),x,y), TimeDependentCircle(r,Point(x,y))
{
  
}

TimeDependentEnrichmentInclusion::TimeDependentEnrichmentInclusion( Function & r, double x, double y) : EnrichmentInclusion( VirtualMachine().eval(r, Point(0,0,0,0)),x,y), TimeDependentCircle(r,Point(x,y))
{
  
}

TimeDependentEnrichmentInclusion::~TimeDependentEnrichmentInclusion() 
{
  
}
	
void TimeDependentEnrichmentInclusion::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	freeIds[dtree].clear();
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
		{
			std::vector<size_t> idpts = cache[i]->clearEnrichment(getPrimitive()) ;
			for(size_t j = 0 ; j < idpts.size() ; j++)
				freeIds[dtree].insert(idpts[j]) ;
		}
		cache[i]->enrichmentUpdated = true ;
	}
	cache = dtree->getConflictingElements(getPrimitive()) ;
	if(cache.empty())
	{
		std::vector<DelaunayTriangle *> candidates = dtree->getConflictingElements(&getCenter()) ;
		for(size_t i = 0 ; i < candidates.size() ; i++)
		{
			if(candidates[i]->isTriangle && candidates[i]->in(getCenter()))
			{
				cache.push_back(static_cast<DelaunayTriangle *>(candidates[i])) ;
				break ;
			}
		}
	}
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
		{
			std::vector<size_t> idpts = cache[i]->clearEnrichment(getPrimitive()) ;
			for(size_t j = 0 ; j < idpts.size() ; j++)
				freeIds[dtree].insert(idpts[j]) ;
		}
		cache[i]->enrichmentUpdated = true ;

	}
// 	std::cout << "update !!!! " << cache.size() << std::endl ;
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

Function getTimeDependentBlendingFunction(const std::map<const Point *, int> & dofIds, const DelaunayTriangle * t)
{
	return Function("1") ;
//	TriElement father(LINEAR_TIME_LINEAR) ;
	
/*	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return father.getShapeFunction(2) ;
	}
	
	if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return 1-father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
	{
		return 1-father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
	{
		return 1-father.getShapeFunction(2) ;
	}*/
	
}


void TimeDependentEnrichmentInclusion::enrich(size_t & lastId, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
//	std::cout << "are you called or what ???" << std::endl ;
  
//	std::cout << lastId << std::endl ;
  
	freeIds.clear() ;
	if(updated)
	{
		update(dtree) ;
	}
 	updated = false ;
	const std::vector<DelaunayTriangle *> & disc  = cache;
	
	if(disc.size() < 2)
	{
		return ;
	}
	
//	std::cout << "hello" << std::endl ;

	std::vector<DelaunayTriangle *> ring ;
	
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
			
		for(size_t j = 0 ; j < nodesIterator.size() && !added ; j++)
		{
			bool bin = false ;
			bool bout = false ;
			if(added)
				break ;
			
			
			if( in( disc[i]->getBoundingPoint( nodesIterator[j][0] ) ) )
				bin = true ;
			else
				bout = true ;

			
			if( in( disc[i]->getBoundingPoint( nodesIterator[j][1] ) ) )
				bin = true ;
			else
				bout = true ;

			if( in( disc[i]->getBoundingPoint( nodesIterator[j][2] ) ) )
				bin = true ;
			else
				bout = true ;
			
			if(bin && bout)
			{
				added = true ;
				ring.push_back(disc[i]) ;
			}
		}
	  
	}
	
//	std::cout << "\n" << ring.size() << "/" << disc.size() << std::endl ;
	
	//then we build a list of points to enrich
	std::set<Point *> points ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			points.insert(&ring[i]->getBoundingPoint(j)) ;
		}
		ring[i]->enrichmentUpdated = false ;
		for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
		{
			ring[i]->getNeighbourhood(j)->enrichmentUpdated = false ;
		}
	}

	//we build a map of the points and corresponding enrichment ids
	std::map<const Point *, int> dofId ;

	for(auto i = points.begin() ; i != points.end() ; ++i)
	{
			dofId[*i] = lastId++ ;
	}
	dofIdCurrent = dofId ;
	
	std::set<std::pair<DelaunayTriangle *, Point *> > enriched ;
	enrichedElem.clear() ;
	//then we iterate on every element
	
	std::map<Point*, size_t> extradofs ;
	
	std::cout << ring.size() << std::endl ;
	
	for(size_t i = 0 ; i < ring.size() ; i++)
	{	  
//		std::cout << "enriching triangle = " << i << "/" << ring.size() << std::endl ;
		enrichedElem.insert(ring[i]) ;
		std::vector<Point> hint ;
		std::vector<Point> inter = intersection(ring[i]->getPrimitive()) ;
		
		if(inter.size() == 4)
		{
			for(double j = 0.1 ; j < 0.9 ; j += 0.1)
			{
				for(double k = 0.1 ; k < 0.9 ; k += 0.1)
				{
					Point h0 = inter[0]*j*k+inter[1]*(1.-j)*k + inter[2]*j*(1.-k)+inter[3]*(1.-j)*(1.-k) ;
					project(&h0);
					hint.push_back(ring[i]->inLocalCoordinates(h0)) ;
				}
			}
			
			hint.push_back(ring[i]->inLocalCoordinates(inter[0])) ;
			hint.push_back(ring[i]->inLocalCoordinates(inter[1])) ;
			hint.push_back(ring[i]->inLocalCoordinates(inter[2])) ;
			hint.push_back(ring[i]->inLocalCoordinates(inter[3])) ;
		}
		

		//we build the enrichment function, first, we get the transforms from the triangle
		//this function returns the distance to the centre
		Function x("x") ; Function y("y") ;
// 		Function zero("0") ; Function one("1") ; Function z("0") ;
// 		zero.setNumberOfDerivatives(4);
// 		zero.setDerivative(XI, z);
// 		zero.setDerivative(ETA, z);
// 		zero.setDerivative(ZETA, z);
// 		zero.setDerivative(TIME_VARIABLE, z);
// 		one.setNumberOfDerivatives(4);
// 		one.setDerivative(XI, z);
// 		one.setDerivative(ETA, z);
// 		one.setDerivative(ZETA, z);
// 		one.setDerivative(TIME_VARIABLE, z);
// 		x.setNumberOfDerivatives(4);
// 		x.setDerivative(XI, one);
// 		x.setDerivative(ETA, zero);
// 		x.setDerivative(ZETA, zero);
// 		x.setDerivative(TIME_VARIABLE, zero);
// 		y.setNumberOfDerivatives(4);
// 		y.setDerivative(ETA, one);
// 		y.setDerivative(XI, zero);
// 		y.setDerivative(ZETA, zero);
// 		y.setDerivative(TIME_VARIABLE, zero);
		
		Function position = f_sqrt((x-getCenter().x)*(x-getCenter().x) +
		                           (y-getCenter().y)*(y-getCenter().y)) ;
					   
// 		Function dfx( "cst x -", getCenter().x) ;
// 		dfx /= position ;
// 		Function dfy( "cst y -", getCenter().y) ;
// 		dfy /= position ;
// 		
// 		dfx.setNumberOfDerivatives(4) ;
// 		dfx.setDerivative(XI, zero);
// 		dfx.setDerivative(ETA, zero);
// 		dfx.setDerivative(ZETA, zero);
// 		dfx.setDerivative(TIME_VARIABLE, zero);
// 		dfy.setNumberOfDerivatives(4) ;
// 		dfy.setDerivative(XI, zero);
// 		dfy.setDerivative(ETA, zero);
// 		dfy.setDerivative(ZETA, zero);
// 		dfy.setDerivative(TIME_VARIABLE, zero);
// 
// 		position.setNumberOfDerivatives(4);
// 		position.setDerivative(XI, dfx);
// 		position.setDerivative(ETA, dfy);
// 		position.setDerivative(ZETA, zero);
// 		position.setDerivative(TIME_VARIABLE, zero);
		
	
		Function hatPrev = 1.-f_abs(position-radiusAtTime( ring[i]->getBoundingPoint(0) ))/radiusAtTime(ring[i]->getBoundingPoint(0));
		Function hatNext = 1.-f_abs(position-radiusAtTime( ring[i]->getBoundingPoint(5) ))/radiusAtTime(ring[i]->getBoundingPoint(5));
		Function dx = ring[i]->getXTransform() ;
		Function dy = ring[i]->getYTransform() ;
		Function dt = ring[i]->getTTransform() ;
		hatPrev.setVariableTransform( XI, dx );
		hatPrev.setVariableTransform( ETA, dy );
 		hatPrev.setVariableTransform( TIME_VARIABLE, dt );
		hatNext.setVariableTransform( XI, dx );
		hatNext.setVariableTransform( ETA, dy );
		hatNext.setVariableTransform( TIME_VARIABLE, dt );
		hatPrev.setNumberOfDerivatives(0);
		hatNext.setNumberOfDerivatives(0);

		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			std::pair<DelaunayTriangle *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
			if(enriched.find(that) == enriched.end())
			{
				enriched.insert(that) ;
				Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;
				
				Function f ;
				if(j < ring[i]->getBoundingPoints().size() / ring[i]->timePlanes())
				{
				 f =  ring[i]->getShapeFunction(j)*(hatPrev - VirtualMachine().eval(hatPrev, p.x, p.y,p.z,p.t)) ;
				}
				else
				{
				 f =  ring[i]->getShapeFunction(j)*(hatNext - VirtualMachine().eval(hatNext, p.x, p.y,p.z,p.t)) ;
				}
				f.setIntegrationHint(hint) ;
				f.setPoint(&ring[i]->getBoundingPoint(j)) ;
				f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
		
				
				ring[i]->setEnrichment( f, getPrimitive()) ;
			}
		}
//		ring[i]->enrichmentUpdated = true ;
		
		hint.clear(); hint.push_back(Point(1./3., 1./3.,0,0));
		for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
		{
			DelaunayTriangle * t = ring[i]->getNeighbourhood(j) ;
			enrichedElem.insert(t) ;
			if(std::binary_search(ring.begin(), ring.end(), t) )
				continue ;
			
//			Function blend = getTimeDependentBlendingFunction(dofId, t) ;
			
			if(!t->enrichmentUpdated)
				t->clearEnrichment( getPrimitive()) ;
			
//			t->enrichmentUpdated = true ;
			bool hinted = false ;
			Function position = f_sqrt((x-getCenter().x)*(x-getCenter().x) +
		                           (y-getCenter().y)*(y-getCenter().y)) ;
			Function hatPrev = 1.-f_abs(position-radiusAtTime(t->getBoundingPoint(0)))/radiusAtTime(t->getBoundingPoint(0));
			Function hatNext = 1.-f_abs(position-radiusAtTime(t->getBoundingPoint(5)))/radiusAtTime(t->getBoundingPoint(5));
			Function dx = t->getXTransform() ;
			Function dy = t->getYTransform() ;
			Function dt = t->getTTransform() ;
			hatPrev.setVariableTransform(XI, dx);
			hatPrev.setVariableTransform(ETA, dy);
 			hatPrev.setVariableTransform(TIME_VARIABLE, dt);
			hatPrev.makeVariableTransformDerivative() ;
			hatNext.setVariableTransform(XI, dx);
			hatNext.setVariableTransform(ETA, dy);
 			hatNext.setVariableTransform(TIME_VARIABLE, dt);
			hatNext.makeVariableTransformDerivative() ;
			hatPrev.setNumberOfDerivatives(0);
			hatNext.setNumberOfDerivatives(0);
			
// 			Function hat = 1./(f_abs(position-getRadius())*0.2+2.*getRadius()) ;
			
			for(size_t k = 0 ; k< t->getBoundingPoints().size() ; k++)
			{
				std::pair<DelaunayTriangle *, Point *> that(t, &t->getBoundingPoint(k) ) ;
				if(enriched.find(that) == enriched.end())
				{
					if(dofId.find(&t->getBoundingPoint(k)) != dofId.end() )
					{
						enriched.insert(that) ;
						Point p = t->inLocalCoordinates(t->getBoundingPoint(k)) ;
						Function f ;
						if(k < t->getBoundingPoints().size() / t->timePlanes())
						{
						f =  t->getShapeFunction(k)*(hatPrev - VirtualMachine().eval(hatPrev, p.x, p.y,p.z,p.t)) ;
						}
						else
						{
						f =  t->getShapeFunction(k)*(hatNext - VirtualMachine().eval(hatNext, p.x, p.y,p.z,p.t)) ;
						}
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(&t->getBoundingPoint(k)) ;
						f.setDofID(dofId[&t->getBoundingPoint(k)]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
				}
			}
		}
	}
	
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(enrichedElem.find(disc[i]) == enrichedElem.end())
			disc[i]->clearAllEnrichment() ;
	}
	
	std::cout << enrichedElem.size() << std::endl ;
	
// 	std::cout << lastId << std::endl ;
	
//	updated = true ;*/
}
	
void TimeDependentEnrichmentInclusion::step(double dt, std::valarray< double >*, const Mu::Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
{
	if(dt > POINT_TOLERANCE_2D)
	{
		if( VirtualMachine().deval(TimeDependentCircle::getRadiusFunction(), TIME_VARIABLE, 0,0,0,dt) > POINT_TOLERANCE_3D)
		{
	  
		changed = true ;
		updated = true ;
		
//		std::cout << "update you bastard!" << std::endl ;
		}
		else
		{

		  
		}
	}
	else
	{
		changed = false ;
		updated = false ;
	}
  
}
	

