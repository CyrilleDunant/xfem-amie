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

using namespace Amie ;

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
	double dt = 1. ;
	if(cache.size())
		dt = cache[0]->getState().getNodalDeltaTime() * 2. ;

	Function time("t") ;
	time += dt ;
	Function radius = getRadiusFunction(time) ;
	Point c = getCenter() ;
	TimeDependentCircle dummy( radius, c) ;
	
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
	cache = dtree->getConflictingElements(& dummy) ;
	
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
	TriElement father(LINEAR) ;
	
	if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
	{
		return father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
	{
		return father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
	{
		return father.getShapeFunction(2) ;
	}
	
	if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
	{
		return 1-father.getShapeFunction(0) ;
	}
	
	if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
	{
		return 1-father.getShapeFunction(1) ;
	}
	
	if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
	{
		return 1-father.getShapeFunction(2) ;
	}
	return Function("1") ;
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
	
// 	if(disc.size() < 2)
// 	{
// 		return ;
// 	}
	
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
		double dt = disc[i]->getState().getNodalDeltaTime() ;
			
		bool bin = false ;
		bool bout = false ;
		for(size_t j = 0 ; j < nodesIterator.size() && !added ; j++)
		{
			if(added)
				break ;

			double ddt = 0. ;
/*			if(j == 0)
				ddt = -dt ;
			if(j == nodesIterator.size()-1)
				ddt = +dt ;*/
			
			Point A = disc[i]->getBoundingPoint( nodesIterator[j][0] ) ;
			A.getT() += ddt ;
			if(in(A))
				bin= true ;
			else
				bout = true ;
			Point B = disc[i]->getBoundingPoint( nodesIterator[j][1] ) ;
			B.getT() += ddt ;
			if(in(B))
				bin= true ;
			else
				bout = true ;
			Point C = disc[i]->getBoundingPoint( nodesIterator[j][2] ) ;
			C.getT() += ddt ;
			if(in(C))
				bin= true ;
			else
				bout = true ;

	  		if(bin && bout)
	  		{
				added = true ;
				ring.push_back(disc[i]) ;
				break ;
 			}
		}
	}
	std::sort(ring.begin(), ring.end()) ;
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
        std::vector<DelaunayTriangle *> neighbourhood = dtree->getNeighbourhood(ring[i]) ;
		for(auto & n : neighbourhood)
		{
			n->enrichmentUpdated = false ;
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
	
	for(size_t i = 0 ; i < ring.size() ; i++)
	{	  
//		std::cout << "enriching triangle = " << i << "/" << ring.size() << std::endl ;
		enrichedElem.insert(ring[i]) ;
		std::vector<Point> hint ;		

		//we build the enrichment function, first, we get the transforms from the triangle
		//this function returns the distance to the centre
		Function x("x") ; Function y("y") ;

		Function position = f_sqrt((x-getCenter().getX())*(x-getCenter().getX()) +
		                           (y-getCenter().getY())*(y-getCenter().getY())) ;
					   
		
		Function hat = 1.-f_abs(position-getRadiusFunction())/getRadiusFunction(); //radiusAtTime( ring[i]->getBoundingPoint(0) )-f_abs(position-radiusAtTime( ring[i]->getBoundingPoint(0) ))^2;
		Function dx = ring[i]->getXTransform() ;
		Function dy = ring[i]->getYTransform() ;
		Function dt = ring[i]->getTTransform() ;
		hat.setVariableTransform( XI, dx );
		hat.setVariableTransform( ETA, dy );
 		hat.setVariableTransform( TIME_VARIABLE, dt );
		

		hat.setNumberOfDerivatives(0);

		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			std::pair<DelaunayTriangle *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
			if(enriched.find(that) == enriched.end())
			{
				enriched.insert(that) ;
				Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;
				
				Function f =  ring[i]->getShapeFunction(j)*(hat - VirtualMachine().eval(hat, p.getX(), p.getY(),p.getZ(),p.getT())) ;

				f.setIntegrationHint(hint) ;
				f.setPoint(&ring[i]->getBoundingPoint(j)) ;
				f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
		
				
				ring[i]->setEnrichment( f, getPrimitive()) ;
			}
		}
//		ring[i]->enrichmentUpdated = true ;
		
		hint.clear(); hint.push_back(Point(1./3., 1./3.,0,0));
        std::vector<DelaunayTriangle *> neighbourhood = dtree->getNeighbourhood(ring[i]) ;
		for(auto & t : neighbourhood)
		{
			enrichedElem.insert(t) ;
			if(std::binary_search(ring.begin(), ring.end(), t) )
				continue ;
			
			Function blend = getTimeDependentBlendingFunction(dofId, t) ;
			
			if(!t->enrichmentUpdated)
				t->clearEnrichment( getPrimitive()) ;
			
//			t->enrichmentUpdated = true ;
			bool hinted = false ;
			Function position = f_sqrt((x-getCenter().getX())*(x-getCenter().getX()) +
		                           (y-getCenter().getY())*(y-getCenter().getY())) ;
			Function hat= 1.-f_abs(position-getRadiusFunction())/getRadiusFunction(); //radiusAtTime(t->getBoundingPoint(0))-f_abs(position-radiusAtTime(t->getBoundingPoint(0)))^2;
			Function dx = t->getXTransform() ;
			Function dy = t->getYTransform() ;
			Function dt = t->getTTransform() ;
			hat.setVariableTransform(XI, dx);
			hat.setVariableTransform(ETA, dy);
 			hat.setVariableTransform(TIME_VARIABLE, dt);
// 			hat.makeVariableTransformDerivative() ;
			hat.setNumberOfDerivatives(0);
			
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
						f =  t->getShapeFunction(k)*(hat - VirtualMachine().eval(hat, p.getX(), p.getY(),p.getZ(),p.getT())) ;
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
	
//	std::cout << enrichedElem.size() << std::endl ;
	
// 	std::cout << lastId << std::endl ;
	
//	updated = true ;*/
}
	
void TimeDependentEnrichmentInclusion::step(double dt, std::valarray< double >*, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
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
//		changed = false ;
//		updated = false ;

		  
		}
	}
	else
	{
		changed = false ;
		updated = false ;
	}
  
}
	



TimeDependentHomogenisingInclusion::TimeDependentHomogenisingInclusion(Feature * father, Function & r, double x, double y, ViscoelasticityAndImposedDeformation * imp) : EnrichmentInclusion( father, VirtualMachine().eval(r, Point(0,0,0,0)),x,y), TimeDependentCircle(r,Point(x,y)), imposed(imp)
{
  
}

TimeDependentHomogenisingInclusion::TimeDependentHomogenisingInclusion( Function & r, double x, double y, ViscoelasticityAndImposedDeformation * imp) : EnrichmentInclusion( VirtualMachine().eval(r, Point(0,0,0,0)),x,y), TimeDependentCircle(r,Point(x,y)), imposed(imp)
{
  
}

TimeDependentHomogenisingInclusion::~TimeDependentHomogenisingInclusion() 
{
  
}
	
void TimeDependentHomogenisingInclusion::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	double dt = 1. ;
	if(cache.size())
		dt = cache[0]->getState().getNodalDeltaTime() * 2. ;

	Function time("t") ;
	time += dt ;
	Function radius = getRadiusFunction() ;
	Point c = getCenter() ;
	TimeDependentCircle dummy( radius, c) ;

	cache = dtree->getConflictingElements(& dummy) ;

	
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
		if(homogeneised.find( cache[i] ) == homogeneised.end() )
		{
			homogeneised[cache[i]] = cache[i]->getBehaviour()->getCopy() ;
		}
	}
// 	std::cout << "update !!!! " << cache.size() << std::endl ;

	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void TimeDependentHomogenisingInclusion::enrich(size_t & lastId, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	if(updated)
	{
		update(dtree) ;
	}
 	updated = false ;
	

	const std::vector<DelaunayTriangle *> & disc  = cache;
	std::vector<DelaunayTriangle *> ring ;
	std::vector<DelaunayTriangle *> inside ;

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
		double t0 = disc[i]->getBoundingPoint(0).getT() ;
			
		bool bin = false ;
		bool bout = false ;
		for(size_t j = 0 ; j < nodesIterator.size() && !added ; j++)
		{
			if(added)
				break ;

			double ddt = 0. ;
/*			if(j == 0)
				ddt = -dt ;
			if(j == nodesIterator.size()-1)
				ddt = +dt ;*/
			
			Point A = disc[i]->getBoundingPoint( nodesIterator[j][0] ) ;
			A.getT() += ddt ;

			Point B = disc[i]->getBoundingPoint( nodesIterator[j][1] ) ;
			B.getT() += ddt ;

			Point C = disc[i]->getBoundingPoint( nodesIterator[j][2] ) ;
			C.getT() += ddt ;

			Triangle dummy(A,B,C) ;
	  		if(circleAtTime(A).intersects(&dummy))
	  		{
				added = true ;
				ring.push_back(disc[i]) ;
				break ;
 			}
		}
		if(!added)
		{
			Point c = disc[i]->getCenter() ; c.getT() = t0+dt/2 ;
			Point a = disc[i]->getBoundingPoint(0) ; a.getT() = t0+dt/2 ;
			if(in(c) && in(a))
				inside.push_back(disc[i]) ;
			else if(in(c))
				ring.push_back(disc[i]) ;
		}
	}

//	std::cout << disc.size() << "\t" << inside.size() << "\t" << ring.size() << std::endl ;

	if(inside.size() == 0 && ring.size() == 0)
	{
		for(size_t i = 0 ; i < disc.size() ; i++)
			ring.push_back(disc[i]) ;
	}


	for(size_t i = 0 ; i < inside.size() ; i++)
	{
		inside[i]->setBehaviour(dtree,imposed->getCopy()) ;
		if(zoneElem.find(inside[i]) == zoneElem.end())
		{
			zoneElem.insert(inside[i]) ;
			inside[i]->clearElementaryMatrix() ;
		}
	}

	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		ring[i]->clearElementaryMatrix() ;
		Point A = ring[i]->getBoundingPoint(3) ;
		Point B = ring[i]->getBoundingPoint(4) ;
		Point C = ring[i]->getBoundingPoint(5) ;
		int k = 0 ;
		int areIn = 0 ;
		while(k < 1000)
		{
			double x = ((double) rand())/RAND_MAX ;
			double y = ((double) rand())/RAND_MAX ;
			if(x+y < 1)
			{
				Point D = A ;
				D += (B-A)*x ;
				D += (C-A)*y ;
				D.getT() = A.getT() ;
				if(in(D))
					areIn++ ;
				k++ ;
			}
		}
		
		Matrix stiff = (imposed->param)*((double) areIn/(double)k) ;
		stiff += (homogeneised[ring[i]]->getTensor(A))*((double) (k-areIn)/(double)k) ;
		size_t b = imposed->blocks ;
		Matrix realStiff(stiff.numRows()/b, stiff.numRows()/b) ;
		for(size_t r = 0 ; r < realStiff.numRows() ; r++)
		{
			for(size_t c = 0 ; c < realStiff.numRows() ; c++)
			{
				realStiff[r][c] = stiff[r][c] ;
			}		
		}

//		(homogeneised[ring[i]]->param).print() ;
		Point c = ring[i]->getCenter() ; c.getT() = A.getT() ;
//		std::cout << (int) dynamic_cast<Viscoelasticity *>(homogeneised[ring[i]])->model ;
		double factor = ((double) areIn/(double)k) * imposed->param[0][0]/realStiff[0][0] ;
		Vector alpha = (imposed->getImposedStrain(ring[i]->getCenter())) * factor;
		ring[i]->setBehaviour(dtree, new ViscoelasticityAndImposedDeformation(PURE_ELASTICITY,realStiff, alpha, b-1)) ;

		
	}

	
}


void TimeDependentHomogenisingInclusion::step(double dt, std::valarray< double >*, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
{
	if(dt > POINT_TOLERANCE_2D)
	{
		if( VirtualMachine().deval(TimeDependentCircle::getRadiusFunction(), TIME_VARIABLE, 0,0,0,dt) > POINT_TOLERANCE_3D)
		{
	  
//		changed = true ;
//		updated = true ;
		
//		std::cout << "update you bastard!" << std::endl ;
		}
		else
		{
//		changed = false ;
//		updated = false ;

		  
		}
	}
	else
	{
		changed = false ;
		updated = false ;
	}
  
}
