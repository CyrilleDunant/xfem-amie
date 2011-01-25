// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion.h"
#include "inclusion.h"
#include "feature_base.h"
#include "sample.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/stiffness.h"
#include "../physics/void_form.h"

using namespace Mu ;

EnrichmentInclusion::EnrichmentInclusion(Feature *father, double radius, double x, double y) : EnrichmentFeature(father), Circle(radius, x, y)
{
	updated = true ;
}

EnrichmentInclusion::EnrichmentInclusion(double radius, double x, double y) :EnrichmentFeature(NULL), Circle(radius, x, y)
{
	updated = true ;
}

EnrichmentInclusion::~EnrichmentInclusion() {}
	
bool EnrichmentInclusion::enrichmentTarget(DelaunayTriangle * t)
{
	Segment s0(*t->first, *t->second) ;
	Segment s1(*t->second, *t->third) ;
	Segment s2(*t->third, *t->first) ;
	return (s0.intersects(getPrimitive()) || s1.intersects(getPrimitive()) || s2.intersects(getPrimitive())) ;
}

void EnrichmentInclusion::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
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
	
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void EnrichmentInclusion::enrich(size_t & lastId, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	freeIds.clear() ;
	if(updated)
	{
		update(dtree) ;
	}
	updated = false ;
		
	const std::vector<DelaunayTriangle *> & disc  = cache;
// 	if(disc.size() < 6)
// 	{
// 		std::cout << "cowardly discarding" << std::endl ;
// 		return ;
// 	}
	TriElement father(disc[0]->getOrder()) ;
	
	if(disc.size() == 1) // special case for really small inclusions
	{
	
		std::vector<Feature *> brother = this->getFather()->getChildren() ;
		std::vector<Feature *> feat ;
		for(size_t i = 0 ; i < brother.size() ; i++)
		{
			if(disc[0]->in(brother[i]->getCenter()))
				feat.push_back(brother[i]) ;
		}
			

		HomogeneisedBehaviour * hom = new HomogeneisedBehaviour(feat, disc[0]) ;
	
		disc[0]->setBehaviour(hom) ;
		
		return ;


	}

//	for(size_t i = 0 ; i < disc.size() ; i++)
//	{
//		disc[i]->setBehaviour(this->getFather()->getBehaviour()->getCopy()) ;	
//	}

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		Segment s0(*disc[i]->first, *disc[i]->second) ;
		Segment s1(*disc[i]->second, *disc[i]->third) ;
		Segment s2(*disc[i]->third, *disc[i]->first) ;
		if(!(getPrimitive()->in(*disc[i]->first) && getPrimitive()->in(*disc[i]->second) && getPrimitive()->in(*disc[i]->third)) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			ring.push_back(disc[i]) ;
		}
	}
	//then we build a list of points to enrich
	std::set<Point *> points ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		for(size_t j = 0 ; j < ring[i]->getBoundingPoints().size() ; j++)
		{
			points.insert(&ring[i]->getBoundingPoint(j)) ;

		}
	}

	
	//we build a map of the points and corresponding enrichment ids
	std::map<Point *, int> dofId ;
	
	for(std::set<Point * >::const_iterator i = points.begin() ; i != points.end() ; ++i)
	{
			dofId[*i] = lastId++ ;
	}
	
	//now, we will start the enrichment itself
	

	std::set<std::pair<DelaunayTriangle *, Point *> > done ;
	//then we iterate on every element
	for(size_t i = 0 ; i < ring.size() ; i++)
	{

		std::vector<Point> hint ;
		hint.push_back(Point(1./3., 1./3.)) ;
		//if the number of intersection points is not 2, we need not do anything

		//we build the enrichment function, first, we get the transforms from the triangle
		Function x = ring[i]->getXTransform() ;
		Function y = ring[i]->getYTransform() ;
		
		//this function returns the distance to the centre
		Function position(getCenter(), x, y) ;
		
		Function hat = f_abs((getRadius() - position)) ;
		for(int j = 0 ; j < ring[i]->getBoundingPoints().size() ; j++)
		{
			std::pair<DelaunayTriangle *, Point *> that(ring[i],&ring[i]->getBoundingPoint(j)) ;
			
			if(done.find(that) != done.end() )
				continue ;
			done.insert(that) ;
			Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;
			Function f = ring[i]->getShapeFunction(j)*(hat - VirtualMachine().eval(hat,p.x, p.y )) ;
			
// 			for(double k = 0 ; k < 1 ; k += 0.01)
// 			{
// 				for(double l = 0 ; l < 1 ; l += 0.01)
// 				{
// 					if(k+l < 1)
// 						std::cout << VirtualMachine().eval(f, k, l) << "  "<<std::flush ;
// 					else
// 						std::cout << "0  "<<std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// 			exit(0) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(&ring[i]->getBoundingPoint(j)) ;
			f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
			ring[i]->setEnrichment( f, getPrimitive()) ;
		}
		
		for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
		{
			DelaunayTriangle * t = ring[i]->getNeighbourhood(j) ;

			if(!t->enrichmentUpdated)
				t->clearEnrichment( getPrimitive()) ;
			t->enrichmentUpdated = true ;
			bool hinted = false ;
			Function hat = f_abs((getRadius() - Function(getCenter(), 
																				t->getXTransform(), t->getYTransform()))) ;
			
			for(int k = 0 ; k < t->getBoundingPoints().size() ; k++)
			{
				std::map<Point*, int>::const_iterator pt = dofId.find(&t->getBoundingPoint(k)) ;
				if(pt != dofId.end())
				{
					std::pair<DelaunayTriangle *, Point *> that(t,&t->getBoundingPoint(k)) ;
					if(done.find(that) != done.end())
						continue ;
					
					Point p = t->inLocalCoordinates(t->getBoundingPoint(k)) ;
					done.insert(that) ;
					Function f = t->getShapeFunction(k)*(hat - VirtualMachine().eval(hat, p.x, p.y )) ;
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
	
bool EnrichmentInclusion::interacts(Feature * f, double d) const { return false ;}
void EnrichmentInclusion::snap(DelaunayTree * dtree) {}
	
bool EnrichmentInclusion::inBoundary(const Point v) const {return false ; }
bool EnrichmentInclusion::inBoundary(const Point *v) const { return false ;}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getElements2D( FeatureTree * dt) 
{ 
	return dt->getElements2D(getPrimitive()) ;
}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getIntersectingTriangles( FeatureTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dt->getElements2D(getPrimitive()) ;
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		Segment s0(*disc[i]->first, *disc[i]->second) ;
		Segment s1(*disc[i]->second, *disc[i]->third) ;
		Segment s2(*disc[i]->third, *disc[i]->first) ;
		if(s0.intersects(getPrimitive()) || s1.intersects(getPrimitive()) || s2.intersects(getPrimitive()))
			ring.push_back(disc[i]) ;
	}
	
	return ring ;
}
	
void EnrichmentInclusion::setInfluenceRadius(double r) { }
	
std::vector<Geometry *> EnrichmentInclusion::getRefinementZones( size_t level) const 
{
	return std::vector<Geometry *>(0) ;
}
	
void EnrichmentInclusion::step(double dt, std::valarray< double >*, const Mu::Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) {}
	
bool EnrichmentInclusion::moved() const {  return updated ;}

void EnrichmentInclusion::setRadius(double newR)
{
	Circle::setRadius(newR) ;
	updated = true ;
//	std::cout << "updating radius" << std::endl ;
}
