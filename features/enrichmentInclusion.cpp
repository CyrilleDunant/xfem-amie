// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion.h"
#include "inclusion.h"
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
	freeIds.clear();
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
		{
			std::vector<size_t> idpts = cache[i]->clearEnrichment(getPrimitive()) ;
			for(size_t j = 0 ; j < idpts.size() ; j++)
				freeIds.insert(idpts[j]) ;
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
				freeIds.insert(idpts[j]) ;
		}
		cache[i]->enrichmentUpdated = true ;

	}
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void EnrichmentInclusion::enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	if(updated)
		update(dtree) ;
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

		Inclusion * circle = new Inclusion(this->getRadius(), this->getCenter()) ;
		circle->setBehaviour(this->getBehaviour()->getCopy()) ;
		std::vector<Feature *> feat ;
		feat.push_back(circle) ;

		HomogeneisedBehaviour hom(feat, disc[0]) ;

		disc[0]->setBehaviour(hom.getCopy()) ;
		return ;


	}

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
	double minrad = getRadius()*2 ;
	double maxrad = 0 ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		for(size_t j = 0 ; j < ring[i]->getBoundingPoints().size() ; j++)
		{
			points.insert(&ring[i]->getBoundingPoint(j)) ;
			minrad = std::min(minrad, dist(getCenter(), ring[i]->getBoundingPoint(j))) ;
			maxrad = std::max(maxrad, dist(getCenter(), ring[i]->getBoundingPoint(j))) ;
		}
	}
	double del = std::min(std::abs(minrad-radius), std::abs(maxrad-radius)) ;
	minrad = radius-del/4 ;
	maxrad = radius+del/4 ;
	
	//we build a map of the points and corresponding enrichment ids
	std::map<Point *, int> dofId ;
	
	for(std::set<Point * >::const_iterator i = points.begin() ; i != points.end() ; ++i)
	{
		if(freeIds.empty())
			dofId[*i] = dtree->getLastNodeId()++ ;
		else
		{
			dofId[*i] = *freeIds.begin() ;
			freeIds.erase(freeIds.begin()) ;
		}
	}
	
	//now, we will start the enrichment itself
	

	
	//then we iterate on every element
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		std::vector<Point> triCircleIntersectionPoints = getPrimitive()->intersection(static_cast<Triangle *>(ring[i])) ;
		std::vector<Segment> insegs ;
		std::vector<Segment> outsegs ;
		
		if(getPrimitive()->in(*ring[i]->first) && !getPrimitive()->in(*ring[i]->second))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->first, *ring[i]->second).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->first));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->second));
		}
		if(!getPrimitive()->in(*ring[i]->first) && getPrimitive()->in(*ring[i]->second))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->first, *ring[i]->second).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->second));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->first));
		}
		if(getPrimitive()->in(*ring[i]->second) && !getPrimitive()->in(*ring[i]->third))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->second, *ring[i]->third).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->second));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->third));
		}
		if(!getPrimitive()->in(*ring[i]->second) && getPrimitive()->in(*ring[i]->third))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->second, *ring[i]->third).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->third));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->second));
		}
		if(getPrimitive()->in(*ring[i]->third) && !getPrimitive()->in(*ring[i]->first))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->third, *ring[i]->first).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->third));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->first));
		}
		if(!getPrimitive()->in(*ring[i]->third) && getPrimitive()->in(*ring[i]->first))
		{
			std::vector<Point> intersectors = Segment(*ring[i]->third, *ring[i]->first).intersection(getPrimitive()) ;
			insegs.push_back(Segment(intersectors[0], *ring[i]->first));
			outsegs.push_back(Segment(intersectors[0], *ring[i]->third));
		}
		
			
		//if there are no intersection points we need not do anything
		if(!triCircleIntersectionPoints.empty())
		{
			std::vector<Point> hint ;
			//if the number of intersection points is not 2, we need not do anything
			if(triCircleIntersectionPoints.size() == 2)
			{
				double d = dist(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[0]), ring[i]->inLocalCoordinates(triCircleIntersectionPoints[1])) ;
				int npts = (int)std::max(4.,round(8.*d)) ;
				std::vector<Point> pts = getPrimitive()->getSamplingBoundingPointsOnArc(npts,triCircleIntersectionPoints[0],triCircleIntersectionPoints[1]  ) ;
				
				hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[0])) ;
				for(size_t k = 0 ; k < pts.size() ;k++)
				{
					if(ring[i]->in(pts[k]))
					{
						hint.push_back(ring[i]->inLocalCoordinates(pts[k])) ;
					}
				}
				hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[1])) ;
			}

			//we build the enrichment function, first, we get the transforms from the triangle
			Function x = ring[i]->getXTransform() ;
			Function y = ring[i]->getYTransform() ;
			
			//this function returns the distance to the centre
			Function position(getCenter(), x, y) ;
			
			Function hat = 1-f_abs((radius - position)/(maxrad-minrad)) ;
			for(int j = 0 ; j < ring[i]->getBoundingPoints().size() ; j++)
			{
				Function f = father.getShapeFunction(j)*(hat - VirtualMachine().eval(hat, ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)))) ;
				f.setIntegrationHint(hint) ;
				f.setPoint(&ring[i]->getBoundingPoint(j)) ;
				f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
				ring[i]->setEnrichment( f, getPrimitive()) ;

				
			}
			
			for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
			{
				DelaunayTriangle * t = ring[i]->getNeighbourhood(j) ;
				if(std::find(ring.begin(), ring.end(), t) == ring.end())
				{
					if(!t->enrichmentUpdated)
						t->clearEnrichment( getPrimitive()) ;
					t->enrichmentUpdated = true ;
					bool hinted = false ;
					Function hat = 1-f_abs((radius - Function(getCenter(), 
					                                 t->getXTransform(), t->getYTransform()))/(maxrad-minrad)) ;

					std::vector<Point> hint;
					hint.push_back(Point(1./3., 1./3.)) ;
					for(int k = 0 ; k < t->getBoundingPoints().size() ; k++)
					{
						std::map<Point*, int>::const_iterator pt = dofId.find(&t->getBoundingPoint(k)) ;
						if(pt != dofId.end())
						{
							int delta = 0 ;
							while( &t->getBoundingPoint(delta) != pt->first)
								delta++ ;
							Function f = father.getShapeFunction(delta)*(hat - VirtualMachine().eval(hat, t->inLocalCoordinates(t->getBoundingPoint(delta)))) ;
							if(!hinted)
							{
								f.setIntegrationHint(hint) ;
								hinted = true ;
							}
							f.setPoint(pt->first) ;
							f.setDofID(dofId[pt->first]) ;
							t->setEnrichment(f, getPrimitive()) ;
						}
					}
				}
			}
// 		
		}
	}

}
	
bool EnrichmentInclusion::interacts(Feature * f, double d) const { return false ;}
void EnrichmentInclusion::snap(DelaunayTree * dtree) {}
	
bool EnrichmentInclusion::inBoundary(const Point v) const {return false ; }
bool EnrichmentInclusion::inBoundary(const Point *v) const { return false ;}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getTriangles( DelaunayTree * dt) 
{ 
	return dt->conflicts(getPrimitive()) ;
}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getIntersectingTriangles( DelaunayTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dt->conflicts(getPrimitive()) ;
	
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
	
bool EnrichmentInclusion::moved() const { return updated ;}

void EnrichmentInclusion::setRadius(double newR)
{
	Circle::setRadius(newR) ;
	updated = true ;
}
