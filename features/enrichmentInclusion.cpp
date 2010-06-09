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

void EnrichmentInclusion::enrich(size_t & counter, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	counter++ ;
	if(updated)
		update(dtree) ;
	updated = false ;
	
	const std::vector<DelaunayTriangle *> & disc  = cache;
// 	if(disc.size() < 6)
// 	{
// 		std::cout << "cowardly discarding" << std::endl ;
// 		return ;
// 	}

	std::valarray<Function> shapefunc = TriElement(LINEAR).getShapeFunctions() ;
	
	if(disc.size() == 1) // special case for really small inclusions
	{

		Inclusion * circle = new Inclusion(this->getRadius(), this->getCenter()) ;
		circle->setBehaviour(this->getBehaviour()->getCopy()) ;
		std::vector<Feature *> feat ;
		feat.push_back(circle) ;

		HomogeneisedBehaviour * hom = new HomogeneisedBehaviour(feat, disc[0]) ;

		disc[0]->setBehaviour(hom->getCopy()) ;
		return ;

/*		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			std::vector<Point> samplingPoints = getPrimitive()->getSamplingBoundingPoints(8) ;
	
			std::map<Point *, int> dofId ;
			
			dofId[disc[i]->first] = counter++ ;
			dofId[disc[i]->second] = counter++ ;
			dofId[disc[i]->third] = counter++ ;
	
			std::vector<Point> hint ;
			for(size_t j = 0 ; j < 8 ; j++)
			{
				hint.push_back(disc[i]->inLocalCoordinates(samplingPoints[j])) ;
			}
			
						//we build the enrichment function, first, we get the transforms from the triangle
			Function x = disc[i]->getXTransform() ;
			Function y = disc[i]->getYTransform() ;
			
				//this function returns the distance to the centre
			Function position(getCenter(), x, y) ;
			
				//finaly, we have the enrichment function
			Function hat = 1-f_abs(position -radius)/radius;
			
				//enriching the first point
			Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->first) ;
			f.setDofID(dofId[disc[i]->first]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
			
				//enriching the second point
			f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->second) ;
			f.setDofID(dofId[disc[i]->second]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
			
				//enriching the third point
			f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->third) ;
			f.setDofID(dofId[disc[i]->third]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
			
			for(size_t j = 0 ; j < disc[i]->neighbourhood.size() ; j++)
			{
				DelaunayTriangle * t = disc[i]->getNeighbourhood(j) ;
				if(!enrichmentTarget(t))
				{
					if(!t->enrichmentUpdated)
						t->clearEnrichment( getPrimitive()) ;
					t->enrichmentUpdated = true ;
					bool hinted = false ;
					Function hat = 1- f_abs(Function(getCenter(), 
					                                 t->getXTransform(), t->getYTransform()) -radius)/radius ;
					std::vector<Point> hint;
					hint.push_back(Point(1./3., 1./3.)) ;
					
					if(dofId.find(t->first) != dofId.end())
					{
						Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->first) ;
						f.setDofID(dofId[t->first]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
					
					if(dofId.find(t->second) != dofId.end())
					{
						Function f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->second) ;
						f.setDofID(dofId[t->second]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
					
					if(dofId.find(t->third) != dofId.end())
					{
						Function f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->third) ;
						f.setDofID(dofId[t->third]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
				}
			}
			
		}
		return ;*/
	}

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		Segment s0(*disc[i]->first, *disc[i]->second) ;
		Segment s1(*disc[i]->second, *disc[i]->third) ;
		Segment s2(*disc[i]->third, *disc[i]->first) ;
		if(!(getPrimitive()->in(*disc[i]->first) && getPrimitive()->in(*disc[i]->second) && getPrimitive()->in(*disc[i]->third)))
		{
			ring.push_back(disc[i]) ;
		}
	}
	std::cout << "ring.size() " << ring.size() << " disk.size()" << disc.size() << std::endl ;
	//then we build a list of points to enrich
	std::vector<Point *> points ;

	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		points.push_back(ring[i]->first) ;
		points.push_back(ring[i]->second) ;
		points.push_back(ring[i]->third) ;
	}
	
	//we make the points in the list unique
// 	std::stable_sort(points.begin(), points.end()) ;
	std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	//we build a map of the points and corresponding enrichment ids
	std::map<Point *, int> dofId ;
	
	for(size_t i = 0 ; i< points.size() ; i++)
	{
		if(freeIds.empty())
			dofId[points[i]] = counter++ ;
		else
		{
			dofId[points[i]] = *freeIds.begin() ;
			freeIds.erase(freeIds.begin()) ;
		}
	}
	
	//now, we will start the enrichment itself
	

	
	//then we iterate on every element
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		std::vector<Point> triCircleIntersectionPoints = getPrimitive()->intersection(static_cast<Triangle *>(ring[i])) ;

		//for convenience
		Point *a = ring[i]->first ;
		Point *b = ring[i]->second ;
		Point *c = ring[i]->third ;

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
			
			//finaly, we have the enrichment function
			Function hat = 1- f_abs(position -radius)/radius;
			
			//enriching the first point
			Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(a) ;
			f.setDofID(dofId[a]) ;
			ring[i]->setEnrichment( f, getPrimitive()) ;
			
			//enriching the second point
			f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
			f.setPoint(b) ;
			f.setDofID(dofId[b]) ;
			ring[i]->setEnrichment( f, getPrimitive()) ;
			
			//enriching the third point
			f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
			f.setPoint(c) ;
			f.setDofID(dofId[c]) ;
			ring[i]->setEnrichment(f, getPrimitive()) ;
			for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
			{
				DelaunayTriangle * t = ring[i]->getNeighbourhood(j) ;
				if(std::find(ring.begin(), ring.end(), t) == ring.end())
				{
					if(!t->enrichmentUpdated)
						t->clearEnrichment( getPrimitive()) ;
					t->enrichmentUpdated = true ;
					bool hinted = false ;
					Function hat = 1- f_abs(Function(getCenter(), 
					                                 t->getXTransform(), t->getYTransform()) -radius)/radius ;
					std::vector<Point> hint;
					hint.push_back(Point(1./3., 1./3.)) ;
					
					if(dofId.find(t->first) != dofId.end())
					{
						Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->first) ;
						f.setDofID(dofId[t->first]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
					
					if(dofId.find(t->second) != dofId.end())
					{
						Function f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->second) ;
						f.setDofID(dofId[t->second]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
					
					if(dofId.find(t->third) != dofId.end())
					{
						Function f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(t->third) ;
						f.setDofID(dofId[t->third]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
				}
			}
		
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
	
void EnrichmentInclusion::step(double dt, std::valarray<double> *, const DelaunayTree * dtree) {}
	
bool EnrichmentInclusion::moved() const { return updated ;}

void EnrichmentInclusion::setRadius(double newR)
{
	Circle::setRadius(newR) ;
	updated = true ;
}
