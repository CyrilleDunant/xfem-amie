// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentRing.h"
#include "../polynomial/vm_function_extra.h"
using namespace Amie ;

EnrichmentRing::EnrichmentRing(Feature *father, double radius, double inradius, double x, double y) : EnrichmentFeature(father), Circle(radius, x, y), self(inradius, x, y)
{
	updated = true ;
}

EnrichmentRing::EnrichmentRing(double radius, double inradius, double x, double y) :EnrichmentFeature(nullptr), Circle(radius, x, y), self(inradius, x, y)
{
	updated = true ;
}

EnrichmentRing::~EnrichmentRing() {}
	
bool EnrichmentRing::enrichmentTarget(DelaunayTriangle * t)
{
	return (!(getPrimitive()->in(*t->first) && getPrimitive()->in(*t->second) && getPrimitive()->in(*t->third)) || !(self.in(*t->first) && self.in(*t->second) && self.in(*t->third)) )&& t->getBehaviour()->type != VOID_BEHAVIOUR ;
	return t->intersects(getPrimitive())||t->intersects(&self)  ;
}

void EnrichmentRing::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
		{
			cache[i]->clearEnrichment(getPrimitive()) ;
			cache[i]->clearEnrichment(&self) ;
		}
		cache[i]->enrichmentUpdated = true ;
	}
	cache = dtree->getConflictingElements(getPrimitive()) ;
	std::vector<DelaunayTriangle *> inElements = dtree->getConflictingElements(&self) ;
	cache.insert(cache.end(), inElements.begin(), inElements.end()) ;
	std::sort(cache.begin(), cache.end()) ;
	auto e = std::unique(cache.begin(), cache.end()) ;
	cache.erase(e, cache.end()) ;
	if(cache.empty())
	{
		std::vector<DelaunayTriangle *> candidates = dtree->getConflictingElements(&getCenter()) ;
		for(size_t i = 0 ; i < candidates.size() ; i++)
		{
			if(candidates[i]->isTriangle && candidates[i]->in(getCenter()))
			{
				cache.push_back(candidates[i]) ;
				break ;
			}
		}
	}
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
		{
			cache[i]->clearEnrichment(getPrimitive()) ;
			cache[i]->clearEnrichment(&self) ;
		}
		cache[i]->enrichmentUpdated = true ;
	}
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void EnrichmentRing::enrich(size_t & lastId,  Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
	if(updated)
		update(dtree) ;
	updated = false ;

	
	const std::vector<DelaunayTriangle *> & disc  = cache;

	std::valarray<Function> shapefunc = TriElement(LINEAR).getShapeFunctions() ;
	
	if(disc.size() < 6) // special case for really small inclusions
	{
		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			std::vector<Point> samplingPoints = getSamplingBoundingPoints(8) ;
	
			std::map<Point *, int> dofId ;
			
			dofId[disc[i]->first] = lastId++ ;
			dofId[disc[i]->second] = lastId++ ;
			dofId[disc[i]->third] = lastId++ ;
	
			std::vector<Point> hint ;
			for(size_t j = 0 ; j < 8 ; j++)
			{
				hint.push_back(disc[i]->inLocalCoordinates(samplingPoints[j])) ;
			}
			
			//we build the enrichment function, first, we get the transforms from the triangle
			//this function returns the distance to the centre
			Function position(getCenter(),disc[i]) ;
			
				//finaly, we have the enrichment function
			Function hat = Function("1")-f_abs(position -radius)/radius;
			
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
		}
		return ;
	}


	//we build a list of the two circles:

	std::vector<const Circle *> circles ;
	circles.push_back(getPrimitive()) ;
	circles.push_back(&self) ;
	for(size_t h = 0 ; h < 2 ; h++)
	{
		//then we select those that are cut by the circle
		std::vector<DelaunayTriangle *> ring ;
		
		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			if(!(!circles[h]->in(*disc[i]->first) && !circles[h]->in(*disc[i]->second) && !circles[h]->in(*disc[i]->third)) && !(circles[h]->in(*disc[i]->first) && circles[h]->in(*disc[i]->second) && circles[h]->in(*disc[i]->third)) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				ring.push_back(disc[i]) ;
		}
		
		//then we build a list of points to enrich
		std::vector<Point *> points ;
	
		for(size_t i = 0 ; i < ring.size() ; i++)
		{
			points.push_back(ring[i]->first) ;
			points.push_back(ring[i]->second) ;
			points.push_back(ring[i]->third) ;
		}
		
		//we make the points in the list unique
		std::stable_sort(points.begin(), points.end()) ;
		auto e = std::unique(points.begin(), points.end()) ;
		points.erase(e, points.end()) ;
		
		//we build a map of the points and corresponding enrichment ids
		std::map<Point *, int> dofId ;
		
		for(size_t i = 0 ; i< points.size() ; i++)
		{
			dofId[points[i]] =lastId++ ;
		}
		//now, we will start the enrichment itself
		
	
		
		//then we iterate on every element
		for(size_t i = 0 ; i < ring.size() ; i++)
		{
			std::vector<Point> triCircleIntersectionPoints;
			if(circles[h]->intersects(static_cast<Triangle *>(ring[i])))
				triCircleIntersectionPoints = circles[h]->intersection(static_cast<Triangle *>(ring[i])) ;
			
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
					double d = dist(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[0]),triCircleIntersectionPoints[1]) ;
					double dr = std::abs(getRadius()-self.getRadius()) ;
					int numPoints = 2.*round(d/dr) ;
					if(dr > POINT_TOLERANCE && numPoints < 32)
					{
						hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[0])) ;
						std::vector<Point> pts = circles[h]->getSamplingBoundingPointsOnArc(numPoints,triCircleIntersectionPoints[0],triCircleIntersectionPoints[1]  ) ;
						for(size_t k = 0 ; k < pts.size() ;k++)
						{
							if(ring[i]->in(pts[k]))
							{
								hint.push_back(ring[i]->inLocalCoordinates(pts[k])) ;
							}
						}
						hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[1])) ;
					}
					//else the ring is too thin. There is no point to try and mesh it

				}

				//this function returns the distance to the centre
				Function position(getCenter(), ring[i]) ;
				
				//finaly, we have the enrichment function
				Function hatOut = Function("1")-f_abs(position -radius)/radius;
				Function hatIn = Function("1")-f_abs(position - getInRadius())/getInRadius();
				Function hat ;

				if(h == 0)
				{
					hat = hatOut  ;
				}
				else
				{
					hat = hatIn ;
				}
				//enriching the first point
				Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
				f.setIntegrationHint(hint) ;
				f.setPoint(a) ;
				f.setDofID(dofId[a]) ;
				ring[i]->setEnrichment( f, circles[h]) ;
				
				//enriching the second point
				f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
				f.setPoint(b) ;
				f.setDofID(dofId[b]) ;
				ring[i]->setEnrichment( f, circles[h]) ;
				
				//enriching the third point
				f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
				f.setPoint(c) ;
				f.setDofID(dofId[c]) ;
				ring[i]->setEnrichment(f, circles[h]) ;
                std::vector<DelaunayTriangle * > neighbourhood = dtree->getNeighbourhood(ring[i]) ;
				for(auto & t : neighbourhood)
				{
					if(((circles[h]->in(*disc[i]->first) && circles[h]->in(*disc[i]->second) && circles[h]->in(*disc[i]->third)) || (!circles[h]->in(*disc[i]->first) && !circles[h]->in(*disc[i]->second) && !circles[h]->in(*disc[i]->third)))  && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						if(!t->enrichmentUpdated)
							t->clearEnrichment(circles[h]) ;
						t->enrichmentUpdated = true ;
						
						Function hatOut = Function("1")-f_abs(Function(getCenter(),t) -radius)/radius;
						Function hatIn = Function("1")-f_abs(Function(getCenter(),t) - getInRadius())/getInRadius();
						Function hat ;

						if(h == 0)
						{
							hat = hatOut  ;
						}
						else
						{
							hat = hatIn ;
						}
						
						std::vector<Point> hint;
						hint.push_back(Point(1./3., 1./3.)) ;
						
						bool hinted = false ;
						
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
							t->setEnrichment(f, circles[h]) ;
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
							t->setEnrichment(f, circles[h]) ;
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
							t->setEnrichment(f, circles[h]) ;
						}
					}
				}
			
			}
		}
	}
}
	
bool EnrichmentRing::interacts(Feature * f, double d) const { return false ;}
void EnrichmentRing::snap(DelaunayTree * dtree) {}
	
bool EnrichmentRing::inBoundary(const Point v) const {return false ; }
bool EnrichmentRing::inBoundary(const Point *v) const { return false ;}
	
std::vector<DelaunayTriangle *> EnrichmentRing::getElements2D( FeatureTree * dt) 
{ 
	return dt->get2DMesh()->getConflictingElements(getPrimitive()) ;
}
	
std::vector<DelaunayTriangle *> EnrichmentRing::getBoundingElements2D( FeatureTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dt->get2DMesh()->getConflictingElements(static_cast<Circle *>(this)) ;
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(disc[i]->getPrimitive()) || self.intersects(disc[i]->getPrimitive()))
			ring.push_back(disc[i]) ;
	}
	
	return ring ;
}
	
void EnrichmentRing::setInfluenceRadius(double r) { }
	
std::vector<Geometry *> EnrichmentRing::getRefinementZones( size_t level) const 
{
	return std::vector<Geometry *>(0) ;
}
	
void EnrichmentRing::step(double dt, std::valarray<double> *, const Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree) {}
	
bool EnrichmentRing::moved() const { return updated ;}

void EnrichmentRing::setRadius(double newR)
{
	Circle::setRadius(newR) ;
	updated = true ;
}

void EnrichmentRing::setInRadius(double newR)
{
	self.setRadius(newR) ;
	updated = true ;
}

double EnrichmentRing::getInRadius() const
{
	return self.getRadius() ;
}
