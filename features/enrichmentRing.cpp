// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentRing.h"

using namespace Mu ;

EnrichmentRing::EnrichmentRing(Feature *father, double radius, double inradius, double x, double y) : EnrichmentFeature(father), Circle(radius, x, y), self(inradius, x, y)
{
	updated = true ;
}

EnrichmentRing::EnrichmentRing(double radius, double inradius, double x, double y) :EnrichmentFeature(NULL), Circle(radius, x, y), self(inradius, x, y)
{
	updated = true ;
}

EnrichmentRing::~EnrichmentRing() {}
	
bool EnrichmentRing::enrichmentTarget(DelaunayTriangle * t)
{
	return t->intersects(getPrimitive())||t->intersects(&self)  ;
}

void EnrichmentRing::update(DelaunayTree * dtree)
{
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
			cache[i]->clearEnrichment(static_cast<const Circle *>(this)) ;
		cache[i]->enrichmentUpdated = true ;
	}
	cache = dtree->conflicts(getPrimitive()) ;
	std::vector<DelaunayTriangle *> inElements = dtree->conflicts(&self) ;
	cache.insert(cache.end(), inElements.begin(), inElements.end()) ;
	if(cache.empty())
	{
		std::vector<DelaunayTreeItem *> candidates = dtree->conflicts(&getCenter()) ;
		for(size_t i = 0 ; i < candidates.size() ; i++)
		{
			if(candidates[i]->isTriangle && static_cast<DelaunayTriangle *>(candidates[i])->in(getCenter()))
			{
				cache.push_back(static_cast<DelaunayTriangle *>(candidates[i])) ;
				break ;
			}
		}
	}
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
			cache[i]->clearEnrichment(static_cast<const Circle *>(this)) ;
		cache[i]->enrichmentUpdated = true ;
	}
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void EnrichmentRing::enrich(size_t & counter,  DelaunayTree * dtree)
{
// 	counter++ ;
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
	
	if(disc.size() < 6) // special case for really small inclusions
	{
		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			std::vector<Point> samplingPoints = static_cast<Circle *>(this)->getSamplingBoundingPoints(8) ;
	
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
			disc[i]->setEnrichment(f, static_cast<Circle *>(this)) ;
			
				//enriching the second point
			f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->second) ;
			f.setDofID(dofId[disc[i]->second]) ;
			disc[i]->setEnrichment(f, static_cast<Circle *>(this)) ;
			
				//enriching the third point
			f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->third) ;
			f.setDofID(dofId[disc[i]->third]) ;
			disc[i]->setEnrichment(f, static_cast<Circle *>(this)) ;
		}
		return ;
	}


	//we build a list of the two circles:

	std::vector<const Circle *> circles ;
	circles.push_back(this->getPrimitive()) ;
	circles.push_back(&self) ;
	for(size_t h = 0 ; h < 2 ; h++)
	{
		//then we select those that are cut by the circle
		std::vector<DelaunayTriangle *> ring ;
		
		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			if(circles[h]->intersects(static_cast<Triangle *>(disc[i])))
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
	// 	std::stable_sort(points.begin(), points.end()) ;
		std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
		points.erase(e, points.end()) ;
		
		//we build a map of the points and corresponding enrichment ids
		std::map<Point *, int> dofId ;
		
		for(size_t i = 0 ; i< points.size() ; i++)
		{
			dofId[points[i]] = counter++ ;
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
					if(numPoints > 8)
						numPoints = 8 ;
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
				//we build the enrichment function, first, we get the transforms from the triangle
				Function x = ring[i]->getXTransform() ;
				Function y = ring[i]->getYTransform() ;
				//this function returns the distance to the centre
				Function position(getCenter(), x, y) ;
				
				//finaly, we have the enrichment function
				Function hatOut = (1-f_abs(position -radius)/radius);
				Function hatIn = (1-f_abs(position - getInRadius())/getInRadius());
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
				ring[i]->setEnrichment( f, static_cast<Circle *>(this)) ;
				
				//enriching the second point
				f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
				f.setPoint(b) ;
				f.setDofID(dofId[b]) ;
				ring[i]->setEnrichment( f, static_cast<Circle *>(this)) ;
				
				//enriching the third point
				f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
				f.setPoint(c) ;
				f.setDofID(dofId[c]) ;
				ring[i]->setEnrichment(f, static_cast<Circle *>(this)) ;
				for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
				{
					DelaunayTriangle * t = ring[i]->getNeighbourhood(j) ;
					if(!enrichmentTarget(t))
					{
						if(!t->enrichmentUpdated)
							t->clearEnrichment(static_cast<Circle *>(this)) ;
						t->enrichmentUpdated = true ;
						Function hat = 1- f_abs(Function(getCenter(), 
										t->getXTransform(), t->getYTransform()) -radius)/radius ;
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
							t->setEnrichment(f, static_cast<Circle *>(this)) ;
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
							t->setEnrichment(f, static_cast<Circle *>(this)) ;
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
							t->setEnrichment(f, static_cast<Circle *>(this)) ;
						}
					}
				}
			
			}
		}
	}
}
	
bool EnrichmentRing::interacts(Feature * f) const { return false ;}
void EnrichmentRing::snap(DelaunayTree * dtree) {}
	
bool EnrichmentRing::inBoundary(const Point v) const {return false ; }
bool EnrichmentRing::inBoundary(const Point *v) const { return false ;}
	
std::vector<DelaunayTriangle *> EnrichmentRing::getTriangles( DelaunayTree * dt) 
{ 
	return dt->conflicts(static_cast<Circle *>(this)) ;
}
	
std::vector<DelaunayTriangle *> EnrichmentRing::getIntersectingTriangles( DelaunayTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dt->conflicts(static_cast<Circle *>(this)) ;
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Triangle *>(disc[i])) || self.intersects(static_cast<Triangle *>(disc[i])))
			ring.push_back(disc[i]) ;
	}
	
	return ring ;
}
	
void EnrichmentRing::setInfluenceRadius(double r) { }
	
std::vector<Geometry *> EnrichmentRing::getRefinementZones( size_t level) const 
{
	return std::vector<Geometry *>(0) ;
}
	
void EnrichmentRing::step(double dt, std::valarray<double> *, const DelaunayTree * dtree) {}
	
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
