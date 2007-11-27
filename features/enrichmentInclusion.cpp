// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion.h"

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
	return static_cast<Triangle *>(t)->intersection(static_cast<Circle *>(this)).size() == 2 ;
}

void EnrichmentInclusion::update(DelaunayTree * dtree)
{
	if(cache.empty())
		cache = dtree->conflicts(static_cast<Circle *>(this)) ;
	else
	{
		std::vector<DelaunayTriangle *> temp ;
		std::vector<DelaunayTriangle *> cleanup ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			cache[i]->visited = true ;
			if(cache[i]->isConflicting(static_cast<Circle *>(this)))
				temp.push_back(cache[i]) ;
			
			cleanup.push_back(cache[i]) ;
		}
		
		
		std::vector<DelaunayTriangle *> toCheck ;
		
		for(size_t i = 0 ; i < temp.size() ; i++)
		{

			for(size_t j = 0 ; j < temp[i]->neighbour.size() ; j++)
			{
				if(!temp[i]->neighbour[j]->visited && temp[i]->neighbour[j]->isTriangle)
				{
					cleanup.push_back(static_cast<DelaunayTriangle *>(temp[i]->neighbour[j])) ;
					toCheck.push_back(static_cast<DelaunayTriangle *>(temp[i]->neighbour[j])) ;
				}
			}
		}
		
		std::sort(toCheck.begin(), toCheck.end()) ;
		std::vector<DelaunayTriangle *>::iterator uend = std::unique(toCheck.begin(), toCheck.end()) ;
		
		
		while(!toCheck.empty())
		{
			std::vector<DelaunayTriangle *> newSet ;
			
			for(std::vector<DelaunayTriangle *>::iterator i = toCheck.begin() ; i != uend ; ++i)
			{
				(*i)->visited = true ;
				cleanup.push_back(*i) ;
				
				if((*i)->isConflicting(static_cast<Circle *>(this)))
				{
					temp.push_back(*i) ;
				}
			}
			
			for(std::vector<DelaunayTriangle *>::iterator i = toCheck.begin() ; i != uend ; ++i)
			{
				
				for(size_t j = 0 ; j< (*i)->neighbour.size() ; j++)
				{
					if((*i)->neighbour[j]->isTriangle && !(*i)->neighbour[j]->visited )
					{
						newSet.push_back(static_cast<DelaunayTriangle *>((*i)->neighbour[j])) ;
					}
				}
			}
			
			toCheck.clear() ;
			std::sort(newSet.begin(), newSet.end()) ;
			uend = std::unique(newSet.begin(), newSet.end()) ;
			toCheck.insert(toCheck.end(),newSet.begin(),uend ) ;
			uend = toCheck.end() ;
			
		}
		
		for(size_t i = 0 ; i < cleanup.size() ; i++)
			cleanup[i]->visited = false ;
		
		std::sort(temp.begin(), temp.end()) ;
		std::vector<DelaunayTriangle *>::iterator en = std::unique(temp.begin(), temp.end()) ;
		temp.erase(en, temp.end()) ;
		
		cache = temp ;
	}
}

void EnrichmentInclusion::enrich(size_t & counter,  DelaunayTree * dtree)
{
	updated = false ;
// 	dtree->getTriangles() ;
	//first we get All the triangles affected
	update(dtree) ;
	const std::vector<DelaunayTriangle *> & disc  = cache;
// 	if(cache.empty())
// 		disc = dtree->conflicts(static_cast<Circle *>(this)) ;
// 	cache = disc ;
	//first, we get the basis functions for an archetypal element
	std::valarray<Function> shapefunc = TriElement(LINEAR).getShapeFunctions() ;
	
	if(disc.size() == 1) // special case for really small inclusions
	{
		static_cast<Circle *>(this)->sampleBoundingSurface(8) ;

		std::map<Point *, int> dofId ;
		
		dofId[disc[0]->first] = counter++ ;
		dofId[disc[0]->second] = counter++ ;
		dofId[disc[0]->third] = counter++ ;

		std::vector<Point> hint ;
		for(size_t i = 0 ; i < 8 ; i++)
		{
			hint.push_back(disc[0]->inLocalCoordinates(getBoundingPoint(i))) ;
		}
		
					//we build the enrichment function, first, we get the transforms from the triangle
		Function x = disc[0]->getXTransform() ;
		Function y = disc[0]->getYTransform() ;
		
			//this function returns the distance to the centre
		Function position(getCenter(), x, y) ;
		
			//finaly, we have the enrichment function
		Function hat = 1- f_abs(position -this->getRadius());
		
			//enriching the first point
		Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
		f.setIntegrationHint(hint) ;
		f.setPointID(disc[0]->first->id) ;
		f.setDofID(dofId[disc[0]->first]) ;
		disc[0]->setEnrichment(std::pair<size_t, Function>(dofId[disc[0]->first], f)) ;
		
			//enriching the second point
		f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
		f.setIntegrationHint(hint) ;
		f.setPointID(disc[0]->second->id) ;
		f.setDofID(dofId[disc[0]->second]) ;
		disc[0]->setEnrichment(std::pair<size_t, Function>(dofId[disc[0]->second], f)) ;
		
			//enriching the third point
		f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
		f.setIntegrationHint(hint) ;
		f.setPointID(disc[0]->third->id) ;
		f.setDofID(dofId[disc[0]->third]) ;
		disc[0]->setEnrichment(std::pair<size_t, Function>(dofId[disc[0]->third], f)) ;
		return ;
	}

	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Triangle *>(disc[i])))
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
		std::vector<Point> triCircleIntersectionPoints = static_cast<Circle *>(this)->intersection(static_cast<Triangle *>(ring[i])) ;
		
		//for convenience
		Point *a = ring[i]->first ;
		Point *b = ring[i]->second ;
		Point *c = ring[i]->third ;

		//if the number of intersection points is not 2, we need not do anything
		if(triCircleIntersectionPoints.size() == 2)
		{
			//our baseline
			Line frontier(triCircleIntersectionPoints[0],
			              triCircleIntersectionPoints[1]-triCircleIntersectionPoints[0]) ;
			Point mid = (triCircleIntersectionPoints[1]+triCircleIntersectionPoints[0])*.5 ;
			Point q1 = (triCircleIntersectionPoints[1]+mid)*.5 ;
			Point q4 = (mid+triCircleIntersectionPoints[0])*.5 ;
			this->project(&mid) ;
			this->project(&q1) ;
			this->project(&q4) ;
			
			//for the enrichment, we want to specify points which will form the subtriangulation
			std::vector<Point> hint ;
			hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[0])) ;
			hint.push_back(ring[i]->inLocalCoordinates(triCircleIntersectionPoints[1])) ;
			hint.push_back(ring[i]->inLocalCoordinates(mid)) ;
			hint.push_back(ring[i]->inLocalCoordinates(q1)) ;
			hint.push_back(ring[i]->inLocalCoordinates(q4)) ;

			//we build the enrichment function, first, we get the transforms from the triangle
			Function x = ring[i]->getXTransform() ;
			Function y = ring[i]->getYTransform() ;
			
			//this function returns the distance to the centre
			Function position(getCenter(), x, y) ;
			
			//finaly, we have the enrichment function
			Function hat = 1- f_abs(position -this->getRadius());
			
			//enriching the first point
			Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
			f.setIntegrationHint(hint) ;
			f.setPointID(a->id) ;
			f.setDofID(dofId[a]) ;
			ring[i]->setEnrichment(std::pair<size_t, Function>(dofId[a], f)) ;
			
			//enriching the second point
			f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPointID(b->id) ;
			f.setDofID(dofId[b]) ;
			ring[i]->setEnrichment(std::pair<size_t, Function>(dofId[b], f)) ;
			
			//enriching the third point
			f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPointID(c->id) ;
			f.setDofID(dofId[c]) ;
			ring[i]->setEnrichment(std::pair<size_t, Function>(dofId[c], f)) ;
			for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
			{
				DelaunayTriangle * t = ring[i]->neighbourhood[j] ;
				if(!enrichmentTarget(t))
				{
					Function hat = 1- f_abs(Function(getCenter(), 
					                                 t->getXTransform(), t->getYTransform()) -this->getRadius()) ;
					std::vector<Point> hint;
					hint.push_back(Point(1./3., 1./3.)) ;
					
					if(dofId.find(t->first) != dofId.end())
					{
						Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(0,1))) ;
						f.setIntegrationHint(hint) ;
						f.setPointID(t->first->id) ;
						f.setDofID(dofId[t->first]) ;
						t->setEnrichment(std::pair<size_t, Function>(dofId[t->first], f)) ;
					}
					
					if(dofId.find(t->second) != dofId.end())
					{
						Function f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,0))) ;
						f.setIntegrationHint(hint) ;
						f.setPointID(t->second->id) ;
						f.setDofID(dofId[t->second]) ;
						t->setEnrichment(std::pair<size_t, Function>(dofId[t->second], f)) ;
					}
					
					if(dofId.find(t->third) != dofId.end())
					{
						Function f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(1,0))) ;
						f.setIntegrationHint(hint) ;
						f.setPointID(t->third->id) ;
						f.setDofID(dofId[t->third]) ;
						t->setEnrichment(std::pair<size_t, Function>(dofId[t->third], f)) ;
					}
				}
			}
		
		}
	}

}
	
bool EnrichmentInclusion::interacts(Feature * f) const { return false ;}
void EnrichmentInclusion::snap(DelaunayTree * dtree) {}
	
bool EnrichmentInclusion::inBoundary(const Point v) const {return false ; }
bool EnrichmentInclusion::inBoundary(const Point *v) const { return false ;}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getTriangles( DelaunayTree * dt) 
{ 
	return dt->conflicts(static_cast<Circle *>(this)) ;
}
	
std::vector<DelaunayTriangle *> EnrichmentInclusion::getIntersectingTriangles( DelaunayTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dt->conflicts(static_cast<Circle *>(this)) ;
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTriangle *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Triangle *>(disc[i])))
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
