// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion.h"

using namespace Mu ;

EnrichmentInclusion::EnrichmentInclusion(Feature *father, double radius, double x, double y) : EnrichmentFeature(father), Circle(radius, x, y)
{
	
}

EnrichmentInclusion::EnrichmentInclusion(double radius, double x, double y) :EnrichmentFeature(NULL), Circle(radius, x, y)
{
	
}

EnrichmentInclusion::~EnrichmentInclusion() {}
	
bool EnrichmentInclusion::enrichmentTarget(DelaunayTriangle * t)
{
	return static_cast<Triangle *>(t)->intersection(static_cast<Circle *>(this)).size() == 2 ;
}

void EnrichmentInclusion::enrich(size_t & counter,  DelaunayTree * dtree)
{
	dtree->getTriangles() ;
	//first we get All the triangles affected
	std::vector<DelaunayTriangle *> disc = dtree->conflicts(static_cast<Circle *>(this)) ;

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
	
	//first, we get the basis functions for an archetypal element
	std::valarray<Function> shapefunc = TriElement(LINEAR).getShapeFunctions() ;
	
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
			Point mid = (triCircleIntersectionPoints[1]+triCircleIntersectionPoints[0])/2 ;
			Point q1 = (triCircleIntersectionPoints[1]+mid)/2 ;
			Point q4 = (mid+triCircleIntersectionPoints[0])/2 ;
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
			
			if(this->in(*a) && this->in(*b))
			{

				hint.push_back((ring[i]->inLocalCoordinates(*a) +
				                ring[i]->inLocalCoordinates(*b) +
				               hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*c)
			                   )/3
				              ) ;
			}
			if(this->in(*a) && this->in(*c))
			{
				hint.push_back((ring[i]->inLocalCoordinates(*a) +
				                ring[i]->inLocalCoordinates(*c) +
				                hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*b)
				               )/3
				              ) ;
			}
			if(this->in(*b) && this->in(*c))
			{
				hint.push_back((ring[i]->inLocalCoordinates(*b) +
				                ring[i]->inLocalCoordinates(*c) +
				                hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*a)
				               )/3
				              ) ;
			}
			if(!this->in(*a) && !this->in(*b))
			{
				hint.push_back((ring[i]->inLocalCoordinates(*b) +
				                ring[i]->inLocalCoordinates(*a) +
				                hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*c)
				               )/3
				              ) ;
			}
			if(!this->in(*a) && !this->in(*c))
			{
				hint.push_back((ring[i]->inLocalCoordinates(*c) +
				                ring[i]->inLocalCoordinates(*a) +
				                hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*b)
				               )/3
				              ) ;

			}
			if(!this->in(*b) && !this->in(*c))
			{
				hint.push_back((ring[i]->inLocalCoordinates(*c) +
				                ring[i]->inLocalCoordinates(*b) +
				                hint[0] +
				                hint[1])/4
				              ) ;
				hint.push_back((hint[0] +
				                hint[1] +
				                ring[i]->inLocalCoordinates(*a)
				               )/3
				              ) ;
			}
			
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
					hint.push_back(Point(0, .5)) ;
					hint.push_back(Point(.5, .5)) ;
					hint.push_back(Point(.5, 0)) ;
					
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
	
bool EnrichmentInclusion::moved() const { return false ;}
