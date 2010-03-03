// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion3d.h"

using namespace Mu ;

EnrichmentInclusion3D::EnrichmentInclusion3D(Feature *father, double radius, double x, double y, double z) : EnrichmentFeature(father), Sphere(radius, x, y, z)
{
	updated = true ;
}

EnrichmentInclusion3D::EnrichmentInclusion3D(double radius, double x, double y, double z) :EnrichmentFeature(NULL), Sphere(radius, x, y, z)
{
	updated = true ;
}

EnrichmentInclusion3D::~EnrichmentInclusion3D() {}
	
bool EnrichmentInclusion3D::enrichmentTarget(DelaunayTetrahedron * t)
{
	int pointsin = in(*t->first) + in(*t->second) + in(*t->third) + in(*t->fourth) ;
		if(intersects(t->getPrimitive()) && pointsin != 4)
			return true ;
		
		return false ;
}

void EnrichmentInclusion3D::update(DelaunayTree3D * dtree)
{
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
			cache[i]->clearEnrichment(getPrimitive()) ;
		cache[i]->enrichmentUpdated = true ;
	}
	cache = dtree->conflicts(getPrimitive()) ;
	if(cache.empty())
	{
		std::vector<DelaunayTreeItem3D *> candidates = dtree->conflicts(&getCenter()) ;
		for(size_t i = 0 ; i < candidates.size() ; i++)
		{
			DelaunayTetrahedron * t = dynamic_cast<DelaunayTetrahedron *>(candidates[i]) ;
			if(t && (t->in(getCenter()) || this->intersects(t->getPrimitive())) )
			{
				cache.push_back(t) ;
				break ;
			}
		}
	}
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		if(!cache[i]->enrichmentUpdated)
			cache[i]->clearEnrichment(getPrimitive()) ;
		cache[i]->enrichmentUpdated = true ;

	}
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

void EnrichmentInclusion3D::enrich(size_t & counter,  DelaunayTree3D * dtree)
{
	counter++ ;
	if(updated)
		update(dtree) ;
	updated = false ;

	VirtualMachine vm ;
	const std::vector<DelaunayTetrahedron *> & disc  = cache;

	std::valarray<Function> shapefunc = TetrahedralElement(LINEAR).getShapeFunctions() ;
	
	if(disc.size() < 6) // special case for really small inclusions
	{
		for(size_t i = 0 ; i < disc.size() ; i++)
		{
			std::vector<Point> samplingPoints = getPrimitive()->getSamplingBoundingPoints(8) ;
	
			std::map<Point *, int> dofId ;
			
			dofId[disc[i]->first] = counter++ ;
			dofId[disc[i]->second] = counter++ ;
			dofId[disc[i]->third] = counter++ ;
			dofId[disc[i]->fourth] = counter++ ;
	
			std::vector<Point> hint ;
			for(size_t j = 0 ; j < 8 ; j++)
			{
				hint.push_back(disc[i]->inLocalCoordinates(samplingPoints[j])) ;
			}
			
						//we build the enrichment function, first, we get the transforms from the triangle
			Function x = disc[i]->getXTransform() ;
			Function y = disc[i]->getYTransform() ;
			Function z = disc[i]->getZTransform() ;
			
				//this function returns the distance to the centre
			Function position(getCenter(), x, y, z) ;
			
				//finaly, we have the enrichment function
			Function hat = 1-f_abs(position -radius)/radius;
			
				//enriching the first point
			Function f = shapefunc[0]*(hat - vm.eval(hat, Point(1,0,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->first) ;
			f.setDofID(dofId[disc[i]->first]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
			
				//enriching the second point
			f = shapefunc[1]*(hat - vm.eval(hat, Point(0,1,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->second) ;
			f.setDofID(dofId[disc[i]->second]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
			
				//enriching the third point
			f = shapefunc[2]*(hat - vm.eval(hat, Point(0,0,1))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->third) ;
			f.setDofID(dofId[disc[i]->third]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;

				//enriching the fourth point
			f = shapefunc[3]*(hat - vm.eval(hat, Point(0,0,0))) ;
			f.setIntegrationHint(hint) ;
			f.setPoint(disc[i]->third) ;
			f.setDofID(dofId[disc[i]->third]) ;
			disc[i]->setEnrichment(f, getPrimitive()) ;
		}
		return ;
	}

	//then we select those that are cut by the circle
	std::vector<DelaunayTetrahedron *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(enrichmentTarget(static_cast<DelaunayTetrahedron *>(disc[i])) )
		{
			ring.push_back(disc[i]) ;
		}
	}
	
	//then we build a list of points to enrich
	std::vector<Point *> points ;
	double rmin = radius ;
	double rmax = radius ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		double d = dist(getCenter(), *ring[i]->first) ;
		if( d > rmax)
			rmax = d ;
		if(d < rmin)
			rmin = d ;
		d = dist(getCenter(), *ring[i]->second) ;
		if( d > rmax)
			rmax = d ;
		if(d < rmin)
			rmin = d ;
		d = dist(getCenter(), *ring[i]->third) ;
		if( d > rmax)
			rmax = d ;
		if(d < rmin)
			rmin = d ;
		d = dist(getCenter(), *ring[i]->fourth) ;
		if( d > rmax)
			rmax = d ;
		if(d < rmin)
			rmin = d ;
		points.push_back(ring[i]->first) ;
		points.push_back(ring[i]->second) ;
		points.push_back(ring[i]->third) ;
		points.push_back(ring[i]->fourth) ;
	}
	
	//we make the points in the list unique
	std::stable_sort(points.begin(), points.end()) ;
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
		std::vector<Point> tetSphereIntersectionPoints = getPrimitive()->intersection(static_cast<Tetrahedron *>(ring[i])) ;

		//for convenience
		Point *a = ring[i]->first ;
		Point *b = ring[i]->second ;
		Point *c = ring[i]->third ;
		Point *d = ring[i]->fourth ;
		std::vector<Point> hint ;
		//if there are no intersection points we need not do anything
		if(!tetSphereIntersectionPoints.empty())
		{
			
			//if the number of intersection points is not 3, we need not do anything
			if(tetSphereIntersectionPoints.size() >= 3)
			{
				Point bary ;
				for(size_t j = 0 ; j < tetSphereIntersectionPoints.size() ; j++)
				{
					hint.push_back(ring[i]->inLocalCoordinates(tetSphereIntersectionPoints[j])) ;
					bary += tetSphereIntersectionPoints[j]/tetSphereIntersectionPoints.size() ;
				}
				hint.push_back(ring[i]->inLocalCoordinates(bary)) ;
			}
		}

		//we build the enrichment function, first, we get the transforms from the triangle
		Function x = ring[i]->getXTransform() ;
		Function y = ring[i]->getYTransform() ;
		Function z = ring[i]->getZTransform() ;
		
		//this function returns the distance to the centre
		Function position(getCenter(), x, y, z) ;
		position -= radius ;
		Function hatNorm = 1.-f_abs(position) ;/*vm.eval(f_abs(position), Point(1,0,0))*shapefunc[0] +
												vm.eval(f_abs(position), Point(0,1,0))*shapefunc[1] +
												vm.eval(f_abs(position), Point(0,0,1))*shapefunc[2] +
												vm.eval(f_abs(position), Point(0,0,0))*shapefunc[3] -
												f_abs( vm.eval(position, Point(1,0,0))*shapefunc[0] +
															vm.eval(position, Point(0,1,0))*shapefunc[1] +
															vm.eval(position, Point(0,0,1))*shapefunc[2] +
															vm.eval(position, Point(0,0,0))*shapefunc[3]);*/
		//finaly, we have the enrichment function
		Function hat = hatNorm ; //1-f_abs((rmax-position)/(rmin-radius) + ((rmax-rmin)/(radius-rmin)-1));

		//enriching the first point
		Function f = shapefunc[0]*(hat-vm.eval(hat, Point(1,0,0))) ;
		f.setIntegrationHint(hint) ;
		f.setPoint(c) ;
		f.setDofID(dofId[c]) ;
		ring[i]->setEnrichment( f, getPrimitive()) ;

		//enriching the second point
		f = shapefunc[1]*(hat - vm.eval(hat, Point(0,1,0))) ;
		f.setIntegrationHint(hint) ;
		f.setPoint(d) ;
		f.setDofID(dofId[d]) ;
		ring[i]->setEnrichment( f, getPrimitive()) ;
		
		//enriching the third point
		f = shapefunc[2]*(hat - vm.eval(hat, Point(0,0,1))) ;
		f.setIntegrationHint(hint) ;
		f.setPoint(a) ;
		f.setDofID(dofId[a]) ;
		ring[i]->setEnrichment(f, getPrimitive()) ;
		
		//enriching the fourth point
		f = shapefunc[3]*(hat - vm.eval(hat, Point(0,0,0))) ;
		f.setIntegrationHint(hint) ;
		f.setPoint(b) ;
		f.setDofID(dofId[b]) ;
		ring[i]->setEnrichment(f, getPrimitive()) ;

		for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
		{
			DelaunayTetrahedron * t = ring[i]->getNeighbourhood(j) ;
			if(!enrichmentTarget(t))
			{
				if(!t->enrichmentUpdated)
					t->clearEnrichment(getPrimitive()) ;
				t->enrichmentUpdated = true ;
				bool hinted = false ;

				//this function returns the distance to the centre
				Function position(getCenter(), t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
				position -= radius ;
				//finaly, we have the enrichment function
				Function hat = 1-f_abs(position)/*vm.eval(f_abs(position), Point(1,0,0))*shapefunc[0] +
												vm.eval(f_abs(position), Point(0,1,0))*shapefunc[1] +
												vm.eval(f_abs(position), Point(0,0,1))*shapefunc[2] +
												vm.eval(f_abs(position), Point(0,0,0))*shapefunc[3] -
												f_abs( vm.eval(position, Point(1,0,0))*shapefunc[0] +
																vm.eval(position, Point(0,1,0))*shapefunc[1] +
																vm.eval(position, Point(0,0,1))*shapefunc[2] +
																vm.eval(position, Point(0,0,0))*shapefunc[3])*/;
				hint.push_back(Point(1./4., 1./4., 1./4.)) ;
				
				if(dofId.find(t->third) != dofId.end())
				{
					Function f = shapefunc[0]*(hat - VirtualMachine().eval(hat, Point(1,0,0))) ;
					if(!hinted)
					{
						f.setIntegrationHint(hint) ;
						hinted = true ;
					}
					f.setPoint(t->third) ;
					f.setDofID(dofId[t->third]) ;
					t->setEnrichment(f, getPrimitive()) ;
				}
				
				if(dofId.find(t->fourth) != dofId.end())
				{
					Function f = shapefunc[1]*(hat - VirtualMachine().eval(hat, Point(0,1,0))) ;
					if(!hinted)
					{
						f.setIntegrationHint(hint) ;
						hinted = true ;
					}
					f.setPoint(t->fourth) ;
					f.setDofID(dofId[t->fourth]) ;
					t->setEnrichment(f, getPrimitive()) ;
				}
				
				if(dofId.find(t->first) != dofId.end())
				{
					Function f = shapefunc[2]*(hat - VirtualMachine().eval(hat, Point(0,0,1))) ;
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
					Function f = shapefunc[3]*(hat - VirtualMachine().eval(hat, Point(0,0,0))) ;
					if(!hinted)
					{
						f.setIntegrationHint(hint) ;
						hinted = true ;
					}
					f.setPoint(t->second) ;
					f.setDofID(dofId[t->second]) ;
					t->setEnrichment(f, getPrimitive()) ;
				}
			}
		}
	}

}
	
bool EnrichmentInclusion3D::interacts(Feature * f, double d) const { return false ;}
	
bool EnrichmentInclusion3D::inBoundary(const Point & v) const {return false ; }
	
std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getTetrahedrons( const DelaunayTree3D * dt) 
{ 
	return dt->conflicts(getPrimitive()) ;
}
	
std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getIntersectingTetrahedrons( DelaunayTree3D * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTetrahedron *> disc = dt->conflicts(getPrimitive()) ;
	
	//then we select those that are cut by the circle
	std::vector<DelaunayTetrahedron *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(this->intersects(static_cast<Tetrahedron *>(disc[i])))
			ring.push_back(disc[i]) ;
	}
	
	return ring ;
}
	
void EnrichmentInclusion3D::setInfluenceRadius(double r) { }
	
std::vector<Geometry *> EnrichmentInclusion3D::getRefinementZones( size_t level) const 
{
	return std::vector<Geometry *>(0) ;
}
	
void EnrichmentInclusion3D::step(double dt, std::valarray<double> *, const DelaunayTree * dtree) {}
	
bool EnrichmentInclusion3D::moved() const { return updated ;}

void EnrichmentInclusion3D::setRadius(double newR)
{
	Sphere::setRadius(newR) ;
	updated = true ;
}
