// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion3d.h"

using namespace Mu ;

EnrichmentInclusion3D::EnrichmentInclusion3D(Feature *father, double radius, double x, double y, double z) : EnrichmentFeature(father), Sphere(radius, x, y, z)
{
	updated = true ;
}

EnrichmentInclusion3D::EnrichmentInclusion3D(double radius, double x, double y, double z) :EnrichmentFeature(nullptr), Sphere(radius, x, y, z)
{
	updated = true ;
}

EnrichmentInclusion3D::~EnrichmentInclusion3D() {}
	
bool EnrichmentInclusion3D::enrichmentTarget(DelaunayTetrahedron * t)
{
	int pointsin = in(*t->first) + in(*t->second) + in(*t->third) + in(*t->fourth) ;
	if(pointsin && pointsin != 4 || intersects(t->getPrimitive()) || t->in(getCenter()))
		return true ;
		
	return false ;
}

void EnrichmentInclusion3D::update(Mesh<DelaunayTetrahedron,DelaunayTreeItem3D> * dtree)
{
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		cache[i]->enrichmentUpdated = true ;
	}
	cache = dtree->getConflictingElements(getPrimitive()) ;
	if(cache.empty())
	{
		std::vector<DelaunayTetrahedron *> candidates = dtree->getConflictingElements(&getCenter()) ;
		for(size_t i = 0 ; i < candidates.size() ; i++)
		{
			if(candidates[i]->isTetrahedron() && candidates[i]->in(getCenter()))
			{
				cache.push_back(static_cast<DelaunayTetrahedron *>(candidates[i])) ;
				break ;
			}
		}
	}
	for(size_t i = 0 ; i < cache.size() ; i++)
	{
		cache[i]->enrichmentUpdated = true ;

	}
	
	if(cache.empty())
		std::cout << "cache empty !" << std::endl ;
}

Function getBlendingFunction(const std::map<const Point *, int> & dofIds, const DelaunayTetrahedron * t)
{
// 	return Function("1") ;
	
	if(t->getOrder() == LINEAR)
	{
		TetrahedralElement father(LINEAR) ;
		
// 		std::cout << std::endl ;
// 	for(size_t i = 0 ; i < 4 ; i++)
// 		{
// 			for(size_t j = 0 ; j < 4 ; j++)
// 			{
// 				std::cout << i << "  " << j << "  " << VirtualMachine().eval(father.getShapeFunction(i), t->inLocalCoordinates(t->getBoundingPoint(j))) << std::endl ;
// 			}
// 		}
// 		std::cout << std::endl ;
		/*exit(0) */;
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() &&
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(1) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(1)+t->getShapeFunction(2)+t->getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(0)+t->getShapeFunction(2)+t->getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(1)+t->getShapeFunction(0)+t->getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(1)+t->getShapeFunction(2)+t->getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(2) + t->getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(1) + t->getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(1) + t->getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return t->getShapeFunction(0) + t->getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(0) + t->getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return t->getShapeFunction(0) + t->getShapeFunction(1) ;
		}
	}
	else
	{
		TetrahedralElement father(LINEAR) ;
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(1) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			
			return father.getShapeFunction(1)+father.getShapeFunction(0)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return  father.getShapeFunction(0)+father.getShapeFunction(2)+father.getShapeFunction(1) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(2)+father.getShapeFunction(1)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0)+father.getShapeFunction(3)+father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(1);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(2) + father.getShapeFunction(0);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(2) + father.getShapeFunction(3) ;
		}
	}
	
	
	return Function("1") ;
}

void EnrichmentInclusion3D::enrich(size_t & lastId,  Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree)
{
// 	lastId = 40000;
	if(updated)
		update(dtree) ;
	updated = false ;

	VirtualMachine vm ;
	const std::vector<DelaunayTetrahedron *> & disc  = cache;

	//then we select those that are cut by the circle
	std::vector<DelaunayTetrahedron *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(enrichmentTarget(disc[i]) )
		{
			ring.push_back(disc[i]) ;
		}
	}
	
	if(ring.size() < 2)
		return ;
	
	Order order = LINEAR ;
	if(!ring.empty())
		order = ring[0]->getOrder() ;
	TetrahedralElement father(order) ;
	//sorting the element for later usage of std::find
	std::sort(ring.begin(), ring.end()) ;
	
	//then we build a list of points to enrich
	std::set<Point *> points ;
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			points.insert(&ring[i]->getBoundingPoint(j)) ;
		}
	}
	
	//we build a map of the points and corresponding enrichment ids
	std::map<const Point *, int> dofId ;
	std::map<Point*, size_t> extradofs ;
	for(auto i = points.begin() ; i != points.end() ; ++i)
	{
			dofId[*i] = lastId++ ;
	}
	
	//now, we will start the enrichment itself
	std::set<std::pair<DelaunayTetrahedron *, Point *> > enriched ;
	
	//then we iterate on every element
	
	for(size_t i = 0 ; i < ring.size() ; i++)
	{
		std::vector<Point> tetSphereIntersectionPoints = getPrimitive()->intersection(ring[i]->getPrimitive()) ;

		std::vector<Point> hint ;
// 		if there are no intersection points we need not do anything
		size_t isize = tetSphereIntersectionPoints.size() ;
		
// 		for(size_t j = 0 ; j < isize ; j++)
// 		{
// 			for(size_t k = j+1 ; k < isize ; k++)
// 			{
// 				tetSphereIntersectionPoints.push_back((tetSphereIntersectionPoints[j]*.6666666666666666666667+tetSphereIntersectionPoints[k])*.3333333333333333333);
// 			}
// 		}
		
		for(size_t j = 0 ; j < tetSphereIntersectionPoints.size() ; j++)
		{
			Point localintersect = ring[i]->inLocalCoordinates(tetSphereIntersectionPoints[j]) ;
			bool add = true ;
			for(size_t k = 0 ; k < father.getBoundingPoints().size() ; k++)
			{
				if(squareDist3D(father.getBoundingPoint(k), localintersect) < .0001 || !ring[i]->in(tetSphereIntersectionPoints[j]))
				{
					add = false ;
					break ;
				}
			}
			
			for(size_t k = 0 ; k < hint.size() ; k++)
			{
				if(squareDist3D(hint[k], localintersect) < .0001)
				{
					add = false ;
					break ;
				}
			}
			
			if(add)
			{
				hint.push_back(localintersect) ;
			}
		}

// 		if(ring[i]->getOrder() == QUADRATIC)
// 		{
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(1))) ;
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(3))) ;
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(5))) ;
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(7))) ;
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(8))) ;
// 			hint.push_back(ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(9))) ;
// 		}
		
		
		//this function returns the distance to the centre
		Function position(getCenter(),ring[i]) ;
		Function hat =  getRadius()-f_abs(position-getRadius());
		
		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			std::pair<DelaunayTetrahedron *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
			if(enriched.find(that) == enriched.end())
			{
				enriched.insert(that) ;
				Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;
				Function f =  ring[i]->getShapeFunction(j)*(hat - vm.eval(hat, p.x, p.y, p.z)) ;
				f.setIntegrationHint(hint) ;
				f.setPoint(&ring[i]->getBoundingPoint(j)) ;
				f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
				ring[i]->setEnrichment( f, getPrimitive()) ;
			}
		}
		hint.clear();
		hint.push_back(Point(.25, .25, .25));
		for(size_t j = 0 ; j < ring[i]->neighbourhood.size() ; j++)
		{
			DelaunayTetrahedron * t = ring[i]->getNeighbourhood(j) ;
			if(std::binary_search(ring.begin(), ring.end(), t))
				continue ;
			
			Function blend = getBlendingFunction(dofId, t) ;
			
			if(!t->enrichmentUpdated)
				t->clearEnrichment( getPrimitive()) ;
			
			t->enrichmentUpdated = true ;
			bool hinted = false ;
			Function position(getCenter(), t) ;
			Function hat = getRadius()-f_abs(position-getRadius()) ;

			
			for(size_t k = 0 ; k< t->getBoundingPoints().size() ; k++)
			{
				std::pair<DelaunayTetrahedron *, Point *> that(t, &t->getBoundingPoint(k) ) ;
				
				if( enriched.find(that) == enriched.end())
				{
					if(dofId.find(&t->getBoundingPoint(k)) != dofId.end() )
					{
						enriched.insert(that) ;
						Point p = t->inLocalCoordinates(t->getBoundingPoint(k)) ;
						Function f = t->getShapeFunction(k)*(hat - VirtualMachine().eval(hat, p.x, p.y, p.z)) ;
						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(&t->getBoundingPoint(k)) ;
						f.setDofID(dofId[&t->getBoundingPoint(k)]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
					else 
					{
						enriched.insert(that) ;
						Point p = t->inLocalCoordinates(t->getBoundingPoint(k)) ;
						Function f = t->getShapeFunction(k)*(hat*blend - VirtualMachine().eval(hat*blend, p.x, p.y, p.z)) ;

						if(!hinted)
						{
							f.setIntegrationHint(hint) ;
							hinted = true ;
						}
						f.setPoint(&t->getBoundingPoint(k)) ;
						if(extradofs.find(&t->getBoundingPoint(k)) == extradofs.end())
							extradofs[&t->getBoundingPoint(k)] = lastId++ ;
						
						f.setDofID(extradofs[&t->getBoundingPoint(k)]) ;
						t->setEnrichment(f, getPrimitive()) ;
					}
				}
			}
		}
	}

}
	
bool EnrichmentInclusion3D::interacts(Feature * f, double d) const { return false ;}
	
bool EnrichmentInclusion3D::inBoundary(const Point & v) const {return false ; }
	
std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getElements3D( FeatureTree * dt) 
{ 
	return dt->getElements3D(getPrimitive()) ;
}
	
std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getBoundingElements3D( FeatureTree * dt)
{
	//first we get All the triangles affected
	std::vector<DelaunayTetrahedron *> disc = dt->getElements3D(getPrimitive()) ;
	
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

void EnrichmentInclusion3D::step(double dt, std::valarray<double> *, const Mesh<Mu::DelaunayTriangle, Mu::DelaunayTreeItem> * dtree) {}
	
bool EnrichmentInclusion3D::moved() const { return updated ;}

void EnrichmentInclusion3D::setRadius(double newR)
{
	Sphere::setRadius(newR) ;
	updated = true ;
}
