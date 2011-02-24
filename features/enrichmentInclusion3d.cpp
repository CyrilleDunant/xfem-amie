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
	if(pointsin != 0 && pointsin != 4)
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
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() &&
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(1) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(1)+father.getShapeFunction(2)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end()&& dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(0)+father.getShapeFunction(2)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(1)+father.getShapeFunction(0)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(1)+father.getShapeFunction(2)+father.getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(2) + father.getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end()&& dofIds.find(&(t->getBoundingPoint(3))) != dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(1))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(3))) == dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(1) ;
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
			return father.getShapeFunction(0) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(1) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(1)+father.getShapeFunction(2)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return  father.getShapeFunction(0)+father.getShapeFunction(2)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(0)+father.getShapeFunction(1)+father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0)+father.getShapeFunction(1)+father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(2) + father.getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(3);
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(1) + father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) != dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(3) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(2) ;
		}
		
		if(dofIds.find(&(t->getBoundingPoint(0))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(2))) != dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(4))) == dofIds.end() && 
			dofIds.find(&(t->getBoundingPoint(6))) == dofIds.end())
		{
			return father.getShapeFunction(0) + father.getShapeFunction(1) ;
		}
	}
	
	
	return Function("1") ;
}

void EnrichmentInclusion3D::enrich(size_t & lastId,  Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree)
{
	if(updated)
		update(dtree) ;
	updated = false ;

	VirtualMachine vm ;
	const std::vector<DelaunayTetrahedron *> & disc  = cache;

	std::valarray<Function> shapefunc = TetrahedralElement(LINEAR).getShapeFunctions() ;

	//then we select those that are cut by the circle
	std::vector<DelaunayTetrahedron *> ring ;
	
	for(size_t i = 0 ; i < disc.size() ; i++)
	{
		if(enrichmentTarget(static_cast<DelaunayTetrahedron *>(disc[i])) )
		{
			ring.push_back(disc[i]) ;
		}
	}
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
// 		std::vector<Point> tetSphereIntersectionPoints = getPrimitive()->intersection(ring[i]->getPrimitive()) ;

		std::vector<Point> hint ;
// 		if there are no intersection points we need not do anything
// 		if(!tetSphereIntersectionPoints.empty())
// 		{
// 			
// 			//if the number of intersection points is not 3, we need not do anything
// 			if(tetSphereIntersectionPoints.size() >= 3)
// 			{
// // 				Point bary ;
// 				for(size_t j = 0 ; j < tetSphereIntersectionPoints.size() ; j++)
// 				{
// 					hint.push_back(ring[i]->inLocalCoordinates(tetSphereIntersectionPoints[j])) ;
// // 					bary += tetSphereIntersectionPoints[j]/tetSphereIntersectionPoints.size() ;
// 				}
// // 				project(&bary);
// // 				hint.push_back(ring[i]->inLocalCoordinates(bary)) ;
// 			}
// 		}
// 		hint.clear();
		//we build the enrichment function, first, we get the transforms from the triangle
		Function x = ring[i]->getXTransform() ;
		Function y = ring[i]->getYTransform() ;
		Function z = ring[i]->getZTransform() ;
		
		//this function returns the distance to the centre
		Function position(getCenter(), x, y, z) ;
		Function hat = f_abs(position-getRadius());
		
		for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
		{
			std::pair<DelaunayTetrahedron *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
			if(enriched.find(that) == enriched.end())
			{
				enriched.insert(that) ;
				Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;
				Function f =  ring[i]->getShapeFunction(j)*(hat - VirtualMachine().eval(hat, p.x, p.y, p.z)) ;
				f.setIntegrationHint(hint) ;
				f.setPoint(&ring[i]->getBoundingPoint(j)) ;
				f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
				ring[i]->setEnrichment( f, getPrimitive()) ;
			}
		}
		hint.clear();
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
			Function position(getCenter(), t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
			Function hat = f_abs(position-getRadius()) ;

			
			for(size_t k = 0 ; k< t->getBoundingPoints().size() ; k++)
			{
				std::pair<DelaunayTetrahedron *, Point *> that(t, &t->getBoundingPoint(k) ) ;
				
				if(dofId.find(&t->getBoundingPoint(k)) != dofId.end() && enriched.find(that) == enriched.end())
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
				else if(dofId.find(&t->getBoundingPoint(k)) == dofId.end() && enriched.find(that) == enriched.end())
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
	
void EnrichmentInclusion3D::step(double dt, std::valarray<double> *, const DelaunayTree * dtree) {}
	
bool EnrichmentInclusion3D::moved() const { return updated ;}

void EnrichmentInclusion3D::setRadius(double newR)
{
	Sphere::setRadius(newR) ;
	updated = true ;
}
