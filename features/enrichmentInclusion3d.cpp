// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion3d.h"
#include "../polynomial/vm_function_extra.h"
// #include "../polynomial/vm_function_base.h"



using namespace Amie ;



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
    
        bool in0 = in(t->getBoundingPoint(0)) ;
        bool in1 = in(t->getBoundingPoint(1)) ;
        bool in2 = in(t->getBoundingPoint(2)) ; 
        bool in3 = in(t->getBoundingPoint(3)) ;
        
         if( in0 == in1 && in0 == in2  && in0 != in3 )
         {
            return true ;
         }
         else if ( in0 == in1 && in0 == in3 && in0 != in2)
         {
            return true ;
         }
         else if ( in0 == in2 && in0 == in3 && in0 != in1)
         {
            return true ;
         }
         else if ( in1 == in2 && in1 == in3 && in0 != in1)
         {
            return true ;
         }

        

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
            if(candidates[i]->isTetrahedron && candidates[i]->in(getCenter()))
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
    TetrahedralElement father(LINEAR) ;

    size_t idxFirst = 0 ;
    size_t idxSecond = 0 ;
    size_t idxThird = 0 ;
    size_t idxFourth = 0 ;

    while(idxFirst < t->getBoundingPoints().size()&&t->inLocalCoordinates(t->getBoundingPoint(idxFirst)) != father.getBoundingPoint(0)  )
    {
        idxFirst++ ;
    }
    while(idxSecond < t->getBoundingPoints().size()&&t->inLocalCoordinates(t->getBoundingPoint(idxSecond)) != father.getBoundingPoint(1) )
    {
        idxSecond++ ;
    }
    while(idxThird < t->getBoundingPoints().size()&&t->inLocalCoordinates(t->getBoundingPoint(idxThird)) != father.getBoundingPoint(2) )
    {
        idxThird++ ;
    }
    while(idxFourth < t->getBoundingPoints().size()&&t->inLocalCoordinates(t->getBoundingPoint(idxFourth)) != father.getBoundingPoint(3) )
    {
        idxFourth++ ;
    }

    size_t a = 0 ;
    size_t b = 1 ;
    size_t c = 2 ;
    size_t d = 3 ;

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(a) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(b) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(c) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(d) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(b)+father.getShapeFunction(c)+father.getShapeFunction(d) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(a)+father.getShapeFunction(c)+father.getShapeFunction(d) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(b)+father.getShapeFunction(a)+father.getShapeFunction(d) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(b)+father.getShapeFunction(c)+father.getShapeFunction(a) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(c) + father.getShapeFunction(d);
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(b) + father.getShapeFunction(d);
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(b) + father.getShapeFunction(c) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) != dofIds.end())
    {
        return father.getShapeFunction(a) + father.getShapeFunction(d) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(a) + father.getShapeFunction(c) ;
    }

    if(dofIds.find(&t->getBoundingPoint(idxFirst)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxSecond)) != dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxThird)) == dofIds.end() &&
            dofIds.find(&t->getBoundingPoint(idxFourth)) == dofIds.end())
    {
        return father.getShapeFunction(a) + father.getShapeFunction(b) ;
    }

    Function one("1") ;
    Function zero("0") ;
    one.setNumberOfDerivatives(3);
    one.setDerivative(XI, zero);
    one.setDerivative(ETA, zero);
    one.setDerivative(ZETA, zero);
    return one ;
}

void EnrichmentInclusion3D::enrich(size_t & lastId,  Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree)
{
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

    TetrahedralElement father(LINEAR) ;
    //sorting the element for later usage of std::find
    std::sort(ring.begin(), ring.end()) ;

    //then we build a list of points to enrich
    std::set<Point *> points ;
    double maxd = 0 ;
    for(size_t i = 0 ; i < ring.size() ; i++)
    {
        maxd = std::max(maxd, ring[i]->getRadius()) ;
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


        Function hat ;
        bool in0 = in(ring[i]->getBoundingPoint(0)) ;
        bool in1 = in(ring[i]->getBoundingPoint(1)) ;
        bool in2 = in(ring[i]->getBoundingPoint(2)) ; 
        bool in3 = in(ring[i]->getBoundingPoint(3)) ;
        
         if( in0 == in1 && in0 == in2 )
         {
             hat = Function(getPrimitive(), 
                            ring[i]->getBoundingPoint(3), 
                            TriPoint(ring[i]->getBoundingPoint(0), 
                                     ring[i]->getBoundingPoint(1), 
                                     ring[i]->getBoundingPoint(2)), 
                            ring[i]) ;
         }
         else if ( in0 == in1 && in0 == in3 )
         {
             hat = Function(getPrimitive(), 
                            ring[i]->getBoundingPoint(2), 
                            TriPoint(ring[i]->getBoundingPoint(0), 
                                     ring[i]->getBoundingPoint(1), 
                                     ring[i]->getBoundingPoint(3)), 
                            ring[i]) ;
         }
         else if ( in0 == in2 && in0 == in3 )
         {
             hat = Function(getPrimitive(), 
                            ring[i]->getBoundingPoint(1),  
                            TriPoint(ring[i]->getBoundingPoint(0), 
                                     ring[i]->getBoundingPoint(2), 
                                     ring[i]->getBoundingPoint(3)), 
                            ring[i]) ;
         }
         else if ( in1 == in2 && in1 == in3 )
         {
             hat = Function(getPrimitive(), 
                            ring[i]->getBoundingPoint(0),  
                            TriPoint(ring[i]->getBoundingPoint(1), 
                                     ring[i]->getBoundingPoint(2), 
                                     ring[i]->getBoundingPoint(3)), 
                            ring[i]) ;
         }
         else
         {
             std::cout << "oops ?" << std::endl ;
             hat = Function("1") ;
         }

        for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
        {
            std::pair<DelaunayTetrahedron *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
            if(enriched.find(that) == enriched.end())
            {
                enriched.insert(that) ;
                Function f =  father.getShapeFunction(j)*hat  ;

//                 for(double l = 0 ; l < 1 ; l += .01)
//                 {
//                     for(double k = 0 ; k < 1 ; k += .01)
//                     {
//                             if(l+k +.1< 1)
//                                 std::cout << vm.eval(f, l, k, .1)<< "  " << std::flush ;
//                             else
//                                 std::cout << 0<< "  " << std::flush ;
//                     }
//                     std::cout << std::endl ;
//                 }

                f.setIntegrationHint(hint) ;
                f.setPoint(&ring[i]->getBoundingPoint(j)) ;
                f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
                ring[i]->setEnrichment( f, getPrimitive()) ;
            }
        }
//         exit(0) ;
    }

}

bool EnrichmentInclusion3D::interacts(Feature * f, double d) const {
    return false ;
}

bool EnrichmentInclusion3D::inBoundary(const Point & v) const {
    return false ;
}

std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getElements3D( FeatureTree * dt)
{
    return dt->get3DMesh()->getConflictingElements(getPrimitive()) ;
}

std::vector<DelaunayTetrahedron *> EnrichmentInclusion3D::getBoundingElements3D( FeatureTree * dt)
{
    //first we get All the triangles affected
    std::vector<DelaunayTetrahedron *> disc = dt->get3DMesh()->getConflictingElements(getPrimitive()) ;

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

void EnrichmentInclusion3D::step(double dt, std::valarray<double> *, const Mesh<Amie::DelaunayTriangle, Amie::DelaunayTreeItem> * dtree) {}

bool EnrichmentInclusion3D::moved() const {
    return updated ;
}

void EnrichmentInclusion3D::setRadius(double newR)
{
    Sphere::setRadius(newR) ;
    updated = true ;
}
