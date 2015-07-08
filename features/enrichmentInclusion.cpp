// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "enrichmentInclusion.h"
#include "inclusion.h"
#include "feature_base.h"
#include "sample.h"
#include "../physics/stiffness.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/void_form.h"
#include "../polynomial/vm_function_extra.h"

using namespace Amie ;

EnrichmentInclusion::EnrichmentInclusion(Feature *father, double radius, double x, double y) : EnrichmentFeature(father), Circle(radius, x, y)
{
    updated = true ;
}

EnrichmentInclusion::EnrichmentInclusion(double radius, double x, double y) :EnrichmentFeature(nullptr), Circle(radius, x, y)
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

Function EnrichmentInclusion::functionTozero(const DelaunayTriangle * t)
{
    TriElement father(LINEAR) ;
    if(in(t->getBoundingPoint(0)) && !in(t->getBoundingPoint(1))  && !in(t->getBoundingPoint(2)))
    {
        return father.getShapeFunction(0)*(Function("1")-father.getShapeFunction(0)) ;
    }

    if(!in(t->getBoundingPoint(0)) && in(t->getBoundingPoint(1))  && !in(t->getBoundingPoint(2)) )
    {
        return father.getShapeFunction(1)*(Function("1")-father.getShapeFunction(1)) ;
    }

    if(!in(t->getBoundingPoint(0)) && !in(t->getBoundingPoint(1))  && in(t->getBoundingPoint(2)) )
    {
        return father.getShapeFunction(2)*(Function("1")-father.getShapeFunction(2)) ;
    }
    return Function("0") ;

}

Function EnrichmentInclusion::getBlendingFunction(const std::map<const Point *, int> & dofIds, const DelaunayTriangle * t)
{
    if(t->getOrder() == QUADRATIC)
        return Function("1") ;
//      return Function("1") ;

    // if(t->getOrder() == QUADRATIC)
    // {
    //  TriElement father(QUADRATIC) ;
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(0) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
    //  }
    //
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(2) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(3);
    //  }
    //
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(4) + 0.25*father.getShapeFunction(3)+ 0.25*father.getShapeFunction(5);
    //  }
    //
    //  if(dofIds.find(t->first) == dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(2)+father.getShapeFunction(3)+father.getShapeFunction(4) + 0.25*father.getShapeFunction(1)+ 0.25*father.getShapeFunction(5);
    //  }
    //
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) == dofIds.end() && dofIds.find(t->third) != dofIds.end())
    //  {
    //      return father.getShapeFunction(0) + father.getShapeFunction(5) + father.getShapeFunction(4) + 0.25*father.getShapeFunction(1) +0.25*father.getShapeFunction(3);
    //  }
    //
    //  if(dofIds.find(t->first) != dofIds.end() && dofIds.find(t->second) != dofIds.end() && dofIds.find(t->third) == dofIds.end())
    //  {
    //      return father.getShapeFunction(1)+father.getShapeFunction(0)+father.getShapeFunction(2) + 0.25*father.getShapeFunction(3) + 0.25*father.getShapeFunction(5);
    //  }
    // }


    TriElement father(LINEAR) ;
    //  Function f ;
    //  for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
    //  {
    //      if(dofIds.find(&(t->getBoundingPoint(i))) != dofIds.end())
    //          f += father.getShapeFunction(i) ;
    //  }
    //  return f ;

    if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
    {
        return father.getShapeFunction(0) ;
    }

    if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
    {
        return father.getShapeFunction(1) ;
    }

    if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
    {
        return father.getShapeFunction(2) ;
    }

    if(dofIds.find(&t->getBoundingPoint(0)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
    {
        return Function("1")-father.getShapeFunction(0) ;
    }

    if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) == dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) != dofIds.end())
    {
        return Function("1")-father.getShapeFunction(1) ;
    }

    if(dofIds.find(&t->getBoundingPoint(0)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(1)) != dofIds.end() && dofIds.find(&t->getBoundingPoint(2)) == dofIds.end())
    {
        return Function("1")-father.getShapeFunction(2) ;
    }

    return Amie::Function("0") ;
} 

void EnrichmentInclusion::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
    freeIds[dtree].clear();
    for(size_t i = 0 ; i < cache.size() ; i++)
    {
        if(!cache[i]->enrichmentUpdated)
        {
            std::vector<size_t> idpts = cache[i]->clearEnrichment(getPrimitive()) ;
            for(size_t j = 0 ; j < idpts.size() ; j++)
                freeIds[dtree].insert(idpts[j]) ;
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
                freeIds[dtree].insert(idpts[j]) ;
        }
        cache[i]->enrichmentUpdated = true ;

    }
    if(cache.empty())
        std::cout << "cache empty !" << std::endl ;
}

void EnrichmentInclusion::enrich(size_t & lastId, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
    freeIds.clear() ;
    if(updated)
    {
        update(dtree) ;
    }
    updated = false ;
    std::vector<DelaunayTriangle *> & disc  = cache;

    if(disc.size() < 6)
    {
        DelaunayTriangle * toHomogenise = disc[0];
        for(size_t i = 0 ; i < disc.size() ; i++)
        {
            if(disc[i]->in(Circle::getCenter()))
            {
                toHomogenise = disc[i] ;
                break ;
            }
        }
        disc.clear() ;
        disc.push_back(toHomogenise) ;
        cache = disc ;
        
        for(size_t i = 0 ; i < disc.size() ; i++)
        {
            HomogeneisedBehaviour * hom = dynamic_cast<HomogeneisedBehaviour *>(disc[i]->getBehaviour());
            if(hom)
            {
                if(disc.size() < 6)
                {
                    std::vector<Feature *> brother ;
                    if(getFather())
                        brother = getFather()->getChildren() ;
                    std::vector<Feature *> feat ;
                    for(size_t j = 0 ; j < brother.size() ; j++)
                    {
                        if(disc[i]->in(brother[j]->getCenter()))
                            feat.push_back(brother[j]) ;
                    }
                    hom->updateEquivalentBehaviour(feat, disc[i]) ;

                }
            
                
            }
            else
            {
                std::vector< Feature *> brother ;
                if(getFather())
                    brother = getFather()->getChildren() ;
                std::vector< Feature *> feat ;
                for(size_t j= 0 ; j < brother.size() ; j++)
                {
                    if(disc[i]->in(brother[j]->getCenter()))
                        feat.push_back(brother[j]) ;
                }
                HomogeneisedBehaviour * hom2 = new HomogeneisedBehaviour(feat, disc[i]) ;
                disc[i]->setBehaviour(dtree,hom2) ;
                disc[i]->getBehaviour()->setSource(getPrimitive()) ;
                
            }
        } 
        return ;
    }
    
    if(disc.size() < 6)
        return ;

    //then we select those that are cut by the circle
    std::vector<DelaunayTriangle *> ring ;

    for(size_t i = 0 ; i < disc.size() ; i++)
    {

        Segment s0(*disc[i]->first, *disc[i]->second) ;
        Segment s1(*disc[i]->second, *disc[i]->third) ;
        Segment s2(*disc[i]->third, *disc[i]->first) ;
        if(!(s0.intersection(getPrimitive()).empty() && s1.intersection(getPrimitive()).empty() && s2.intersection(getPrimitive()).empty())&& disc[i]->getBehaviour()->type != VOID_BEHAVIOUR)
            ring.push_back(disc[i]) ;
    }
    //then we build a list of points to enrich
    std::set<Point *> points ;
    for(size_t i = 0 ; i < ring.size() ; i++)
    {
        for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
        {
            points.insert(&ring[i]->getBoundingPoint(j)) ;
        }
        ring[i]->enrichmentUpdated = false ;
        std::vector<DelaunayTriangle *> neighbourhood = dtree->getNeighbourhood(ring[i]) ;
        for(size_t j = 0 ; j < neighbourhood.size() ; j++)
        {
            neighbourhood[j]->enrichmentUpdated = false ;
        }
    }
    std::sort(ring.begin(), ring.end()) ;
    //we build a map of the points and corresponding enrichment ids
    std::map<const Point *, int> dofId ;

    for(auto i = points.begin() ; i != points.end() ; ++i)
    {
        dofId[*i] = lastId++ ;
    }

    std::set<std::pair<DelaunayTriangle *, Point *> > enriched ;
    std::set<DelaunayTriangle *> enrichedElem;
    //then we iterate on every element

    std::map<Point*, size_t> extradofs ;
    TriElement father(ring[0]->getOrder()) ;
    for(size_t i = 0 ; i < ring.size() ; i++)
    {
        enrichedElem.insert(ring[i]) ;
        std::vector<Point> hint ;
        std::vector<Point> inter = intersection(ring[i]->getPrimitive()) ;

        if(inter.size() == 2)
        {
            for(double j = 0.1 ; j < 0.9 ; j += 0.1)
            {
                Point h0 = inter[0]*j+inter[1]*(1.-j) ;
                project(&h0);
                hint.push_back(ring[i]->inLocalCoordinates(h0)) ;
            }

            hint.push_back(ring[i]->inLocalCoordinates(inter[0])) ;
            hint.push_back(ring[i]->inLocalCoordinates(inter[1])) ;
        }


        //we build the enrichment function, first, we get the transforms from the triangle
        //this function returns the distance to the centre
        Function x = XTransform(ring[i]->getBoundingPoints(), father.getShapeFunctions()) ;
        Function y = YTransform(ring[i]->getBoundingPoints(), father.getShapeFunctions()) ;
        
        Function position = f_sqrt((x-getCenter().getX())*(x-getCenter().getX()) +
                                   (y-getCenter().getY())*(y-getCenter().getY())) ;
        Function hat = getRadius()-f_abs(position-getRadius(), false) ;
        
         if(in(ring[i]->getBoundingPoint(0)) == in(ring[i]->getBoundingPoint(1)))
         {
             hat = Function(getPrimitive(), ring[i]->getBoundingPoint(2), Segment(ring[i]->getBoundingPoint(0),ring[i]->getBoundingPoint(1)), ring[i]) ;
         }
         else if(in(ring[i]->getBoundingPoint(0)) == in(ring[i]->getBoundingPoint(2)))
         {
             hat = Function(getPrimitive(), ring[i]->getBoundingPoint(1), Segment(ring[i]->getBoundingPoint(0),ring[i]->getBoundingPoint(2)), ring[i]) ;
         }
         else if(in(ring[i]->getBoundingPoint(1)) == in(ring[i]->getBoundingPoint(2)))
         {
             hat = Function(getPrimitive(), ring[i]->getBoundingPoint(0), Segment(ring[i]->getBoundingPoint(1),ring[i]->getBoundingPoint(2)), ring[i]) ;
         }
         else
         {
             std::cout << "oops ?" << std::endl ;
             hat = Function("1") ;
         }
         
//         for(double n = 0 ; n < 1 ; n+= 0.01)
//         {
//             for (double m = 0 ;  m < 1 ; m += 0.01)
//             {
//                 if(m+n < 1)
//                     std::cout << VirtualMachine().eval(hat, n, m) << "   "<<std::flush;
//                 else
//                     std::cout << 0 << "   "<<std::flush;
//             }
//             std::cout << std::endl;    
//         }
       

        for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j++)
        {
            std::pair<DelaunayTriangle *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
            if(enriched.find(that) == enriched.end())
            {
                enriched.insert(that) ;
                Point p = ring[i]->inLocalCoordinates(ring[i]->getBoundingPoint(j)) ;

                Function f = father.getShapeFunction(j)*hat /*- VirtualMachine().eval(father.getShapeFunction(j)*hat, p.getX(), p.getY())*//**functionTozero(ring[i])*/;

//                 for(double n = 0 ; n < 1 ; n+= 0.01)
//                 {
//                     for (double m = 0 ;  m < 1 ; m += 0.01)
//                     {
//                         if(m+n < 1)
//                             std::cout << VirtualMachine().eval(f, n, m) << "   "<<std::flush;
//                         else
//                             std::cout << 0 << "   "<<std::flush;
//                     }
//                     std::cout << std::endl;    
//                 }
                
                f.setIntegrationHint(hint) ;
                f.setPoint(&ring[i]->getBoundingPoint(j)) ;
                f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;
                ring[i]->setEnrichment( -f, getPrimitive()) ;
            }
        } 
        
//         if(i == 8)
//             exit (0) ;
// //         exit(0) ;
//         hint.clear();
//         hint.push_back(Point(1./3., 1./3.));
// 
//         std::vector<DelaunayTriangle *> neighbourhood = dtree->getNeighbourhood(ring[i]) ;
//         for(auto & t : neighbourhood)
//         {
//             enrichedElem.insert(t) ;
//             if(std::binary_search(ring.begin(), ring.end(), t) )
//                 continue ;
// 
//             Function blend = getBlendingFunction(dofId, t) ;
// 
//             if(!t->enrichmentUpdated)
//                 t->clearEnrichment( getPrimitive()) ;
// 
//             t->enrichmentUpdated = true ;
//             bool hinted = false ;
//             Function x = XTransform(t->getBoundingPoints(), father.getShapeFunctions()) ;
//             Function y = YTransform(t->getBoundingPoints(), father.getShapeFunctions()) ;
//             
//             Function position = f_sqrt((x-getCenter().getX())*(x-getCenter().getX()) +
//                                        (y-getCenter().getY())*(y-getCenter().getY())) ;
//             Function hat = getRadius()-f_abs(position-getRadius(), false) ;
// 
//             for(size_t k = 0 ; k< t->getBoundingPoints().size() ; k++)
//             {
//                 std::pair<DelaunayTriangle *, Point *> that(t, &t->getBoundingPoint(k) ) ;
//                 if(enriched.find(that) == enriched.end())
//                 {   
//                     enriched.insert(that) ;           
//                     Point p = t->inLocalCoordinates(t->getBoundingPoint(k)) ;
//                     Function f = (father.getShapeFunction(k)*hat- VirtualMachine().eval(father.getShapeFunction(k)*hat, p.getX(), p.getY()))*blend;
//                     
// //                     for(double n = 0 ; n < 1 ; n+= 0.01)
// //                     {
// //                         for (double m = 0 ;  m < 1 ; m += 0.01)
// //                         {
// //                             if(m+n < 1)
// //                                 std::cout << VirtualMachine().eval(f, n, m) << "   "<<std::flush;
// //                             else
// //                                 std::cout << 0 << "   "<<std::flush;
// //                         }
// //                         std::cout << std::endl;    
// //                     }
// //                     exit (0) ;
//                     
//                     if(!hinted)
//                     {
//                         f.setIntegrationHint(hint) ;
//                         hinted = true ;
//                     }
//                     f.setPoint(&t->getBoundingPoint(k)) ;
//                     
//                     if(dofId.find(&t->getBoundingPoint(k)) != dofId.end() )
//                     {
//                         f.setDofID(dofId[&t->getBoundingPoint(k)]) ;
//                         
//                     }
//                     else
//                     {
//                         if(extradofs.find(&t->getBoundingPoint(k)) == extradofs.end())
//                             extradofs[&t->getBoundingPoint(k)] = lastId++ ;
//                         f.setDofID(extradofs[&t->getBoundingPoint(k)]) ; 
//                     }
//                     t->setEnrichment(f, getPrimitive()) ;
//                     
//                 }
//             }
//         }
//     
        
    }

    for(size_t i = 0 ; i < disc.size() ; i++)
    {
        if(enrichedElem.find(disc[i]) == enrichedElem.end())
            disc[i]->clearAllEnrichment() ;
    }
//	updated = true ;
}




bool EnrichmentInclusion::interacts(Feature * f, double d) const {
    return false ;
}
bool EnrichmentInclusion::inBoundary(const Point v) const {
    return false ;
}
bool EnrichmentInclusion::inBoundary(const Point *v) const {
    return false ;
}

std::vector<DelaunayTriangle *> EnrichmentInclusion::getElements2D( FeatureTree * dt)
{
    return dt->get2DMesh()->getConflictingElements(getPrimitive()) ;
}

std::vector<DelaunayTriangle *> EnrichmentInclusion::getIntersectingTriangles( FeatureTree * dt)
{
    //first we get All the triangles affected
    std::vector<DelaunayTriangle *> disc = dt->get2DMesh()->getConflictingElements(getPrimitive()) ;

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

void EnrichmentInclusion::step(double dt, std::valarray< double >*, const Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) {}

bool EnrichmentInclusion::moved() const {
    return updated ;
}

void EnrichmentInclusion::setRadius(double newR)
{
    Circle::setRadius(newR) ;
    updated = true ;
}

