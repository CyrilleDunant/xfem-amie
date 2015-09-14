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
    bool in0 = in(*t->first)   ;
    bool in1 = in(*t->second)  ;
    bool in2 = in(*t->third)   ;
    
    if(!t->intersects(getPrimitive()))
        return false ;
//     if(t->intersection(getPrimitive()).size() != 2)
//         return false ;
    if((in0 && in1 && !in2) ||
       (in0 && !in1 && in2) || 
       (!in0 && in1 && in2) ||
       (in0 && !in1 && !in2) ||
       (!in0 && !in1 && in2) || 
       (!in0 && in1 && !in2)
    )
        return true ;
        
    return false ;
    
//     Segment s0(*t->first, *t->second) ;
//     Segment s1(*t->second, *t->third) ;
//     Segment s2(*t->third, *t->first) ;
//     return (s0.intersects(getPrimitive()) || s1.intersects(getPrimitive()) || s2.intersects(getPrimitive())) ;
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
     if(cache.empty()) // assume that this is the first step
     {
         DelaunayTriangle * tri = dtree->getUniqueConflictingElement(&getCenter()) ;
         //if(getRadius() < tri->getRadius())
         //{
             std::map<double, Point> pd;
             pd[dist(*tri->first, getCenter())] = *tri->first ;
             pd[dist(*tri->second, getCenter())] = *tri->second ;
             pd[dist(*tri->third, getCenter())] = *tri->third ;
             setCenter(pd.begin()->second);
         //}
     }
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

    if(disc.size() < 4)
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
                if(disc.size() < 4)
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

    //then we select those that are cut by the circle
    std::vector<DelaunayTriangle *> ring ;

    for(size_t i = 0 ; i < disc.size() ; i++)
    {
        if(enrichmentTarget(disc[i]))
            ring.push_back(disc[i]) ;
    }
    
    if(ring.empty())
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
                if(disc.size() < 4)
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
    //then we build a list of points to enrich
    std::set<Point *> points ;
    
    int factor  = ring[0]->getOrder() ;
    for(size_t i = 0 ; i < ring.size() ; i++)
    {
        for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j+=factor)
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

    TriElement father(LINEAR) ;
    for(size_t i = 0 ; i < ring.size() ; i++)
    {
        enrichedElem.insert(ring[i]) ;
        std::vector<Point> hint ;
        std::vector<Point> inter = intersection(ring[i]->getPrimitive()) ;

        if(inter.size() == 2)
        {
            Point linter0 = ring[i]->inLocalCoordinates(inter[0]) ;
            Point linter1 = ring[i]->inLocalCoordinates(inter[1]) ;

            double n  = 3 ;
            for(double j= 1./n ; j < .9999  ; j+=1./n)
            {
                Point h0 = inter[0]*j+inter[1]*(1.-j) ;
                project(&h0);
                hint.push_back(ring[i]->inLocalCoordinates(h0)) ;
            }

            hint.push_back(linter0) ;
            hint.push_back(linter1) ;     

        }

        Function hat ;
        Function hatdx ;
        Function hatdy ;
        
        int factor  = ring[i]->getOrder() ;
        
        if(in(ring[i]->getBoundingPoint(0*factor)) == in(ring[i]->getBoundingPoint(1*factor)))
        {
            hat   = Function(getPrimitive(), ring[i]->getBoundingPoint(2*factor), Segment(ring[i]->getBoundingPoint(0*factor),ring[i]->getBoundingPoint(1*factor)), ring[i]) ;
//            hatdx = Function(getPrimitive(), ring[i]->getBoundingPoint(2*factor), Segment(ring[i]->getBoundingPoint(0*factor),ring[i]->getBoundingPoint(1*factor)), ring[i], XI) ;
//            hatdy = Function(getPrimitive(), ring[i]->getBoundingPoint(2*factor), Segment(ring[i]->getBoundingPoint(0*factor),ring[i]->getBoundingPoint(1*factor)), ring[i], ETA) ;
        }
        else if(in(ring[i]->getBoundingPoint(0*factor)) == in(ring[i]->getBoundingPoint(2*factor)))
        {
            hat   = Function(getPrimitive(), ring[i]->getBoundingPoint(1*factor), Segment(ring[i]->getBoundingPoint(2*factor),ring[i]->getBoundingPoint(0*factor)), ring[i]) ;
//            hatdx = Function(getPrimitive(), ring[i]->getBoundingPoint(1*factor), Segment(ring[i]->getBoundingPoint(2*factor),ring[i]->getBoundingPoint(0*factor)), ring[i], XI) ;
//            hatdy = Function(getPrimitive(), ring[i]->getBoundingPoint(1*factor), Segment(ring[i]->getBoundingPoint(2*factor),ring[i]->getBoundingPoint(0*factor)), ring[i], ETA) ;
        }
        else 
        {
            hat   = Function(getPrimitive(), ring[i]->getBoundingPoint(0*factor), Segment(ring[i]->getBoundingPoint(1*factor),ring[i]->getBoundingPoint(2*factor)), ring[i]) ;
//            hatdx = Function(getPrimitive(), ring[i]->getBoundingPoint(0*factor), Segment(ring[i]->getBoundingPoint(1*factor),ring[i]->getBoundingPoint(2*factor)), ring[i], XI) ;
//            hatdy = Function(getPrimitive(), ring[i]->getBoundingPoint(0*factor), Segment(ring[i]->getBoundingPoint(1*factor),ring[i]->getBoundingPoint(2*factor)), ring[i], ETA) ;
        }
         
        for(size_t j = 0 ; j< ring[i]->getBoundingPoints().size() ; j+= factor)
        {
            bool hinted = false ;
            std::pair<DelaunayTriangle *, Point *> that(ring[i], &ring[i]->getBoundingPoint(j) ) ;
            if(enriched.find(that) == enriched.end())
            {
                enriched.insert(that) ;

                Function f = father.getShapeFunction(j)*hat ;

                if(!hinted)
                {
                    f.setIntegrationHint(hint) ;
                    hinted = true ;
                }
                f.setPoint(&ring[i]->getBoundingPoint(j)) ;
                f.setDofID(dofId[&ring[i]->getBoundingPoint(j)]) ;

//                 f.setNumberOfDerivatives(2);
//                 Function fdx = (father.getShapeFunction(j).d(XI) *hat + father.getShapeFunction(j)*hatdx) ;
//                 Function fdy = (father.getShapeFunction(j).d(ETA)*hat + father.getShapeFunction(j)*hatdy) ;
//                 f.setDerivative(XI , fdx);
//                 f.setDerivative(ETA, fdy);
    
                ring[i]->setEnrichment( f, getPrimitive()) ;
                
                
                
                
//                 for(double j = 0 ; j < 1 ; j+=0.001)
//                 {
//                     for(double k = 0 ; k < 1 ; k+=0.001)
//                     {
//                         if(j+k <= 1 && j >= 0 && k >= 0)
//                             std::cout << VirtualMachine().deval(f,XI,  j,k) << "  " << std::flush ;
//                         else
//                             std::cout << 0 << "  " << std::flush ;
//                     }
//                     std::cout << std::endl ;
//                 }
//                 for(double j = 0 ; j < 1 ; j+=0.001)
//                 {
//                     for(double k = 0 ; k < 1 ; k+=0.001)
//                     {
//                         if(j+k <= 1 && j >= 0 && k >= 0)
//                             std::cout << VirtualMachine().deval(f,ETA,  j,k) << "  " << std::flush ;
//                         else
//                             std::cout << 0 << "  " << std::flush ;
//                     }
//                     std::cout << std::endl ;
//                 }
//                 for(double j = 0 ; j < 1 ; j+=0.05)
//                 {
//                     for(double k =  0 ; k < 1 ; k+=0.05)
//                     {
//                         if(j+k <= 1 && j >= 0 && k >= 0)
//                             std::cout << VirtualMachine().eval(f, j,k) << "  " << std::flush ;
//                         else
//                             std::cout << 0 << "  " << std::flush ;
//                     }
//                     std::cout << std::endl ;
//                 }
                
                
                
            }
        }
//         exit(0) ;
        
    }
//      exit(0) ;

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

