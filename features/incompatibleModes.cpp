// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "incompatibleModes.h"
#include "inclusion.h"
#include "feature_base.h"
#include "sample.h"

namespace Amie {

IncompatibleModes::IncompatibleModes(Feature *father) : EnrichmentFeature(father)
{
    updated = true ;
}

IncompatibleModes::IncompatibleModes() :EnrichmentFeature(nullptr)
{
    updated = true ;
}

IncompatibleModes::~IncompatibleModes() {}

bool IncompatibleModes::enrichmentTarget(DelaunayTriangle * t)
{
  return true ;
}



void IncompatibleModes::update(Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
     if(cache.empty()) // assume that this is the first step
     {
        for(auto e = dtree->begin() ; e != dtree->end() ; e++)
        {
            if(e->getBehaviour() && e->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;
            e->enrichmentUpdated = true ;
            cache.push_back(e) ;
        }
     }
     else
     {
         return ;
     }

}

void IncompatibleModes::enrich(size_t & lastId, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
    freeIds.clear() ;
    if(updated)
    {
        update(dtree) ;
    }
    updated = false ;

    //then we build a list of points to enrich
    std::set<Point *> points ;
    
    int firstID = lastId ;
    for(size_t i = 0 ; i < cache.size() ; i++)
    {
        cache[i]->enrichmentUpdated = false ;
    }

    //then we iterate on every element

    TriElement father(LINEAR) ;
    Function enrfunc0("1 9 / x x * -") ;
    Function enrfunc1("1 9 / y y * -") ;
    
    for(size_t i = 0 ; i < cache.size() ; i++)
    {

        Function f = enrfunc0 ;

        f.setPoint(nullptr) ;
        f.setDofID(firstID+i*2) ;
        lastId++ ;

        cache[i]->setEnrichment( f, nullptr) ;
                    
        f = enrfunc1 ;

        f.setPoint(nullptr) ;
        f.setDofID(firstID+i*2+1) ;
        lastId++ ;

        cache[i]->setEnrichment( f, nullptr) ;

        
        
    }

//	updated = true ;
}




bool IncompatibleModes::interacts(Feature * f, double d) const {
    return false ;
}
bool IncompatibleModes::inBoundary(const Point & v, double d) const {
    return false ;
}
bool IncompatibleModes::inBoundary(const Point *v, double d) const {
    return false ;
}

std::vector<DelaunayTriangle *> IncompatibleModes::getElements2D( FeatureTree * dt)
{
    return cache ;
}

std::vector<DelaunayTriangle *> IncompatibleModes::getIntersectingTriangles( FeatureTree * dt)
{
    return cache ;
}

void IncompatibleModes::setInfluenceRadius(double r) { }

std::vector<Geometry *> IncompatibleModes::getRefinementZones( size_t level) const
{
    return std::vector<Geometry *>(0) ;
}

void IncompatibleModes::step(double dt, std::valarray< double >*, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) {}

bool IncompatibleModes::moved() const {
    return updated ;
}


}
