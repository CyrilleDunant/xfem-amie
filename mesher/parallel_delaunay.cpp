//
// C++ Implementation: delaunay
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parallel_delaunay.h"
#include <omp.h>
#include <limits>
#include <algorithm>

// #define DEBUG
// #undef DEBUG

using namespace Amie ;

bool ParallelDelaunayTree::isSame(const DelaunayTreeItem * i0, const DelaunayTreeItem * i1) const 
{
    if(i0->isTriangle && !i1->isTriangle)
        return false ;
    if(i0->isPlane && !i1->isPlane)
        return false ;
    return i0->isVertex( i1->first ) && i0->isVertex( i1->second ) && i0->isVertex( i1->third )  ;
}

bool ParallelDelaunayTree::inDomain(int domain_index, const DelaunayTriangle * tet) const 
{
    if(domains[domain_index]->in(tet->getCenter()))
    {
        bool unique  = true ;
        for(size_t k = domain_index+1 ; k < domains.size() ; k++)
        {
            if(domains[k]->in(tet->getCenter()))
            {
                unique = false ;
                break ;
            }
        }
        if(unique)
            return true ;
    }
    return false ;
}

int ParallelDelaunayTree::getMesh(const DelaunayTreeItem * self) const 
{
    for(size_t i = 0 ; i < meshes.size() ; i++)
        if(self->tree == meshes[i])
            return i ;
    return -1 ;
}

int ParallelDelaunayTree::getDomain(const DelaunayTriangle * tet) const 
{
    int mesh = getMesh(tet) ;
    std::valarray<bool> inDomain(false, domains.size()) ;
    #pragma omp parallel for
    for(size_t domain_index = 0 ; domain_index < domains.size() ;  domain_index++)
    {
        inDomain[domain_index] = domains[domain_index]->in(tet->getCenter()) ;
    }
    
    for(size_t domain_index = 0 ; domain_index < domains.size() ;  domain_index++)
    {
        if(inDomain[domain_index])
        {
            bool unique  = true ;
            for(size_t k = domain_index+1 ; k < domains.size() ; k++)
            {
                if(inDomain[k])
                {
                    unique = false ;
                    break ;
                }
            }
            if( unique && mesh == domain_index )
                return domain_index ;
        }
    }
    return -1 ;
}

int ParallelDelaunayTree::getDomain(const Point & center) const 
{
    for(size_t domain_index = 0 ; domain_index < domains.size() ;  domain_index++)
    {
        if(domains[domain_index]->in(center))
        {
            bool unique  = true ;
            for(size_t k = domain_index+1 ; k < domains.size() ; k++)
            {
                if(domains[k]->in(center))
                {
                    unique = false ;
                    break ;
                }
            }
            if(unique)
                return domain_index ;
        }
    }
    return -1 ;
}

ParallelDelaunayTree::ParallelDelaunayTree(Point * p0,  Point *p1,  Point *p2, const std::vector<Geometry *> & domains) : domains(domains)
{
    std::vector<std::vector<DelaunayTreeItem *> > newElements(domains.size()) ;

    for(size_t i = 0 ; i < domains.size() ; i++)
    {
        meshes.push_back(new DelaunayTree(p0, p1, p2));

    }

    global_counter = meshes[0]->global_counter ;
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getConflictingElements(const Point  * p)
{
    std::valarray<std::vector<DelaunayTriangle *> > conflicts(meshes.size()) ;
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTriangle *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ; j++)
        {
            if(getDomain(tmpConflicts[i]) != 1)
                conflicts[i].push_back(tmpConflicts[j]) ;
        }
    }
    std::vector<DelaunayTriangle *> ret = conflicts[0] ;
    for(size_t i = 1 ; i < conflicts.size() ; i++)
        ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());

    return ret ;
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getConflictingElements(const Geometry  * p)
{
    std::valarray<std::vector<DelaunayTriangle *> > conflicts(meshes.size()) ;
    
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTriangle *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ;  j++)
        {
            if(inDomain(i, tmpConflicts[j]))
                conflicts[i].push_back(tmpConflicts[j]);
        }
    }

    std::vector<DelaunayTriangle *> ret = conflicts[0] ;
    for(size_t i = 1 ; i < conflicts.size() ; i++)
        ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());

    return ret ;
}

bool interacts(const Geometry * g , DelaunayTreeItem * item)
{

//     checkfather = false ;
    if( item->isPlane )
    {
        return false ;
    }
    else if( item->isTriangle)
    {
        if(g->in(*item->first)  || 
           g->in(*item->second) || 
           g->in(*item->third)   
          )
        {
            return true ;
        }
            
        for(size_t j = 0 ; j < item->neighbour.size() ;  j++)
        {
            if(item->getNeighbour(j)->isPlane)
                continue ;
            
            DelaunayTreeItem * n = item->getNeighbour(j) ;
            
            if(g->in(*(n->first))  || 
               g->in(*(n->second)) || 
               g->in(*(n->third))   )
            {
                return true ;
            }
        }

    }    

    return false ;
}

void ParallelDelaunayTree::insert(Point * p)
{   
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {   
//         #pragma omp task
//         {
            bool isVertex = false ;
            bool noInteraction = true ;
            std::vector<DelaunayTreeItem *> cons = meshes[i]->conflicts(p) ;
            std::vector<DelaunayTreeItem *> newElems ;
            if(cons.empty())
            {
                std::cout << "Failed insertion : in nothing !" << std::endl ;
            }
            
            for(size_t j = 0 ; j < cons.size() ; j++)
            {
                if(cons[j]->isVertex(p))
                {
                    isVertex = true ;
                    break ;
                }
            }
            if(!isVertex)
            {
            
                if(domains[i]->in(*p))
                    noInteraction = false ;
                else
                {
                    for(const auto & c : cons )
                    {
                        if(interacts(domains[i], c))
                        {
                            noInteraction = false ;
                            break ;
                        }
                    }
                }
                
                if(!noInteraction || global_counter < 5 )
                {
                    meshes[i]->neighbourhood = false ;
                    newElems = meshes[i]->addElements(cons, p) ;
                }
                
                unsigned int maxIdx = 0 ;
                for(size_t j = 0 ; j < newElems.size() ; j++)
                    maxIdx = std::max(maxIdx, newElems[j]->index) ;


                if(i == meshes.size()-1)
                    p->getId()= global_counter++ ;
            }
//         }
    }
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getElements()
{
    std::vector<std::vector<DelaunayTriangle *>> tris(meshes.size()) ;

    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTriangle *> tmp = meshes[i]->getElements() ;
        for(const auto & t : tmp)
        {

            if(getDomain(t) == i)
                tris[i].push_back(t);
        }
    }

    std::vector<DelaunayTriangle *> ret ;
    for(size_t i = 0 ; i < tris.size() ;  i++)
        ret.insert(ret.end(), tris[i].begin(), tris[i].end());

    return ret ;
}

std::vector<DelaunayTriangle *> ParallelDelaunayTree::getNeighbourhood(DelaunayTriangle * element) const
{
    std::vector<DelaunayTriangle *> ret ;
    int domain = getDomain(element) ;
    if(domain == -1)
    {
        for(const auto & idx : element->neighbourhood)
        {
            ret.push_back((DelaunayTriangle *)meshes[domain]->tree[idx]);
        }
        return ret ;
    }
    else
    {
        for(const auto & idx : element->neighbourhood)
        {
            DelaunayTriangle * n = (DelaunayTriangle *)meshes[domain]->tree[idx] ;
            if(getDomain(n) == domain)
                ret.push_back(n);
            else
            {
                int altDomain = getDomain(n->getCenter()) ;
                DelaunayTriangle * an = meshes[altDomain]->getUniqueConflictingElement(&n->getCenter()) ;
                ret.push_back(an);
            }
        }
        return ret ;
    }
};

void ParallelDelaunayTree::addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep)
{   
    std::vector<DelaunayTriangle *> tri = getElements() ;
    std::valarray<bool> visited(false, size()) ;
    double timeSlice = timestep ;
    
    for(auto & i :tri)
    {
        visited[i->index] = true ;
            
        size_t nodes_per_plane = nodes_per_side*3+3 ;
        
        std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
        std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
        
        for(size_t plane = 0 ; plane < time_planes ; plane++)
        {
            for(size_t side = 0 ; side < 3 ; side++)
            {
                Point a(i->getBoundingPoint(side)) ;
                Point b(i->getBoundingPoint((side+1)%3)) ;
                
                if(time_planes> 1)
                {
                    a.getT() = -timestep/2 + ((double) plane / (double) (time_planes-1))*timeSlice ;
                    b.getT() = -timestep/2 + ((double) plane / (double) (time_planes-1))*timeSlice ;
                }
                for(size_t node = 0 ; node < nodes_per_side+1 ; node++)
                {
                    double fraction = (double)(node)/((double)nodes_per_side+1) ;
                    Point proto = a*(1.-fraction) + b*fraction ;
                    Point * foundPoint = nullptr ;
                    
                    for(size_t j = 0 ; j< i->getBoundingPoints().size() ; j++)
                    {
                        if(i->getBoundingPoint(j) == proto)
                        {
                            foundPoint = &i->getBoundingPoint(j) ;
                            break ;
                        }
                    }
                    
                    if(!foundPoint)
                    {
                        std::vector<DelaunayTriangle *> neighbourhood = getNeighbourhood(i) ;
                        for(size_t j = 0 ; j < neighbourhood.size() ; j++)
                        {
                            if(visited[neighbourhood[j]->index])
                            {
                                DelaunayTriangle * n = neighbourhood[j] ;
                                for(size_t k = 0 ; k < n->getBoundingPoints().size() ; k++)
                                {
                                    if(n->getBoundingPoint(k) == proto)
                                    {
                                        foundPoint = &n->getBoundingPoint(k) ;
                                        break ;
                                    }
                                }
                                
                                if(foundPoint)
                                    break ;
                            }
                        }
                    }
                    
                    if(!done[nodes_per_plane*plane+side*(nodes_per_side+1)+node])
                    {
                        if(foundPoint)
                        {
                            newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = foundPoint ;
                        }
                        else
                        {
                            additionalPoints.push_back(new Point(proto) ) ;
                            newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = additionalPoints.back();
                            newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]->getId() = global_counter++ ;
                        }
                        
                        done[nodes_per_plane*plane+side*(nodes_per_side+1)+node] = true ;
                    }
                }
            }
        }
        
        i->setBoundingPoints(newPoints) ;
    }


}

void ParallelDelaunayTree::setElementOrder( Order elemOrder, double dt )
{
    switch( elemOrder )
    {
        case CONSTANT:
        {
            break ;
        }
        case LINEAR:
        {
            break ;
        }
        case QUADRATIC:
        {
            addSharedNodes( 1, 1, 0 ) ;
            break ;
        }
        case CUBIC:
        {
            addSharedNodes( 2, 1, 0 ) ;
            break ;
        }
        case QUADRIC:
        {
            addSharedNodes( 3, 1, 0 ) ;
            break ;
        }
        case QUINTIC:
        {
            addSharedNodes( 3, 1, 0 ) ;
            break ;
        }
        case CONSTANT_TIME_LINEAR:
        {
            addSharedNodes( 0, 2, dt ) ;
            break ;
        }
        case CONSTANT_TIME_QUADRATIC:
        {
            addSharedNodes( 0, 3, dt ) ;
            break ;
        }
        case LINEAR_TIME_LINEAR:
        {
            addSharedNodes( 0, 2, dt ) ;
            break ;
        }
        case LINEAR_TIME_QUADRATIC:
        {
            addSharedNodes( 0, 3, dt ) ;
            break ;
        }
        case QUADRATIC_TIME_LINEAR:
        {
            addSharedNodes( 1, 2, dt ) ;
            break ;
        }
        case QUADRATIC_TIME_QUADRATIC:
        {
            addSharedNodes( 1, 3, dt ) ;
            break ;
        }
        case CUBIC_TIME_LINEAR:
        {
            addSharedNodes( 2, 2, dt ) ;
            break ;
        }
        case CUBIC_TIME_QUADRATIC:
        {
            addSharedNodes( 2, 3, dt ) ;
            break ;
        }
        case QUADRIC_TIME_LINEAR:
        {
            addSharedNodes( 3, 2, dt ) ;
            break ;
        }
        case QUADRIC_TIME_QUADRATIC:
        {
            addSharedNodes( 3, 3, dt ) ;
            break ;
        }
        case QUINTIC_TIME_LINEAR:
        {
            addSharedNodes( 3, 2, dt ) ;
            break ;
        }
        case QUINTIC_TIME_QUADRATIC:
        {
            addSharedNodes( 3, 3, dt ) ;
            break ;
        }
        default:
            break ;

    }
}

void ParallelDelaunayTree::extrude(double dt)
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->extrude(dt) ;
    }

    int initial_global_counter = global_counter ;
    

    std::vector<DelaunayTriangle *> tmp = getElements() ;
    for(size_t j = 0 ; j < tmp.size() ;  j++)
    {

        for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
        {
            if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
        }
    }

}

void ParallelDelaunayTree::extrude(const Vector & dt)
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->extrude(dt) ;
    }

    int initial_global_counter = global_counter ;

    std::vector<DelaunayTriangle *> tmp = getElements() ;
    for(size_t j = 0 ; j < tmp.size() ;  j++)
    {
        for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
        {
            if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
        }
    }
}

std::vector<Point * > & ParallelDelaunayTree::getAdditionalPoints()
{
    return additionalPoints ;
}

const std::vector<Point * > & ParallelDelaunayTree::getAdditionalPoints() const
{
    return additionalPoints ;
}

std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree()
{
    return tree ;
}

const std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree() const
{
    return tree ;
}

