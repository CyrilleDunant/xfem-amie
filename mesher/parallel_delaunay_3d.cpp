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

#include "parallel_delaunay_3d.h"
#include <omp.h>
#include <limits>
#include <algorithm>

// #define DEBUG
// #undef DEBUG

using namespace Amie ;


ParallelDelaunayTree3D::ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2,  Point *p3, const std::vector<Geometry *> & domains) : domains(domains)
{
    
    for(size_t i = 0 ; i < domains.size() ; i++)
    {
        meshes.push_back(new DelaunayTree3D(p0, p1, p2, p3));
        elementMap.push_back(std::vector<int>());
        outElements.push_back(std::vector<bool>());
    }
    global_counter = meshes[0]->global_counter ;
    //there needs to be one more mesh than domains: it contains all the boundary elements
//     meshes.push_back(new DelaunayTree3D(p0, p1, p2, p3));
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getConflictingElements(const Point  * p)
{
    std::valarray<std::vector<DelaunayTetrahedron *> > conflicts(meshes.size()) ;
    std::valarray<std::vector<DelaunayTetrahedron *> > boundaries(meshes.size()) ;
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ; j++)
        {
            if(elementMap[i][tmpConflicts[j]->index] >= 0)
                conflicts[i].push_back(tmpConflicts[j]) ;
        }
    }
    std::vector<DelaunayTetrahedron *> ret = conflicts[0] ;
    for(size_t i = 1 ; i < conflicts.size() ; i++)
        ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());

    return ret ;
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getConflictingElements(const Geometry  * p)
{
    std::valarray<std::vector<DelaunayTetrahedron *> > conflicts(meshes.size()) ;
    std::valarray<std::vector<DelaunayTetrahedron *> > boundaries(meshes.size()) ;
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ; j++)
        {
            if(elementMap[i][tmpConflicts[j]->index] >= 0)
                conflicts[i].push_back(tmpConflicts[j]) ;
        }
    }

    std::vector<DelaunayTetrahedron *> ret = conflicts[0] ;
    for(size_t i = 1 ; i < conflicts.size() ; i++)
        ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());

    return ret ;
}

bool interacts(const Geometry * g , const DelaunayTreeItem3D * item)
{
    if( item->isSpace() )
    {
        std::vector<Point> bbox = g->getBoundingBox() ;
        
        for(size_t i = 0 ; i < bbox.size() ; i++)
            if(item->inCircumSphere(bbox[i]) || item->onCircumSphere(bbox[i]))
                return true ;
        return false ;
    }
    else if( item->isTetrahedron() && item->isAlive())
    {
        if(g->in(*item->first)  || 
           g->in(*item->second) || 
           g->in(*item->third)  || 
           g->in(*item->fourth) 
          )
            return true ;
            
        const DelaunayTetrahedron * t = static_cast<const DelaunayTetrahedron *>(item) ;
        Sphere test(t->getRadius(), t->getCircumCenter() ) ;
        if( g->intersects(&test) )
            return true ;
        
        for(size_t j = 0 ; j < t->neighbourhood.size() ;  j++)
        {
           if(g->in(*(t->getNeighbourhood(j)->first))  || 
                g->in(*(t->getNeighbourhood(j))->second) || 
                g->in(*(t->getNeighbourhood(j))->third)  || 
                g->in(*(t->getNeighbourhood(j))->fourth) 
                )
            return true ;
            
          Sphere test(t->getNeighbourhood(j)->getRadius(), t->getNeighbourhood(j)->getCenter() ) ;
          if( g->intersects(&test) )
            return true ;
          
          for(size_t k = 0 ; k < t->getNeighbourhood(j)->neighbourhood.size() ;  k++)
          {
            if(g->in(*(t->getNeighbourhood(j)->getNeighbourhood(k))->first)  || 
                g->in(*(t->getNeighbourhood(j)->getNeighbourhood(k))->second) || 
                g->in(*(t->getNeighbourhood(j)->getNeighbourhood(k))->third)  || 
                g->in(*(t->getNeighbourhood(j)->getNeighbourhood(k))->fourth) 
                )
                return true ;
                
            Sphere test(t->getNeighbourhood(j)->getNeighbourhood(k)->getRadius(), t->getNeighbourhood(j)->getNeighbourhood(k)->getCenter() ) ;
            if( g->intersects(&test) )
                return true ;
          }
        }
        
        return false ;

    }
    else if(item->isDeadTetrahedron())
    {
        if(g->in(*item->first)  || 
           g->in(*item->second) || 
           g->in(*item->third)  || 
           g->in(*item->fourth) 
          )
            return true ;
            
        Sphere test(item->getRadius(), static_cast<const DelaunayDeadTetrahedron *>(item)->center ) ;
        if( g->intersects(&test) )
            return true ;
        return false ;
    }
}

void ParallelDelaunayTree3D::insert(Point * p)
{
    std::valarray<std::vector<DelaunayTreeItem3D *>> newElems(meshes.size()) ;
    std::valarray<bool> insertion(false, meshes.size()) ;
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {   
        bool isVertex = false ;
        bool allout = true ;
        bool inMesh = true ;
        std::vector<DelaunayTreeItem3D *> cons = meshes[i]->conflicts(p) ;
        if(cons.empty())
        {
            inMesh = false ;
            std::cout << "Failed insertion : in nothing !" << std::endl ;
        }

        meshes[i]->neighbourhood = false ;


        if(!domains[i]->in(*p))
        {
            for(size_t j = 0 ; j < cons.size() ; j++)
            {
                if(cons[j]->isVertex(p))
                    isVertex = true ; 
            }
            
            for(size_t j = 0 ; j < cons.size() && allout; j++)
            {
                
                for(size_t k = 0 ; k < outElements[i].size() && allout ; k++)
                {
                    if(outElements[i][cons[j]->index])
                        allout = false ;
                }
                
            }
            
        }
        else
        {
            allout = false ;
        }
        
        if((!allout && !isVertex && inMesh) || global_counter < 8)
        {
            insertion[i] = true ;
            newElems[i] = meshes[i]->addElements(cons, p) ;
        }
        
    }

    p->getId()= global_counter++ ;
    
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        if(insertion[i])
        {
            for( size_t j = 0 ; j < newElems[i].size() ; j++ )
            {
                if(!newElems[i][j]->isTetrahedron())
                    continue ;
                
                static_cast<DelaunayTetrahedron *>(newElems[i][j])->neighbourhood.resize( 0 ) ;
                newElems[i][j]->clearVisited() ;
            }



            for( size_t j = 0 ; j < newElems[i].size() ; j++ )
            {
                if(!newElems[i][j]->isTetrahedron())
                    continue ;
                
                DelaunayTetrahedron * tri = static_cast<DelaunayTetrahedron *>(newElems[i][j]) ;
    //          if( i % 10000 == 0 )
    //              std::cerr << "\r building neighbourhood... element " << i << "/" << ret.size() << std::flush ;

                std::vector<DelaunayTetrahedron *> tocheck ;
                std::vector<DelaunayTetrahedron *> toclean ;

                for( size_t k = 0 ; k < tri->neighbour.size() ; k++ )
                {
                    if( !tri->getNeighbour( k )->visited() && tri->getNeighbour( k )->isTetrahedron() )
                    {
                        tocheck.push_back( static_cast<DelaunayTetrahedron *>( tri->getNeighbour( k ) ) );
                        tri->getNeighbour( k )->visited() = true ;
                        toclean.push_back( tocheck.back() ) ;
                        tri->addNeighbourhood( tocheck.back() ) ;
                    }
                }

                for( size_t k = 0 ; k < tocheck.size() ; k++ )
                {
                    if( tocheck[k]->numberOfCommonVertices( tri ) > 0 )
                        tri->addNeighbourhood( tocheck[k] ) ;
                }

                while( !tocheck.empty() )
                {
                    std::vector<DelaunayTetrahedron *> tocheck_temp ;

                    for( size_t k = 0 ; k < tocheck.size() ; k++ )
                    {
                        for( size_t l = 0 ; l < tocheck[k]->neighbour.size() ; l++ )
                        {
                            if( tocheck[k]->getNeighbour( l )->isTetrahedron()
                                && tocheck[k]->getNeighbour( l )->isAlive()
                                && !tocheck[k]->getNeighbour( l )->visited()
                                && tocheck[k]->getNeighbour( l ) != tri
                                && static_cast<DelaunayTetrahedron *>( tocheck[k]->getNeighbour( l ) )->numberOfCommonVertices( tri ) > 0
                            )
                            {
                                tocheck_temp.push_back( static_cast<DelaunayTetrahedron *>( tocheck[k]->getNeighbour( l ) ) );
                                tocheck[k]->getNeighbour( l )->visited() = true ;
                                toclean.push_back( tocheck_temp.back() ) ;
                                tri->addNeighbourhood( tocheck_temp.back() ) ;
                            }
                        }
                    }

                    tocheck = tocheck_temp ;
                }

                for( size_t k = 0 ; k < toclean.size() ; k++ )
                {
                    toclean[k]->clearVisited() ;
                }

            }
        }
    }     
         
    #pragma omp parallel for
    for(size_t i = 0 ; i < newElems.size() ; i++)
    {
//         if(!insertion[i])
//         {
//             std::cout << "pop !" << std::endl ;
//             continue ;
//         }
//         else
//             std::cout << "plouf !" << std::endl ;
        
        int maxIdx = 0 ;
        for(size_t j = 0 ; j < newElems[i].size() ; j++)
            maxIdx = std::max(maxIdx, newElems[i][j]->index) ;

        for(int j = elementMap[i].size() ; j < maxIdx+1 ; j++)
        {
            elementMap[i].push_back(0) ;
            outElements[i].push_back(true) ;
        }

        for(size_t j = 0 ; j < newElems[i].size() ; j++)
        {
            elementMap[i][newElems[i][j]->index] = newElems[i][j]->index ;
            outElements[i][newElems[i][j]->index] = interacts(domains[i], newElems[i][j]) ;
        }
    }
    
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        if(!insertion[i])
            continue ;
        for(size_t j = 0 ; j < newElems[i].size() ; j++)
        {
            for(size_t k = i+1 ; k < meshes.size() ; k++)
            {
                if(!insertion[k])
                    continue ;
                for(size_t l = 0 ; l < newElems[k].size() ; l++)
                {
                    if(newElems[i][j]->isVertex( newElems[k][l]->first)  &&
                       newElems[i][j]->isVertex( newElems[k][l]->second) &&
                       newElems[i][j]->isVertex( newElems[k][l]->third)  &&
                       newElems[i][j]->isVertex( newElems[k][l]->fourth) &&
                       elementMap[i][newElems[i][j]->index] >= 0)
                            elementMap[i][newElems[i][j]->index] = -newElems[i][j]->index ;
                }
            }
        }
    }
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getElements()
{
    std::vector<std::vector<DelaunayTetrahedron *>> tris(meshes.size()) ;

    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
        for(size_t j = 0 ; j < tmp.size() ;  j++)
        {
//             std::cout  << elementMap[i].size() << " vs " << tmp[j]->index << std::endl ;
            if(elementMap[i][tmp[j]->index] >= 0)
                tris[i].push_back(tmp[j]);
        }
    }

    std::vector<DelaunayTetrahedron *> ret ;
    for(size_t i = 0 ; i < tris.size() ;  i++)
        ret.insert(ret.end(), tris[i].begin(), tris[i].end());

    return ret ;
}

void ParallelDelaunayTree3D::setElementOrder(Order o, double dt )
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->setElementOrder(o, dt) ;
    }
    int initial_global_counter = global_counter ;
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
        for(size_t j = 0 ; j < tmp.size() ;  j++)
        {
            if(elementMap[i][tmp[j]->index] >= 0)
            {
                for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
                {
                    if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                        tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
                }
            }
        }
    }
}

void ParallelDelaunayTree3D::extrude(double dt)
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->extrude(dt) ;
    }

    int initial_global_counter = global_counter ;
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
        for(size_t j = 0 ; j < tmp.size() ;  j++)
        {
            if(elementMap[i][tmp[j]->index] >= 0)
            {
                for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
                {
                    if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                        tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
                }
            }
        }
    }
}

void ParallelDelaunayTree3D::extrude(const Vector & dt)
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->extrude(dt) ;
    }

    int initial_global_counter = global_counter ;
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
        for(size_t j = 0 ; j < tmp.size() ;  j++)
        {
            if(elementMap[i][tmp[j]->index] >= 0)
            {
                for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
                {
                    if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                        tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
                }
            }
        }
    }
}

std::vector<Point * > & ParallelDelaunayTree3D::getAdditionalPoints()
{
    return additionalPoints ;
}

const std::vector<Point * > & ParallelDelaunayTree3D::getAdditionalPoints() const
{
    return additionalPoints ;
}

std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree()
{
    return tree ;
}

const std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree() const
{
    return tree ;
}

