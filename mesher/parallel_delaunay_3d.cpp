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

bool ParallelDelaunayTree3D::isSame(const DelaunayTreeItem3D * i0, const DelaunayTreeItem3D * i1) const
{
    if(i0->isTetrahedron() && !i1->isTetrahedron())
        return false ;
    if(i0->isSpace() && !i1->isSpace())
        return false ;
    return i0->isVertex( i1->first ) && i0->isVertex( i1->second ) && i0->isVertex( i1->third )  && i0->isVertex( i1->fourth ) ;
}

bool ParallelDelaunayTree3D::inDomain( int domain_index, const DelaunayTetrahedron * tet) const
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

int ParallelDelaunayTree3D::getMesh(const DelaunayTreeItem3D * self) const
{
    for(size_t i = 0 ; i < meshes.size() ; i++)
        if(self->tree == meshes[i])
            return i ;
    return -1 ;
}

int ParallelDelaunayTree3D::getDomain(const DelaunayTetrahedron * tet) const
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

int ParallelDelaunayTree3D::getDomain(const Point & center) const
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

ParallelDelaunayTree3D::ParallelDelaunayTree3D(Point * p0,  Point *p1,  Point *p2,  Point *p3, const std::vector<const Geometry *> & domains) : Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D>(SPACE_THREE_DIMENSIONAL), domains(domains)
{
    std::vector<std::vector<DelaunayTreeItem3D *> > newElements(domains.size()) ;

    for(size_t i = 0 ; i < domains.size() ; i++)
    {
        meshes.push_back(new DelaunayTree3D(p0, p1, p2, p3));
    }

    global_counter = meshes[0]->global_counter ;
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getConflictingElements(const Point  * p)
{
    std::valarray<std::vector<DelaunayTetrahedron *> > conflicts(meshes.size()) ;
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ; j++)
        {
            if(getDomain(tmpConflicts[i]) != 1)
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
    
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {
        std::vector<DelaunayTetrahedron *> tmpConflicts = meshes[i]->getConflictingElements(p) ;
        for(size_t j = 0 ; j < tmpConflicts.size() ;  j++)
        {
            if(inDomain( i, tmpConflicts[j]))
                conflicts[i].push_back(tmpConflicts[j]);
        }
    }

    std::vector<DelaunayTetrahedron *> ret = conflicts[0] ;
    for(size_t i = 1 ; i < conflicts.size() ; i++)
        ret.insert(ret.end(), conflicts[i].begin(), conflicts[i].end());

    return ret ;
}

bool interacts(const Geometry * g , DelaunayTreeItem3D * item)
{

//     checkfather = false ;
    if( item->isSpace() )
    {
        return false ;
    }
    else if( item->isTetrahedron())
    {
        if(g->in(*item->first)  || 
           g->in(*item->second) || 
           g->in(*item->third)  || 
           g->in(*item->fourth) 
          )
        {
            return true ;
        }
            
        for(size_t j = 0 ; j < item->neighbour.size() ;  j++)
        {
            if(item->getNeighbour(j)->isSpace())
                continue ;
            
            DelaunayTreeItem3D * n = item->getNeighbour(j) ;
            
            if(g->in(*(n->first))  || 
               g->in(*(n->second)) || 
               g->in(*(n->third))  || 
               g->in(*(n->fourth))  )
            {
                return true ;
            }
        }

    }    

    return false ;
}

void ParallelDelaunayTree3D::insert(Point * p)
{   
    #pragma omp parallel for schedule(static,1)
    for(size_t i = 0 ; i < meshes.size() ; i++)
    {   
        bool isVertex = false ;
        bool noInteraction = true ;
        std::vector<DelaunayTreeItem3D *> cons = meshes[i]->conflicts(p) ;
        std::vector<DelaunayTreeItem3D *> newElems ;
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
        if(isVertex)
            continue ;
        
        if(domains[i]->in(*p))
             noInteraction = false ;
        else
        {
            for(auto & c : cons )
            {
                if(interacts(domains[i], c))
                {
                    noInteraction = false ;
                    break ;
                }
            }
        }

        if(!noInteraction || global_counter < 9 )
        {
            newElems = meshes[i]->addElements(cons, p) ;
        }
        
        int maxIdx = 0 ;
        for(size_t j = 0 ; j < newElems.size() ; j++)
            maxIdx = std::max(maxIdx, newElems[j]->index) ;


        if(!isVertex && i == meshes.size()-1)
            p->getId()= global_counter++ ;
    }
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getElements() const
{
    std::vector<std::vector<DelaunayTetrahedron *>> tris(meshes.size()) ;

    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        std::vector<DelaunayTetrahedron *> tmp = meshes[i]->getElements() ;
        for(const auto & t : tmp)
        {

            if(getDomain( t) == i)
                tris[i].push_back(t);
        }
    }

    std::vector<DelaunayTetrahedron *> ret ;
    for(size_t i = 0 ; i < tris.size() ;  i++)
        ret.insert(ret.end(), tris[i].begin(), tris[i].end());

    return ret ;
}

std::vector<DelaunayTetrahedron *> ParallelDelaunayTree3D::getNeighbourhood(DelaunayTetrahedron * element) const
{
    std::vector<DelaunayTetrahedron *> ret ;
    int domain = getDomain(element) ;
    if(domain == -1)
    {
        for(const auto & idx : element->neighbourhood)
        {
            ret.push_back((DelaunayTetrahedron *)meshes[domain]->tree[idx]);
        }
        return ret ;
    }
    else
    {
        for(const auto & idx : element->neighbourhood)
        {
            DelaunayTetrahedron * n = (DelaunayTetrahedron *)meshes[domain]->tree[idx] ;
            if(getDomain(n) == domain)
                ret.push_back(n);
            else
            {
                int altDomain = getDomain(n->getCenter()) ;
                DelaunayTetrahedron * an = meshes[altDomain]->getUniqueConflictingElement(&n->getCenter()) ;
                ret.push_back(an);
            }
        }
        return ret ;
    }
};

void ParallelDelaunayTree3D::addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep)
{
    std::vector<DelaunayTetrahedron *> tet = getElements() ;
    std::map<DelaunayTetrahedron *, bool> visited ;
    for(const auto & t : tet)
        visited[t] = false ;

    if( nodes_per_side > 1 )
        nodes_per_side = 1 ;

    std::vector<size_t> positions ;

    if( nodes_per_side )
    {
        positions.push_back( 0 ) ;
        positions.push_back( 1 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 2 ) ;
        positions.push_back( 3 ) ;
        positions.push_back( 4 ) ;

        positions.push_back( 4 ) ;
        positions.push_back( 5 ) ;
        positions.push_back( 6 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 7 ) ;
        positions.push_back( 6 ) ;

        positions.push_back( 6 ) ;
        positions.push_back( 8 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 9 ) ;
        positions.push_back( 4 ) ;
    }
    else
    {
        positions.push_back( 0 ) ;
        positions.push_back( 1 ) ;

        positions.push_back( 1 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 2 ) ;
        positions.push_back( 3 ) ;

        positions.push_back( 3 ) ;
        positions.push_back( 0 ) ;

        positions.push_back( 3 ) ;
        positions.push_back( 1 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 2 ) ;
    }
    
    std::cerr << "setting order... elements 0/" << tet.size() << std::flush ;
    for( size_t i = 0 ; i < tet.size() ; i++ )
    {
        if(i%1000 == 0)
            std::cerr << "\rsetting order... elements " << i+1 << "/" << tet.size() << std::flush ;

        std::vector<std::pair<Point, Point> > sides ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 1 ), tet[i]->getBoundingPoint( 2 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 2 ), tet[i]->getBoundingPoint( 3 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 0 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 2 ) ) ) ;

        
        visited[tet[i]] = true ;

        size_t nodes_per_plane = nodes_per_side * 6 + 4 ;

        std::valarray<Point *> newPoints( ( Point * )nullptr, nodes_per_plane * time_planes ) ;
        std::valarray<bool> done( false, nodes_per_plane * time_planes ) ;

        for( size_t plane = 0 ; plane < time_planes ; plane++ )
        {
            size_t current = 0 ;

            for( size_t side = 0 ; side < 6 ; side++ )
            {
                Point a( sides[side].first ) ;
                Point b( sides[side].second ) ;

                if( time_planes > 1 )
                {
                    a.getT() = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
                    b.getT() = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
                }

                for( size_t node = 0 ; node < nodes_per_side + 2 ; node++ )
                {
                    double fraction = ( double )( node ) / ( ( double )nodes_per_side + 1 ) ;
                    Point proto = a * ( 1. - fraction ) + b * fraction ;
                    Point *foundPoint = nullptr ;

                    for( size_t j = 0 ; j < tet[i]->getBoundingPoints().size() ; j++ )
                    {
                        if( tet[i]->getBoundingPoint( j ) == proto )
                        {
                            foundPoint = &tet[i]->getBoundingPoint( j ) ;
                            break ;
                        }
                    }

                    if( !foundPoint )
                    {
                        std::vector<DelaunayTetrahedron *> toTest =  getNeighbourhood(tet[i]);

                        for( size_t j = 0 ; j < toTest.size() ; j++ )
                        {
                            if( visited[toTest[j]] && tet[i] != toTest[j])
                            {
                                for( size_t k = 0 ; k < toTest[j]->getBoundingPoints().size() ; k++ )
                                {
                                    if( toTest[j]->getBoundingPoint( k ) == proto )
                                    {
                                        foundPoint = &toTest[j]->getBoundingPoint( k ) ;
                                        break ;
                                    }
                                }

                                if( foundPoint )
                                {
                                    break ;
                                }
                            }
                        }

                    }

                    if( !done[positions[current] + plane * nodes_per_plane] )
                    {
                        if( foundPoint )
                        {
                            newPoints[positions[current] + plane * nodes_per_plane]  = foundPoint ;
                        }
                        else
                        {
                            additionalPoints.push_back( new Point( proto ) ) ;
                            newPoints[positions[current] + plane * nodes_per_plane]  = additionalPoints.back() ;
                            newPoints[positions[current] + plane * nodes_per_plane]->getId() = global_counter++ ;
                        }

                        done[positions[current] + plane * nodes_per_plane] = true ;
                    }

                    current++ ;
                }
            }
        }

        for( size_t k = 0 ; k < newPoints.size() ; k++ )
        {
            if( !newPoints[k] )
            {
                std::cout << "ouch !" << std::endl ;

                for( size_t k = 0 ; k < newPoints.size() ; k++ )
                    if( newPoints[k] )
                        newPoints[k]->print() ;

                exit( 0 ) ;
            }
        }

        tet[i]->setBoundingPoints( newPoints ) ;

    }
    
    std::cerr << "setting order... elements " << tet.size() << "/" << tet.size() << " ...done."<< std::endl ;
}

void ParallelDelaunayTree3D::setElementOrder( Order elemOrder, double dt )
{
    if(allElementsCacheID != -1)
    {
        caches[allElementsCacheID].clear() ;
        coefs[allElementsCacheID].clear() ;
        elementMap[allElementsCacheID].clear() ;
        allElementsCacheID = -1 ;
    }
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

void ParallelDelaunayTree3D::extrude(double dt)
{
    #pragma omp parallel for
    for(size_t i = 0 ; i < meshes.size() ;  i++)
    {
        meshes[i]->global_counter = global_counter ;
        meshes[i]->extrude(dt) ;
    }

    int initial_global_counter = global_counter ;
    

    std::vector<DelaunayTetrahedron *> tmp = getElements() ;
    for(size_t j = 0 ; j < tmp.size() ;  j++)
    {

        for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
        {
            if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
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

    std::vector<DelaunayTetrahedron *> tmp = getElements() ;
    for(size_t j = 0 ; j < tmp.size() ;  j++)
    {
        for(size_t k = 0 ; k < tmp[j]->getBoundingPoints().size() ; k++)
        {
            if( tmp[j]->getBoundingPoint(k).getId() >= initial_global_counter )
                tmp[j]->getBoundingPoint(k).getId() = global_counter++ ;
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

// std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree()
// {
//     return tree ;
// }
// 
// const std::vector<DelaunayTreeItem3D *> & ParallelDelaunayTree3D::getTree() const
// {
//     return tree ;
// }


unsigned int ParallelDelaunayTree3D::generateCache(const Geometry * locus, const Geometry * source, Function smoothing )
{
    VirtualMachine vm ;
    std::vector<double> co ;
    std::vector<DelaunayTetrahedron *> elems = getConflictingElements(locus) ;
    //search for first empty cache slot ;
    if(caches.empty())
    {
        caches.push_back(std::vector<int>());
        coefs.push_back(std::vector<std::vector<double>>());
        elementMap.push_back(std::vector<int>());
    }
    size_t position = 0;
    for( ; position < caches.size() ; position++)
    {
        if(caches[position].empty())
            break ;
    }
    if(position == caches.size())
    {
        caches.push_back(std::vector<int>());
        coefs.push_back(std::vector<std::vector<double>>());
        elementMap.push_back(std::vector<int>());
    }

    for(auto & element : elems)
    {
        if(source && element->getBehaviour()->getSource() != source)
            continue ;
        
        if(locus->in(element->getCenter()))
        {
            caches[position].push_back(element->index) ;
            elementMap[position].push_back(getDomain(element));
            coefs[position].push_back(std::vector<double>()) ;
            Function x = element->getXTransform() ;
            Function y = element->getYTransform() ;
            Function z = element->getZTransform() ;
            Function t = element->getTTransform() ;
            for(size_t i = 0 ; i < element->getGaussPoints().gaussPoints.size() ; i++)
            {
                double xx = vm.eval(x, element->getGaussPoints().gaussPoints[i].first) ;
                double xy = vm.eval(y, element->getGaussPoints().gaussPoints[i].first) ;
                double xz = vm.eval(z, element->getGaussPoints().gaussPoints[i].first) ;
                double xt = vm.eval(t, element->getGaussPoints().gaussPoints[i].first) ;

                coefs[position].back().push_back(vm.eval(smoothing, xx, xy, xz, xt));
            }
        }
    }
    
    return position ;
}

 unsigned int ParallelDelaunayTree3D::generateCache ()
{
    std::vector<DelaunayTetrahedron *> elems = getElements() ;
    //search for first empty cache slot ;
    if(caches.empty())
    {
        caches.push_back(std::vector<int>());
        coefs.push_back(std::vector<std::vector<double>>());
        elementMap.push_back(std::vector<int>());
    }
    size_t position = 0;
    for( ; position < caches.size() ; position++)
    {
        if(caches[position].empty())
            break ;
    }
    if(position == caches.size())
    {
        caches.push_back(std::vector<int>());
        coefs.push_back(std::vector<std::vector<double>>());
        elementMap.push_back(std::vector<int>());
    }

    for(auto & element : elems)
    {
        caches[position].push_back(element->index) ;
        elementMap[position].push_back(getDomain(element));
        coefs[position].push_back(std::vector<double>()) ;

        for(size_t i = 0 ; i < element->getGaussPoints().gaussPoints.size() ; i++)
        {
            coefs[position].back().push_back(1);
        }

    }
    allElementsCacheID = position ;
    return position ;
};

Vector ParallelDelaunayTree3D::getField( FieldType f, unsigned int cacheID, int dummy, double t) const
{
    VirtualMachine vm ;
    size_t blocks = 0 ;
    std::vector<int> treeStarts ;
    treeStarts.push_back(caches[cacheID][0]) ;
    for(size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++)
    {
        blocks = static_cast<DelaunayTetrahedron *>(meshes[elementMap[cacheID][i]]->getInTree(caches[cacheID][i]))->getBehaviour()->getNumberOfDegreesOfFreedom()/3 ;
    }
    Vector ret(0., fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ;
    Vector buffer(ret) ;
    double w = 0 ;
    for(size_t i = 0 ; i < caches[cacheID].size() ; i++)
    {
        double v = static_cast<DelaunayTetrahedron *>(meshes[elementMap[cacheID][i]]->getInTree(caches[cacheID][i]))->getState().getAverageField(f, buffer, &vm, dummy, t, coefs[cacheID][i]) ;
        ret += buffer * v ;
        w +=v ;
    }
    return ret/w ;
}

Vector ParallelDelaunayTree3D::getField( FieldType f, int dummy, double t) const
{
    VirtualMachine vm ;
    size_t blocks = 0 ;
    
    std::vector<DelaunayTetrahedron *> elems = getElements() ;
    for(size_t i = 0 ; i < elems.size() && !blocks; i++)
    {
        blocks = elems[i]->getBehaviour()->getNumberOfDegreesOfFreedom()/3 ;
    }
    Vector ret(0., fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ;
    Vector buffer(ret) ;
    double w = 0 ;
    for(size_t i = 0 ; i < elems.size() ; i++)
    { 
        double v = elems[i]->getState().getAverageField(f, buffer, &vm, dummy, t) ;
        ret += buffer * v ;
        w +=v ;
    }
    return ret/w ;
}

    Vector ParallelDelaunayTree3D::getSmoothedField (  FieldType f0, unsigned int cacheID, IntegrableEntity * e,int dummy , double t ) const {
        Vector first ;
        Vector strain ;
        Vector stress ;
        Vector strainrate ;
        Vector buffer ;
        int tsize = 6 ;
        int psize = 3 ;

        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        VirtualMachine vm ;
        if ( f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
            //we first need to compute the strain field
            if ( !spaceTime ) {
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

                    double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, &vm, 0, t, coefs[cacheID][i] );
                    if ( !strain.size() ) {
                        strain.resize ( 0., buffer.size() );
                    }
                    strain += buffer*v ;
                    sumFactors += v ;
                }
                strain /= sumFactors ;
            } else {
                double sumFactors ( 0 ) ;
                Vector tmpstrain ;
                Vector tmpstrainrate ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;


                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, &vm, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrain.size() ) {
                        tmpstrain.resize ( 0., buffer.size() );
                    }
                    tmpstrain += buffer*v ;
                    sumFactors += v ;
                }
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, &vm, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrainrate.size() ) {
                        tmpstrainrate.resize ( 0., buffer.size() );
                    }
                    tmpstrainrate += buffer*v ;
                }
                tmpstrain /= sumFactors ;
                tmpstrainrate /=sumFactors ;

                Vector tmpstress = tmpstrain*e->getBehaviour()->getTensor ( Point() ) + ( Vector ) ( tmpstrainrate*e->getBehaviour()->getViscousTensor ( Point() ) ) ;
                stress.resize ( tsize, 0. ) ;
                strain.resize ( tsize, 0. ) ;
                for ( size_t i = 0 ; i < tsize ; i++ ) {
                    stress[i] = tmpstress[i] ;
                    strain[i] = tmpstrain[i] ;
                }
            }

            if ( f0 == PRINCIPAL_STRAIN_FIELD ) {
                first.resize ( psize );
                first = toPrincipal ( strain ) ;
            }
            if ( f0 == REAL_STRESS_FIELD ) {
                first.resize ( tsize );
                if ( !spaceTime ) {
                    first = strain*e->getBehaviour()->getTensor ( e->getCenter() ) ;
                } else {
                    first = stress ;
                }
            }
            if ( f0 == EFFECTIVE_STRESS_FIELD ) {
                first.resize ( tsize );
                first = strain*e->getBehaviour()->param ;
            }
            if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                first.resize ( psize );
                first = toPrincipal ( strain*e->getBehaviour()->param ) ;
            }
            if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
                first.resize ( psize );
                if ( !spaceTime ) {
                    first = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() ) ) ;
                } else {
                    first = toPrincipal ( stress ) ;
                }
            }

        } else {
            double sumFactors ( 0 ) ;
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                double v = ci->getState().getAverageField ( f0, buffer, &vm, dummy, t, coefs[cacheID][i] );
                if ( first.size() != buffer.size()) {
                    first.resize ( buffer.size(), 0. );
                }
                
                first += buffer*v ;
                sumFactors += v ;
            }

            first /= sumFactors ;
        }

        return first ;
    }

    std::pair<Vector, Vector> ParallelDelaunayTree3D::getSmoothedFields ( FieldType f0, FieldType f1, unsigned int cacheID, IntegrableEntity * e ,int dummy, double t ) const {
        Vector first ;
        Vector second ;
        Vector strain ;
        Vector stress ;
        Vector strainrate ;
        Vector buffer ;
        int tsize = 6 ;
        int psize = 3 ;

        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        VirtualMachine vm ;
        if ( f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f0 == STRAIN_FIELD ||
                f1 == PRINCIPAL_STRAIN_FIELD || f1 == REAL_STRESS_FIELD || f1 == EFFECTIVE_STRESS_FIELD || f1 == PRINCIPAL_REAL_STRESS_FIELD || f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f1 == STRAIN_FIELD
           ) {
            //we first need to compute the strain field
            if ( !spaceTime ) {
                buffer.resize ( tsize, 0. );
                strain.resize ( tsize, 0. );
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                    if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                        double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, &vm, dummy, t, coefs[cacheID][i] );
                        strain += buffer*v ;
                        sumFactors += v ;
                    }
                }
                strain /= sumFactors ;
            } else {
                size_t blocks = 0 ;
                for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
                    DelaunayTetrahedron *ci  = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                    blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions  ;
                }
                Vector tmpstrain ;
                Vector tmpstrainrate ;

                tmpstrain.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. ) ;
                buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. );
                tmpstrainrate.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, spaceDimensions, blocks ), 0. ) ;
                double sumFactors ( 0 ) ;


                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;


                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, &vm, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrain.size() ) {
                        tmpstrain.resize ( buffer.size(), 0. );
                    }
                    tmpstrain += buffer*v ;
                    sumFactors += v ;
                }
                buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD,spaceDimensions, blocks ), 0. );
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, &vm, 0, t, coefs[cacheID][i] );
                    if ( !tmpstrainrate.size() ) {
                        tmpstrainrate.resize ( buffer.size(), 0. );
                    }
                    tmpstrainrate += buffer*v ;
                }
                tmpstrain /= sumFactors ;
                tmpstrainrate /=sumFactors ;

                Vector tmpstress = tmpstrain*e->getBehaviour()->getTensor ( Point() ) + ( Vector ) ( tmpstrainrate*e->getBehaviour()->getViscousTensor ( Point() ) ) ;
                stress.resize ( tsize, 0. ) ;
                strain.resize ( tsize, 0. ) ;
                for ( size_t i = 0 ; i < tsize ; i++ ) {
                    stress[i] = tmpstress[i] ;
                    strain[i] = tmpstrain[i] ;
                }
            }

            if ( f0 == PRINCIPAL_STRAIN_FIELD ) {
                first.resize ( psize );
                first = toPrincipal ( strain ) ;
            }
            if ( f1 == PRINCIPAL_STRAIN_FIELD ) {
                second.resize ( psize );
                second = toPrincipal ( strain ) ;
            }
            if ( f0 == REAL_STRESS_FIELD ) {
                first.resize ( tsize );
                if ( !spaceTime ) {
                    first = strain*e->getBehaviour()->getTensor ( e->getCenter() ) ;
                } else {
                    first = stress ;
                }
            }
            if ( f1 == REAL_STRESS_FIELD ) {
                second.resize ( tsize );
                if ( !spaceTime ) {
                    second = strain*e->getBehaviour()->getTensor ( e->getCenter() ) ;
                } else {
                    second = stress ;
                }
            }
            if ( f0 == EFFECTIVE_STRESS_FIELD ) {
                first.resize ( tsize );
                first = strain*e->getBehaviour()->param ;
            }
            if ( f1 == EFFECTIVE_STRESS_FIELD ) {
                second.resize ( tsize );
                second = strain*e->getBehaviour()->param ;
            }
            if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                first.resize ( psize );
                first = toPrincipal ( strain*e->getBehaviour()->param ) ;
            }
            if ( f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                second.resize ( psize );
                second = toPrincipal ( strain*e->getBehaviour()->param ) ;
            }
            if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
                first.resize ( psize );
                if ( !spaceTime ) {
                    first = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() ) ) ;
                } else {
                    first = toPrincipal ( stress ) ;
                }
            }
            if ( f1 == PRINCIPAL_REAL_STRESS_FIELD ) {
                second.resize ( psize );
                if ( !spaceTime ) {
                    second = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() ) ) ;
                } else {
                    second = toPrincipal ( stress ) ;
                }
            }

            if ( f0 == STRAIN_FIELD ) {
                first.resize ( tsize ) ;
                first = strain ;
            }
            if ( f1 == STRAIN_FIELD ) {
                second.resize ( tsize ) ;
                second = strain ;
            }

        } else {
            double sumFactors ( 0 ) ;

            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

                double v = ci->getState().getAverageField ( f0, buffer, &vm, dummy, t, coefs[cacheID][i] );
                if ( !first.size() ) {
                    first.resize ( 0., buffer.size() );
                }
                first += buffer*v ;
                sumFactors += v ;
            }
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                    double v = ci->getState().getAverageField ( f1, buffer, &vm, dummy, t,coefs[cacheID][i] );
                    if ( !second.size() ) {
                        second.resize ( 0., buffer.size() );
                    }
                    second += buffer*v ;
                }
            }
            first /= sumFactors ;
            second /= sumFactors ;
        }


        return std::make_pair ( first, second ) ;
    }
    