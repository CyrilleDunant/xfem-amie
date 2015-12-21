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
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <limits>
#include <algorithm>

// #define DEBUG
// #undef DEBUG

namespace Amie
{

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
            if( unique && mesh == (int)domain_index )
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

ParallelDelaunayTree::ParallelDelaunayTree(Point * p0,  Point *p1,  Point *p2, const std::vector<const Geometry *> & domains) : Mesh< Amie::DelaunayTriangle, Amie::DelaunayTreeItem >(SPACE_TWO_DIMENSIONAL), domains(domains)
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

            if(getDomain(t) == (int)i)
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
}

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

// std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree()
// {
//     return tree ;
// }
// 
// const std::vector<DelaunayTreeItem *> & ParallelDelaunayTree::getTree() const
// {
//     return tree ;
// }

unsigned int ParallelDelaunayTree::generateCache(const Geometry * locus, const Geometry * source, Function smoothing )
{
    VirtualMachine vm ;
    std::vector<double> co ;
    std::vector<DelaunayTriangle *> elems = getConflictingElements(locus) ;
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
            if(!source)
               continue ;
            Function x = element->getXTransform() ;
            Function y = element->getYTransform() ;
            Function t = element->getTTransform() ;
            for(size_t i = 0 ; i < element->getGaussPoints().gaussPoints.size() ; i++)
            {
                double xx = vm.eval(x, element->getGaussPoints().gaussPoints[i].first) ;
                double xy = vm.eval(y, element->getGaussPoints().gaussPoints[i].first) ;
                double xz = 0 ;
                double xt = vm.eval(t, element->getGaussPoints().gaussPoints[i].first) ;

                coefs[position].back().push_back(vm.eval(smoothing, xx, xy, xz, xt));
            }
        }
    }
    
    return position ;
}

unsigned int ParallelDelaunayTree::generateCache ()
{
    std::vector<DelaunayTriangle *> elems = getElements() ;
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
}

Vector ParallelDelaunayTree::getField( FieldType f, int cacheID, double t, int index) 
{
    VirtualMachine vm ;
    size_t blocks = 0 ;
    std::vector<int> treeStarts ;
    treeStarts.push_back(caches[cacheID][0]) ;
    for(size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++)
    {
        blocks = static_cast<DelaunayTriangle *>(meshes[elementMap[cacheID][i]]->getInTree(caches[cacheID][i]))->getBehaviour()->getNumberOfDegreesOfFreedom()/2 ;
    }
    Vector ret(0., fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ;
    Vector buffer(ret) ;
    double w = 0 ;
    for(size_t i = 0 ; i < caches[cacheID].size() ; i++)
    {
        double v = static_cast<DelaunayTriangle *>(meshes[elementMap[cacheID][i]]->getInTree(caches[cacheID][i]))->getState().getAverageField(f, buffer, nullptr, t, coefs[cacheID][i], index) ;
        ret += buffer * v ;
        w +=v ;
    }
    return ret/w ;
}

Vector ParallelDelaunayTree::getField( FieldType f, double t, int index) 
{
    VirtualMachine vm ;
    size_t blocks = 0 ;
    
    std::vector<DelaunayTriangle *> elems = getElements() ;
    for(size_t i = 0 ; i < elems.size() && !blocks; i++)
    {
        blocks = elems[i]->getBehaviour()->getNumberOfDegreesOfFreedom()/2 ;
    }
    Vector ret(0., fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ;
    Vector buffer(ret) ;
    double w = 0 ;
    std::vector<double> coeffs ;
    for(size_t i = 0 ; i < elems.size() ; i++)
    { 
        double v = elems[i]->getState().getAverageField(f, buffer, nullptr, t, coeffs, index) ;
        ret += buffer * v ;
        w +=v ;
    }
    return ret/w ;
}

    Vector ParallelDelaunayTree::getSmoothedField (  FieldType f0, int cacheID, IntegrableEntity * e, double t  , const std::vector<bool> & restrict, int index ) {
    Vector first ;
    Vector strain ;
    Vector stress ;
    Vector strainrate ;
    Vector buffer ;
    int tsize = 3 ;
    int psize = 2 ;

    bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
    VirtualMachine vm ;
    if ( f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
        //we first need to compute the strain field
        if ( !spaceTime ) {
            double sumFactors ( 0 ) ;

            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                IntegrableEntity *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

                double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
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
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;


                double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
                if ( !tmpstrain.size() ) {
                    tmpstrain.resize ( 0., buffer.size() );
                }
                tmpstrain += buffer*v ;
                sumFactors += v ;
            }
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
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
            for ( int i = 0 ; i < tsize ; i++ ) {
                stress[i] = tmpstress[i] ;
                strain[i] = tmpstrain[i] ;
            }
        }

        if ( f0 == PRINCIPAL_STRAIN_FIELD ) {
            first.resize ( psize );
            first = toPrincipal ( strain, DOUBLE_OFF_DIAGONAL_VALUES ) ;
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
            first = toPrincipal ( strain*e->getBehaviour()->param, SINGLE_OFF_DIAGONAL_VALUES ) ;
        }
        if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
            first.resize ( psize );
            if ( !spaceTime ) {
                first = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() ), SINGLE_OFF_DIAGONAL_VALUES ) ;
            } else {
                first = toPrincipal ( stress , SINGLE_OFF_DIAGONAL_VALUES) ;
            }
        }

    } else {
        double sumFactors ( 0 ) ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
            double v = ci->getState().getAverageField ( f0, buffer, nullptr, t, coefs[cacheID][i], index );
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

std::pair<Vector, Vector> ParallelDelaunayTree::getSmoothedFields ( FieldType f0, FieldType f1, int cacheID, IntegrableEntity * e , double t ,  const std::vector<bool> & restrict, int index0, int index1 ) {
    Vector first ;
    Vector second ;
    Vector strain ;
    Vector stress ;
    Vector strainrate ;
    Vector buffer ;
    int tsize = 3 ;
    int psize = 2 ;

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
                IntegrableEntity *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                    double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
                    strain += buffer*v ;
                    sumFactors += v ;
                }
            }
            strain /= sumFactors ;
        } else {
            size_t blocks = 0 ;
            for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
                DelaunayTriangle *ci  = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
                blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions  ;
            }
            Vector tmpstrain ;
            Vector tmpstrainrate ;

            tmpstrain.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. ) ;
            buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. );
            tmpstrainrate.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, spaceDimensions, blocks ), 0. ) ;
            double sumFactors ( 0 ) ;


            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;


                double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
                if ( !tmpstrain.size() ) {
                    tmpstrain.resize ( buffer.size(), 0. );
                }
                tmpstrain += buffer*v ;
                sumFactors += v ;
            }
            buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD,spaceDimensions, blocks ), 0. );
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

                double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, nullptr, t, coefs[cacheID][i] );
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
            for ( int i = 0 ; i < tsize ; i++ ) {
                stress[i] = tmpstress[i] ;
                strain[i] = tmpstrain[i] ;
            }
        }

        if ( f0 == PRINCIPAL_STRAIN_FIELD ) {
            first.resize ( psize );
            first = toPrincipal ( strain , DOUBLE_OFF_DIAGONAL_VALUES) ;
        }
        if ( f1 == PRINCIPAL_STRAIN_FIELD ) {
            second.resize ( psize );
            second = toPrincipal ( strain , DOUBLE_OFF_DIAGONAL_VALUES) ;
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
                second = strain*e->getBehaviour()->getTensor ( e->getCenter()  ) ;
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
            first = toPrincipal ( strain*e->getBehaviour()->param  , SINGLE_OFF_DIAGONAL_VALUES) ;
        }
        if ( f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
            second.resize ( psize );
            second = toPrincipal ( strain*e->getBehaviour()->param  , SINGLE_OFF_DIAGONAL_VALUES) ;
        }
        if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
            first.resize ( psize );
            if ( !spaceTime ) {
                first = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() )  , SINGLE_OFF_DIAGONAL_VALUES) ;
            } else {
                first = toPrincipal ( stress  , SINGLE_OFF_DIAGONAL_VALUES) ;
            }
        }
        if ( f1 == PRINCIPAL_REAL_STRESS_FIELD ) {
            second.resize ( psize );
            if ( !spaceTime ) {
                second = toPrincipal ( strain*e->getBehaviour()->getTensor ( e->getCenter() )  , SINGLE_OFF_DIAGONAL_VALUES) ;
            } else {
                second = toPrincipal ( stress  , SINGLE_OFF_DIAGONAL_VALUES) ;
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
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;

            double v = ci->getState().getAverageField ( f0, buffer, nullptr, t, coefs[cacheID][i], index0 );
            if ( !first.size() ) {
                first.resize ( 0., buffer.size() );
            }
            first += buffer*v ;
            sumFactors += v ;
        }
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
            DelaunayTriangle *ci = static_cast<DelaunayTriangle *> ( meshes[elementMap[cacheID][i]]->getInTree ( caches[cacheID][i] ) ) ;
            if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                double v = ci->getState().getAverageField ( f1, buffer, nullptr,  t,coefs[cacheID][i], index1 );
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

}
