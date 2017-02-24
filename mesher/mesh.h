/*
    mesh abstract implementation for AMIE
    Copyright (C) 2010  Cyrille Dunant

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#ifndef MESH_H
#define MESH_H

#include <vector>
#include <limits>
#include <set>
#include "../utilities/matrixops.h"
#include "element_checker.h"
#include "../elements/integrable_entity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../physics/dual_behaviour.h"
#include "../physics/viscoelasticity.h"

//     inline Vector operator*(const Vector &v , const Amie::Matrix &m ) ;
//     inline Amie::MtV operator*(const Amie::Matrix& mm, const Vector& v) ;
//     inline Amie::MtM operator*(const Amie::Matrix& mm, const Amie::Matrix& mmm) ;
namespace Amie
{


template <class ETYPE, class EABSTRACTTYPE>
class Mesh
{
protected:
    SpaceDimensionality spaceDimensions ;
    std::map<const Mesh<ETYPE, EABSTRACTTYPE> *, std::map<Point *, std::pair<ETYPE *, std::vector<double> > > > cache ;
    std::map<const Mesh<ETYPE, EABSTRACTTYPE> *, std::map<Point *, Point * > > pointcache ;
    std::vector<std::vector<int>> caches ;
    std::vector<std::vector<std::vector<double>>> coefs ;
    int allElementsCacheID ;
    virtual std::vector<ETYPE *> getElements() = 0;

    std::vector<Vector> strainCache ;
    std::vector<Vector> strainGenViscoCache ;
    std::vector<Vector> strainRateGenViscoCache ;
    std::vector<Vector> stressCache ;
    std::vector<Vector> stressGenViscoCache ;
    std::vector<Vector> stressImposedCache ;
    std::vector<Vector> bufferCache ;
    std::vector<Vector> firstResultCache ;
    std::vector<Vector> secondResultCache ;

public:

    virtual std::vector<int> & getCache ( int cacheID ) {
        return caches[cacheID] ;
    } ;
    int getAllElementsCacheID() const { return allElementsCacheID ; }
    virtual int addToTree ( EABSTRACTTYPE * toAdd ) = 0 ;
    virtual EABSTRACTTYPE * getInTree ( int index ) const = 0 ;
    virtual std::vector<Point * > & getAdditionalPoints() = 0 ;
    virtual const std::vector<Point * > & getAdditionalPoints() const = 0 ;
    virtual void extrude ( double dt ) = 0 ;
    virtual void extrude ( const Vector & dt ) = 0 ;
    virtual double getInternalScale() const {
        return 1. ;
    } ;
    virtual bool valid() const { return true ; }

public:
    Mesh(SpaceDimensionality spaceDimensions) : spaceDimensions(spaceDimensions), allElementsCacheID(-1) {} ;
    virtual ~Mesh() {} ;

    virtual std::vector<ETYPE *> getConflictingElements ( const Point  * p )  = 0;
    virtual std::vector<ETYPE *> getConflictingElements ( const Geometry * g ) = 0;
    virtual std::vector<ETYPE *> getNeighbourhood ( ETYPE * element ) const = 0 ;
    virtual std::vector<ETYPE *> getConflictingElements ( const Segment  * s, double width = 0.001 )  
    {
        Point A = s->first() ;
        Point B = s->second() ;
        Point n = s->normal()*width ;
        OrientedRectangle rect( A-n, B-n, B+n, A+n) ;
        return getConflictingElements( &rect) ;
    }


    virtual bool step(double dt) {
        return false ;
    }

    virtual std::vector<ETYPE *> getNeighbouringElementsInGeometry ( ETYPE * start , const Geometry * g ) {
        if ( !start ) {
            return std::vector<ETYPE *>() ;
        }

        std::set<ETYPE *> to_test ;
        std::set<ETYPE *> found ;
        found.insert ( start ) ;
        std::vector<ETYPE *> neighbourhood = getNeighbourhood ( start ) ;
        for ( const auto & neighbour : neighbourhood ) {
            if ( neighbour->timePlanes() > 1 ) {
//                std::cout << "\n!!!!!!!!!! wahh?" << std::endl ;
                if ( g->in ( neighbour->getCenter() ) || neighbour->in ( g->getCenter() ) || g->intersects ( neighbour->getPrimitive() ) ) {
                    to_test.insert ( neighbour ) ;
                } else if ( g->getGeometryType() == TIME_DEPENDENT_CIRCLE ) {

                    for ( size_t j = 0 ; j < neighbour->timePlanes() ; j++ ) {
                        /*		            Point c = neighbour->getCenter() ;
                        		            c.getT() = neighbour->getBoundingPoint ( neighbour->getBoundingPoints().size() * j / neighbour->timePlanes() ).getT() ;
                        		            if ( g->in ( c ) ) {
                        		                to_test.insert ( neighbour ) ;
                        				break ;
                        		            }*/

                        for ( size_t k = 0 ; k <  neighbour->getBoundingPoints().size() / neighbour->timePlanes()-1 ; k++ ) {
                            Point A = neighbour->getBoundingPoint ( neighbour->getBoundingPoints().size() * j / neighbour->timePlanes() + k ) ;
                            Point B = neighbour->getBoundingPoint ( neighbour->getBoundingPoints().size() * j / neighbour->timePlanes() + k + 1 ) ;
                            Segment s ( A,B ) ;
                            if ( s.intersects ( g ) || g->in ( A ) || g->in ( B ) ) {
                                to_test.insert ( neighbour ) ;
                            }
                        }

                    }
                }
            } else if ( g->in ( neighbour->getCenter() ) || neighbour->in ( g->getCenter() ) || g->intersects ( neighbour->getPrimitive() ) ) {
                to_test.insert ( neighbour ) ;
            }
        }
        found.insert ( to_test.begin(), to_test.end() ) ;

        while ( !to_test.empty() ) {
            std::set<ETYPE *> new_test ;
            for ( auto j = to_test.begin() ; j != to_test.end() ; j++ ) {
                neighbourhood = getNeighbourhood ( *j ) ;
                for ( const auto & neighbour : neighbourhood ) {
                    if ( to_test.find ( neighbour ) == to_test.end()
                            && found.find ( neighbour ) == found.end() ) {

                        if ( neighbour->timePlanes() > 1 ) {
                            if ( g->in ( neighbour->getCenter() ) || neighbour->in ( g->getCenter() ) || g->intersects ( neighbour->getPrimitive() ) ) {
                                new_test.insert ( neighbour ) ;
                            } else if ( g->getGeometryType() == TIME_DEPENDENT_CIRCLE ) {

                                for ( size_t k = 0 ; k < neighbour->getBoundingPoints().size()-1 ; k++ ) {
                                    Point A = neighbour->getBoundingPoint ( k ) ;
                                    Point B = neighbour->getBoundingPoint ( k+1 ) ;
                                    Segment s ( A,B ) ;
                                    if ( g->in ( A ) || g->in ( B ) || s.intersects ( g ) ) {
                                        new_test.insert ( neighbour ) ;
                                        break ;
                                    }

                                }
                            }
                        } else if ( g->in ( neighbour->getCenter() ) || neighbour->in ( g->getCenter() ) || g->intersects ( neighbour->getPrimitive() ) ) {
                            new_test.insert ( neighbour ) ;
                        }

                    }
                }
            }
            to_test = new_test ;
            found.insert ( new_test.begin(), new_test.end() ) ;
        }

        std::vector<ETYPE *> ret ( found.begin(),found.end() ) ;
        return ret ;
    }

    virtual ETYPE * getUniqueConflictingElement ( const Point  * p ) {
        std::vector<ETYPE *> elements = getConflictingElements ( p ) ;
        for ( const auto & element : elements ) {
            if (  *p == element->getCenter() ) {
                return element ;
            }
        }
        for ( const auto & element : elements ) {
            if ( element->in ( *p ) ) {
                return element ;
            }
        }

        return nullptr ;
    }

    virtual void setElementOrder ( Order o, double dt = 0. ) = 0;
    virtual void insert ( Point * , double ) = 0 ;
    template <class ETARGETTYPE>
    /** \brief Return the displacements in source mesh projected on current mesh.
    */
    void project ( Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false ) {
        if ( cache.find ( mesh ) == cache.end() ) {
            std::map<Point *, std::pair<ETYPE *, std::vector<double> > > projectionCache ;
            std::map<Point *, Point * > projectionPointCache ;
            std::vector<ETYPE *> selfElements = getElements() ;
            size_t numDofs = 0 ;
            double rav = 0 ;
            double ecount  = 0;
            std::set<Point *> points ;
            for ( size_t i = 0 ; i < selfElements.size() ; i++ ) {
                if ( selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR ) {
                    numDofs = std::max ( selfElements[i]->getBehaviour()->getNumberOfDegreesOfFreedom(),numDofs ) ;
                    for ( size_t j = 0 ; j < selfElements[i]->getBoundingPoints().size() ; j++ ) {
                        points.insert ( &selfElements[i]->getBoundingPoint ( j ) ) ;
                    }
                    rav += selfElements[i]->getRadius() ;
                    ecount++ ;
                }

            }
            rav /= ecount ;
            if ( projection.size() == 0 ) {
                projection.resize ( numDofs*getLastNodeId() ) ;
                projection = 0 ;
            }

            for ( auto i = points.begin() ; i != points.end() ; i++ ) {
                std::vector<ETARGETTYPE *> targets ;

                std::vector<ETARGETTYPE *> temp = mesh->getConflictingElements ( *i ) ;
                for ( size_t k = 0 ; k < temp.size() ; k++ ) {
                    if ( temp[k]->getBehaviour()->type != VOID_BEHAVIOUR ) {
                        targets.push_back ( temp[k] ) ;
                    }
                }

                if ( !targets.empty() ) {
                    std::sort ( targets.begin(), targets.end() ) ;
                    auto e = std::unique ( targets.begin(), targets.end() ) ;
                    targets.erase ( e,targets.end() ) ;
                } else {
                    Circle c ( rav*2., * ( *i ) ) ;
                    targets =  mesh->getConflictingElements ( &c ) ;
                }


                if ( targets.empty() ) {
                    std::cout << "failed projection, empty mesh" << std::endl ;
                    ( *i )->print() ;
                    exit ( 0 ) ;
                }
                projectionPointCache[ ( *i )] = nullptr ;
                std::map<double, ETARGETTYPE *> coincidentElements;
                for ( size_t k = 0 ; k < targets.size() ; k++ ) {
                    Point proj ( * ( *i ) ) ;
                    targets[k]->project ( &proj ) ;
                    if ( targets[k]->in ( * ( *i ) ) || dist ( proj, * ( *i ) ) < 128.*POINT_TOLERANCE*POINT_TOLERANCE ) {
                        coincidentElements[dist ( * ( *i ), targets[k]->getCenter() )] = targets[k] ;
                    }
                }

                if ( !coincidentElements.empty() ) {
                    Vector disps ( 0., numDofs ) ;
                    projectionCache[ ( *i )] = std::make_pair ( coincidentElements.begin()->second, coincidentElements.begin()->second->getState().getInterpolatingFactors ( * ( *i ), false ) ) ;

                    for ( size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size(); j++ ) {
                        if ( coincidentElements.begin()->second->getBoundingPoint ( j ) == * ( *i ) ) {
                            projectionPointCache[ ( *i )] = &coincidentElements.begin()->second->getBoundingPoint ( j ) ;
                            break ;
                        }
                    }

                    if ( projectionPointCache[ ( *i )] && source.size() ) {
                        projectionCache.erase ( projectionCache.find ( *i ) ) ;
                        for ( size_t k = 0 ; k < numDofs ; k++ ) {
                            disps[k] = source[projectionPointCache[ ( *i )]->getId() *numDofs+k] ;
                        }
                    } else {
                        projectionPointCache.erase ( projectionPointCache.find ( *i ) ) ;
                        if ( coincidentElements.begin()->second->getBehaviour()->type != VOID_BEHAVIOUR ) {
                            for ( size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size() ; j++ ) {
                                int id = coincidentElements.begin()->second->getBoundingPoint ( j ).getId() ;
                                double d = projectionCache[ ( *i )].second[j] ;
                                for ( size_t k = 0 ; k < numDofs ; k++ ) {
                                    disps[k] += source[id*numDofs+k]*d ;
                                }
                            }
                        }
                    }

                    for ( size_t k = 0 ; k < numDofs ; k++ ) {
                        projection[ ( *i )->getId() *numDofs+k] = disps[k] ;
                    }
                } else {
                    for ( size_t k = 0 ; k < targets.size() ; k++ ) {
                        Point proj ( * ( *i ) ) ;
                        targets[k]->project ( &proj ) ;
                        coincidentElements[dist ( proj, * ( *i ) )] = targets[k] ;
                    }
                    Vector disps ( 0., numDofs );
                    projectionCache[ ( *i )] = std::make_pair ( coincidentElements.begin()->second, coincidentElements.begin()->second->getState().getInterpolatingFactors ( * ( *i ), false ) ) ;

                    for ( size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size(); j++ ) {
                        if ( coincidentElements.begin()->second->getBoundingPoint ( j ) == * ( *i ) ) {
                            projectionPointCache[ ( *i )] = &coincidentElements.begin()->second->getBoundingPoint ( j ) ;
                            break ;
                        }
                    }

                    if ( projectionPointCache[ ( *i )] && source.size() ) {
                        projectionCache.erase ( projectionCache.find ( *i ) ) ;
                        int id = projectionPointCache[ ( *i )]->getId() ;
                        for ( size_t k = 0 ; k < numDofs ; k++ ) {
                            disps[k] = source[id*numDofs+k] ;
                        }
                    } else {
                        projectionPointCache.erase ( projectionPointCache.find ( *i ) ) ;
                        if ( coincidentElements.begin()->second->getBehaviour()->type != VOID_BEHAVIOUR ) {
                            for ( size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size() ; j++ ) {
                                double d = projectionCache[ ( *i )].second[j] ;
                                int id = coincidentElements.begin()->second->getBoundingPoint ( j ).getId() ;
                                for ( size_t k = 0 ; k < numDofs ; k++ ) {
                                    disps[k] += source[id*numDofs+k]*d ;
                                }
                            }
                        }
                    }

                    for ( size_t k = 0 ; k < numDofs ; k++ ) {
                        projection[ ( *i )->getId() *numDofs+k] = disps[k] ;
                    }
                }

            }
            pointcache[mesh] = projectionPointCache ;
            cache[mesh] = projectionCache ;
        } else {
            std::vector<ETYPE *> selfElements = getElements() ;
            size_t numDofs = 0 ;
            for ( size_t i = 0 ; i < selfElements.size() ; i++ ) {
                if ( selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR ) {
                    numDofs = std::max ( selfElements[i]->getBehaviour()->getNumberOfDegreesOfFreedom(),numDofs ) ;
                    if ( numDofs ) {
                        break ;
                    }
                }
            }

            projection = 0 ;
            Vector disps ( 0.,numDofs ) ;

            for ( auto i = cache[mesh].begin() ; i != cache[mesh].end() ; ++i ) {
                disps = 0 ;
                if ( i->second.first->getBehaviour()->type != VOID_BEHAVIOUR ) {
                    for ( size_t j = 0 ; j < i->second.first->getBoundingPoints().size() ; j++ ) {
                        double d = i->second.second[j] ;
                        int id = i->second.first->getBoundingPoint ( j ).getId() ;
                        for ( size_t k = 0 ; k < numDofs ; k++ ) {
                            disps[k] += source[id*numDofs+k]*d ;
                        }
                    }
                }
                for ( size_t k = 0 ; k < numDofs ; k++ ) {
                    projection[i->first->getId() *numDofs+k] = disps[k] ;
                }

            }

            for ( auto j = pointcache[mesh].begin() ; j != pointcache[mesh].end() ; ++j ) {

                for ( size_t k = 0 ; k < numDofs ; k++ ) {
                    projection[j->first->getId() *numDofs+k] = source[j->second->getId() *numDofs+k] ;
                }

            }

        }

    } 

    template <class ETARGETTYPE>
    /** \brief Return the displacements in source mesh projected on current mesh.
    */
    void leastSquareProject ( const Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false ) {

    }

    virtual size_t getLastNodeId() const = 0;
    virtual size_t size() const = 0 ;

    virtual void deleteCache ( int cacheID ) {
        caches[cacheID].clear() ;
        coefs[cacheID].clear() ;
    }

    virtual unsigned int generateCache( Segment * source, double width = 0.001, bool side = false) {
        //search for first empty cache slot ;
        getElements() ;
        if ( caches.empty() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }
        size_t position = 0;
        for ( ; position < caches.size() ; position++ ) {
            if ( caches[position].empty() ) {
                break ;
            }
        }
        if ( position == caches.size() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }

        std::vector<ETYPE *> elems = getConflictingElements ( source, width ) ;
        if(elems.empty())
            elems = getConflictingElements ( &source->first() ) ;

        for ( auto & element : elems ) {
            if(source->intersects(dynamic_cast<IntegrableEntity *>(element)))
            {
                bool on = true ;
                if(side)
                {
                     on = false ;
                     size_t count = 0 ;
                     for(size_t i = 0 ; i < element->getBoundingPoints().size()/element->timePlanes() ; i++)
                     {
                         if( source->on(element->getBoundingPoint(i)) )
                             count++ ;
                     }
                     if(count >= 2)
                         on = true ;
                }
                if(on)
                {
                    caches[position].push_back ( element->index ) ;
                    coefs[position].push_back ( std::vector<double>() ) ;
                }
            }
        }
        return position ;

    }

    virtual unsigned int generateCacheOut( std::vector<unsigned int> ids )
    {
        getElements() ;
        if ( caches.empty() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }
        size_t position = 0;
        for ( ; position < caches.size() ; position++ ) {
            if ( caches[position].empty() ) {
                break ;
            }
        }
        if ( position == caches.size() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }

        if(allElementsCacheID == -1)
            return position ;

        std::valarray<bool> out ;
        std::valarray<bool> check ;
        int maxIndex = 0 ;
        for(size_t i = 0 ; i < caches[allElementsCacheID].size() ; i++)
        {
            if(caches[allElementsCacheID][i] > maxIndex)
                maxIndex = caches[allElementsCacheID][i] ;
        }
        out.resize(maxIndex+1) ;
        check.resize(maxIndex+1) ;
        out = false ;
        check = false ;

        for(size_t i = 0 ; i < caches[allElementsCacheID].size() ; i++)
            check[ caches[allElementsCacheID][i] ] = true ;

        for(size_t j = 0 ; j < ids.size() ; j++)
        {
            for(size_t i = 0 ; i < caches[ids[j]].size() ; i++)
                out[ caches[ids[j]][i] ] = true ;
        }

        for(size_t i = 0 ; i < check.size() ; i++)
        {
            if(check[i] && !out[i])
            {
                caches[position].push_back( i ) ;
                coefs[position].push_back ( std::vector<double>() ) ;
            }
        }
        return position ;
    }

    virtual unsigned int generateCache( std::vector<Geometry *> source) {
        //search for first empty cache slot ;
        getElements() ;
        if ( caches.empty() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }
        size_t position = 0;
        for ( ; position < caches.size() ; position++ ) {
            if ( caches[position].empty() ) {
                break ;
            }
        }
        if ( position == caches.size() ) {
            caches.push_back ( std::vector<int>() );
            coefs.push_back ( std::vector<std::vector<double>>() );
        }

        for(size_t i = 0 ; i < source.size() ; i++)
        {
            std::cerr << "\rgenerating cache for feature " << i << "/" << source.size() << std::flush ;
            std::vector<ETYPE *> elems = getConflictingElements ( source[i] ) ;
            if(elems.empty())
                elems = getConflictingElements ( &source[i]->getCenter() ) ;

            for ( auto & element : elems ) {
                if(source[i]->in(element->getCenter()) && element->getBehaviour() && element->getBehaviour()->getSource() == source[i])
                {
                    caches[position].push_back ( element->index ) ;
                    coefs[position].push_back ( std::vector<double>() ) ;
                }
            }
        }
        std::cerr << "\rgenerating cache for feature " << source.size() << "/" << source.size() << " ... done" << std::endl ;

        return position ;

    }

    virtual unsigned int generateCache ( const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function ( "1" ) ) {
        size_t position = 0;
        getElements() ;
        #pragma omp critical
        {
            VirtualMachine vm ;
            std::vector<double> co ;
            std::vector<ETYPE *> elems = getConflictingElements ( locus ) ;
            if(elems.empty())
                elems = getConflictingElements ( &locus->getCenter() ) ;
            //search for first empty cache slot ;
            if ( caches.empty() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }
            else
            {
                for ( ; position < caches.size() ; position++ ) {
                    if ( caches[position].empty() ) {
                        break ;
                    }
                }
            }
            if ( position == caches.size() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }
            for ( auto & element : elems ) {
                if ( (source && element->getBehaviour() && element->getBehaviour()->getSource() != source) || (element->getBehaviour() && element->getBehaviour()->type == VOID_BEHAVIOUR)) {
                    continue ;
                }
                if ( locus->in ( element->getCenter() ) ) {
                    caches[position].push_back ( element->index ) ;
                    coefs[position].push_back ( std::vector<double>() ) ;
                    if(element->getOrder() >= CONSTANT_TIME_LINEAR)
                    {
		      double a = (spaceDimensions == SPACE_TWO_DIMENSIONAL) ?element->area() : element->volume() ;
                        if(source)
                        {
                            Function x = element->getXTransformAtCentralNodalTime() ;
                            Function y = element->getYTransformAtCentralNodalTime() ;
                            Function z = element->getZTransformAtCentralNodalTime() ;
    //		            Function t = element->getTTransform() ;
			    
                            GaussPointArray gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray( element, 0. ) ;
                            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                                double xx = vm.eval ( x, gp.gaussPoints[i].first ) ;
                                double xy = vm.eval ( y, gp.gaussPoints[i].first ) ;
                                double xz = vm.eval ( z, gp.gaussPoints[i].first ) ;

                                coefs[position].back().push_back ( vm.eval ( smoothing, xx, xy, xz, 0. )*a/gp.gaussPoints.size() );
                            }
                        }
                        else
                        {
                            GaussPointArray gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray( element, 0. ) ;
                            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                                coefs[position].back().push_back ( a/gp.gaussPoints.size() );
                            }
                        }
                    }
                    else
                    {
		      double a = (element->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ?element->area() : element->volume() ;
                        if(source)
                        {
                            Function x = element->getXTransform() ;
                            Function y = element->getYTransform() ;
                            Function z = element->getZTransform() ;
                            GaussPointArray gp = element->getGaussPoints() ;
			    

                            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                                double xx = vm.eval ( x, gp.gaussPoints[i].first ) ;
                                double xy = vm.eval ( y, gp.gaussPoints[i].first ) ;
                                double xz = vm.eval ( z, gp.gaussPoints[i].first ) ;

                                coefs[position].back().push_back ( vm.eval ( smoothing, xx, xy, xz, 0. )*a/gp.gaussPoints.size() );
                            }
                            
                            if(gp.gaussPoints.size() == 1)
			    {
			      coefs[position].back().back() = 0 ;
			      for(size_t i = 0 ; i <  element->getBoundingPoints().size() ; i++)
				coefs[position].back().back() +=  vm.eval ( smoothing, element->getBoundingPoint(i).x,  element->getBoundingPoint(i).y,  element->getBoundingPoint(i).z, 0. )*a ;
			      coefs[position].back().back() /= element->getBoundingPoints().size() ;
			    }
                        }
                        else
                        {
                           GaussPointArray gp = element->getGaussPoints() ;

                            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                                coefs[position].back().push_back ( a/gp.gaussPoints.size() );
                            }
                        }
                    }
                }
            }
        }
        double tot = 0 ;
	for(size_t i = 0 ;  i < coefs[position].size() ; i++)
	  for(size_t j = 0 ;  j < coefs[position][i].size() ; j++)
	    tot += coefs[position][i][j] ; 
	for(size_t i = 0 ;  i < coefs[position].size() ; i++)
	  for(size_t j = 0 ;  j < coefs[position][i].size() ; j++) 
	    coefs[position][i][j] /= tot ;
        return position ;
    } ;

    virtual void updateCache ( size_t position, const Geometry * source = nullptr, Function smoothing = Function ( "1" ) ) {

    #pragma omp critical
    {
        VirtualMachine vm ;

        size_t iter = 0 ;
        for ( auto element = begin(position) ; element != end(position) ; element++ ) {


            coefs[position][iter].clear() ;
            if ( (source && element->getBehaviour() && element->getBehaviour()->getSource() != source) || (element->getBehaviour() && element->getBehaviour()->type == VOID_BEHAVIOUR)) {
                continue ;
            }
            if(element->getOrder() >= CONSTANT_TIME_LINEAR)
            {

                Function x = element->getXTransformAtCentralNodalTime() ;
                Function y = element->getYTransformAtCentralNodalTime() ;
                Function z = element->getZTransformAtCentralNodalTime() ;
//                  Function t = element->getTTransform() ;
                GaussPointArray gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray( element, 0. ) ;
		double a = (spaceDimensions == SPACE_TWO_DIMENSIONAL) ?element->area() : element->volume() ;
                for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                    double xx = vm.eval ( x, gp.gaussPoints[i].first ) ;
                    double xy = vm.eval ( y, gp.gaussPoints[i].first ) ;
                    double xz = vm.eval ( z, gp.gaussPoints[i].first ) ;

                    coefs[position][iter].push_back ( vm.eval ( smoothing, xx, xy, xz, 0. )*a/gp.gaussPoints.size() );
                }
            }
            else
            {
                Function x = element->getXTransform() ;
                Function y = element->getYTransform() ;
                Function z = element->getZTransform() ;
                GaussPointArray gp = element->getGaussPoints() ;
		double a = (element->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ?element->area() : element->volume() ;

                for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
                    double xx = vm.eval ( x, gp.gaussPoints[i].first ) ;
                    double xy = vm.eval ( y, gp.gaussPoints[i].first ) ;
                    double xz = vm.eval ( z, gp.gaussPoints[i].first ) ;

                    if(element->getBehaviour()->fractured())
                        coefs[position][iter].push_back ( vm.eval ( smoothing, xx, xy, xz, 0. )*a/gp.gaussPoints.size() );
                    else
                        coefs[position][iter].push_back ( 0. );
                }
                
		if(gp.gaussPoints.size() == 1)
		{
		  coefs[position][iter].back() = 0 ;
		  for(size_t i = 0 ; i <  element->getBoundingPoints().size() ; i++)
		    coefs[position][iter].back() +=  vm.eval ( smoothing, element->getBoundingPoint(i).x,  element->getBoundingPoint(i).y,  element->getBoundingPoint(i).z, 0. )*a ;
		  coefs[position][iter].back() /= element->getBoundingPoints().size() ;
		}

            }
            iter++ ;
        }
        
        double tot = 0 ;
        for(size_t i = 0 ;  i < coefs[position].size() ; i++)
	  for(size_t j = 0 ;  j < coefs[position][i].size() ; j++)
	    tot += coefs[position][i][j] ; 
	for(size_t i = 0 ;  i < coefs[position].size() ; i++)
	  for(size_t j = 0 ;  j < coefs[position][i].size() ; j++) 
	    coefs[position][i][j] /= tot ;
    }

    } ;

    
    virtual unsigned int generateCache ()
    {
//        std::cout << omp_get_max_threads() << std::endl ;
#ifdef HAVE_OPENMP
      int maxthread = std::max(omp_get_max_threads(), 8) ;
#else
      int maxthread = 1 ;
#endif
      
        if(strainCache.size() == 0) { strainCache.resize( maxthread, Vector(0) ) ; }
        if(strainGenViscoCache.size() == 0) { strainGenViscoCache.resize( maxthread, Vector(0) ) ; }
        if(strainRateGenViscoCache.size() == 0) { strainRateGenViscoCache.resize( maxthread, Vector(0) ) ; }
        if(stressCache.size() == 0) { stressCache.resize( maxthread, Vector(0) ) ; }
        if(stressGenViscoCache.size() == 0) { stressGenViscoCache.resize( maxthread, Vector(0) ) ; }
        if(stressImposedCache.size() == 0) { stressImposedCache.resize( maxthread, Vector(0) ) ; }
        if(bufferCache.size() == 0) { bufferCache.resize( maxthread, Vector(0) ) ; }
        if(firstResultCache.size() == 0) { firstResultCache.resize( maxthread, Vector(0) ) ; }
        if(secondResultCache.size() == 0) { secondResultCache.resize( maxthread, Vector(0) ) ; }

        getElements() ;
        #pragma omp critical
        {
            //search for first empty cache slot ;
            if ( caches.empty() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }
            size_t position = 0;
            for ( ; position < caches.size() ; position++ ) {
                if ( caches[position].empty() ) {
                    break ;
                }
            }
            if ( position == caches.size() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }

            for ( size_t i = 0 ;  i < size() ; i++ ) {
                ETYPE * elem = dynamic_cast<ETYPE *>(getInTree(i)) ;

                if(!elem)
                    continue ;
                caches[position].push_back ( getInTree(i)->index ) ;
                coefs[position].push_back ( std::vector<double>() ) ;

//             for ( size_t j = 0 ; j < elem->getGaussPoints().gaussPoints.size() ; j++ )
//                 coefs[position].back().push_back ( 1 ) ;
            }
            allElementsCacheID = position ;
        }
        return allElementsCacheID ;
    };

    virtual unsigned int generateCache (ElementChecker * c)
    {
        getElements() ;
        size_t position = 0;
        #pragma omp critical
        {

            //search for first empty cache slot ;
            if ( caches.empty() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }

            for ( ; position < caches.size() ; position++ ) {
                if ( caches[position].empty() ) {
                    break ;
                }
            }
            if ( position == caches.size() ) {
                caches.push_back ( std::vector<int>() );
                coefs.push_back ( std::vector<std::vector<double>>() );
            }

            for ( size_t i = 0 ;  i < size() ; i++ ) {
                ETYPE * elem = dynamic_cast<ETYPE *>(getInTree(i)) ;

                if(!elem)
                    continue ;
                if(c->checkElement(elem))
                {
                    caches[position].push_back ( getInTree(i)->index ) ;
                    coefs[position].push_back ( std::vector<double>() ) ;
                }

            }
        }
        return position ;
    };

    virtual void clearCaches()
    {
        caches.clear();
        coefs.clear();
        allElementsCacheID = -1 ;
    }

    virtual double getField ( std::string f, int cacheID, std::string correction = std::string("correction_factor") ) {
        double ret = 0. ;
        unsigned int realID = cacheID ;
        if(cacheID == -1)
            realID = allElementsCacheID ;
        if( caches[realID].size() == 0)
            return 0. ;
        double w = 0. ;
        std::map<std::string, double> dummy ;
        dummy[correction] = 1 ;
        bool is2d = (static_cast<ETYPE *> ( getInTree ( caches[realID][0] ) )->spaceDimensions() == SPACE_TWO_DIMENSIONAL) ;
        for ( size_t i = 0 ; i < caches[realID].size() ; i++ ) {
            if( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * >(&static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->getState()))
            {
                if(dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * >(&static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->getState())->has(f))
                {
                    double a = ((is2d) ? static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->area() : static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->volume())*dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & >(static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->getState()).get(correction, dummy) ;
                    ret += dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & >(static_cast<ETYPE *> ( getInTree ( caches[realID][i] ) )->getState()).get(f, dummy)*a ;
                    w += a ;
                }
            }
        }
        if(w > 0)
            return ret/w ;
        return 0. ;
    }

    //virtual void getAverageField( Amie::FieldType f, Vector& ret, Amie::VirtualMachine* vm = nullptr, int dummy = 0, double t = 0, std::vector< double > weights = std::vector<double>()) ;
    virtual Vector getField ( FieldType f, int cacheID, double t, int index = 0 ) {
        if(cacheID == -1 && allElementsCacheID == -1)
        {
            VirtualMachine vm ;
            size_t blocks = 0 ;

            std::vector<ETYPE *> elems = getElements() ;
            std::vector<double> coefs ;
            for ( size_t i = 0 ; i < elems.size() && !blocks; i++ ) {
                if(elems[i]->getBehaviour()->type == VOID_BEHAVIOUR)
                    continue ;
                
                blocks = elems[i]->getBehaviour()->getNumberOfDegreesOfFreedom() /spaceDimensions ;
            }
            Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
            Vector buffer ( ret ) ;
	    Vector retcomp( ret ) ;
	    Vector temp( ret ) ;
	    Vector tempsum( ret ) ;
            double w = 0 ;
            for ( size_t i = 0 ; i < elems.size() ; i++ ) {

                double v = elems[i]->getState().getAverageField ( f, buffer, &vm, t, coefs, index ) ;
		temp = buffer * v - retcomp ;
		tempsum = ret+temp ;
		retcomp = tempsum-ret-temp ;
		ret = tempsum ;
		w +=v ;
		
		if(std::isnan(v) || std::isnan(buffer[0]) || std::isnan(buffer[1]))
		{
		  std::cout << buffer[0] << "  " << buffer[1] << std::endl ;
		  exit(0) ;
		}
            }
            if(w > POINT_TOLERANCE*POINT_TOLERANCE)
	      ret /=w ;
	    
	    return ret ;
        }
        if(cacheID == -1)
            cacheID = allElementsCacheID ;
        VirtualMachine vm ;
        size_t blocks = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
            blocks = static_cast<ETYPE *>(getInTree ( caches[cacheID][i] ))->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions ;
        }
        Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
        if( caches[cacheID].size() == 0)
            return ret ;
        Vector buffer ( ret ) ;
	Vector retcomp( ret ) ;
	Vector temp( ret ) ;
	Vector tempsum( ret ) ;
        double w = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
            double v = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->getState().getAverageField ( f, buffer, &vm, t, coefs[cacheID][i], index ) ;
	    temp = buffer * v - retcomp ;
	    tempsum = ret+temp ;
	    retcomp = tempsum-ret-temp ;
            ret = tempsum ;
            w +=v ;
        }
        if(w > POINT_TOLERANCE*POINT_TOLERANCE)
	  ret /=w ;
	
	return ret ;
    }
    
    virtual Vector getField ( FieldType f, const Point & p, Vector & ret, double t = 1, int index = 0 ) 
    {
        ETYPE * e = getUniqueConflictingElement(&p) ;        
        Point p_ = p ;
        p_.t = t ;
        if(e && e->state)
            e->getState().getField(f, p_, ret, false, nullptr, index) ;
        else
            ret = 0 ;
        return ret ;
    }

    virtual Vector getField( FieldType f, Segment * source, int cacheID , double t = 0, int index = 0) {
        int realCacheID = cacheID ;
        if(cacheID == -1)
            realCacheID = allElementsCacheID ;

        VirtualMachine vm ;
        size_t blocks = 0 ;
        for ( size_t i = 0 ; i < caches[realCacheID].size() && !blocks; i++ ) {
            blocks = static_cast<ETYPE *>(getInTree ( caches[realCacheID][i] ))->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions ;
        }
        Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
        Vector buffer ( ret ) ;
        double w = 0 ;
        if( caches[realCacheID].size() == 0)
           return ret ;

        double realT = static_cast<ETYPE *>( getInTree ( caches[realCacheID][0] ) )->getBoundingPoint(0).getT()*(1.-t)/2. + static_cast<ETYPE *>( getInTree ( caches[realCacheID][0] ) )->getBoundingPoint( static_cast<ETYPE *>(getInTree ( caches[realCacheID][0] ))->getBoundingPoints().size()-1 ).getT()*(1.+t)/2. ;

        for ( size_t i = 0 ; i < caches[realCacheID].size() ; i++ ) {
            std::vector<Point> inter = source->intersection( dynamic_cast<IntegrableEntity *>( getInTree ( caches[realCacheID][i] ) ) ) ;
            if(source->intersects( dynamic_cast<IntegrableEntity *>( getInTree ( caches[realCacheID][i] ) ) ) && inter.size() == 2)
            {
                Segment local(inter[0], inter[1]) ;
                Point mid = local.midPoint() ;
                mid.setT( realT ) ;
                static_cast<ETYPE *> ( getInTree ( caches[realCacheID][i] ) )->getState().getField ( f, mid, buffer, false, nullptr, index) ;
                ret += buffer * local.norm() ;
                w += local.norm() ;
            }
        }
        return ret/w ;
    }

    virtual double getField( std::string f, Segment * source, int cacheID = -1, double t = 0) {
        int realCacheID = cacheID ;
        if(cacheID == -1)
            realCacheID = allElementsCacheID ;

        VirtualMachine vm ;
        double ret = 0. ;
        double w = 0 ;
        if( caches[realCacheID].size() == 0)
           return ret ;

        std::map<std::string, double> emptymap ;

        for ( size_t i = 0 ; i < caches[realCacheID].size() ; i++ ) {
            std::vector<Point> inter = source->intersection( dynamic_cast<IntegrableEntity *>( getInTree ( caches[realCacheID][i] ) ) ) ;
            if(source->intersects( dynamic_cast<IntegrableEntity *>( getInTree ( caches[realCacheID][i] ) ) ) && inter.size() == 2)
            {
                Segment local(inter[0], inter[1]) ;
                if( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * >(&static_cast<ETYPE *> ( getInTree ( caches[realCacheID][i] ) )->getState()))
                {
                    if(dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * >(&static_cast<ETYPE *> ( getInTree ( caches[realCacheID][i] ) )->getState())->has(f))
                    {
                        ret += dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & >(static_cast<ETYPE *> ( getInTree ( caches[realCacheID][i] ) )->getState()).get(f, emptymap) * local.norm() ;
                        w += local.norm() ;
                    }
                }
            }
        }
        if(w > POINT_TOLERANCE*POINT_TOLERANCE)
            return ret/w ;
        return 0. ;
    }

    virtual double getArea( int cacheID)
    {
        double a = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ )
            a += static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->area() ;
        return a ;
    }

    virtual double getVolume( int cacheID)
    {
        double v = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ )
            v += static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->volume() ;
        return v ;
    }

    virtual Vector getSmoothedField (  FieldType f0, int cacheID, IntegrableEntity * e, double t, const std::vector<bool> & restrict = std::vector<bool>(), int index = 0) {
        unsigned int tsize = 3 ;
        unsigned int psize = 2 ;
#ifdef HAVE_OPENMP
        int thread = omp_get_thread_num() ;
#else
	int thread = 0 ;
#endif

        if ( spaceDimensions == SPACE_THREE_DIMENSIONAL ) {
            tsize = 6 ;
            psize = 3 ;

        }

        firstResultCache[thread] = 0 ;
        strainCache[thread] = 0 ;
        stressCache[thread] = 0 ;
        bufferCache[thread] = 0 ;
//        Vector first ;
//        Vector strain ;
//        Vector stress ;
//        Vector buffer ;

        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        if(spaceTime) {
            strainGenViscoCache[thread] = 0 ;
            strainRateGenViscoCache[thread] = 0 ;
            stressGenViscoCache[thread] = 0 ;
            stressImposedCache[thread] = 0 ;
        }
        VirtualMachine vm ;
        if (
            f0 == REAL_STRESS_FIELD || 
            f0 == EFFECTIVE_STRESS_FIELD || 
            f0 == PRINCIPAL_REAL_STRESS_FIELD || 
            f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
            //we first need to compute the strain field 
            if(strainCache[thread].size() != tsize ) { strainCache[thread].resize(tsize, 0.); }
            if(bufferCache[thread].size() != tsize ) { bufferCache[thread].resize(tsize, 0.); }
            if ( !spaceTime ) {
                double sumFactors ( 0 ) ;
               
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    
                    if(ci->getBehaviour()->fractured() || ci->getBehaviour()->getSource() != e->getBehaviour()->getSource())
                        continue ;
//                     std::cout << "x" << std::endl ;
                    double v = 0; 

                    if(e != ci  || restrict.empty()|| ( e == ci && (restrict.size() !=  coefs[cacheID][i].size())))
                        v = ci->getState().getAverageField ( MECHANICAL_STRAIN_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                    else if(e == ci && restrict.size() ==  coefs[cacheID][i].size())
                    {
                        std::vector<double> effCoef = coefs[cacheID][i] ;

                        for(size_t j = 0 ; j < restrict.size() ; j++)
                            effCoef[j] *= !restrict[j] ;
                        
                        v = ci->getState().getAverageField ( MECHANICAL_STRAIN_FIELD, bufferCache[thread], &vm, t, effCoef ); 
                    }                   
                    strainCache[thread] += bufferCache[thread]*v ;
                    sumFactors += v ;
// 		    std::cout << "y" << std::endl ;

                }
                if(sumFactors > POINT_TOLERANCE*POINT_TOLERANCE)
                    strainCache[thread] /= sumFactors ;
            } else {
                size_t blocks = 0 ;
                for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
                    ETYPE *ci  = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions  ;
                }

//                Vector tmpstrain ;
//                Vector tmpstrainrate ;

                size_t bsize = fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ) ;
                if( bufferCache[thread].size() != bsize ) { bufferCache[thread].resize( bsize, 0. ) ; }
                if( strainGenViscoCache[thread].size() != bsize ) { strainGenViscoCache[thread].resize( bsize, 0. ) ; }
                if( strainRateGenViscoCache[thread].size() != bsize ) { strainRateGenViscoCache[thread].resize( bsize, 0. ) ; }
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                    if(ci->getBehaviour()->fractured() || ci->getBehaviour()->getSource() != e->getBehaviour()->getSource())
                        continue ;

                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                    strainGenViscoCache[thread] += bufferCache[thread]*v ;
                    sumFactors += v ;
                }
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                    strainRateGenViscoCache[thread] += bufferCache[thread]*v ;
                }
                strainGenViscoCache[thread] /= sumFactors ;
                strainRateGenViscoCache[thread] /=sumFactors ;

                double sum = 0 ; 
                if(strainCache[thread].size() != tsize ) { strainCache[thread].resize(tsize, 0.); }
                if(stressCache[thread].size() != tsize ) { stressCache[thread].resize(tsize, 0.); }
                if(stressImposedCache[thread].size() != tsize ) { stressImposedCache[thread].resize(tsize, 0.); }
                Viscoelasticity * v = dynamic_cast<Viscoelasticity *>(e->getBehaviour()) ;
                for ( unsigned int i = 0 ; i < tsize ; i++ ) {
                     strainCache[thread][i] = strainGenViscoCache[thread][i] ;
                     if( v->model >= MAXWELL && (f0 == MECHANICAL_STRAIN_FIELD || f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ) )
                          {
                              for(size_t j = 1 ; j < strainGenViscoCache[thread].size()/strainCache[thread].size() ; j++)
                                  strainCache[thread][i] -= strainGenViscoCache[thread][ j*tsize + i ] ;
                          } 
                     }
                if(f0 != PRINCIPAL_MECHANICAL_STRAIN_FIELD) { 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty()&& restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        stressGenViscoCache[thread] = strainGenViscoCache[thread]*e->getBehaviour()->getTensor ( p )  + ( Vector ) ( strainRateGenViscoCache[thread]*e->getBehaviour()->getViscousTensor ( p ) ) ;
                        stressImposedCache[thread] = e->getBehaviour()->getImposedStress( p ) ;
                        for ( unsigned int i = 0 ; i < tsize ; i++ ) {
                            stressCache[thread][i] += (stressGenViscoCache[thread][i]-stressImposedCache[thread][i])*e->getGaussPoints().gaussPoints[j].second ;
                        }
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    stressCache[thread] /= sum ;
               }
            }

            if ( f0 == REAL_STRESS_FIELD ) {
                firstResultCache[thread].resize ( tsize, 0. );
                if ( !spaceTime ) {
                    double sum = 0 ; 
                    Point p ;
// 		    std::cout << "z" << std::endl ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty() && restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        firstResultCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p ))*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    if(sum > POINT_TOLERANCE*POINT_TOLERANCE)
                        firstResultCache[thread] /= sum ;
		    
// 		    std::cout << "w" << std::endl ;
                    
                } else {
                    firstResultCache[thread] = stressCache[thread] ;
                }
            }
            else if ( f0 == EFFECTIVE_STRESS_FIELD ) {
                firstResultCache[thread].resize ( tsize );
                if ( !spaceTime ) {
                    firstResultCache[thread] = (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param ;
                } else {
                    firstResultCache[thread] = stressCache[thread] ;
                }
            }
            else if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                firstResultCache[thread].resize ( psize );
                if ( !spaceTime ) {
                    firstResultCache[thread] = toPrincipal( (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param , SINGLE_OFF_DIAGONAL_VALUES) ;
                } else {
                    firstResultCache[thread] = toPrincipal( stressCache[thread] , SINGLE_OFF_DIAGONAL_VALUES) ;
                }
            }
            else if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
                firstResultCache[thread].resize ( psize, 0. );
                if ( !spaceTime ) {
                    if(stressCache[thread].size() != tsize ) { stressCache[thread].resize(tsize, 0.); }
                    double sum = 0 ; 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty() && restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        stressCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p ))*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    if(sum > POINT_TOLERANCE*POINT_TOLERANCE)
                        stressCache[thread] /= sum ;
                    firstResultCache[thread] = toPrincipal( stressCache[thread] , SINGLE_OFF_DIAGONAL_VALUES) ;
                } else {
                    firstResultCache[thread] = toPrincipal( stressCache[thread] , SINGLE_OFF_DIAGONAL_VALUES) ;
                }
            }

        } else {
            double sumFactors ( 0 ) ;
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                
                double v = 0 ;
                               
                if(e != ci  || restrict.empty() || ( e == ci && (restrict.size() !=  coefs[cacheID][i].size())))
                    v = ci->getState().getAverageField ( f0, bufferCache[thread], &vm, t, coefs[cacheID][i], index );
                else if(e == ci && restrict.size() ==  coefs[cacheID][i].size())
                {
                    std::vector<double> effCoef = coefs[cacheID][i] ;

                    for(size_t j = 0 ; j < restrict.size() ; j++)
                        effCoef[j] *= !restrict[j] ;
                    
                    v = ci->getState().getAverageField ( f0, bufferCache[thread], &vm, t, effCoef, index ); 
                }
                    
                if ( firstResultCache[thread].size() != bufferCache[thread].size() ) {
                    firstResultCache[thread].resize ( bufferCache[thread].size(), 0. );
                }
                firstResultCache[thread] += bufferCache[thread]*v ;
                sumFactors += v ;
                
            }
            if(sumFactors > POINT_TOLERANCE*POINT_TOLERANCE)
                firstResultCache[thread] /= sumFactors ;
        }

        return firstResultCache[thread] ;
    }

    virtual std::pair<Vector, Vector> getSmoothedFields ( FieldType f0, FieldType f1, int cacheID, IntegrableEntity * e , double t, const std::vector<bool> & restrict = std::vector<bool>(), int index0 = 0, int index1 = 1  ) {

#ifdef HAVE_OPENMP
        int thread = omp_get_thread_num() ;
#else
	int thread = 0 ;
#endif

        firstResultCache[thread] = 0 ;
        secondResultCache[thread] = 0 ;
        strainCache[thread] = 0 ;
        stressCache[thread] = 0 ;
        bufferCache[thread] = 0 ;
//        Vector first ;
//        Vector strain ;
//        Vector stress ;
//        Vector buffer ;

        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        if(spaceTime) {
            strainGenViscoCache[thread] = 0 ;
            strainRateGenViscoCache[thread] = 0 ;
            stressGenViscoCache[thread] = 0 ;
            stressImposedCache[thread] = 0 ;
        }
//        Vector first ;
//        Vector second ;
//        Vector strain ;
//        Vector stress ;
//        Vector buffer ;
        unsigned int tsize = 3 ;
        unsigned int psize = 2 ;
        if ( spaceDimensions == SPACE_THREE_DIMENSIONAL ) {
            tsize = 6 ;
            psize = 3 ;
        }
        VirtualMachine vm ;
        if (
             f0 == REAL_STRESS_FIELD || 
             f0 == EFFECTIVE_STRESS_FIELD || 
             f0 == PRINCIPAL_REAL_STRESS_FIELD || 
             f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || 
             f0 == MECHANICAL_STRAIN_FIELD || 
             f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ||
               
                f1 == REAL_STRESS_FIELD || 
                f1 == EFFECTIVE_STRESS_FIELD || 
                f1 == PRINCIPAL_REAL_STRESS_FIELD || 
                f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || 
                f1 == MECHANICAL_STRAIN_FIELD || 
                f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD
           ) {
            //we first need to compute the strain field
            if ( !spaceTime ) {
                bufferCache[thread].resize ( tsize, 0. );
                if(strainCache[thread].size() != tsize ) { strainCache[thread].resize ( tsize, 0. ); }
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    if(ci->getBehaviour()->fractured())
                        continue ;
                    if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                       
                        double v = 0; 

                        if(e != ci  || restrict.empty()|| ( e == ci && (restrict.size() !=  coefs[cacheID][i].size())))
                            v = ci->getState().getAverageField ( MECHANICAL_STRAIN_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                        else if(e == ci && restrict.size() ==  coefs[cacheID][i].size())
                        {
                            std::vector<double> effCoef = coefs[cacheID][i] ;

                            for(size_t j = 0 ; j < restrict.size() ; j++)
                                effCoef[j] *= !restrict[j] ;
                            
                            v = ci->getState().getAverageField ( MECHANICAL_STRAIN_FIELD, bufferCache[thread], &vm, t, effCoef ); 
                        } 

                        strainCache[thread] += bufferCache[thread]*v ;
                        sumFactors += v ;
                    }
                }
                if(sumFactors > POINT_TOLERANCE*POINT_TOLERANCE)
                    strainCache[thread] /= sumFactors ;
            } else {
                size_t blocks = 0 ;
                for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
                    ETYPE *ci  = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions  ;
                }

                size_t bsize = fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ) ;
                if( strainGenViscoCache[thread].size() != bsize ) { strainGenViscoCache[thread].resize ( bsize, 0. ) ; }
                if( bufferCache[thread].size() != bsize ) { bufferCache[thread].resize ( bsize, 0. ) ; }
                if( strainRateGenViscoCache[thread].size() != bsize ) { strainRateGenViscoCache[thread].resize ( bsize, 0. ) ; }
                double sumFactors ( 0 ) ;


                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                    if(ci->getBehaviour()->fractured() || ci->getBehaviour()->getSource() != e->getBehaviour()->getSource())
                        continue ;

                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                    strainGenViscoCache[thread] += bufferCache[thread]*v ;
                    sumFactors += v ;
                }
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, bufferCache[thread], &vm, t, coefs[cacheID][i] );
                    strainRateGenViscoCache[thread] += bufferCache[thread]*v ;
                }
                strainGenViscoCache[thread] /= sumFactors ;
                strainRateGenViscoCache[thread] /=sumFactors ;

                double sum = 0 ; 
                if(strainCache[thread].size() != tsize ) { strainCache[thread].resize ( tsize, 0. ); }
                if(stressCache[thread].size() != tsize ) { stressCache[thread].resize ( tsize, 0. ); }
                Viscoelasticity * v = dynamic_cast<Viscoelasticity *>(e->getBehaviour()) ;
                for ( unsigned int i = 0 ; i < tsize ; i++ ) {
                     strainCache[thread][i] = strainGenViscoCache[thread][i] ;
                     if( v->model >= MAXWELL && (f0 == MECHANICAL_STRAIN_FIELD || f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD || f1 == MECHANICAL_STRAIN_FIELD || f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD) )
                          {
                              for(size_t j = 1 ; j < strainGenViscoCache[thread].size()/strainCache[thread].size() ; j++)
                                  strainCache[thread][i] -= strainGenViscoCache[thread][ j*tsize + i ] ;
                          } 
                     }
                Point p ;
                for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                {
                    if(!restrict.empty()&& restrict.size() == e->getGaussPoints().gaussPoints.size())
                        if(restrict[j])
                            continue ;
                    p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                    Vector tmpstress = strainGenViscoCache[thread]*e->getBehaviour()->getTensor ( p ) + ( Vector ) ( strainGenViscoCache[thread]*e->getBehaviour()->getViscousTensor ( p ) ) ;
                    Vector imposed = e->getBehaviour()->getImposedStress( p ) ;
                    for ( unsigned int i = 0 ; i < tsize ; i++ ) {
                        stressCache[thread][i] += (tmpstress[i]-imposed[i])*e->getGaussPoints().gaussPoints[j].second ;
                    }
                    sum += e->getGaussPoints().gaussPoints[j].second ;
                }
                stressCache[thread] /= sum ;
            }


            if ( f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ) {
                firstResultCache[thread].resize ( psize );
                firstResultCache[thread] = toPrincipal( strainCache[thread], DOUBLE_OFF_DIAGONAL_VALUES ) ;
            }
            if ( f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ) {
                secondResultCache[thread].resize ( psize );
                secondResultCache[thread] = toPrincipal ( strainCache[thread] , DOUBLE_OFF_DIAGONAL_VALUES ) ;
            }
            if ( f0 == REAL_STRESS_FIELD ) {
                firstResultCache[thread].resize ( tsize, 0. );
                if ( !spaceTime ) {
                    double sum = 0 ; 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty()&& restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        firstResultCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p ))*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    if(sum > POINT_TOLERANCE*POINT_TOLERANCE)
                        firstResultCache[thread] /= sum ;
                } else {
                    firstResultCache[thread] = stressCache[thread] ;
                }
            }
            if ( f1 == REAL_STRESS_FIELD ) {
                secondResultCache[thread].resize ( tsize );
                if ( !spaceTime ) {
                    double sum = 0 ; 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty()&& restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        secondResultCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p ))*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    secondResultCache[thread] /= sum ;
                } else {
                    secondResultCache[thread] = stressCache[thread] ;
                }
            }
            if ( f0 == EFFECTIVE_STRESS_FIELD ) {
                firstResultCache[thread].resize ( tsize );
                if ( !spaceTime ) {
                    firstResultCache[thread] = (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param ;
                } else {
                    firstResultCache[thread] = stressCache[thread] ;
                }
            }
            if ( f1 == EFFECTIVE_STRESS_FIELD ) {
                secondResultCache[thread].resize ( tsize );
                if ( !spaceTime ) {
                    secondResultCache[thread] = (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param ;
                } else {
                    secondResultCache[thread] = stressCache[thread] ;
                }
            }
            if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                firstResultCache[thread].resize ( psize );
                if ( !spaceTime ) {
                    firstResultCache[thread] = toPrincipal ( (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param, SINGLE_OFF_DIAGONAL_VALUES  ) ;
                } else {
                    firstResultCache[thread] = toPrincipal( stressCache[thread], SINGLE_OFF_DIAGONAL_VALUES  ) ;
                }
            }
            if ( f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                secondResultCache[thread].resize ( psize );
                if ( !spaceTime ) {
                    secondResultCache[thread] = toPrincipal( (strainCache[thread]-e->getBehaviour()->getImposedStrain( Point() ))*e->getBehaviour()->param , SINGLE_OFF_DIAGONAL_VALUES ) ;
                } else {
                    secondResultCache[thread] = toPrincipal( stressCache[thread], SINGLE_OFF_DIAGONAL_VALUES  ) ;
                }
            }
            if ( f0 == PRINCIPAL_REAL_STRESS_FIELD ) {
                firstResultCache[thread].resize ( psize );
                if ( !spaceTime ) {
                    stressCache[thread].resize ( tsize, 0. );
                    double sum = 0 ; 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty() && restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        stressCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p) )*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    if(sum > POINT_TOLERANCE*POINT_TOLERANCE)
                        stressCache[thread] /= sum ;
                    firstResultCache[thread] = toPrincipal ( stressCache[thread] , SINGLE_OFF_DIAGONAL_VALUES  ) ;
                } else {
                    firstResultCache[thread] = toPrincipal ( stressCache[thread] , SINGLE_OFF_DIAGONAL_VALUES  ) ;
                }
            }
            if ( f1 == PRINCIPAL_REAL_STRESS_FIELD ) {
                secondResultCache[thread].resize ( psize );
                if ( !spaceTime ) {
                    stressCache[thread].resize ( tsize, 0. );
                    double sum = 0 ; 
                    Point p ;
                    for(size_t j = 0 ; j < e->getGaussPoints().gaussPoints.size() ; j++)
                    {
                        if(!restrict.empty()&& restrict.size() == e->getGaussPoints().gaussPoints.size())
                            if(restrict[j])
                                continue ;
                        p.set(e->getGaussPoints().gaussPoints[j].first.x,e->getGaussPoints().gaussPoints[j].first.y,e->getGaussPoints().gaussPoints[j].first.z,t) ;
                        stressCache[thread] += ((strainCache[thread]-e->getBehaviour()->getImposedStrain( p))*e->getBehaviour()->getTensor ( p )-e->getBehaviour()->getImposedStress( p ))*e->getGaussPoints().gaussPoints[j].second;
                        sum += e->getGaussPoints().gaussPoints[j].second ;
                    }
                    if(sum > POINT_TOLERANCE*POINT_TOLERANCE)
                        stressCache[thread] /= sum ;
                    secondResultCache[thread] = toPrincipal (  stressCache[thread], SINGLE_OFF_DIAGONAL_VALUES  ) ;
                } else {
                    secondResultCache[thread] = toPrincipal ( stressCache[thread], SINGLE_OFF_DIAGONAL_VALUES  ) ;
                }
            }

            if ( f0 == MECHANICAL_STRAIN_FIELD ) {
                firstResultCache[thread].resize( tsize ) ;
                firstResultCache[thread] = strainCache[thread] ;
            }

            if ( f1 == MECHANICAL_STRAIN_FIELD ) {
                secondResultCache[thread].resize( tsize ) ;
                secondResultCache[thread] = strainCache[thread] ;
            }

        } else {
            double sumFactors ( 0 ) ;

            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

               if(ci->getBehaviour()->fractured())
                        continue ;
                double v = ci->getState().getAverageField ( f0, bufferCache[thread], &vm, t, coefs[cacheID][i], index0 );
                if ( !firstResultCache[thread].size() ) {
                    firstResultCache[thread].resize ( 0., bufferCache[thread].size() );
                }
                firstResultCache[thread] += bufferCache[thread]*v ;
                sumFactors += v ;
            }
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
		
                if(ci->getBehaviour()->fractured())
                        continue ;    
                if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                    double v = ci->getState().getAverageField ( f1, bufferCache[thread], &vm, t,coefs[cacheID][i], index1 );
                    if ( !secondResultCache[thread].size() ) {
                        secondResultCache[thread].resize ( 0., bufferCache[thread].size() );
                    }
                    secondResultCache[thread] += bufferCache[thread]*v ;
                }
            }
            if(sumFactors > POINT_TOLERANCE*POINT_TOLERANCE)
            {
                firstResultCache[thread] /= sumFactors ;
                secondResultCache[thread] /= sumFactors ;
            }
        }


        return std::make_pair ( firstResultCache[thread], secondResultCache[thread] ) ;
    }

    class iterator
    {
    private:
        Mesh<ETYPE, EABSTRACTTYPE> * msh ;
        size_t cacheID ;
        size_t position ;

    public:
        iterator( Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { /*msh->getElements() ;*/ } 
        iterator( Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t position) : msh(msh), position(position)
        {
//             msh->getElements() ;
            if(msh->allElementsCacheID == -1)
                msh->allElementsCacheID = msh->generateCache() ;
            cacheID = msh->allElementsCacheID ;
        } ;

        bool operator ==(const iterator & i) const
        {
            return i.position == position ;
        }
        bool operator !=(const iterator & i) const
        {
            return i.position != position ;
        }

        bool operator <=(const iterator & i)const
        {
            return position <= i.position ;
        }
        bool operator <(const iterator & i) const
        {
            return position < i.position ;
        }
        bool operator >=(const iterator & i)const
        {
            return position >= i.position ;
        }
        bool operator >(const iterator & i) const
        {
            return position > i.position ;
        }

        iterator& operator++() {
            position++ ;
            // actual increment takes place here
            return *this;
        }

        iterator operator++(int) {
            iterator tmp(*this); // copy
            operator++(); // pre-increment
            return tmp;   // return old value
        }

        iterator& operator+=(int i) {
            position +=i ;
            return *this;
        }

        friend iterator operator+(iterator lhs,  int i)
        {
            return lhs += i;
        }

        iterator& operator--() {
            position-- ;
            return *this;
        }

        iterator operator--(int) {
            iterator tmp(*this); // copy
            operator--(); // pre-increment
            return tmp;   // return old value
        }

        iterator& operator-=(int i) {
            position -=i ;
            return *this;
        }

        friend iterator operator-(iterator lhs,  int i)
        {
            return lhs -= i;
        }

        ETYPE * operator-> ( ) {
            return msh->getElement(cacheID, position) ;
        }

        operator ETYPE * () {
            return msh->getElement(cacheID, position) ;
        }

        size_t size() const
        {
            return msh->caches[cacheID].size() ;
        }
        size_t getPosition() const
        {
            return position ;
        }

        size_t getId() const {
            return cacheID ;
        }

    } ;

    virtual const ETYPE * getElement(size_t cacheID, size_t position) const 
    {
        return static_cast<const ETYPE *>(getInTree(caches[cacheID][position])) ;
    };
    
    virtual  ETYPE * getElement(size_t cacheID, size_t position)  
    {
        return static_cast<ETYPE *>(getInTree(caches[cacheID][position])) ;
    };
    
    class const_iterator
    {
    private:
        const Mesh<ETYPE, EABSTRACTTYPE> * msh ;
        size_t cacheID ;
        size_t position ;
    public:
        const_iterator( const Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { } ;

        const_iterator( const Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t position) : msh(msh), position(position)
        {
            cacheID = msh->allElementsCacheID ;
        } ;

        bool operator ==(const const_iterator & i) const
        {
            return i.position == position ;
        }
        bool operator !=(const const_iterator & i) const
        {
            return i.position != position ;
        }
        bool operator <=(const const_iterator & i)const
        {
            return position <= i.position ;
        }
        bool operator <(const const_iterator & i) const
        {
            return position < i.position ;
        }
        bool operator >=(const const_iterator & i)const
        {
            return position >= i.position ;
        }
        bool operator >(const const_iterator & i) const
        {
            return position > i.position ;
        }

        const_iterator& operator++() {
            position++ ;
            // actual increment takes place here
            return *this;
        }

        const_iterator operator++(int) {
            iterator tmp(*this); // copy
            operator++(); // pre-increment
            return tmp;   // return old value
        }

        const_iterator& operator+=(int i) {
            position +=i ;
            return *this;
        }

        friend const_iterator operator+(const_iterator lhs,  int i)
        {
            return lhs += i;
        }

        const_iterator& operator--() {
            position-- ;
            return *this;
        }

        const_iterator operator--(int) {
            iterator tmp(*this); // copy
            operator--(); // pre-increment
            return tmp;   // return old value
        }

        const_iterator& operator-=(int i) {
            position -=i ;
            return *this;
        }

        friend const_iterator operator-(const_iterator lhs,  int i)
        {
            return lhs -= i;
        }

        const ETYPE * operator-> ( ) const {
            return msh->getElement(cacheID, position) ;
        }

        operator const ETYPE * () const {
            return msh->getElement(cacheID, position) ;
        }

        size_t size() const
        {
            return msh->caches[cacheID].size() ;
        }
        size_t getPosition() const
        {
            return position ;
        }

        size_t getId() const {
            return cacheID ;
        }

    } ;

    iterator begin()
    {
        return iterator(this, 0) ;
    }
    
    const_iterator cbegin() const
    {
        return const_iterator(this, 0) ;
    }
    
    iterator end()
    {
        if(allElementsCacheID == -1)
            allElementsCacheID = generateCache() ;
        return iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
    }
    
    const_iterator cend() const
    {

        return const_iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
    }
    
    iterator begin( size_t cacheID)
    {
        return iterator(this, cacheID, 0) ;
    }
    
    const_iterator cbegin(size_t cacheID)
    {
        return const_iterator(this, cacheID, 0) ;
    }
    
    iterator end(size_t cacheID)
    {
        return iterator(this,cacheID, caches[cacheID].size()) ;
    }
    
    const_iterator cend(size_t cacheID)
    {
        return const_iterator(this,cacheID, caches[cacheID].size()) ;
    }

} ;



template<class ETYPE, class EABSTRACTTYPE>
class SingleElementMesh : public Mesh<ETYPE, EABSTRACTTYPE>
{
protected:
    ETYPE * element ;
    std::vector<EABSTRACTTYPE *> tree ;
    std::vector<Point *> points ;
    std::map<int *, int> trans ;
    size_t global_counter ;

    virtual std::vector<ETYPE *> getNeighbourhood ( ETYPE * element ) const {
        std::vector<ETYPE *> ret = { element };
        return ret ;
    };
    virtual std::vector<ETYPE *> getElements() {
        std::vector<ETYPE *> ret ;
        ret.push_back ( element ) ;
        return ret ;
    }

    void addSharedNodes ( size_t nodes_per_side, size_t time_planes, double timestep ) {
        for ( auto  i = tree.begin() ; i != tree.end() ; ++i ) {
            delete element->cachedGps ;
            element->cachedGps = nullptr ;
            size_t nodes_per_plane = nodes_per_side*3+3 ;

            std::valarray<Point *> newPoints ( nodes_per_plane*time_planes ) ;
            std::valarray<bool> done ( false, nodes_per_plane*time_planes ) ;

            for ( size_t plane = 0 ; plane < time_planes ; plane++ ) {
                for ( size_t side = 0 ; side < 3 ; side++ ) {
                    Point a ( static_cast<ETYPE *> ( *i )->getBoundingPoint ( side ) ) ;
                    Point b ( static_cast<ETYPE *> ( *i )->getBoundingPoint ( ( side+1 ) %3 ) ) ;

                    if ( time_planes> 1 ) {
                        a.getT() = ( double ) plane* ( timestep/ ( double ) ( time_planes-1 ) )-timestep/2.;
                        b.getT() = ( double ) plane* ( timestep/ ( double ) ( time_planes-1 ) )-timestep/2.;
                    }
                    for ( size_t node = 0 ; node < nodes_per_side+1 ; node++ ) {
                        double fraction = ( double ) ( node ) / ( ( double ) nodes_per_side+1 ) ;
                        Point proto = a* ( 1.-fraction ) + b*fraction ;

                        for ( size_t j = 0 ; j< static_cast<ETYPE *> ( *i )->getBoundingPoints().size() ; j++ ) {
                            if ( static_cast<ETYPE *> ( *i )->getBoundingPoint ( j ) == proto ) {
                                break ;
                            }
                        }


                        if ( !done[nodes_per_plane*plane+side* ( nodes_per_side+1 ) +node] ) {

                            points.push_back ( new Point ( proto ) ) ;
                            newPoints[nodes_per_plane*plane+side* ( nodes_per_side+1 ) +node]  = points.back();
                            newPoints[nodes_per_plane*plane+side* ( nodes_per_side+1 ) +node]->getId() = global_counter++ ;

                            done[nodes_per_plane*plane+side* ( nodes_per_side+1 ) +node] = true ;
                        }
                    }
                }
            }

            static_cast<ETYPE *> ( *i )->setBoundingPoints ( newPoints ) ;
        }

    }

public:

    virtual size_t size() const {
        return 1. ;
    } ;
    int elementLayer ( ETYPE * e ) const {
        return -1 ;
    }
    void addElementToLayer ( const ETYPE * element, int layer ) { } ;

    SingleElementMesh ( ETYPE * element, SpaceDimensionality spaceDimensions ) : Mesh<ETYPE, EABSTRACTTYPE>(spaceDimensions), element ( element ), global_counter ( 0 ) {
        tree.push_back ( element ) ;
        for ( size_t i = 0 ; i < element->getBoundingPoints().size() ; ++i ) {
            if ( element->getBoundingPoint ( i ).getId() < 0 ) {
                element->getBoundingPoint ( i ).getId() = global_counter ;
            }
            trans[& ( element->getBoundingPoint ( i ).getId() )] = global_counter ;
            global_counter++ ;
        }
    } ;

    virtual ~SingleElementMesh() {
        for ( size_t i = 0 ; i < points.size() ; i++ ) {
            delete points[i] ;
        }
    } ;

    virtual std::vector<ETYPE *> getConflictingElements ( const Point  * p ) {
        std::vector<ETYPE *> ret ;
        if ( element->in ( *p ) ) {
            ret.push_back ( element ) ;
        }
        return ret ;
    }

    virtual void extrude ( double dt ) {
        std::cout << "should extrude.." << std::endl  ;
    } ;
    virtual void extrude ( const Vector & dt ) {
        std::cout << "should extrude.." << std::endl  ;
    } ;

    virtual std::vector<Point * > & getAdditionalPoints() {
        return points ;
    };
    virtual const std::vector<Point * > & getAdditionalPoints() const {
        return points ;
    };

    virtual std::vector<ETYPE *> getConflictingElements ( const Geometry * g ) {
        std::vector<ETYPE *> ret ;
        if ( element->intersects ( g ) || g->in ( element->getCenter() ) || element->in ( g->getCenter() ) ) {
            ret.push_back ( element ) ;
        }

        return ret ;
    }
    /** Does nothing as this is a special-purpose mesh*/
    virtual void setElementOrder ( Order o, double dt = 0. ) {
        if(Mesh<ETYPE, EABSTRACTTYPE>::allElementsCacheID != -1)
        {
            Mesh<ETYPE, EABSTRACTTYPE>::caches[Mesh<ETYPE, EABSTRACTTYPE>::allElementsCacheID].clear() ;
            Mesh<ETYPE, EABSTRACTTYPE>::coefs[Mesh<ETYPE, EABSTRACTTYPE>::allElementsCacheID].clear() ;
            Mesh<ETYPE, EABSTRACTTYPE>::allElementsCacheID = -1 ;
        }
        switch ( o ) {
        case CONSTANT: {
            break ;
        }
        case LINEAR: {
            break ;
        }
        case QUADRATIC: {
            addSharedNodes ( 1,1,0 ) ;
            break ;
        }
        case CUBIC: {
            addSharedNodes ( 2,1,0 ) ;
            break ;
        }
        case QUADRIC: {
            addSharedNodes ( 3,1,0 ) ;
            break ;
        }
        case QUINTIC: {
            addSharedNodes ( 3,1,0 ) ;
            break ;
        }
        case CONSTANT_TIME_LINEAR: {
            addSharedNodes ( 0,2,dt ) ;
            break ;
        }
        case CONSTANT_TIME_QUADRATIC: {
            addSharedNodes ( 0,3,dt ) ;
            break ;
        }
        case LINEAR_TIME_LINEAR: {
            addSharedNodes ( 0,2,dt ) ;
            break ;
        }
        case LINEAR_TIME_QUADRATIC: {
            addSharedNodes ( 0,3,dt ) ;
            break ;
        }
        case QUADRATIC_TIME_LINEAR: {
            addSharedNodes ( 1,2,dt ) ;
            break ;
        }
        case QUADRATIC_TIME_QUADRATIC: {
            addSharedNodes ( 1,3,dt ) ;
            break ;
        }
        case CUBIC_TIME_LINEAR: {
            addSharedNodes ( 2,2,dt ) ;
            break ;
        }
        case CUBIC_TIME_QUADRATIC: {
            addSharedNodes ( 2,3,dt ) ;
            break ;
        }
        case QUADRIC_TIME_LINEAR: {
            addSharedNodes ( 3,2,dt ) ;
            break ;
        }
        case QUADRIC_TIME_QUADRATIC: {
            addSharedNodes ( 3,3,dt ) ;
            break ;
        }
        case QUINTIC_TIME_LINEAR: {
            addSharedNodes ( 3,2,dt ) ;
            break ;
        }
        case QUINTIC_TIME_QUADRATIC: {
            addSharedNodes ( 3,3,dt ) ;
            break ;
        }
        default:
            break ;

        }
    }

    /** Does nothing as this is a special-purpose mesh*/
    virtual void insert ( Point * ) {
    }

    virtual size_t getLastNodeId() const {
        return global_counter ;
    }
    virtual int addToTree ( EABSTRACTTYPE * toAdd ) {
        tree.push_back ( toAdd ) ;
        return tree.size() -1 ;
    }

    virtual EABSTRACTTYPE * getInTree ( int index ) const {
        return tree[std::abs ( index )] ;
    }

} ;
}


#endif // MESH_H
