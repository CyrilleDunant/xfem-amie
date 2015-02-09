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
public:

    virtual std::vector<int> & getCache ( unsigned int cacheID ) {
        return caches[cacheID] ;
    } ;
    virtual int addToTree ( EABSTRACTTYPE * toAdd ) = 0 ;
    virtual EABSTRACTTYPE * getInTree ( int index ) const = 0 ;
    virtual std::vector<Point * > & getAdditionalPoints() = 0 ;
    virtual const std::vector<Point * > & getAdditionalPoints() const = 0 ;
    virtual void extrude ( double dt ) = 0 ;
    virtual void extrude ( const Vector & dt ) = 0 ;
    virtual double getInternalScale() const {
        return 1. ;
    } ;
public:
    Mesh(SpaceDimensionality spaceDimensions) : spaceDimensions(spaceDimensions), allElementsCacheID(-1) {} ;
    virtual ~Mesh() {} ;

    virtual std::vector<ETYPE *> getConflictingElements ( const Point  * p )  = 0;
    virtual std::vector<ETYPE *> getConflictingElements ( const Geometry * g ) = 0;
    virtual std::vector<ETYPE *> getNeighbourhood ( ETYPE * element ) const = 0 ;
    
    virtual bool step(double dt) { return false ; }

    virtual std::vector<ETYPE *> getNeighbouringElementsInGeometry ( ETYPE * start , const Geometry * g ) {
        if ( !start ) {
// 					std::cout << "nullptr" << std::endl ;
            return std::vector<ETYPE *>() ;
        }

        std::set<ETYPE *> to_test ;
        std::set<ETYPE *> found ;
        found.insert ( start ) ;
        std::vector<ETYPE *> neighbourhood = getNeighbourhood ( start ) ;
        for ( const auto & neighbour : neighbourhood ) {
            if ( neighbour->timePlanes() > 1 ) {
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
            if ( element->in ( *p ) ) {
                return element ;
            }
        }

        return nullptr ;
    }

    virtual void setElementOrder ( Order o, double dt = 0. ) = 0;
    virtual void insert ( Point * ) = 0 ;
    template <class ETARGETTYPE>
    /** \brief Return the displacements in source mesh projected on current mesh.
    */
    void project ( Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false ) {
        if ( cache.find ( mesh ) == cache.end() ) {
            std::map<Point *, std::pair<ETYPE *, std::vector<double> > > projectionCache ;
            std::map<Point *, Point * > projectionPointCache ;
            std::vector<ETYPE *> selfElements = getElements() ;
            int pointCount = 0 ;
            size_t idCount = 0 ;
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
                    if ( targets[k]->in ( * ( *i ) ) || dist ( proj, * ( *i ) ) < 128.*POINT_TOLERANCE_2D ) {
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

    } ;

    template <class ETARGETTYPE>
    /** \brief Return the displacements in source mesh projected on current mesh.
    */
    void leastSquareProject ( const Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false ) {

    }

    virtual size_t getLastNodeId() const = 0;
    virtual size_t size() const = 0 ;

    virtual void deleteCache ( unsigned int cacheID ) {
        caches[cacheID].clear() ;
        coefs[cacheID].clear() ;
    }

    virtual unsigned int generateCache( Geometry * source) {
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

        std::vector<ETYPE *> elems = getConflictingElements ( source ) ;
        if(elems.empty())
            elems = getConflictingElements ( &source->getCenter() ) ;

        for ( auto & element : elems ) {
            if(source->in(element->getCenter()) && element->getBehaviour() && element->getBehaviour()->getSource() == source)
            {
                caches[position].push_back ( element->index ) ;
                coefs[position].push_back ( std::vector<double>() ) ;
            }
        }
        return position ;

    }

    virtual unsigned int generateCache( std::vector<Geometry *> source) {
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

        for(size_t i = 0 ; i < source.size() ; i++)
        {
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
	std::cout << caches[position].size() << std::endl ;

        return position ;

    }

    virtual unsigned int generateCache ( const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function ( "1" ) ) {
        size_t position = 0;
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
                if ( source && element->getBehaviour()->getSource() != source) {
                    continue ;
                }

                if ( locus->in ( element->getCenter() ) ) {
                    caches[position].push_back ( element->index ) ;
                    coefs[position].push_back ( std::vector<double>() ) ;
                    if(!source)
                        continue ;
                    if(element->getOrder() >= CONSTANT_TIME_LINEAR)
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
		                double xt = 0. ;//vm.eval ( t, gp.gaussPoints[i].first ) ;

		                coefs[position].back().push_back ( vm.eval ( smoothing, xx, xy, xz, xt ) );
			    }
                    }
		    else
		    {
		            Function x = element->getXTransform() ;
		            Function y = element->getYTransform() ;
		            Function z = element->getZTransform() ;
		            Function t = element->getTTransform() ;
		            GaussPointArray gp = element->getGaussPoints() ;
		            if(element->getOrder() >= CONSTANT_TIME_LINEAR)
		                gp = GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray( element, 0. ) ;
		            for ( size_t i = 0 ; i < gp.gaussPoints.size() ; i++ ) {
		                double xx = vm.eval ( x, gp.gaussPoints[i].first ) ;
		                double xy = vm.eval ( y, gp.gaussPoints[i].first ) ;
		                double xz = vm.eval ( z, gp.gaussPoints[i].first ) ;
		                double xt = vm.eval ( t, gp.gaussPoints[i].first ) ;

		                coefs[position].back().push_back ( vm.eval ( smoothing, xx, xy, xz, xt ) );
		            }
		    }
                }
            }
        }

        return position ;
    } ;

    virtual unsigned int generateCache ()
    {
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
    //virtual void getAverageField( Amie::FieldType f, Vector& ret, Amie::VirtualMachine* vm = nullptr, int dummy = 0, double t = 0, std::vector< double > weights = std::vector<double>()) ;
    Vector getField ( FieldType f, unsigned int cacheID, int dummy = 0, double t = 0 ) {
        VirtualMachine vm ;
        size_t blocks = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
            blocks = static_cast<ETYPE *>(getInTree ( caches[cacheID][i] ))->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions ;
        }
        Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
        Vector buffer ( ret ) ;
        double w = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {

            double v = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->getState().getAverageField ( f, buffer, nullptr, dummy, t, coefs[cacheID][i] ) ;
            ret += buffer * v ;
            w +=v ;
        }
        return ret/w ;
    }

    double getArea( unsigned int cacheID)
    {
        double a = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ )
            a += static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->area() ;
        return a ;
    }

    double getVolume( unsigned int cacheID)
    {
        double v = 0 ;
        for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ )
            v += static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) )->volume() ;
        return v ;
    }

    Vector getField ( FieldType f, int dummy = 0, double t = 0 ) {
        if(allElementsCacheID == -1)
        {
            VirtualMachine vm ;
            size_t blocks = 0 ;

            std::vector<ETYPE *> elems = getElements() ;
            for ( size_t i = 0 ; i < elems.size() && !blocks; i++ ) {
                blocks = elems[i]->getBehaviour()->getNumberOfDegreesOfFreedom() /spaceDimensions ;
            }
            Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
            Vector buffer ( ret ) ;
            double w = 0 ;
            for ( size_t i = 0 ; i < elems.size() ; i++ ) {

                double v = elems[i]->getState().getAverageField ( f, buffer, nullptr, dummy, t ) ;
                ret += buffer * v ;
                w +=v ;
            }
            return ret/w ;
        }

        VirtualMachine vm ;
        size_t blocks = 0 ;


        for ( auto i = begin() ; i  != end() && !blocks; i++ ) {
            blocks = i->getBehaviour()->getNumberOfDegreesOfFreedom() /spaceDimensions ;
        }
        Vector ret ( 0., fieldTypeElementarySize ( f, spaceDimensions, blocks ) ) ;
        Vector buffer ( ret ) ;
        double w = 0 ;
        for ( auto i = begin() ; i  != end() ; i++ ) {
            double v = i->getState().getAverageField ( f, buffer, nullptr, dummy, t ) ;
            ret += buffer * v ;
            w +=v ;
        }
        return ret/w ;
    }

    Vector getSmoothedField (  FieldType f0, unsigned int cacheID, IntegrableEntity * e,int dummy = 0, double t = 0 ) const {
        Vector first ;
        Vector strain ;
        Vector stress ;
        Vector strainrate ;
        Vector buffer ;
        int tsize = 3 ;
        int psize = 2 ;
        if (spaceDimensions == SPACE_THREE_DIMENSIONAL ) {
            tsize = 6 ;
            psize = 3 ;
        }
        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        VirtualMachine vm ;
        if ( f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
            //we first need to compute the strain field
            if ( !spaceTime ) {
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                    double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, nullptr, 0, t, coefs[cacheID][i] );
                    if ( !strain.size() ) {
                        strain.resize ( buffer.size(), 0. );
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
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;


                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, nullptr, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrain.size() ) {
                        tmpstrain.resize ( buffer.size(), 0. );
                    }
                    tmpstrain += buffer*v ;
                    sumFactors += v ;
                }
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, nullptr, dummy, t, coefs[cacheID][i] );
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
                Vector imposed = e->getBehaviour()->getImposedStress( Point() ) ;
                for ( size_t i = 0 ; i < tsize ; i++ ) {
                    stress[i] = tmpstress[i]-imposed[i] ;
                    strain[i] = tmpstrain[i] ;
                }
		if( e->getBehaviour()->hasInducedForces() )
			stress -= e->getBehaviour()->getImposedStress(Point(0,0,0,t)) ;
            }
//            std::cout << "here" << strain.size() << std::endl ;

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
                if ( !spaceTime ) {
                    first = strain*e->getBehaviour()->param ;
                } else {
                    first = stress ;
                }
            }
            if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                first.resize ( psize );
                if ( !spaceTime ) {
                    first = toPrincipal( strain*e->getBehaviour()->param ) ;
                } else {
                    first = toPrincipal( stress ) ;
                }
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
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                double v = ci->getState().getAverageField ( f0, buffer, nullptr, dummy, t, coefs[cacheID][i] );
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

    std::pair<Vector, Vector> getSmoothedFields ( FieldType f0, FieldType f1, unsigned int cacheID, IntegrableEntity * e ,int dummy = 0, double t = 0 ) const {
        Vector first ;
        Vector second ;
        Vector strain ;
        Vector stress ;
        Vector strainrate ;
        Vector buffer ;
        int tsize = 3 ;
        int psize = 2 ;
        if ( spaceDimensions == SPACE_THREE_DIMENSIONAL ) {
            tsize = 6 ;
            psize = 3 ;
        }
        bool spaceTime = e->getOrder() >= CONSTANT_TIME_LINEAR ;
        VirtualMachine vm ;
        if ( f0 == PRINCIPAL_STRAIN_FIELD || f0 == REAL_STRESS_FIELD || f0 == EFFECTIVE_STRESS_FIELD || f0 == PRINCIPAL_REAL_STRESS_FIELD || f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f0 == STRAIN_FIELD || f0 == MECHANICAL_STRAIN_FIELD || f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ||
		f1 == PRINCIPAL_STRAIN_FIELD || f1 == REAL_STRESS_FIELD || f1 == EFFECTIVE_STRESS_FIELD || f1 == PRINCIPAL_REAL_STRESS_FIELD || f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f1 == STRAIN_FIELD || f1 == MECHANICAL_STRAIN_FIELD || f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD
           ) {
            //we first need to compute the strain field
            if ( !spaceTime ) {
                buffer.resize ( tsize, 0. );
                strain.resize ( tsize, 0. );
                double sumFactors ( 0 ) ;

                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    IntegrableEntity *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                        double v = ci->getState().getAverageField ( STRAIN_FIELD, buffer, nullptr, dummy, t, coefs[cacheID][i] );
                        strain += buffer*v ;
                        sumFactors += v ;
                    }
                }
                strain /= sumFactors ;
            } else {
                size_t blocks = 0 ;
                for ( size_t i = 0 ; i < caches[cacheID].size() && !blocks; i++ ) {
                    ETYPE *ci  = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                    blocks = ci->getBehaviour()->getNumberOfDegreesOfFreedom() / spaceDimensions  ;
                }
                Vector tmpstrain ;
                Vector tmpstrainrate ;

                tmpstrain.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. ) ;
                buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, spaceDimensions, blocks ), 0. );
                tmpstrainrate.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, spaceDimensions, blocks ), 0. ) ;
                double sumFactors ( 0 ) ;


                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;


                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, buffer, nullptr, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrain.size() ) {
                        tmpstrain.resize ( buffer.size(), 0. );
                    }
                    tmpstrain += buffer*v ;
                    sumFactors += v ;
                }
                buffer.resize ( fieldTypeElementarySize ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD,spaceDimensions, blocks ), 0. );
                for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                    ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                    double v = ci->getState().getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, buffer, nullptr, dummy, t, coefs[cacheID][i] );
                    if ( !tmpstrainrate.size() ) {
                        tmpstrainrate.resize ( buffer.size(), 0. );
                    }
                    tmpstrainrate += buffer*v ;
                }
                tmpstrain /= sumFactors ;
                tmpstrainrate /=sumFactors ;

                Vector tmpstress = tmpstrain*e->getBehaviour()->getTensor ( Point(0,0,0,t) ) + ( Vector ) ( tmpstrainrate*e->getBehaviour()->getViscousTensor ( Point(0,0,0,t) ) ) ;
                stress.resize ( tsize, 0. ) ;
                strain.resize ( tsize, 0. ) ;
                for ( size_t i = 0 ; i < tsize ; i++ ) {
                    stress[i] = tmpstress[i] ;
                    strain[i] = tmpstrain[i] ;
                }
		if( e->getBehaviour()->hasInducedForces() )
			stress -= e->getBehaviour()->getImposedStress(Point(0,0,0,t)) ;
            }

            if ( f0 == PRINCIPAL_STRAIN_FIELD ) {
                first.resize ( psize );
                first = toPrincipal ( strain ) ;
            }
            if ( f1 == PRINCIPAL_STRAIN_FIELD ) {
                second.resize ( psize );
                second = toPrincipal ( strain ) ;
            }
            if ( f0 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ) {
                first.resize ( psize );
		if(e->getBehaviour() && e->getBehaviour()->hasInducedForces())
			strain -= e->getBehaviour()->getImposedStrain(Point(0,0,0,t)) ;
                first = toPrincipal ( strain ) ;
            }
            if ( f1 == PRINCIPAL_MECHANICAL_STRAIN_FIELD ) {
                second.resize ( psize );
		if(e->getBehaviour() && e->getBehaviour()->hasInducedForces())
			strain -= e->getBehaviour()->getImposedStrain(Point(0,0,0,t)) ;
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
                if ( !spaceTime ) {
                    first = strain*e->getBehaviour()->param ;
                } else {
                    first = stress ;
                }
            }
            if ( f1 == EFFECTIVE_STRESS_FIELD ) {
                second.resize ( tsize );
                if ( !spaceTime ) {
                    second = strain*e->getBehaviour()->param ;
                } else {
                    second = stress ;
                }
            }
            if ( f0 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                first.resize ( psize );
                if ( !spaceTime ) {
                    first = toPrincipal ( strain*e->getBehaviour()->param ) ;
                } else {
                    first = toPrincipal( stress ) ;
                }
            }
            if ( f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD ) {
                second.resize ( psize );
                if ( !spaceTime ) {
                    second = toPrincipal ( strain*e->getBehaviour()->param ) ;
                } else {
                    second = toPrincipal( stress ) ;
                }
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
	    if ( f0 == MECHANICAL_STRAIN_FIELD ) {
		first.resize( tsize ) ;
		first = strain ;
		if(e->getBehaviour() && e->getBehaviour()->hasInducedForces())
			first -= e->getBehaviour()->getImposedStrain(Point(0,0,0,t)) ;
	    }

	    if ( f1 == MECHANICAL_STRAIN_FIELD ) {
		second.resize( tsize ) ;
		second = strain ;
		if(e->getBehaviour() && e->getBehaviour()->hasInducedForces())
			second -= e->getBehaviour()->getImposedStrain(Point(0,0,0,t)) ;
	    }

        } else {
            double sumFactors ( 0 ) ;

            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;

                double v = ci->getState().getAverageField ( f0, buffer, nullptr, dummy, t, coefs[cacheID][i] );
                if ( !first.size() ) {
                    first.resize ( 0., buffer.size() );
                }
                first += buffer*v ;
                sumFactors += v ;
            }
            for ( size_t i = 0 ; i < caches[cacheID].size() ; i++ ) {
                ETYPE *ci = static_cast<ETYPE *> ( getInTree ( caches[cacheID][i] ) ) ;
                if ( ci->getBehaviour()->getSource() == e->getBehaviour()->getSource() ) {
                    double v = ci->getState().getAverageField ( f1, buffer, nullptr, dummy, t,coefs[cacheID][i] );
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

    class iterator
    {
    private:
        Mesh<ETYPE, EABSTRACTTYPE> * msh ;
        size_t cacheID ;
        size_t position ;

    public:
        iterator( Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t cacheID, size_t position) : msh(msh), cacheID(cacheID), position(position) { } ;
        iterator( Mesh<ETYPE, EABSTRACTTYPE> * msh, size_t position) : msh(msh), position(position)
        {
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
            return static_cast<ETYPE *>(msh->getInTree(msh->caches[cacheID][position])) ;
        }

        operator ETYPE * () {
            return static_cast<ETYPE *>(msh->getInTree(msh->caches[cacheID][position])) ;
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
            return static_cast<const ETYPE *>(msh->getInTree(msh->caches[cacheID][position])) ;
        }

        operator const ETYPE * () const {
            return static_cast<const ETYPE *>(msh->getInTree(msh->caches[cacheID][position])) ;
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
        if(allElementsCacheID == -1)
            allElementsCacheID = generateCache() ;
        return iterator(this,allElementsCacheID, caches[allElementsCacheID].size()) ;
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
                        Point * foundPoint = nullptr ;

                        for ( size_t j = 0 ; j< static_cast<ETYPE *> ( *i )->getBoundingPoints().size() ; j++ ) {
                            if ( static_cast<ETYPE *> ( *i )->getBoundingPoint ( j ) == proto ) {
                                foundPoint = &static_cast<ETYPE *> ( *i )->getBoundingPoint ( j ) ;
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
        for ( int i = 0 ; i < points.size() ; i++ ) {
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



// 		virtual std::vector<EABSTRACTTYPE *> & getTree() {return tree ; }
// 		virtual const std::vector<EABSTRACTTYPE *> & getTree() const {return tree ; }
} ;
}


#endif // MESH_H
