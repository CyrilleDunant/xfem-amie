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

#include "delaunay_3d.h"
#include "../features/boundarycondition.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <limits>

// #define DEBUG
// #undef DEBUG

using namespace Amie ;

DelaunayTreeItem3D::DelaunayTreeItem3D( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, DelaunayTreeItem3D *father,  const Point *c ) : stepson( 0 ), neighbour( 0 ), son( 0 )
{
    tree = t ;
    this->stepfather = -1 ;

    if( father )
        this->father = father->index;
    else
        this->father = -1 ;

    this->creator  = c ;
    this->dead = false ;
    this->erased = false ;

    this->isSpace = false ;
    this->isTetrahedron = false ;
    this->isDeadTetrahedron = false ;
    this->first = nullptr ;
    this->second = nullptr ;
    this->third = nullptr ;
    this->fourth = nullptr ;
    index = 0 ;

    if( t )
    {
        index = t->addToTree(this) ;
    }
}

DelaunayTreeItem3D::~DelaunayTreeItem3D()
{
// 	tree->tree[index] = nullptr ;
}


size_t DelaunayTreeItem3D::numberOfCommonVertices( const DelaunayTreeItem3D *s ) const
{
    if( this->isTetrahedron && s->isTetrahedron )
    {
        Point *inter[4] ;
        Point **e = std::set_intersection( &first, &first + 4, &s->first, &s->first + 4, &inter[0] )  ;
        return e - &inter[0] ;
    }

    if( this->isTetrahedron && s->isSpace )
    {
        return coplanarCount( &first, 4, *s->first, *s->second, *s->third, tree->getInternalScale()  ) ;
    }

    if( this->isSpace && s->isSpace )
    {
        return coplanarCount( &s->first, 3, *first, *second, *third, tree->getInternalScale()  ) ;
    }

    if( this->isSpace && s->isTetrahedron )
    {
        return coplanarCount( &s->first, 4, *first, *second, *third, tree->getInternalScale()  ) ;
    }

    return 0 ;
}


void Star3D::updateNeighbourhood()
{
    int count = 0 ;
    int soncount = 0 ;

    for( auto i = treeitem.begin() ; i != treeitem.end() ; ++i )
    {
        count += ( *i )->son.size() + ( *i )->stepson.size() + ( *i )->neighbour.size() ;
        soncount += ( *i )->son.size() ;
    }

    std::valarray<DelaunayTreeItem3D *> items( count ) ;
    count = 0 ;

    for( auto i = treeitem.begin() ; i != treeitem.end() ; ++i )
    {
        for( size_t j = 0 ; j < ( *i )->son.size() ; j++ )
        {
            items[count++] = ( *i )->getSon( j ) ;
        }

    }

    for( auto i = treeitem.begin() ; i != treeitem.end() ; ++i )
    {
        for( size_t j = 0 ; j < ( *i )->stepson.size() ; j++ )
        {
            if( !( *i )->getStepson( j )->dead && ( !( *i )->getStepson( j )->inCircumSphere( *creator ) || ( *i )->getStepson( j )->onCircumSphere( *creator ) ) )
                items[count++] = ( *i )->getStepson( j ) ;
        }

        for( size_t j = 0 ; j < ( *i )->neighbour.size() ; j++ )
        {
            if( !( *i )->getNeighbour( j )->dead && ( !( *i )->getNeighbour( j )->inCircumSphere( *creator ) || ( *i )->getNeighbour( j )->onCircumSphere( *creator ) ) ) 
                items[count++] = ( *i )->getNeighbour( j ) ;
        }

    }

    std::sort( &items[0], &items[count] ) ;
    auto e = std::unique( &items[0], &items[count] ) ;

    if( !items.size() )
        return ;

    for( DelaunayTreeItem3D **i = &items[0] ; i != e + 1 && i != &items[count] ; ++i )
    {
        if( !( *i )->isSpace )
        {
            DelaunayTreeItem3D *ii = *i ;
            size_t ins = ii->neighbour.size() ;

            if( ins != 4 )
            {
                for( DelaunayTreeItem3D **j = i + 1 ; j != e + 1 && j != &items[count] ; ++j )
                {

                    DelaunayTreeItem3D *jj = *j ;

                    size_t jns = jj->neighbour.size() ;

                    if( !jj->isSpace && jns != 4 && ii->numberOfCommonVertices( jj ) == 3 )
                    {

                        if( std::find( &ii->neighbour[0], &ii->neighbour[ins], jj->index ) ==  &ii->neighbour[ins] )
                        {
                            std::valarray<unsigned int>  newneighbours( ins + 1 ) ;
                            std::copy( &ii->neighbour[0], &ii->neighbour[ins], &newneighbours[0] ) ;
                            newneighbours[ins] = jj->index ;
                            ii->neighbour.resize( ins + 1 ) ;
                            ii->neighbour = newneighbours ;
                        }

                        if( std::find( &jj->neighbour[0], &jj->neighbour[jns], ii->index ) ==  &jj->neighbour[jns] )
                        {
                            std::valarray<unsigned int>  newneighbours( jns + 1 ) ;
                            std::copy( &jj->neighbour[0], &jj->neighbour[jns], &newneighbours[0] ) ;
                            newneighbours[jns] = ii->index ;
                            jj->neighbour.resize( jns + 1 ) ;
                            jj->neighbour = newneighbours ;
                        }

                        ins = ii->neighbour.size() ;

                        if( ins == 4 )
                            break ;
                    }
                }
            }
        }
    }

}

void DelaunayTree3D::extrude(const Vector & dt)
{
    std::map<Point *, std::vector<Point *> > points ;
    std::map<Point *, Point *> pointsInTetrahedron ;
    std::map<Point *, std::vector<Point *> >::iterator finder ;

    std::vector<DelaunayTetrahedron *> tri = getTetrahedrons() ;
    double beginning = tri[0]->getBoundingPoint(0).getT() ;
    double end = tri[0]->getBoundingPoint(0).getT() ;
    for(size_t i = 1 ; i < tri[0]->getBoundingPoints().size() ; i++)
    {
        if(tri[0]->getBoundingPoint(i).getT() < beginning)
            beginning = tri[0]->getBoundingPoint(i).getT() ;
        if(tri[0]->getBoundingPoint(i).getT() > end)
            end = tri[0]->getBoundingPoint(i).getT() ;
    }

    int indexOfLastTimePlane = (tri[0]->timePlanes()-1)*tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;
    int pointsPerTimePlane = tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
        {
            Point * current = &tri[i]->getBoundingPoint(j) ;
            current->getT() = dt[0]+(dt[1]-dt[0])*(current->getT() - beginning)/(end-beginning) ;
//			tri[i]->setBoundingPoint(j, current) ;
        }
    }


    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
        {
            Point * current = &tri[i]->getBoundingPoint(j) ;
            finder = points.find(current) ;
            if(finder == points.end())
            {
                Point * current = &tri[i]->getBoundingPoint(j) ;
                current->getT() = dt[0]+(dt[1]-dt[0])*(current->getT() - beginning)/(end-beginning) ;

                std::vector<Point *> newPoints ;
                if( (int)j < indexOfLastTimePlane)
                {
                    if( (int)j < pointsPerTimePlane )
                    {
                        newPoints.push_back(&tri[i]->getBoundingPoint(indexOfLastTimePlane+j)) ;
                    }
                    size_t ifirst = newPoints.size() ;
                    for(size_t k = ifirst+1 ; k < dt.size()-1 ; k++)
                    {
                        Point * next = new Point(current->getX(), current->getY()) ;
                        next->getT() = dt[k]+(dt[k+1]-dt[k])*(current->getT()-beginning)/(end-beginning) ;
                        next->setId(global_counter++) ;
                        newPoints.push_back(next) ;
                    }
                }
                else
                {
                    size_t jp = j - indexOfLastTimePlane ;
                    Point * previous = &tri[i]->getBoundingPoint(jp) ;
                    std::vector<Point *> previousPoints = points.find(previous)->second ;
                    for(size_t k = 1 ; k < previousPoints.size() ; k++)
                        newPoints.push_back(previousPoints[k]) ;
                    Point * next = new Point(current->getX(), current->getY()) ;
                    size_t k = dt.size()-2 ;
                    next->getT() = dt[k]+(dt[k+1]-dt[k])*(current->getT()-beginning)/(end-beginning) ;
                    next->setId(global_counter++) ;
                    newPoints.push_back(next) ;
                }

                points.insert(std::pair<Point *, std::vector<Point *> >(current, newPoints)) ;
            }
        }
    }

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        std::valarray<Point *> newPoints(tri[i]->getBoundingPoints().size()) ;
        for(size_t k = 1 ; k < dt.size()-1 ; k++)
        {
            for(size_t j = 0 ; j < newPoints.size() ; j++)
            {
                std::vector<Point *> found = points.find(&tri[i]->getBoundingPoint(j))->second ;
                newPoints[j] = found[k-1] ;
            }

            DelaunayTetrahedron * toInsert = new DelaunayTetrahedron(tri[i]->tree, nullptr, newPoints[0],newPoints[pointsPerTimePlane/4],newPoints[pointsPerTimePlane*2/4],newPoints[pointsPerTimePlane*3/4],newPoints[0]) ;
            toInsert->setOrder(tri[i]->getOrder()) ;
            toInsert->setBoundingPoints(newPoints) ;
            toInsert->setBehaviour(tri[i]->getState().getMesh3D(), tri[i]->getBehaviour()->getCopy()) ;
        }
    }

    std::cout << getTetrahedrons().size() << "\t" << tri.size() << std::endl ;

}

void DelaunayTree3D::addSharedNodes(DelaunayTree3D * dt)
{
    std::vector<DelaunayTetrahedron *> tet = getTetrahedrons() ;
    std::vector<DelaunayTetrahedron *> tets = dt->getTetrahedrons() ;

    for(size_t i = 0 ; i <  tets.size() ; i++)
    {
        std::valarray<Point *> newPoints = tets[i]->getBoundingPoints() ;
        tet[i]->setBoundingPoints(newPoints) ;
    }
}

void DelaunayTree3D::addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep, const TetrahedralElement *father )
{

    std::vector<DelaunayTetrahedron *> tet = getTetrahedrons() ;
    std::valarray<bool> visited(false, tree.size()) ;

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

    for( size_t i = 0 ; i < tet.size() ; i++ )
    {
        delete tet[i]->cachedGps ;
        tet[i]->cachedGps = nullptr ;
        tet[i]->behaviourUpdated = true ;

        std::vector<std::pair<Point, Point> > sides ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 1 ), tet[i]->getBoundingPoint( 2 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 2 ), tet[i]->getBoundingPoint( 3 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 0 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 2 ) ) ) ;


        visited[tet[i]->index] = true ;

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
                        for( size_t j = 0 ; j < tet[i]->neighbourhood.size() ; j++ )
                        {
                            if( visited[tet[i]->neighbourhood[j]] )
                            {
                                for( size_t k = 0 ; k < tet[i]->getNeighbourhood( j )->getBoundingPoints().size() ; k++ )
                                {
                                    if( tet[i]->getNeighbourhood( j )->getBoundingPoint( k ) == proto )
                                    {
                                        foundPoint = &tet[i]->getNeighbourhood( j )->getBoundingPoint( k ) ;
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
            if( !newPoints[k] )
            {
                std::cout << "ouch !" << std::endl ;

                for( size_t k = 0 ; k < newPoints.size() ; k++ )
                    if( newPoints[k] )
                        newPoints[k]->print() ;

                exit( 0 ) ;
            }

        tet[i]->setBoundingPoints( newPoints ) ;

// 		if( tet[i]->volume() < 0 )
// 		{
// 			for( int j = 0 ; j < tet[i]->getBoundingPoints().size() ; j++ )
// 				tet[i]->getBoundingPoint( j ).print() ;
//
// 			exit( 0 ) ;
// 		}
    }


    for( size_t i = 0 ; i < tet.size() ; i++ )
    {
        if( father != nullptr )
            refresh( father ) ;
    }
}

void DelaunayTree3D::refresh( const TetrahedralElement *father )
{

    for( size_t i = 0 ; i < tree.size() ; i++ )
    {
        if( !tree[i]->dead && tree[i]->isTetrahedron )
        {
            static_cast<DelaunayTetrahedron *>( tree[i] )->refresh( father ) ;

        }
    }
}

size_t DelaunayTree3D::numPoints() const
{
    return this->global_counter ;
}


void DelaunayTreeItem3D::flatConflicts(std::valarray<bool> & visited, std::vector<DelaunayTreeItem3D *> & toTest,  std::vector<DelaunayTreeItem3D *>  & ret, const Geometry *g )
{
    if( visited[index] )
    {
        return ;
    }

    visited[index] = true ;


    if(
        !inCircumSphere( g->getCenter() )
        &&
        (
            ( !g->in( *first )
              && !g->in( *second )
              && !g->in( *third )
              && !g->in( *fourth )
              && isTetrahedron )
            ||
            (
                g->in( *fourth )
                && isSpace
            )
        )
    )
    {
        return ;
    }

    for( size_t i  = 0 ;  i < stepson.size() ; i++ )
    {
        if( !visited[stepson[i]]
                &&
                !(
                    !getStepson( i )->inCircumSphere( g->getCenter() )
                    &&
                    (
                        ( !g->in( *getStepson( i )->first )
                          && !g->in( *getStepson( i )->second )
                          && !g->in( *getStepson( i )->third )
                          && !g->in( *getStepson( i )->fourth )
                          && getStepson( i )->isTetrahedron
                        )
                        ||
                        (
                            g->in( *getStepson( i )->fourth )
                            && getStepson( i )->isSpace
                        )
                    )
                )
          )
        {
            toTest.push_back( getStepson( i ) ) ;
        }
    }

    for( size_t i  = 0 ;  i < son.size() ; i++ )
    {
        if(
            !visited[son[i]]
            &&
            !(
                !getSon( i )->inCircumSphere( g->getCenter() )
                &&
                (
                    ( !g->in( *getSon( i )->first )
                      && !g->in( *getSon( i )->second )
                      && !g->in( *getSon( i )->third )
                      && !g->in( *getSon( i )->fourth )
                      && getSon( i )->isTetrahedron
                    )
                    ||
                    (
                        g->in( *getSon( i )->fourth )
                        && getSon( i )->isSpace
                    )
                )
            )
        )
        {
            toTest.push_back( getSon( i ) ) ;
        }
    }

    if( !dead && isTetrahedron )
    {
        ret.push_back( static_cast<DelaunayTetrahedron *>( this ) ) ;
    }

    for( size_t i  = 0 ;  i < neighbour.size() ; i++ )
    {
        if(
            !visited[neighbour[i]]
            &&
            !(
                !getNeighbour( i )->inCircumSphere( g->getCenter() )
                &&
                (
                    (
                        !g->in( *getNeighbour( i )->first )
                        && !g->in( *getNeighbour( i )->second )
                        && !g->in( *getNeighbour( i )->third )
                        && !g->in( *getNeighbour( i )->fourth )
                        && getNeighbour( i )->isTetrahedron
                    )
                    ||
                    (
                        g->in( *getNeighbour( i )->fourth )
                        && getNeighbour( i )->isSpace
                    )
                )
            )
        )
        {
            toTest.push_back( getNeighbour( i ) ) ;
        }
    }
}

void DelaunayTreeItem3D::flatConflicts(std::valarray<bool> & visited, std::vector< DelaunayTreeItem3D * >& toTest, std::vector< DelaunayTreeItem3D * > & ret, const Point *p, int threadid )
{
//     if(threadid >= 0 && index%(threadid+1) != 0)
//         return ;

    if(visited[index])
        return  ;

    visited[index] = true ;

    if(!inCircumSphere(*p))
        return  ;

// 	std::vector< DelaunayTreeItem3D * > toteststepson ;
// 	std::vector< DelaunayTreeItem3D * > totestson ;
// 	std::vector< DelaunayTreeItem3D * > totestneighbour ;

// #pragma omp parallel
// {
// 	#pragma omp sections
// 	{
// 	#pragma omp section
// 	{
    for ( size_t i  = 0 ;  i < stepson.size() ; i++)
    {
//         bool limit = false ;
        if(!visited[stepson[i]])
        {
//             if(getStepson(i)->isTetrahedron && !getStepson(i)->isDeadTetrahedron)
//             {
//                 DelaunayTetrahedron * t = static_cast<DelaunayTetrahedron *>(getStepson(i)) ;
// //                 limit = std::abs(squareDist3D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius())
// //                         < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
//             }
//             if(getStepson(i)->isDeadTetrahedron)
//             {
//                 DelaunayDeadTetrahedron * t = static_cast<DelaunayDeadTetrahedron *>(getStepson(i)) ;
// //                 limit = std::abs(squareDist3D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius())
// //                         < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
//             }

            if( (getStepson(i)->inCircumSphere(*p)) /*|| limit*/)
            {
                toTest.push_back(getStepson(i)) ;
            }
        }
    }
// 	}

// 	#pragma omp section
// 	{
    for (size_t i  = 0 ;  i < son.size() ; i++)
    {
//         bool limit = false ;
        if(!visited[son[i]])
        {

            if( (getSon(i)->inCircumSphere(*p)) /*|| limit*/)
            {
                toTest.push_back(getSon(i)) ;
            }
        }
    }
// 	}

// 	#pragma omp section
// 	{
    for (size_t i  = 0 ;  i < neighbour.size() ; i++)
    {
//         bool limit = false ;
        if(!visited[neighbour[i]])
        {

            if( (getNeighbour(i)->inCircumSphere(*p)) /*|| limit*/)
            {
                toTest.push_back(getNeighbour(i)) ;
            }
        }
    }

    if(!dead)
    {
        ret.push_back(this) ;
    }
}


void DelaunayTreeItem3D::removeNeighbour( DelaunayTreeItem3D *t )
{

    unsigned int* e = std::find(&neighbour[0], &neighbour[neighbour.size()], t->index) ;

    if(e !=  &neighbour[neighbour.size()])
    {
        std::valarray<unsigned int>  newneighbours(neighbour.size()-1) ;
        int iterator = 0 ;
        for(size_t i = 0 ; i < neighbour.size() ; i++)
        {
            if((int)neighbour[i] == t->index)
                continue ;

            newneighbours[iterator] = neighbour[i] ;
            iterator++ ;
        }

        neighbour.resize(neighbour.size()-1) ;
        neighbour = newneighbours ;
    }
}

void DelaunayTetrahedron::removeNeighbourhood( DelaunayTetrahedron *t )
{
    unsigned int *e = std::find( &neighbourhood[0], &neighbourhood[neighbourhood.size()], t->index ) ;

    if( e !=  &neighbourhood[neighbourhood.size()] )
    {
        std::valarray<unsigned int>  newneighbours( neighbourhood.size() - 1 ) ;
        std::copy( &neighbourhood[0], e, &newneighbours[0] ) ;
        std::copy( e + 1, &neighbourhood[neighbourhood.size()], &newneighbours[e - &neighbourhood[0]] ) ;
        neighbourhood.resize( neighbourhood.size() - 1 ) ;
        neighbourhood = newneighbours ;
    }
}

void DelaunayTreeItem3D::addNeighbour( DelaunayTreeItem3D *t )
{
    assert( t != this ) ;

    if( std::find( &neighbour[0], &neighbour[neighbour.size()], t->index ) !=  &neighbour[neighbour.size()] )
    {
        return ;
    }

    if( !t->dead )
    {
        std::valarray<unsigned int>  newneighbours( neighbour.size() + 1 ) ;
        std::copy( &neighbour[0], &neighbour[neighbour.size()], &newneighbours[0] ) ;
        newneighbours[neighbour.size()] = t->index ;
        neighbour.resize( neighbour.size() + 1 ) ;
        neighbour = newneighbours ;
    }
}

void DelaunayTetrahedron::addNeighbourhood( DelaunayTetrahedron *t )
{
    assert( t != this ) ;

    if( std::find( &neighbourhood[0], &neighbourhood[neighbourhood.size()], t->index ) !=  &neighbourhood[neighbourhood.size()] )
    {
        return ;
    }

    if( !t->dead )
    {
        std::valarray<unsigned int>  newneighbours( neighbourhood.size() + 1 ) ;
        std::copy( &neighbourhood[0], &neighbourhood[neighbourhood.size()], &newneighbours[0] ) ;
        newneighbours[neighbourhood.size()] = t->index ;
        neighbourhood.resize( neighbourhood.size() + 1 ) ;
        neighbourhood = newneighbours ;
    }
}

DelaunayTetrahedron *DelaunayTetrahedron::getNeighbourhood( size_t i ) const
{
    return static_cast<DelaunayTetrahedron *>( tree->getInTree(neighbourhood[i]) ) ;
}

DelaunayTreeItem3D *DelaunayTreeItem3D::getNeighbour( size_t i ) const
{
    return tree->getInTree(neighbour[i]) ;
}

DelaunayTreeItem3D *DelaunayTreeItem3D::getSon( size_t i ) const
{
    return tree->getInTree(son[i]) ;
}

DelaunayTreeItem3D *DelaunayTreeItem3D::getStepson( size_t i ) const
{
    if(i < stepson.size())
        return tree->getInTree(stepson[i]) ;
    return nullptr ;
}

DelaunayTreeItem3D * DelaunayTreeItem3D::getFather() const
{
    return tree->getInTree(father) ;
}
DelaunayTreeItem3D * DelaunayTreeItem3D::getStepfather() const
{
    return tree->getInTree(stepfather) ;
}

void DelaunayTreeItem3D::erase( const Point *p )
{
    dead = true ;
    killer  = p ;
}

void DelaunayTreeItem3D::kill( const Point *p )
{
    dead = true ;
    killer  = p ;

    for( size_t i = 0 ; i < this->neighbour.size() ; i++ )
    {
        this->getNeighbour( i )->removeNeighbour( this ) ;
    }
}

void DelaunayTreeItem3D::addStepson( DelaunayTreeItem3D *s )
{
    if( s == this )
        return ;

    std::valarray<unsigned int>  newstepson( stepson.size() + 1 ) ;
    std::copy( &stepson[0], &stepson[stepson.size()], &newstepson[0] ) ;
    newstepson[stepson.size()] = s->index ;
    stepson.resize( stepson.size() + 1 ) ;
    stepson = newstepson ;
    s->setStepfather( this ) ;
    addNeighbour( s ) ;
}

void DelaunayTreeItem3D::addSon( DelaunayTreeItem3D *s )
{
    std::valarray<unsigned int>  newson( son.size() + 1 ) ;
    std::copy( &son[0], &son[son.size()], &newson[0] ) ;
    newson[son.size()] = s->index ;
    son.resize( son.size() + 1 ) ;
    son = newson ;
}

void DelaunayTreeItem3D::removeSon( DelaunayTreeItem3D *t )
{
    if( !son.size() )
        return ;

    unsigned int *e = std::find( &son[0], &son[son.size()], t->index ) ;

    if( e !=  &son[son.size()] )
    {
        std::valarray<unsigned int>  newson( son.size() - 1 ) ;
        std::copy( &son[0], e, &newson[0] ) ;
        std::copy( e + 1, &son[son.size()], &newson[e - &son[0]] ) ;
        son.resize( son.size() - 1 ) ;
        son = newson ;
    }
}

void DelaunayTreeItem3D::removeStepson( DelaunayTreeItem3D *t )
{
    if( stepson.size() == 0 )
        return ;

    unsigned int *e = std::find( &stepson[0], &stepson[stepson.size()], t->index ) ;

    if( e !=  &stepson[stepson.size()] )
    {
        std::valarray<unsigned int>  newstepson( stepson.size() - 1 ) ;
        std::copy( &stepson[0], e, &newstepson[0] ) ;
        std::copy( e + 1, &stepson[stepson.size()], &newstepson[e - &stepson[0]] ) ;
        stepson.resize( stepson.size() - 1 ) ;
        stepson = newstepson ;

        if( newstepson.size() == 0 )
            t->stepfather = -1 ;

    }
}

void DelaunayTreeItem3D::setStepfather( DelaunayTreeItem3D *s )
{
    if( stepfather != -1 )
        tree->getInTree(stepfather)->removeStepson( this ) ;

    stepfather = s->index ;
    addNeighbour( s ) ;
}

void DelaunayTreeItem3D::setFather( DelaunayTreeItem3D *s )
{
    if( father != -1 )
        tree->getInTree(father)->removeSon( this ) ;

    father = s->index ;
}

DelaunayTetrahedron::~DelaunayTetrahedron()
{
}

void DelaunayTetrahedron::kill( const Point *p )
{
    this->DelaunayTreeItem3D::kill( p ) ;
    getBoundingPoints().resize( 0 ) ;
    getInPoints().resize( 0 ) ;
}

DelaunayTetrahedron::DelaunayTetrahedron( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, DelaunayTreeItem3D *father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point *c ) : TetrahedralElement( p0, p1, p2, p3 ), DelaunayTreeItem3D( t, father, c )
{
    first = &getBoundingPoint( 0 ) ;
    second = &getBoundingPoint( 1 ) ;
    third = &getBoundingPoint( 2 ) ;
    fourth = &getBoundingPoint( 3 ) ;

    std::sort( &first, &first + 4 ) ;

    isSpace = false ;
    isTetrahedron = true ;
    isDeadTetrahedron = false ;

}

DelaunayTetrahedron::DelaunayTetrahedron( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t, DelaunayTreeItem3D *father,  Point *p0,  Point *p1, Point *p2, Point *p3,  Point *p4,  Point *p5, Point *p6, Point *p7, Point *c ) : TetrahedralElement( p0, p1, p2, p3, p4, p5, p6, p7 ), DelaunayTreeItem3D( t, father, c )
{

    first = &getBoundingPoint( 0 ) ;
    second = &getBoundingPoint( 2 ) ;
    third = &getBoundingPoint( 4 ) ;
    fourth = &getBoundingPoint( 6 ) ;

    std::sort( &first, &first + 4 ) ;

    isSpace = false ;
    isTetrahedron = true ;
    isDeadTetrahedron = false ;

}

DelaunayTetrahedron::DelaunayTetrahedron() : DelaunayTreeItem3D( nullptr, nullptr, nullptr )
{
    first = &getBoundingPoint( 0 ) ;
    second = &getBoundingPoint( 1 ) ;
    third = &getBoundingPoint( 2 ) ;
    fourth = &getBoundingPoint( 3 ) ;

    isSpace = false ;
    isTetrahedron = true ;
    isDeadTetrahedron = false ;
    index = 0 ;
}

DelaunayDemiSpace::~DelaunayDemiSpace()
{
}

bool DelaunayTetrahedron::isVertex( const Point *p ) const
{
    return squareDist3D( *p, *first ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE
           || squareDist3D( *p, *second ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE
           || squareDist3D( *p, *third ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE
           || squareDist3D( *p, *fourth ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE ;
    return ( dist( p, first ) < 0.0000001 || dist( p, second ) < 0.0000001 || dist( p, third ) < 0.0000001 || dist( p, fourth ) < 0.0000001 ) ;
}

bool DelaunayDemiSpace::isVertex( const Point *p ) const
{
    return squareDist3D( *p, *first ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE
           || squareDist3D( *p, *second ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE
           || squareDist3D( *p, *third ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE ;
}

bool DelaunayTetrahedron::isVertexByID( const Point *p ) const
{
    return ( p->getId() == first->getId() || p->getId() == second->getId() || p->getId() == third->getId() || p->getId() == fourth->getId() ) ;
}


bool DelaunayTetrahedron::hasVertexByID( const std::valarray<Point *> * p ) const
{
    for( size_t i = 0 ; i < p->size() ; i++ )
    {
        if( ( *p )[i]->getId() == first->getId() || ( *p )[i]->getId() == second->getId() || ( *p )[i]->getId() == third->getId() || ( *p )[i]->getId() == fourth->getId() )
            return true ;
    }

    return false ;
}

DelaunayDeadTetrahedron::DelaunayDeadTetrahedron( DelaunayTetrahedron *parent ) : DelaunayTreeItem3D( *parent ),
    center( parent->getCircumCenter() ),
    radius( parent->getRadius() )
{
    index = parent->index ;
    tree = parent->tree ;

    stepson.resize( parent->stepson.size() ) ;
    stepson = parent->stepson ;

    neighbour.resize( parent->neighbour.size() ) ;
    neighbour = parent->neighbour ;

    son.resize( parent->son.size() ) ;
    son = parent->son ;

    father = parent->father ;
    stepfather = parent->stepfather ;

    isTetrahedron = true ;
    isSpace = false ;
    isDeadTetrahedron = true ;

    first = parent->first ;
    second = parent->second ;
    third = parent->third ;
    fourth = parent->fourth ;
    dead = true ;
}

DelaunayDeadTetrahedron::~DelaunayDeadTetrahedron() { } 

const Point * DelaunayDeadTetrahedron::getCircumCenter() const
{
    return &center ;
}

double DelaunayDeadTetrahedron::getRadius() const
{
    return radius ;
}



std::vector< Point *> DelaunayDeadTetrahedron::commonSurface( DelaunayTreeItem3D *t )
{
    std::vector<Point *> ret ;

    if( t == this )
    {
        return ret ;
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr );
    }

    if( t->isTetrahedron )
    {
        if( this->isVertexByID( t->first ) && this->isVertexByID( t->second ) && this->isVertexByID( t->third ) )
        {
            ret.push_back( t->first ) ;
            ret.push_back( t->second ) ;
            ret.push_back( t->third ) ;
        }


        if( this->isVertexByID( t->first ) && this->isVertexByID( t->second ) && this->isVertexByID( t->fourth ) )
        {
            ret.push_back( t->first ) ;
            ret.push_back( t->second ) ;
            ret.push_back( t->fourth ) ;
        }

        if( this->isVertexByID( t->first ) && this->isVertexByID( t->fourth ) && this->isVertexByID( t->third ) )
        {
            ret.push_back( t->first ) ;
            ret.push_back( t->third ) ;
            ret.push_back( t->fourth ) ;
        }

        if( this->isVertexByID( t->fourth ) && this->isVertexByID( t->second ) && this->isVertexByID( t->third ) )
        {
            ret.push_back( t->third ) ;
            ret.push_back( t->second ) ;
            ret.push_back( t->fourth ) ;
        }


    }
    else
    {
        Point A( *first*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        Point B( *third*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        Point C( *third*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        Point D( *third*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        Point E( *third*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;
        double fst = std::abs( triProduct( A, B, C ) ) ;
        fst = std::max( fst, std::abs( triProduct( A, B, D ) ) ) ;
        fst = std::max( fst, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *first*tree->getInternalScale() - *third*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *third*tree->getInternalScale() ) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double ftf = std::abs( triProduct( A, B, C ) ) ;
        ftf = std::max( ftf, std::abs( triProduct( A, B, D ) ) ) ;
        ftf = std::max( ftf, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *first*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double fsf = std::abs( triProduct( A, B, C ) ) ;
        fsf = std::max( fsf, std::abs( triProduct( A, B, D ) ) ) ;
        fsf = std::max( fsf, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *third*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double tsf = std::abs( triProduct( A, B, C ) ) ;
        tsf = std::max( tsf, std::abs( triProduct( A, B, D ) ) ) ;
        tsf = std::max( tsf, std::abs( triProduct( A, B, E ) ) ) ;

        if( fst <= ftf
                && fst <= fsf
                && fst <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( second ) ;
            ret.push_back( third ) ;
        }

        else if( ftf <= fst
                 && ftf <= fsf
                 && ftf <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( third ) ;
            ret.push_back( fourth ) ;
        }

        else if( fsf <= fst
                 && fsf <= ftf
                 && fsf <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( second ) ;
            ret.push_back( fourth ) ;
        }

        else
        {
            ret.push_back( third ) ;
            ret.push_back( second ) ;
            ret.push_back( fourth ) ;
        }

    }

    return ret;
}


bool DelaunayDeadTetrahedron::inCircumSphere( const Point &p ) const
{

    if( p.getX() > center.getX() + 1.001 * radius )
        return false ;

    if( p.getX() < center.getX() - 1.001 * radius )
        return false ;

    if( p.getY() > center.getY() + 1.001 * radius )
        return false ;

    if( p.getY() < center.getY() - 1.001 * radius )
        return false ;

    if( p.getZ() > center.getZ() + 1.001 * radius )
        return false ;

    if( p.getZ() < center.getZ() - 1.001 * radius )
        return false ;

    Point pr = p*tree->getInternalScale()  ;
    double d = sqrt(
                   ( center.getX()*tree->getInternalScale()  - pr.getX() ) * ( center.getX()*tree->getInternalScale() - pr.getX() ) +
                   ( center.getY()*tree->getInternalScale()  - pr.getY() ) * ( center.getY()*tree->getInternalScale() - pr.getY() ) +
                   ( center.getZ()*tree->getInternalScale()  - pr.getZ() ) * ( center.getZ()*tree->getInternalScale() - pr.getZ() ) ) ;
    return  d - radius*tree->getInternalScale()  < POINT_TOLERANCE * radius*tree->getInternalScale()  ;
}

bool DelaunayDeadTetrahedron::onCircumSphere( const Point &p ) const
{
    if( p.getX() > center.getX() + 1.001 * radius )
        return false ;

    if( p.getX() < center.getX() - 1.001 * radius )
        return false ;

    if( p.getY() > center.getY() + 1.001 * radius )
        return false ;

    if( p.getY() < center.getY() - 1.001 * radius )
        return false ;

    if( p.getZ() > center.getZ() + 1.001 * radius )
        return false ;

    if( p.getZ() < center.getZ() - 1.001 * radius )
        return false ;

    Point pr = p*tree->getInternalScale()  ;
    double d = sqrt(
                   ( center.getX()*tree->getInternalScale()  - pr.getX() ) * ( center.getX()*tree->getInternalScale() - pr.getX() ) +
                   ( center.getY()*tree->getInternalScale()  - pr.getY() ) * ( center.getY()*tree->getInternalScale() - pr.getY() ) +
                   ( center.getZ()*tree->getInternalScale()  - pr.getZ() ) * ( center.getZ()*tree->getInternalScale() - pr.getZ() ) ) ;
    return  std::abs( d - radius*tree->getInternalScale()  ) < POINT_TOLERANCE * radius*tree->getInternalScale()  ;
}

bool DelaunayDeadTetrahedron::isNeighbour( const DelaunayTreeItem3D *t ) const
{
    size_t cv = this->numberOfCommonVertices( t ) ;

    return ( cv == 3 );

}

bool DelaunayDeadTetrahedron::isVertex( const Point *p ) const
{
    double sc = tree->getInternalScale() ;
    return squareDist3D( *p*sc, *first*sc ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE*sc*sc
           || squareDist3D( *p*sc, *second*sc ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE*sc*sc
           || squareDist3D( *p*sc, *third*sc ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE*sc*sc
           || squareDist3D( *p*sc, *fourth*sc ) < 128.*POINT_TOLERANCE * POINT_TOLERANCE*sc*sc ;
}

bool DelaunayDeadTetrahedron::isVertexByID( const Point *p )
{
    return p == first || p == second || p == third || p == fourth ;
}

bool DelaunayDeadTetrahedron::normalisedIn( const Point &p ) const
{
    double internalScale = tree->getInternalScale() ;
    Point pf = *first  * internalScale ;
    Point ps = *second * internalScale ;
    Point pt = *third  * internalScale ;
    Point po = *fourth * internalScale ;
    return Tetrahedron( pf, ps, pt, po ).in( p*internalScale ) ;
}

bool DelaunayTetrahedron::normalisedIn( const Point &p ) const
{
    double internalScale = tree->getInternalScale() ;
    Point pf = *first  * internalScale ;
    Point ps = *second * internalScale ;
    Point pt = *third  * internalScale ;
    Point po = *fourth * internalScale ;
    return Tetrahedron( pf, ps, pt, po ).in( p*internalScale ) ;
}

void DelaunayDeadTetrahedron::print() const
{

    std::cout << "(" << first->getX() << ", " << first->getY() << ", " << first->getZ() << ") " ;
    std::cout << "(" << second->getX() << ", " << second->getY() << ", " << second->getZ() << ") " ;
    std::cout << "(" << third->getX() << ", " << third->getY() << ", " << third->getZ() << ") " ;
    std::cout << "(" << fourth->getX() << ", " << fourth->getY() << ", " << fourth->getZ() << ") " ;
    std::cout <<  ":: " << !dead << std::endl ;
}

std::vector<Point *> DelaunayTetrahedron::commonSurface( DelaunayTreeItem3D *t )
{
    std::vector<Point *> ret ;

    if( t == this )
    {
        return ret ;
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr );
    }

    if( t->isTetrahedron )
    {
        std::vector<Point *> inter( 4 ) ;
        auto e = std::set_intersection( &first, &first + 4, &t->first, &t->first + 4, inter.begin() )  ;
        inter.erase( e, inter.end() ) ;
        return inter ;
    }
    else
    {
        Point A( *first*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        Point B( *third*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        Point C( *third*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        Point D( *third*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        Point E( *third*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;
        double fst = std::abs( triProduct( A, B, C ) ) ;
        fst = std::max( fst, std::abs( triProduct( A, B, D ) ) ) ;
        fst = std::max( fst, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *first*tree->getInternalScale()  - *third*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *third*tree->getInternalScale() ) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double ftf = std::abs( triProduct( A, B, C ) ) ;
        ftf = std::max( ftf, std::abs( triProduct( A, B, D ) ) ) ;
        ftf = std::max( ftf, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *first*tree->getInternalScale()  - *second*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *second*tree->getInternalScale()) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double fsf = std::abs( triProduct( A, B, C ) ) ;
        fsf = std::max( fsf, std::abs( triProduct( A, B, D ) ) ) ;
        fsf = std::max( fsf, std::abs( triProduct( A, B, E ) ) ) ;

        A = ( *third*tree->getInternalScale()  - *second*tree->getInternalScale() ) ;
        B = ( *fourth*tree->getInternalScale() - *second*tree->getInternalScale() ) ;
        C = ( *fourth*tree->getInternalScale() - *t->first*tree->getInternalScale() ) ;
        D = ( *fourth*tree->getInternalScale() - *t->second*tree->getInternalScale() ) ;
        E = ( *fourth*tree->getInternalScale() - *t->third*tree->getInternalScale() ) ;

        double tsf = std::abs( triProduct( A, B, C ) ) ;
        tsf = std::max( tsf, std::abs( triProduct( A, B, D ) ) ) ;
        tsf = std::max( tsf, std::abs( triProduct( A, B, E ) ) ) ;

        if( fst <= ftf
                && fst <= fsf
                && fst <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( second ) ;
            ret.push_back( third ) ;
        }

        else if( ftf <= fst
                 && ftf <= fsf
                 && ftf <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( third ) ;
            ret.push_back( fourth ) ;
        }

        else if( fsf <= fst
                 && fsf <= ftf
                 && fsf <= tsf )
        {
            ret.push_back( first ) ;
            ret.push_back( second ) ;
            ret.push_back( fourth ) ;
        }

        else
        {
            ret.push_back( third ) ;
            ret.push_back( second ) ;
            ret.push_back( fourth ) ;
        }

    }

    return ret;
}

inline std::vector<Point *> DelaunayDemiSpace::commonSurface( DelaunayTreeItem3D *t )
{
    std::vector<Point *> ret ;

    if( t == this )
    {
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr ) ;
        ret.push_back( nullptr );
    }

    if( t->isTetrahedron )
    {
        return t->commonSurface( this ) ;
    }

    return ret;
}

std::vector<Point *> DelaunayRoot3D::commonSurface( DelaunayTreeItem3D *t )
{
    std::vector<Point *> ret ;

    ret.push_back( nullptr ) ;
    ret.push_back( nullptr ) ;
    ret.push_back( nullptr );

    return ret;
}

void DelaunayDemiSpace::merge( DelaunayDemiSpace *p )
{
    if( !dead && p != this && !p->dead )
    {
        if( std::binary_search( &first, &first + 3, p->first ) &&
                std::binary_search( &first, &first + 3, p->second ) &&
                std::binary_search( &first, &first + 3, p->third ) )
        {
            for( size_t i = 0 ; i <  p->neighbour.size() ; i++ )
            {
                for( size_t j = 0 ; j <  neighbour.size() ; j++ )
                {
                    if( p->getNeighbour( i )->isNeighbour( getNeighbour( j ) ) )
                        makeNeighbours3d( getNeighbour( j ), p->getNeighbour( i ) ) ;
                }
            }

            p->kill( first ) ;
            kill( first ) ;
            return ;
        }

        if( isCoplanar( p->first, first, second, third, tree->getInternalScale() ) &&
                isCoplanar( p->second, first, second, third, tree->getInternalScale() ) &&
                isCoplanar( p->third, first, second, third, tree->getInternalScale() ) )
        {
            for( size_t i = 0 ; i <  p->neighbour.size() ; i++ )
            {
                this->addNeighbour( p->getNeighbour( i ) ) ;
                p->getNeighbour( i )->addNeighbour( this ) ;
            }

            for( size_t i = 0 ; i <  p->son.size() ; i++ )
            {
                p->getSon( i )->setFather( this ) ;
                addSon( p->getSon( i ) ) ;
            }

            for( size_t i = 0 ; i <  p->stepson.size() ; i++ )
            {
                std::valarray<unsigned int>  newstepson( stepson.size() + 1 ) ;

                if( stepson.size() != 0 )
                    std::copy( &stepson[0], &stepson[stepson.size()], &newstepson[0] ) ;

                newstepson[stepson.size()] = p->stepson[i] ;
                stepson.resize( newstepson.size() ) ;
                stepson = newstepson ;
                p->getStepson( i )->stepfather = index ;

                addNeighbour( p->getStepson( i ) ) ;
                p->getStepson( i )->addNeighbour( this ) ;
            }

            tree->getInTree(p->father)->addSon( this ) ;
            p->addSon( this ) ;
            p->kill( first ) ;
        }
    }
}

bool DelaunayTetrahedron::inCircumSphere( const Point &p ) const
{
    if( p.getX() > circumCenter.getX() + 1.001 * radius )
        return false ;

    if( p.getX() < circumCenter.getX() - 1.001 * radius )
        return false ;

    if( p.getY() > circumCenter.getY() + 1.001 * radius )
        return false ;

    if( p.getY() < circumCenter.getY() - 1.001 * radius )
        return false ;

    if( p.getZ() > circumCenter.getZ() + 1.001 * radius )
        return false ;

    if( p.getZ() < circumCenter.getZ() - 1.001 * radius )
        return false ;

    double d = dist( circumCenter*tree->getInternalScale(), p*tree->getInternalScale() ) ;
    return  d - radius*tree->getInternalScale() < POINT_TOLERANCE * radius*tree->getInternalScale() ;
}

bool DelaunayTetrahedron::onCircumSphere( const Point &p ) const
{
    if( p.getX() > circumCenter.getX() + 1.001 * radius )
        return false ;

    if( p.getX() < circumCenter.getX() - 1.001 * radius )
        return false ;

    if( p.getY() > circumCenter.getY() + 1.001 * radius )
        return false ;

    if( p.getY() < circumCenter.getY() - 1.001 * radius )
        return false ;

    if( p.getZ() > circumCenter.getZ() + 1.001 * radius )
        return false ;

    if( p.getZ() < circumCenter.getZ() - 1.001 * radius )
        return false ;

    double d = dist( circumCenter*tree->getInternalScale(), p*tree->getInternalScale() ) ;
    return  std::abs( d - radius*tree->getInternalScale() ) < POINT_TOLERANCE * radius*tree->getInternalScale() ;
}

bool DelaunayTetrahedron::isNeighbour( const DelaunayTreeItem3D *t ) const
{

    size_t cv = numberOfCommonVertices( t ) ;

    return ( cv == 3 );

}

bool DelaunayTreeItem3D::isDuplicate( const DelaunayTreeItem3D *t ) const
{

    if( t == this )
        return false;

    if( !isTetrahedron || isDeadTetrahedron )
        return false ;

    if( !t->isTetrahedron || t->isDeadTetrahedron )
        return false ;

    if( squareDist3D( static_cast<const DelaunayTetrahedron *>( t )->getCenter()*tree->getInternalScale(),
                      static_cast<const DelaunayTetrahedron *>( this )->getCenter()*tree->getInternalScale() ) > POINT_TOLERANCE*tree->getInternalScale() )
        return false ;

    if( squareDist3D( static_cast<const DelaunayTetrahedron *>( t )->getCircumCenter()*tree->getInternalScale() ,
                      static_cast<const DelaunayTetrahedron *>( this )->getCircumCenter()*tree->getInternalScale() )  > POINT_TOLERANCE*tree->getInternalScale() )
        return false ;

    if( !std::binary_search( &first, &first + 4, t->first ) )
    {
        return false ;
    }

    if( !std::binary_search( &first, &first + 4, t->second ) )
    {
        return false ;
    }

    if( !std::binary_search( &first, &first + 4, t->third ) )
    {
        return false ;
    }

    if( !std::binary_search( &first, &first + 4, t->fourth ) )
    {
        return false ;
    }

    return true ;
}

void DelaunayTetrahedron::insert(std::valarray<bool> & visited, std::vector<DelaunayTreeItem3D *> & ret, Point *p,  Star3D *s )
{
    if( visited[index] )
        return ;

    visited[index] = true ;

    for( size_t i = 0 ; i < 4 ; i++ )
    {
        bool ins =  !visited[neighbour[i]] && ( !getNeighbour( i )->inCircumSphere( *p ) || getNeighbour( i )->onCircumSphere( *p ) ) ;

        if( ins )
        {

            std::vector< Point *> pp = commonSurface( getNeighbour( i ) ) ;

            if( Tetrahedron( *p*tree->getInternalScale(),
                             *pp[0]*tree->getInternalScale(),
                             *pp[1]*tree->getInternalScale(),
                             *pp[2]*tree->getInternalScale() ).volume() > POINT_TOLERANCE*tree->getInternalScale() )
            {
                DelaunayTetrahedron *ss = new DelaunayTetrahedron( tree, this, p, pp[0], pp[1], pp[2], p ) ;

                addSon( ss ) ;
                getNeighbour( i )->addStepson( ss ) ;
                ret.push_back( ss ) ;

            }
        }
    }

// 	s->updateNeighbourhood() ;


}

void DelaunayTetrahedron::print() const
{

    std::cout << "(" << first->getX() << ", " << first->getY() << ", " << first->getZ() << ", " << first->getId() << ") " ;
    std::cout << "(" << second->getX() << ", " << second->getY() << ", " << second->getZ() << ", " << second->getId() << ") " ;
    std::cout << "(" << third->getX() << ", " << third->getY() << ", " << third->getZ() << ", " << third->getId() << ") " ;
    std::cout << "(" << fourth->getX() << ", " << fourth->getY() << ", " << fourth->getZ() << ", " << fourth->getId() << ") " ;
    std::cout <<  ":: " << !dead << std::endl ;
}

DelaunayDemiSpace::DelaunayDemiSpace( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t, DelaunayTreeItem3D *father,  Point   *_one,  Point   *_two, Point   *_three, Point   *p,  Point *c ) : DelaunayTreeItem3D( t, father, c )
{
    second  = _two ;
    first = _one ;
    third = _three ;
    std::sort( &first, &first + 3 ) ;
    fourth = p ;
    dead = false ;
    pseudonormal = ( ( *second*tree->getInternalScale() ) - ( *first*tree->getInternalScale() ) ) ^( ( *third*tree->getInternalScale() ) - ( *first*tree->getInternalScale() ) );
    pseudonormal /= pseudonormal.norm() ;
// 	if(pseudonormal*(*first-*p) < 0)
// 	{
// 		pseudonormal = vector2^vector1 ;
// 	}
    isSpace = true ;
    isTetrahedron = false ;
    isDeadTetrahedron = false;


}

bool DelaunayDemiSpace::inCircumSphere( const Point &p ) const
{
    double planeConst = first->getX() * tree->getInternalScale() * pseudonormal.getX() +
                        first->getY() * tree->getInternalScale() * pseudonormal.getY() +
                        first->getZ() * tree->getInternalScale() * pseudonormal.getZ() ;
    double signedDistP = p.getX() * tree->getInternalScale() * pseudonormal.getX() +
                         p.getY() * tree->getInternalScale() * pseudonormal.getY() +
                         p.getZ() * tree->getInternalScale() * pseudonormal.getZ() - planeConst;
    double signedDistF = fourth->getX() * tree->getInternalScale() * pseudonormal.getX() +
                         fourth->getY() * tree->getInternalScale()* pseudonormal.getY() +
                         fourth->getZ() * tree->getInternalScale() * pseudonormal.getZ() - planeConst;
    return signedDistF * signedDistP < POINT_TOLERANCE * POINT_TOLERANCE || isCoplanar( &p, first, second, third,tree->getInternalScale() );

}

bool DelaunayDemiSpace::onCircumSphere( const Point &p ) const
{
    double planeConst = first->getX()*tree->getInternalScale() * pseudonormal.getX() +
                        first->getY()*tree->getInternalScale() * pseudonormal.getY() +
                        first->getZ()*tree->getInternalScale() * pseudonormal.getZ() ;
    double signedDistP = p.getX() * pseudonormal.getX() * tree->getInternalScale() +
                         p.getY() * tree->getInternalScale() * pseudonormal.getY() +
                         p.getZ() * tree->getInternalScale() * pseudonormal.getZ() - planeConst;

    return std::abs( signedDistP ) < POINT_TOLERANCE*tree->getInternalScale() ;

}


bool  DelaunayDemiSpace::isNeighbour( const DelaunayTreeItem3D *t )  const
{

    if( t->isTetrahedron )
    {
        return numberOfCommonVertices( t ) == 3;
    }

    return numberOfCommonVertices( t ) >= 2 ;
}

void DelaunayDemiSpace::insert(std::valarray<bool> & visited, std::vector<DelaunayTreeItem3D *> & ret, Point *p, Star3D *s )
{

    if( visited[index] )
        return ;

    visited[index] = true ;

    for( size_t i = 0 ; i < neighbour.size() ; i++ )
    {
        std::vector< Point *> pp = getNeighbour( i )->commonSurface( this ) ;

        if( !visited[neighbour[i]] && ( !getNeighbour( i )->inCircumSphere( *p ) || getNeighbour( i )->onCircumSphere( *p ) ) )
        {
            if( !isCoplanar( p, pp[0], pp[1] , pp[2],tree->getInternalScale() ) )
            {

                DelaunayTetrahedron *ss = new DelaunayTetrahedron( tree, this, p, pp[0], pp[1] , pp[2], p ) ;

                addSon( ss ) ;
                getNeighbour( i )->addStepson( ss ) ;

                ret.push_back( ss ) ;

                DelaunayDemiSpace *p0 = new DelaunayDemiSpace( tree, this, pp[0], p, pp[1], pp[2], p ) ;
                DelaunayDemiSpace *p1 = new DelaunayDemiSpace( tree, this, pp[0], p, pp[2], pp[1], p ) ;
                DelaunayDemiSpace *p2 = new DelaunayDemiSpace( tree, this, pp[2], p, pp[1], pp[0], p ) ;

                addSon( p0 ) ;
                addSon( p1 ) ;
                addSon( p2 ) ;

                ret.push_back( p0 ) ;
                ret.push_back( p1 ) ;
                ret.push_back( p2 ) ;
            }
        }
    }

// 	s->updateNeighbourhood() ;


    return  ;
}

void DelaunayDemiSpace::print() const
{
    std::cout << "###############(" << first->getX() << ", " << first->getY() << ", " << first->getZ() << ") (" <<
              second->getX() << ", " << second->getY() << ", " << second->getZ() << ") (" <<
              third->getX() << ", " << third->getY() << ", " << third->getZ() << ")" << "X (" << fourth->getX() << ", " << fourth->getY() << ", " << fourth->getZ() << ") :: " << !dead  << std::endl ;
}

void makeNeighbours3d( DelaunayTreeItem3D *t0, DelaunayTreeItem3D *t1 )
{
    if( t0 == t1 || t0->dead || t1->dead )
        return ;

    if( t0->isSpace && t1->isSpace )
        return ;

    if( !t1->isNeighbour( t0 ) )
    {
        return ;
    }

    if( std::find( &t0->neighbour[0], &t0->neighbour[t0->neighbour.size()], t1->index ) ==  &t0->neighbour[t0->neighbour.size()] )
    {
        std::valarray<unsigned int>  newneighbours( t0->neighbour.size() + 1 ) ;
        std::copy( &t0->neighbour[0], &t0->neighbour[t0->neighbour.size()], &newneighbours[0] ) ;
        newneighbours[t0->neighbour.size()] = t1->index ;
        t0->neighbour.resize( t0->neighbour.size() + 1 ) ;
        t0->neighbour = newneighbours ;
    }

    if( std::find( &t1->neighbour[0], &t1->neighbour[t1->neighbour.size()], t0->index ) ==  &t1->neighbour[t1->neighbour.size()] )
    {
        std::valarray<unsigned int>  newneighbours( t1->neighbour.size() + 1 ) ;
        std::copy( &t1->neighbour[0], &t1->neighbour[t1->neighbour.size()], &newneighbours[0] ) ;
        newneighbours[t1->neighbour.size()] = t0->index ;
        t1->neighbour.resize( t1->neighbour.size() + 1 ) ;
        t1->neighbour = newneighbours ;
    }

}

void updateNeighbours( std::vector<DelaunayTreeItem3D *> * t )
{
    for( size_t i = 0 ; i < t->size() ; i++ )
    {
        for( size_t j = i + 1 ; j < t->size() ; j++ )
        {
            makeNeighbours3d( ( *t )[i], ( *t )[j] ) ;
        }
    }
}

DelaunayRoot3D::DelaunayRoot3D( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *t, Point *p0, Point *p1, Point *p2, Point *p3 ) : DelaunayTreeItem3D( t, nullptr, nullptr )
{
    isSpace = false ;
    isTetrahedron = false ;
    isDeadTetrahedron = false ;
    this->father = -1 ;
    DelaunayTetrahedron *tet = new DelaunayTetrahedron( t, this, p0, p1, p2, p3, nullptr ) ;
    DelaunayDemiSpace *pl0 = new DelaunayDemiSpace( t, this, p0, p1, p2, p3, nullptr );
    DelaunayDemiSpace *pl1 = new DelaunayDemiSpace( t, this, p0, p1, p3, p2, nullptr );
    DelaunayDemiSpace *pl2 = new DelaunayDemiSpace( t, this, p1, p2, p3, p0, nullptr );
    DelaunayDemiSpace *pl3 = new DelaunayDemiSpace( t, this, p2, p3, p0, p1, nullptr );

    makeNeighbours3d( tet, pl0 ) ;
    makeNeighbours3d( tet, pl1 ) ;
    makeNeighbours3d( tet, pl2 ) ;
    makeNeighbours3d( tet, pl3 ) ;
    makeNeighbours3d( pl1, pl0 ) ;
    makeNeighbours3d( pl2, pl1 ) ;
    makeNeighbours3d( pl0, pl2 ) ;
    makeNeighbours3d( pl3, pl0 ) ;
    makeNeighbours3d( pl3, pl1 ) ;
    makeNeighbours3d( pl3, pl2 ) ;

    addSon( tet ) ;
    addSon( pl0 ) ;
    addSon( pl1 ) ;
    addSon( pl2 ) ;
    addSon( pl3 ) ;


    kill( p0 ) ;
}

void DelaunayRoot3D::print() const
{
    std::cout << "I am root !" << std::endl ;
}

bool DelaunayRoot3D::isVertex( const Point *p ) const
{
    return false ;
}

bool DelaunayRoot3D::inCircumSphere( const Point &p ) const
{
    return true ;
}

bool DelaunayRoot3D::onCircumSphere( const Point &p ) const
{
    return true ;
}


void DelaunayRoot3D::insert(std::valarray<bool> & visited, std::vector<DelaunayTreeItem3D *>& ret, Point *p, Star3D *s )
{

    for( size_t i  = 0 ;  i < son.size() ; i++ )
    {

        std::vector<DelaunayTreeItem3D *> temp ;
        getSon( i )->insert(visited, temp, p, s ) ;

        for( size_t j = 0 ; j < temp.size() ; j++ )
        {
            ret.push_back( temp[j] ) ;
        }
    }

    updateNeighbours( &ret ) ;
}

void DelaunayRoot3D::conflicts(std::valarray<bool> & visited,  std::vector< DelaunayTreeItem3D * > &ret, const Amie::Geometry *g )
{
    visited[index] = true ;

    for( size_t i  = 0 ;  i < son.size() ; i++ )
    {
        std::vector<DelaunayTreeItem3D *> toTest ;
        getSon( i )->flatConflicts(visited, toTest, ret, g ) ;

        while( !toTest.empty() )
        {
            std::vector<DelaunayTreeItem3D *> tempToTest ;

            for( size_t j  = 0 ;  j < toTest.size() ; j++ )
            {
                toTest[j]->flatConflicts(visited, tempToTest, ret, g ) ;
            }

            toTest = tempToTest ;
        }

    }


}

void DelaunayRoot3D::conflicts(std::valarray<bool> & visited, std::vector< DelaunayTreeItem3D * > & ret, const Amie::Point *p )
{

    visited[index] = true ;

    for (size_t i  = 0 ;  i < 5 ; i++)
    {
        std::vector<DelaunayTreeItem3D *> toTest ;

        getSon(i)->flatConflicts(visited,toTest,ret,p) ;

        while(!toTest.empty())
        {
            std::vector<DelaunayTreeItem3D *> tempToTest ;
            for(size_t j  = 0 ;  j < toTest.size() ; j++)
            {
                toTest[j]->flatConflicts(visited,tempToTest,ret,p) ;
            }
            toTest = tempToTest ;
        }

    }

}

Star3D::Star3D( std::vector<DelaunayTreeItem3D *> *t, const Point *p ) :  treeitem( *t ), creator( p )
{
    for( size_t i = 0 ; i < t->size() ; i++ )
    {
        if( ( *t )[i]->isTetrahedron )
        {
            this->edge.push_back( ( *t )[i]->first ) ;
            this->edge.push_back( ( *t )[i]->second ) ;
            this->edge.push_back( ( *t )[i]->third ) ;
            this->edge.push_back( ( *t )[i]->fourth );
        }
        else if( ( *t )[i]->isSpace )
        {
            this->edge.push_back( ( *t )[i]->first ) ;
            this->edge.push_back( ( *t )[i]->second ) ;
            this->edge.push_back( ( *t )[i]->third );
        }
    }

    std::sort( edge.begin(), edge.end() ) ;
    auto e = std::unique( edge.begin(), edge.end() ) ;
    edge.erase( e, edge.end() ) ;
}

size_t Star3D::size()
{
    return edge.size() ;
}

const Point *Star3D::getEdge( size_t i ) const
{
    assert( i < edge.size() ) ;
    return this->edge[i] ;
}

DelaunayTree3D::DelaunayTree3D( Point *p0, Point *p1, Point *p2, Point *p3 ) : Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D>(SPACE_THREE_DIMENSIONAL)
{
    fixedScale = false ;
    global_counter = 4;
    neighbourhood = false ;
    double maxDist = std::min(dist(*p0, *p1),dist(*p0, *p2)) ;
    maxDist = std::min(maxDist,dist(*p0, *p3)) ;
    maxDist = std::min(maxDist,dist(*p1, *p2)) ;
    maxDist = std::min(maxDist,dist(*p1, *p3)) ;
    maxDist = std::min(maxDist,dist(*p2, *p3)) ;

    internalScale = 4e12/maxDist ;

    p0->setId(0)  ;
    p1->setId(1)  ;
    p2->setId(2)  ;
    p3->setId(3) ;

// 	*p0 *= internalScale ;
// 	*p1 *= internalScale ;
// 	*p2 *= internalScale ;
// 	*p3 *= internalScale ;

    DelaunayRoot3D *root = new DelaunayRoot3D( this, p0, p1, p2, p3 ) ;

    space.push_back( static_cast<DelaunayDemiSpace *>( root->getSon( 1 ) ) ) ;
    space.push_back( static_cast<DelaunayDemiSpace *>( root->getSon( 2 ) ) ) ;
    space.push_back( static_cast<DelaunayDemiSpace *>( root->getSon( 3 ) ) ) ;
    space.push_back( static_cast<DelaunayDemiSpace *>( root->getSon( 4 ) ) ) ;


}

DelaunayTree3D::~DelaunayTree3D()
{

    std::valarray<Point *> nularray( 0 ) ;

    for( size_t i = 0 ;  i < this->tree.size() ; i++ )
    {

        if( this->tree[i] && this->tree[i]->isTetrahedron && !this->tree[i]->isDeadTetrahedron )
        {
            DelaunayTetrahedron *t = dynamic_cast<DelaunayTetrahedron *>( tree[i] ) ;

            t->setBoundingPoints( nularray ) ;
        }

    }

    for( size_t i = 0 ;  i < this->tree.size() ; i++ )
    {
        delete this->tree[i] ;
    }

    for( size_t i = 0 ; i < additionalPoints.size() ; i++ )
        delete additionalPoints[i] ;

}

std::vector<DelaunayTreeItem3D *> DelaunayTree3D::addElements(std::vector<DelaunayTreeItem3D *> & cons, Point * p)
{

    p->getId() = global_counter++ ;
    std::valarray<bool> visited(false, tree.size()) ;

    Star3D *s = new Star3D( &cons, p ) ;

    std::vector<DelaunayTreeItem3D *> ret ;

    for( auto & item : cons )
    {
        if( !item->onCircumSphere( *p ) )
        {
            std::vector<DelaunayTreeItem3D *> temp ;
            item->insert(visited, temp, p, s ) ;
            ret.insert( ret.end(), temp.begin(), temp.end() ) ;
        }
    }

    s->updateNeighbourhood() ;
    bool weGotSpaces = false ;


    for( const auto & item : ret)
    {
        if( !item->dead && item->isSpace )
        {
            weGotSpaces = true ;
            space.push_back( static_cast<DelaunayDemiSpace *>( item ) ) ;
        }
    }

    if( weGotSpaces )
    {
        for( size_t k = 0 ; k < space.size() - 1 ; k++ )
        {
            for( size_t j = k + 1 ; j < space.size() ; j++ )
            {
                space[j]->merge( space[k] ) ;
            }
        }
    }

    for( auto & item : ret )
    {
        for( auto & s : space )
        {
            makeNeighbours3d( item, s );
        }
    }

    for( size_t i = 0 ; i < cons.size() ; i++ )
    {
        if( !cons[i]->onCircumSphere( *p ) )
            cons[i]->kill( p ) ;
    }

    bool correct = true ;

    for( const auto & item : ret )
    {
        if( !item->erased && ( ( !item->dead && item->isTetrahedron ) || item->isSpace ) )
        {
// #ifndef NDEBUG
            if( item->isTetrahedron && item->neighbour.size() != 4 )
            {

                std::cout << "we have " << item->neighbour.size() << " neighbours" << std::endl ;
                p->print() ;
                tree[item->father]->print() ;
                item->print() ;

                for( size_t k = 0 ; k < item->neighbour.size() ;  k++ )
                {
                    std::cout << "   --> " << std::flush ;
                    item->getNeighbour( k )->print() ;
                }

                correct = false ;
            }

// #endif
        }
        else
        {

            std::valarray<Point *> nullarray( 0 ) ;

            if( item->isTetrahedron )
                static_cast<DelaunayTetrahedron *>( item )->setBoundingPoints( nullarray ) ;

            tree[item->index] = nullptr ;
            delete item ;
        }
    }

    for( auto i = space.begin() ; i != space.end() ; i++ )
    {
        if( ( *i )->dead )
        {
            space.erase( i ) ;
            i-- ;
        }
    }

    if( !correct )
    {
        std::cout << "inconsistent state, will crash soon" << std::endl ;
// 		print();
        exit( 0 ) ;
    }

    for( const auto & item : cons )
    {

        if( item->dead && item->isTetrahedron && !item->isDeadTetrahedron )
        {
            DelaunayDeadTetrahedron *dt = new DelaunayDeadTetrahedron( static_cast<DelaunayTetrahedron *>( item ) ) ;
            tree[item->index] = dt ;
            delete item ;
        }
    }

// //
    delete s ;
    return ret ;
}

void DelaunayTree3D::insert( Point *p , double minDist)
{
    neighbourhood = false ;
    std::vector<DelaunayTreeItem3D *> cons = conflicts( p ) ;
   
    for( size_t i = 0 ; i < cons.size() ; i++ )
    {
        if( cons[i]->isVertex( p ) )
        {

            return ;
        }
    }

    addElements(cons, p) ;

}

Star3D::~Star3D()
{
}

std::vector<DelaunayTreeItem3D *> DelaunayTree3D::conflicts(const Point *p)
{
    std::vector<DelaunayTreeItem3D *> cons;
    std::valarray<bool> visited(false, tree.size()) ;

    static_cast<DelaunayRoot3D *>(tree[0])->conflicts(visited, cons, p ) ;

    if( cons.empty() )
    {
        for( auto & s : space )
        {

	    std::vector<DelaunayTreeItem3D *> toTest ;

	    s->flatConflicts(visited,toTest,cons,p) ;

	    while(!toTest.empty())
	    {
		std::vector<DelaunayTreeItem3D *> tempToTest ;
		for(size_t j  = 0 ;  j < toTest.size() ; j++)
		{
		    toTest[j]->flatConflicts(visited,tempToTest,cons,p) ;
		}
		toTest = tempToTest ;
	    }
        }
    }

    return cons ;
}

void DelaunayTree3D::purge()
{
    int originalSize = tree.size() ;
    std::cerr << "purging..." << std::flush ;
    std::vector<DelaunayTetrahedron *> tets = getTetrahedrons() ;
    std::valarray<int> indexCorrespondance( tree.size() );
    std::vector<DelaunayTreeItem3D *> newTree ;
    std::vector<DelaunayDemiSpace *> newSpace ;
    newTree.push_back( tree[0] ) ;
    indexCorrespondance[tree[0]->index] = 0 ;
    int counter = 1 ;

    for( size_t i = 0 ; i < tets.size() ; i++ )
    {
        if( i % 1000 == 0 )
            std::cerr << "\r getting elements element " << i << "/" << tets.size() + space.size() << std::flush ;

        newTree.push_back( tets[i] ) ;
        indexCorrespondance[tets[i]->index] = counter++ ;
    }

    for( size_t i = 0 ; i < space.size() ; i++ )
    {
        if( i % 1000 == 0 )
            std::cerr << "\r getting elements element " << i + tets.size() << "/" << tets.size() + space.size() << std::flush ;

        if( !space[i]->dead )
        {
            newTree.push_back( space[i] ) ;
            indexCorrespondance[space[i]->index] = counter++ ;
            newSpace.push_back( space[i] ) ;
        }
    }

    std::cerr << std::endl ;

    tree[0]->neighbour.resize( 0 ) ;
    tree[0]->son.resize( 0 ) ;
    tree[0]->stepson.resize( 0 ) ;
    std::valarray<unsigned int>  newson( newTree.size() - 1 ) ;

    for( size_t i = 1 ; i < newTree.size() ; i++ )
    {
        if( i % 1000 == 0 )
            std::cerr << "\r reconstructing element " << i << "/" << newTree.size() << std::flush ;

        for( size_t j = 0 ; j < newTree[i]->neighbour.size() ; j++ )
        {
            newTree[i]->neighbour[j] = indexCorrespondance[newTree[i]->neighbour[j]] ;
        }

        for( size_t j = 0 ; j < newTree[i]->son.size() ; j++ )
        {
            newTree[i]->son[j] = indexCorrespondance[newTree[i]->son[j]] ;
        }

        for( size_t j = 0 ; j < newTree[i]->stepson.size() ; j++ )
        {
            newTree[i]->stepson[j] = indexCorrespondance[newTree[i]->stepson[j]] ;
        }

        newTree[i]->father = 0 ;
        newTree[i]->stepfather = -1 ;
        newTree[i]->index = i ;
        newson[i - 1] = i ;
    }

    tree[0]->son.resize( newson.size() ) ;
    std::copy( &tree[0]->son[0], &tree[0]->son[tree[0]->son.size()], &newson[0] ) ;

    std::cerr << std::endl ;

    for( size_t i = 1 ; i < tree.size() ; i++ )
    {
        if( !!tree[i]->dead )
            delete tree[i] ;
    }

    tree = newTree ;
    space = newSpace ;
    std::cerr << " ...done. " << tree.size() << " elements kept from " << originalSize << "." << std::endl ;
}

std::vector< DelaunayTetrahedron* > DelaunayTree3D::getConflictingElements(const Point* p)
{
    std::vector< DelaunayTreeItem3D* > targets = conflicts(p) ;
    std::vector<DelaunayTetrahedron*> ret ;
    for(size_t i = 0 ; i < targets.size() ; i++)
    {
        if(targets[i]->isTetrahedron && !targets[i]->isDeadTetrahedron )
            ret.push_back(static_cast< DelaunayTetrahedron *> (targets[i])) ;
    }
    return ret ;
}

std::vector<DelaunayTetrahedron *> DelaunayTree3D::getNeighbourhood(DelaunayTetrahedron * element) const
{
    std::vector<DelaunayTetrahedron *> ret ;
    for(const auto & idx : element->neighbourhood)
    {
        ret.push_back((DelaunayTetrahedron *)tree[idx]);
    }
    return ret ;
}

std::vector<DelaunayTetrahedron *> DelaunayTree3D::conflicts(const Geometry *g )
{

    std::vector<DelaunayTreeItem3D *> cons = conflicts( &g->getCenter() ) ;
    DelaunayTetrahedron *origin = nullptr ;

    for( const auto & item : cons )
    {
        if( item->isTetrahedron && static_cast<DelaunayTetrahedron *>( item )->normalisedIn( g->getCenter() ) )
        {
            origin = static_cast<DelaunayTetrahedron *>( item ) ;
            break ;
        }
    }

    if(origin)
    {
        
        std::vector<DelaunayTetrahedron *> ret = getNeighbouringElementsInGeometry( origin, g ) ;

        return ret ;
    }

    return std::vector<DelaunayTetrahedron *>() ;

}


std::vector<DelaunayDemiSpace *>  DelaunayTree3D::getConvexHull()
{
    std::vector<DelaunayDemiSpace *> ret ;

    for( const auto & s : space )
    {
        if( !s->dead )
            ret.push_back( s ) ;
    }

    return ret ;
}

void DelaunayTree3D::buildNeighbourhoods()
{
    std::valarray<bool> visited(false, size()) ;

    if(!neighbourhood)
    {
        neighbourhood = true ;
//      std::cerr << "\r building neighbourhood... element 0/" << ret.size() << std::flush ;
        for( size_t i = 0 ; i < tree.size() ; i++)
        {
            if(tree[i]->isTetrahedron && !tree[i]->dead)
                static_cast<DelaunayTetrahedron *>(tree[i])->neighbourhood.resize(0) ;
        }



        for( size_t i = 0 ; i < tree.size() ; i++)
        {
            if(!tree[i]->isTetrahedron || tree[i]->dead)
                continue ;
//          if(i%100 == 0)
//              std::cerr << "\r building neighbourhood... element " << i <<"/" << ret.size() << std::flush ;

            std::vector<DelaunayTetrahedron *> tocheck ;
            std::vector<DelaunayTetrahedron *> toclean ;
            for(size_t j = 0 ; j < tree[i]->neighbour.size() ; j++)
            {
                if(tree[i]->getNeighbour(j)->isTetrahedron && !visited[tree[i]->neighbour[j]])
                {
                    tocheck.push_back(static_cast<DelaunayTetrahedron *>(tree[i]->getNeighbour(j)));
                    visited[tree[i]->neighbour[j]] = true ;
                    toclean.push_back(*tocheck.rbegin()) ;
                    static_cast<DelaunayTetrahedron *>(tree[i])->addNeighbourhood(*tocheck.rbegin()) ;
                }
            }

            while(!tocheck.empty())
            {
                std::vector<DelaunayTetrahedron *> tocheck_temp ;
                for(size_t k = 0 ; k < tocheck.size() ; k++)
                {
                    for(size_t j = 0 ; j < tocheck[k]->neighbour.size() ; j++)
                    {
                        if(
                            tocheck[k]->getNeighbour(j)->isTetrahedron
                            && !visited[tocheck[k]->neighbour[j]]
                            && tocheck[k]->getNeighbour(j) != tree[i]
                            && static_cast<DelaunayTetrahedron *>(tocheck[k]->getNeighbour(j))->numberOfCommonVertices(tree[i])
                        )
                        {
                            tocheck_temp.push_back(static_cast<DelaunayTetrahedron *>(tocheck[k]->getNeighbour(j)));
                            visited[tocheck[k]->neighbour[j]] = true ;
                            toclean.push_back(*tocheck_temp.rbegin()) ;
                            static_cast<DelaunayTetrahedron *>(tree[i])->addNeighbourhood(*tocheck_temp.rbegin()) ;
                        }
                    }
                }

                tocheck = tocheck_temp ;
            }
            for(size_t j = 0 ; j < toclean.size() ; j++)
            {
                visited[toclean[j]->index] = false ;
            }

        }
    }

}

std::vector<DelaunayTetrahedron *>  DelaunayTree3D::getTetrahedrons( bool buildNeighbourhood ) 
{
    std::vector<DelaunayTetrahedron *> ret;
    std::valarray<bool> visited(false, tree.size()) ;
    for( const auto & item : tree )
    {
        if( !item->dead && item->isTetrahedron && !item->isDeadTetrahedron)
        {
            ret.push_back( ( DelaunayTetrahedron * )( item ) ) ;
        }
    }
    buildNeighbourhoods();


    return ret ;
}

void DelaunayTree3D::print() const
{
    size_t alive = 0 ;
    std::cout << "we have a total of " << tree.size() << " elements" << std::endl ;

    for( size_t i = 0 ; i < tree.size() ; i++ )
    {
        if( !tree[i]->dead )
        {
            alive++ ;

        }
    }

    std::cout << "of which " << alive << "are alive" << std::endl ;
#ifdef DEBUG

    for( size_t i = 0 ; i < tree.size() ; i++ )
    {
        if( !tree[i]->dead )
        {
            tree[i]->print() ;

            for( size_t j = 0 ; j < tree[i]->neighbour.size() ; j++ )
            {
                std::cout << "\t --> " ;
                tree[i]->getNeighbour( j )->print() ;
            }
        }
    }

#endif
}

void DelaunayTetrahedron::clearElementaryMatrix()
{
    for( size_t i = 0 ; i <  cachedElementaryMatrix.size() ; i++ )
    {
        cachedElementaryMatrix[i].resize( 0 ) ;
    }

    cachedElementaryMatrix.resize( 0 ) ;
}

std::valarray<std::valarray<Matrix> > DelaunayTetrahedron::getTangentElementaryMatrix(VirtualMachine * vm)
{
    return cachedElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > & DelaunayTetrahedron::getElementaryMatrix(VirtualMachine * vm )
{
/*    delete cachedGps ;
    cachedGps = nullptr ;
    behaviourUpdated = true */;
    if( !behaviourUpdated && !enrichmentUpdated && cachedElementaryMatrix.size() )
    {
        return cachedElementaryMatrix ;
    }
    needAssembly = true ;
    clearBoundaryConditions() ;

    size_t dofCount = getShapeFunctions().size() + getEnrichmentFunctions().size() ;
    size_t ndofs = getBehaviour()->getNumberOfDegreesOfFreedom() ;

    if(        cachedElementaryMatrix.size() == 0
            || cachedElementaryMatrix.size() != dofCount
            ||  (cachedElementaryMatrix.size() && cachedElementaryMatrix[0].size() != dofCount))
    {

        std::valarray< Matrix > v_j(Matrix(ndofs, ndofs), dofCount) ;
        cachedElementaryMatrix.resize(dofCount,v_j) ;
        getSubTriangulatedGaussPoints() ;
    }


    getState().updateInverseJacobianCache(Point( .25,.25,.25)) ;

    
    std::valarray<Matrix> Jinv((*getState().JinvCache),  getGaussPoints().gaussPoints.size()) ;
    if(moved)
    {
        for(size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++)
        {
            getState().getInverseJacobianMatrix( getGaussPoints().gaussPoints[i].first, Jinv[i]) ;
        }
    }

    size_t start = 0 ;
    size_t startEnriched = 0 ;
    if(timePlanes() > 1)
    {
        start = getShapeFunctions().size() -  getShapeFunctions().size()/timePlanes() ;
        startEnriched = getEnrichmentFunctions().size() -  getEnrichmentFunctions().size()/timePlanes() ;
    }
    bool cleanup = false ;
    if(!vm)
    {
       vm = new VirtualMachine();
       cleanup = true ;
    }

    if(behaviour->isSymmetric())
    {

        for(size_t i = start ; i < getShapeFunctions().size() ; i++)
        {
            behaviour->apply(getShapeFunction(i),getShapeFunction(i), getGaussPoints(), Jinv, cachedElementaryMatrix[i][i], vm) ;

            for(size_t j = 0 ; j < i ; j++)
            {
                behaviour->apply(getShapeFunction(i), getShapeFunction(j), getGaussPoints(), Jinv,cachedElementaryMatrix[i][j], vm) ;
                cachedElementaryMatrix[i][j].transpose(cachedElementaryMatrix[j][i]) ;
            }

            for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j), getGaussPoints(), Jinv,cachedElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                cachedElementaryMatrix[i][j+getShapeFunctions().size()].transpose(cachedElementaryMatrix[j+getShapeFunctions().size()][i]) ;
            }
        }


        for(size_t i = startEnriched ; i < getEnrichmentFunctions().size() ; i++)
        {
            behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i), getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = 0 ; j < i ; j++)
            {
                behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j), getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()].transpose(cachedElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()]) ;
            }
        }
// 		if(getEnrichmentFunctions().size() == 6)
// 			exit(0) ;
    }
    else
    {

        for(size_t i = start; i < getShapeFunctions().size() ; i++)
        {

            for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
            {
                behaviour->apply(getShapeFunction(i), getShapeFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i][j], vm) ;
                behaviour->apply(getShapeFunction(j), getShapeFunction(i),getGaussPoints(), Jinv,cachedElementaryMatrix[j][i], vm) ;
            }

            for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                behaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),getGaussPoints(), Jinv,cachedElementaryMatrix[j+getShapeFunctions().size()][i], vm) ;
            }
        }


        for(size_t i = startEnriched ; i < getEnrichmentFunctions().size() ; i++)
        {
            behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = i+1 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                behaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;
            }
        }
    }



    enrichmentUpdated = false ;
    behaviourUpdated = false ;

    if( behaviour->hasInducedForces() )
        cachedForces.resize( 0 ) ;

// 	if(getEnrichmentFunctions().size())
// 	{
// 		for( size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++ )
// 			getGaussPoints().gaussPoints[i].first.print() ;
//
// 		for( size_t i = 0 ; i < getShapeFunctions().size() ; i++ )
// 			cachedElementaryMatrix[i][i].print() ;
// 		for( size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++ )
// 			cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()].print() ;
// 		exit(0) ;
// 	}

// 	exit(0) ;
    if(cleanup)
        delete vm ;
    return cachedElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > & DelaunayTetrahedron::getViscousElementaryMatrix(VirtualMachine * vm)
{
 size_t dofCount = getShapeFunctions().size()+getEnrichmentFunctions().size() ;


    if( !behaviourUpdated && !enrichmentUpdated && cachedViscousElementaryMatrix.size() && cachedViscousElementaryMatrix[0].size() == dofCount)
    {
        return cachedViscousElementaryMatrix ;
    }
    size_t ndofs = getBehaviour()->getNumberOfDegreesOfFreedom() ;


    if( cachedViscousElementaryMatrix.size() == 0
            || cachedViscousElementaryMatrix.size() != dofCount
            ||  (cachedViscousElementaryMatrix.size() && cachedViscousElementaryMatrix[0].size() != dofCount))
    {

        std::valarray< Matrix > v_j(Matrix(ndofs, ndofs), dofCount) ;
        cachedViscousElementaryMatrix.resize(dofCount,v_j) ;
        getSubTriangulatedGaussPoints() ;
    }

    getState().updateInverseJacobianCache(Point( .25,.25,.25)) ;
    
    std::valarray<Matrix> Jinv((*getState().JinvCache),  getGaussPoints().gaussPoints.size()) ;
    if(moved)
    {
        for(size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++)
        {
            getState().getInverseJacobianMatrix( getGaussPoints().gaussPoints[i].first, Jinv[i]) ;
        }
    }

    size_t start = 0 ;
    size_t startEnriched = 0 ;
    if(timePlanes() > 1)
    {
        start = getShapeFunctions().size() -  getShapeFunctions().size()/timePlanes() ;
        startEnriched = getEnrichmentFunctions().size() -  getEnrichmentFunctions().size()/timePlanes() ;
    }
    
    bool cleanup = false ;
    if(!vm)
    {
        vm = new VirtualMachine() ;
        cleanup = true ;
    }
    
    if(behaviour->isSymmetric())
    {
//      #pragma omp parallel for
        for(size_t i = start ; i < getShapeFunctions().size() ; i++)
        {

            
            behaviour->applyViscous(getShapeFunction(i), getShapeFunction(i),getGaussPoints(), Jinv, cachedViscousElementaryMatrix[i][i], vm) ;

            for(size_t j = 0 ; j < i ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getShapeFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j], vm) ;
                cachedViscousElementaryMatrix[i][j].transpose(cachedViscousElementaryMatrix[j][i]) ;
            }
            for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                cachedViscousElementaryMatrix[i][j+getShapeFunctions().size()].transpose(cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i]) ;
            }
        }
//      #pragma omp parallel for
        for(size_t i = startEnriched ; i < getEnrichmentFunctions().size() ; i++)
        {

            behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = 0 ; j <  i ; j++)
            {
                behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                cachedViscousElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()].transpose(cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()]) ;
            }
        }
    }
    else
    {
//      #pragma omp parallel for
        for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
        {
            behaviour->applyViscous(getShapeFunction(i), getShapeFunction(i),getGaussPoints(), Jinv, cachedViscousElementaryMatrix[i][i], vm) ;

            for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getShapeFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j], vm) ;
                behaviour->applyViscous(getShapeFunction(j), getShapeFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[j][i], vm) ;
            }
            for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                behaviour->applyViscous(getEnrichmentFunction(j), getShapeFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i], vm) ;
            }
        }

//      #pragma omp parallel for
        for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
        {
            behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = i+1 ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                behaviour->applyViscous(getEnrichmentFunction(j), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;
            }
        }
    }

    enrichmentUpdated = false ;
    behaviourUpdated = false ;

    if(behaviour->hasInducedForces())
        cachedForces.resize(0) ;

    if(cleanup)
        delete vm ;
    
    return cachedViscousElementaryMatrix ;
}

void DelaunayTree3D::extrude(double dt)
{
    std::map<Point *, Point *> points ;

    std::vector<DelaunayTetrahedron *> tet = getTetrahedrons() ;
    double beginning = tet[0]->getBoundingPoint(0).getT() ;
    double end = tet[0]->getBoundingPoint(0).getT() ;
    for(size_t i = 1 ; i < tet[0]->getBoundingPoints().size() ; i++)
    {
        if(tet[0]->getBoundingPoint(i).getT() < beginning)
            beginning = tet[0]->getBoundingPoint(i).getT() ;
        if(tet[0]->getBoundingPoint(i).getT() > end)
            end = tet[0]->getBoundingPoint(i).getT() ;
    }

    int indexOfLastTimePlane = (tet[0]->timePlanes()-1)*tet[0]->getBoundingPoints().size()/tet[0]->timePlanes() ;

    for(size_t i = 0 ; i < tet.size() ; i++)
    {
        for(size_t j = 0 ; j < tet[i]->getBoundingPoints().size() ; j++)
        {
            Point * next = new Point(tet[i]->getBoundingPoint(j).getX(), tet[i]->getBoundingPoint(j).getY(), tet[i]->getBoundingPoint(j).getZ()) ;
            next->getT() = tet[i]->getBoundingPoint(j).getT() ;
            next->getT() = end + dt * (next->getT() - beginning) / (end - beginning) ;
            bool increment = true ;
            if(next->getT() == end)
            {
                next = &tet[i]->getBoundingPoint(j+indexOfLastTimePlane) ;
                increment = false ;
            }
            if(increment && !points.find(&tet[i]->getBoundingPoint(j))->second)
            {
                next->setId(global_counter++) ;
            }
            points.insert(std::pair<Point *, Point *>(&tet[i]->getBoundingPoint(j), next)) ;
        }
    }

    std::map<Point *, Point *>::iterator finder ;
    for(size_t i = 0 ; i < tet.size() ; i++)
    {

        std::valarray<Point *> newPoints(tet[i]->getBoundingPoints().size()) ;
        for(size_t j = 0 ; j < newPoints.size() ; j++)
        {
            newPoints[j] = points.find(&tet[i]->getBoundingPoint(j))->second ;
//			newPoints[j]->print() ;
        }

//		DelaunayTetrahedron * toInsert = new DelaunayDeadTetrahedron(tet[i]->tree, nullptr, a,b,c,d, a) ;
// 		toInsert->setOrder(tet[i]->getOrder()) ;
// 		toInsert->setBoundingPoints(newPoints) ;
// 		toInsert->setBehaviour(tet[i]->getBehaviour()->getCopy()) ;
    }

}

std::vector<Point *> DelaunayTetrahedron::getIntegrationHints() const
{
    std::vector<Point *> to_add ;
    to_add.push_back( new Point( 0, 1, 0 ) ) ;
    to_add.push_back( new Point( 0, 0, 0 ) ) ;
    to_add.push_back( new Point( 1, 0.0 ) ) ;
    to_add.push_back( new Point( 0, 0, 1 ) ) ;

    TetrahedralElement f( LINEAR ) ;

    for( size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++ )
    {

        for( size_t j = 0 ; j < getEnrichmentFunction( i ).getIntegrationHint().size() ; j++ )
        {
            bool go = true ;

            for( size_t k = 0 ; k < to_add.size()  ; k++ )
            {
                if( squareDist3D( getEnrichmentFunction( i ).getIntegrationHint( j ), *to_add[k] ) < .005 )
                {
                    go = false ;
                    break ;
                }
            }

            if( go && f.in( getEnrichmentFunction( i ).getIntegrationHint( j ) ) )
            {
                to_add.push_back( new Point( getEnrichmentFunction( i ).getIntegrationHint( j ) ) ) ;
            }
        }
    }

    return to_add ;

}


const GaussPointArray &DelaunayTetrahedron::getSubTriangulatedGaussPoints()
{
    if( !enrichmentUpdated )
        return *getCachedGaussPoints() ;

    GaussPointArray gp = getGaussPoints() ;

    size_t numberOfRefinements = 2;

    VirtualMachine vm ;

    if( getEnrichmentFunctions().size() > 0  )
    {
        if( getCachedGaussPoints()->regularGrid )
            return *getCachedGaussPoints() ;

        std::vector<std::pair<Point, double> > gp_alternative ;
        double originalSum = 0 ;
//        double fsum = 0 ;
        for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
            originalSum+=gp.gaussPoints[i].second ;

        if( false )
        {
            TetrahedralElement father(LINEAR) ;

            size_t target = 1024 ;

            while(gp_alternative.size() < target)
            {

                Point test = Point((double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX);

                if( father.in( test ) )
                {
                    gp_alternative.push_back( std::make_pair( test, 1. / 2. ) ) ;
                }

            }

//                 double fsum = 0 ;
//                 for( size_t i = 0 ; i < gp_alternative.size() ; i++ )
//                     fsum += gp_alternative[i].second ;
            for( size_t i = 0 ; i < gp_alternative.size() ; i++ )
                gp_alternative[i].second = originalSum/gp_alternative.size() ;


            delete cachedGps ;
            cachedGps = new GaussPointArray(gp_alternative) ;
            cachedGps->regularGrid = true ;
            return *getCachedGaussPoints() ;
        }

        VirtualMachine vm ;
        std::vector<Point *> to_add = getIntegrationHints();
        std::vector<Point *> pointsToCleanup = to_add;
        std::vector<DelaunayTetrahedron *> triangleToCleanup;
        std::vector<DelaunayTetrahedron *> tri ;

        DelaunayTree3D *dt = new DelaunayTree3D( to_add[0], to_add[1], to_add[2], to_add[3] ) ;
        TetrahedralElement f( LINEAR ) ;

        std::random_shuffle( to_add.begin() + 8, to_add.end() ) ;

        for( size_t i = 4 ; i < to_add.size() ; i++ )
        {
            dt->insert( to_add[i], 0 ) ;
        }

        for( size_t i = 0 ; i < numberOfRefinements ; i++ )
        {
            tri = dt->getTetrahedrons( false ) ;
            std::vector<Point> newPoints ;

            for( size_t j = 0 ; j < tri.size() ; j++ )
            {
                if(i%4 == 0)
                {
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->second )*.5);
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->fourth )*.5);
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->third  )*.5);
                }
                else if(i%4 == 1)
                {
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->second )*.5);
                    newPoints.push_back( ( *tri[j]->second + *tri[j]->third  )*.5);
                    newPoints.push_back( ( *tri[j]->second + *tri[j]->fourth )*.5);
                }
                else if(i%4 == 2)
                {
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->third  )*.5);
                    newPoints.push_back( ( *tri[j]->second + *tri[j]->third  )*.5);
                    newPoints.push_back( ( *tri[j]->third  + *tri[j]->fourth )*.5);
                }
                else if(i%4 == 3)
                {
                    newPoints.push_back( ( *tri[j]->first  + *tri[j]->fourth )*.5);
                    newPoints.push_back( ( *tri[j]->second + *tri[j]->fourth )*.5);
                    newPoints.push_back( ( *tri[j]->third  + *tri[j]->fourth )*.5);
                }
            }

            std::vector<Point *> uniquePoints ;

            if( !newPoints.empty() )
            {
                for( size_t j = 0 ; j < newPoints.size() ; j++ )
                {
                    bool unique  = true ;

                    for( size_t k = 0 ; k < uniquePoints.size() ; k++ )
                    {
                        if( squareDist3D( newPoints[j], *uniquePoints[k] ) < .00125 )
                        {
                            unique  = false ;
                            break ;
                        }
                    }

                    if( unique )
                    {
                        for( size_t k = 0 ; k < to_add.size() ; k++ )
                        {
                            if( squareDist3D( newPoints[j], *to_add[k] ) < .00125 )
                            {
                                unique  = false ;
                                break ;
                            }
                        }
                    }

                    if( unique )
                    {
                        for( size_t k = 0 ; k < pointsToCleanup.size() ; k++ )
                        {
                            if( squareDist3D( newPoints[j], *pointsToCleanup[k] ) < .00125 )
                            {
                                unique  = false ;
                                break ;
                            }
                        }
                    }

                    if( unique )
                    {
                        uniquePoints.push_back( new Point(newPoints[j]) );
                    }
                }
            }
            to_add.insert(to_add.end(), uniquePoints.begin(), uniquePoints.end()) ;
            delete dt ;

            dt = new DelaunayTree3D( to_add[0], to_add[1], to_add[2], to_add[3] ) ;

            std::random_shuffle( to_add.begin() + 8, to_add.end() ) ;

            for( size_t i = 4 ; i < to_add.size() ; i++ )
            {
                dt->insert( to_add[i], 0 ) ;
            }

            for( size_t k = 0 ; k < uniquePoints.size() ; k++ )
            {
                pointsToCleanup.push_back( uniquePoints[k] ) ;
            }
        }

        tri = dt->getTetrahedrons( false ) ;
        double jac = volume() * 6. ;

        Matrix J(3,3) ;
        for( size_t i = 0 ; i < tri.size() ; i++ )
        {
            GaussPointArray gp_temp = tri[i]->getGaussPoints() ;

            if( f.in( tri[i]->getCenter() ) )
            {

                Function x = XTransform( tri[i]->getBoundingPoints(), f.getShapeFunctions() ) ;
                Function y = YTransform( tri[i]->getBoundingPoints(), f.getShapeFunctions() ) ;
                Function z = ZTransform( tri[i]->getBoundingPoints(), f.getShapeFunctions() ) ;

                for( size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++ )
                {

                    gp_temp.gaussPoints[j].first.set( vm.eval( x, gp_temp.gaussPoints[j].first ), vm.eval( y, gp_temp.gaussPoints[j].first ), vm.eval( z, gp_temp.gaussPoints[j].first ) ) ;
                    if( !isFather && isMoved() )
                    {
                        getState().getInverseJacobianMatrix(gp_temp.gaussPoints[j].first, J) ;
                        gp_temp.gaussPoints[j].second *= det(J) ;
                    }
                    else
                        gp_temp.gaussPoints[j].second *= jac ;

                    gp_alternative.push_back( gp_temp.gaussPoints[j] ) ;
                }
            }
        }

        delete dt ;

        for( size_t i = 0 ; i < pointsToCleanup.size() ; i++ )
            delete pointsToCleanup[i] ;

        if( gp.gaussPoints.size() < gp_alternative.size() )
        {
            gp.gaussPoints.resize( gp_alternative.size() ) ;
            std::copy( gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0] );
        }
    }

    setCachedGaussPoints( new GaussPointArray( gp ) ) ;
    return *getCachedGaussPoints();
}

void DelaunayTree3D::setElementOrder( Order elemOrder, double dt )
{
    if(allElementsCacheID != -1)
    {
        caches[allElementsCacheID].clear() ;
        coefs[allElementsCacheID].clear() ;
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



void Amie::DelaunayTreeItem3D::print() const
{
}


std::pair<std::vector<DelaunayTetrahedron *>, std::vector<Point *> > Amie::quad( const DelaunayTetrahedron *t )
{
    std::vector<DelaunayTetrahedron * > tris ;
    std::vector<Point *> points ;

    points.push_back( new Point( *t->first + ( *t->second - *t->first )*.5 ) ) ;
    points.push_back( new Point( *t->first + ( *t->third - *t->first )*.5 ) ) ;
    points.push_back( new Point( *t->first + ( *t->fourth - *t->first )*.5 ) ) ;
    points.push_back( new Point( *t->second + ( *t->third - *t->second )*.5 ) ) ;
    points.push_back( new Point( *t->second + ( *t->fourth - *t->second )*.5 ) ) ;
    points.push_back( new Point( *t->third + ( *t->fourth - *t->third )*.5 ) ) ;

    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[0], points[1], points[2], t->first, nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[0], points[3], points[4], t->second, nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[1], points[3], points[5], t->third, nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[2], points[4], points[5], t->fourth, nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[0], points[1], points[2], points[3], nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[0], points[3], points[4], points[1], nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[1], points[3], points[5], points[0], nullptr ) ) ;
    tris.push_back( new DelaunayTetrahedron( nullptr, nullptr, points[2], points[4], points[5], points[1], nullptr ) ) ;
    return std::make_pair( tris, points ) ;
}
