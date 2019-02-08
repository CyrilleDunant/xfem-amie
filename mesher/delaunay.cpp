
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


#include "delaunay.h"
#include <limits>
#include "../features/crack.h"
#include "../utilities/samplingcriterion.h"
#include "../physics/dual_behaviour.h"
#include "../features/inclusion.h"
#include "../features/boundarycondition.h"

#define DEBUG
#undef DEBUG

namespace Amie{

DelaunayTreeItem::DelaunayTreeItem( Mesh<DelaunayTriangle, DelaunayTreeItem> * t, DelaunayTreeItem * father,  const Point * c) : dead(false), m_c(c), index(0),tree(t), father(father), stepfather(nullptr),first(nullptr) ,
    second(nullptr) ,
    third(nullptr),
    isPlane( false ),
    isTriangle( false ),
    isDeadTriangle( false),
    erased(false),
    stepson(0), neighbour(0),  son(0) 
{

    if(t)
    {
        index = t->addToTree(this) ;
    }
}

bool DelaunayTriangle::isConflicting(const Geometry * g) const
{
// 	Circle c(getRadius(), getCircumCenter()) ;
    return g->in(*first)
           || g->in(*second)
           || g->in(*third)
           || in(g->getCenter())
           || g->intersects(getPrimitive()) ;

}

bool DelaunayDeadTriangle::isConflicting(const Geometry * g) const
{
    Triangle t(*first, *second, *third) ;
    return g->in(*first)
           || g->in(*second)
           || g->in(*third)
           || inCircumCircle(g->getCenter())
           || g->intersects(&t) ;

}

bool DelaunayDemiPlane::isConflicting(const Geometry * g) const
{
    if(g->in(*first))
        return true ;
    if(g->in(*second))
        return true ;
    if(inCircumCircle(g->getCenter()))
        return true ;

    return false ;
}

bool DelaunayTreeItem::isConflicting(const Geometry * g) const
{
    return inCircumCircle(g->getCenter()) || g->in(*first) || g->in(*second) || (isTriangle && g->in(*third)) ;
}



DelaunayTreeItem::~DelaunayTreeItem()
{
}

const Point * DelaunayTreeItem::killer() const
{
    return m_k ;
}

const Point * DelaunayTreeItem::creator() const
{
    return m_c ;
}

void  DelaunayTreeItem::setCreator(const Point * p)
{
    m_c = p ;
}

size_t DelaunayTreeItem::numberOfCommonVertices(const DelaunayTreeItem * s) const
{
    if(this->isTriangle && s->isTriangle)
    {
        size_t ret = 0 ;

        if(s->first == this->first || s->first == this->second || s->first == this->third)
            ret++ ;
        if(s->second == this->first || s->second == this->second || s->second == this->third)
            ret++ ;
        if(s->third == this->first || s->third == this->second || s->third == this->third)
            ret++ ;

        assert(ret < 4) ;

        return ret ;
    }
    if(this->isTriangle && s->isPlane)
    {
        size_t ret = 0 ;

        if(isAligned(s->first, s->second, this->first))
            ret++ ;
        if(isAligned(s->first, s->second, this->second))
            ret++ ;
        if(isAligned(s->first, s->second, this->third))
            ret++ ;

        assert(ret < 4) ;

        return ret ;
    }
    if(this->isPlane && s->isPlane)
    {
        size_t ret = 0 ;

        if(isAligned(s->first, s->second, this->first))
            ret++ ;
        if(isAligned(s->first, s->second, this->second))
            ret++ ;

        assert(ret < 4) ;

        return ret ;
    }

    if(this->isPlane && s->isTriangle)
    {
        size_t ret = 0 ;

        if(isAligned(this->first, this->second, s->first))
            ret++ ;
        if(isAligned(this->first, this->second, s->second))
            ret++ ;
        if(isAligned(this->first, this->second, s->third))
            ret++ ;

        assert(ret < 4) ;

        return ret ;
    }

    return 0 ;
}

void Star::updateNeighbourhood()
{
    int count = 0 ;
    int soncount = 0 ;
    for( auto i = treeitem.begin() ; i != treeitem.end() ; ++i)
    {
        count += (*i)->son.size() + (*i)->stepson.size() + (*i)->neighbour.size() ;
        soncount += (*i)->son.size() ;
    }
    std::valarray<DelaunayTreeItem *> items(count) ;
    count = 0 ;
    for(auto i = treeitem.begin() ; i != treeitem.end() ; ++i)
    {
        for(size_t j = 0 ; j < (*i)->son.size() ; j++)
        {
            items[count++] = (*i)->getSon(j) ;
        }

    }
    for( auto i = treeitem.begin() ; i != treeitem.end() ; ++i)
    {
        for(size_t j = 0 ; j < (*i)->stepson.size() ; j++)
        {
            if((*i)->getStepson(j)->isAlive() && ((*i)->getStepson(j)->onCircumCircle(*creator) || !(*i)->getStepson(j)->inCircumCircle(*creator)))
                items[count++] = (*i)->getStepson(j) ;
        }
        for(size_t j = 0 ; j < (*i)->neighbour.size() ; j++)
        {
            if((*i)->getNeighbour(j)->isAlive() && ((*i)->getNeighbour(j)->onCircumCircle(*creator) || !(*i)->getNeighbour(j)->inCircumCircle(*creator)))
                items[count++] = (*i)->getNeighbour(j) ;
        }

    }

    std::sort(&items[soncount], &items[count]) ;
    auto e = std::unique(&items[0], &items[count]) ;

    if(!items.size())
        return ;

    for(DelaunayTreeItem ** i = &items[0] ; i != e+1 && i != &items[count] ; ++i)
    {
// 		if(!(*i)->isPlane)
// 		{
        DelaunayTreeItem * ii = *i ;
        size_t ins = ii->neighbour.size() ;
// 			if(ins != 3 || (*i)->isPlane)
// 			{
        for(DelaunayTreeItem ** j = i+1 ; j != e+1 && j != &items[count] ; ++j)
        {

            DelaunayTreeItem * jj = *j ;

            size_t jns = jj->neighbour.size() ;
            if(ii->numberOfCommonVertices(jj) == 2 /*&& (jns != 3 || (*i)->isPlane)*/)
            {

                if(std::find(&ii->neighbour[0], &ii->neighbour[ins],jj->index) ==  &ii->neighbour[ins])
                {
                    std::valarray<unsigned int>  newneighbours(ins+1) ;
                    std::copy(&ii->neighbour[0], &ii->neighbour[ins], &newneighbours[0]) ;
                    newneighbours[ins] = jj->index ;
                    ii->neighbour.resize(ins+1) ;
                    ii->neighbour = newneighbours ;
                }

                if(std::find(&jj->neighbour[0], &jj->neighbour[jns],ii->index) ==  &jj->neighbour[jns])
                {
                    std::valarray<unsigned int>  newneighbours(jns+1) ;
                    std::copy(&jj->neighbour[0], &jj->neighbour[jns], &newneighbours[0]) ;
                    newneighbours[jns] = ii->index ;
                    jj->neighbour.resize(jns+1) ;
                    jj->neighbour = newneighbours ;
                }
                ins = ii->neighbour.size() ;
// 						if(ins == 3 && !(*i)->isPlane)
// 							break ;
            }
        }
// 			}
// 		}
    }

    return ;
}

bool DelaunayTreeItem::isDuplicate(const DelaunayTreeItem * t) const
{
    return t!=this && this->isTriangle && t->isTriangle && this->numberOfCommonVertices(t) == 3 ;
}

void DelaunayTree::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep, const TriElement * father)
{
    std::vector<DelaunayTriangle *> tri = getTriangles() ;
    if(visited.size() != size())
      visited.resize(size(), false);
    else
      visited = false ;

    double timeSlice = timestep ;

    for(auto & i : tri)
    {
        delete i->cachedGps ;
        i->cachedGps = nullptr ;
        i->behaviourUpdated = true ;
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
                        for(size_t j = 0 ; j < i->neighbourhood.size() ; j++)
                        {
                            if(visited[i->neighbourhood[j]])
                            {
                                DelaunayTriangle * n = i->getNeighbourhood(j) ;
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
// 	i->genGaussPoints() ;
    }  

}

void DelaunayTree::addSharedNodes(DelaunayTree * dt)
{
    std::vector<DelaunayTriangle *> tri = getTriangles() ;
    std::vector<DelaunayTriangle *> tris = dt->getTriangles() ;

    for(size_t i = 0 ; i <  tris.size() ; i++)
    {
        std::valarray<Point *> newPoints = tris[i]->getBoundingPoints() ;
        tri[i]->setBoundingPoints(newPoints) ;
    }
}

void DelaunayTree::extrude(double dt)
{
    std::map<Point *, Point *> points ;

    std::vector<DelaunayTriangle *> tri = getTriangles() ;
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
            Point * next = new Point(tri[i]->getBoundingPoint(j).getX(), tri[i]->getBoundingPoint(j).getY()) ;
            additionalPoints.push_back(next);
            next->getT() = tri[i]->getBoundingPoint(j).getT() ;
            next->getT() = end + dt * (next->getT() - beginning) / (end - beginning) ;
            bool increment = true ;
            if(std::abs(next->getT() - end) < POINT_TOLERANCE)
            {
                next = &tri[i]->getBoundingPoint(j+indexOfLastTimePlane) ;
                increment = false ;
// 				next->print() ;
            }
            if(increment && !points.find(&tri[i]->getBoundingPoint(j))->second)
            {
                next->setId(global_counter++) ;
            }
            points.insert(std::pair<Point *, Point *>(&tri[i]->getBoundingPoint(j), next)) ;
        }
    }

    std::map<Point *, Point *>::iterator finder ;
    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        Point * a = points.find(&tri[i]->getBoundingPoint(0))->second ;
        Point * b = points.find(&tri[i]->getBoundingPoint(pointsPerTimePlane/3))->second ;
        Point * c = points.find(&tri[i]->getBoundingPoint(2*pointsPerTimePlane/3))->second ;

        std::valarray<Point *> newPoints(tri[i]->getBoundingPoints().size()) ;
        for(size_t j = 0 ; j < newPoints.size() ; j++)
        {
            newPoints[j] = points.find(&tri[i]->getBoundingPoint(j))->second ;
//			newPoints[j]->print() ;
        }

        DelaunayTriangle * toInsert = new DelaunayTriangle(tri[i]->tree, nullptr, a,b,c, a) ;
        toInsert->setOrder(tri[i]->getOrder()) ;
        toInsert->setBoundingPoints(newPoints) ;
        toInsert->setBehaviour(tri[i]->getState().getMesh2D(), tri[i]->getBehaviour()->getCopy()) ;
    }
}


void DelaunayTree::extrude(const Vector & dt)
{
    std::map<Point *, std::vector<Point *> > points ;
    std::map<Point *, Point *> pointsInTriangle ;
    std::map<Point *, std::vector<Point *> >::iterator finder ;

    std::vector<DelaunayTriangle *> tri = getTriangles() ;
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
                        additionalPoints.push_back(next);
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

            DelaunayTriangle * toInsert = new DelaunayTriangle(tri[i]->tree, nullptr, newPoints[0],newPoints[pointsPerTimePlane/3],newPoints[pointsPerTimePlane*2/3],newPoints[0]) ;
            toInsert->setOrder(tri[i]->getOrder()) ;
            toInsert->setBoundingPoints(newPoints) ;
            toInsert->setBehaviour(tri[i]->getState().getMesh2D(), tri[i]->getBehaviour()->getCopy()) ;
        }
    }

    std::cout << getTriangles().size() << "\t" << tri.size() << std::endl ;

}

void DelaunayTree::setElementOrder(Order elemOrder, double dt)
{
    if(allElementsCacheID != -1)
    {
        caches[allElementsCacheID].clear() ;
        coefs[allElementsCacheID].clear() ;
        allElementsCacheID = -1 ;
    }
    switch(elemOrder)
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
        addSharedNodes(1,1,0) ;
        break ;
    }
    case CUBIC:
    {
        addSharedNodes(2,1,0) ;
        break ;
    }
    case QUADRIC:
    {
        addSharedNodes(3,1,0) ;
        break ;
    }
    case QUINTIC:
    {
        addSharedNodes(3,1,0) ;
        break ;
    }
    case CONSTANT_TIME_LINEAR:
    {
        addSharedNodes(0,2,dt) ;
        break ;
    }
    case CONSTANT_TIME_QUADRATIC:
    {
        addSharedNodes(0,3,dt) ;
        break ;
    }
    case LINEAR_TIME_LINEAR:
    {
        addSharedNodes(0,2,dt) ;
        break ;
    }
    case LINEAR_TIME_QUADRATIC:
    {
        addSharedNodes(0,3,dt) ;
        break ;
    }
    case QUADRATIC_TIME_LINEAR:
    {
        addSharedNodes(1,2,dt) ;
        break ;
    }
    case QUADRATIC_TIME_QUADRATIC:
    {
        addSharedNodes(1,3,dt) ;
        break ;
    }
    case CUBIC_TIME_LINEAR:
    {
        addSharedNodes(2,2,dt) ;
        break ;
    }
    case CUBIC_TIME_QUADRATIC:
    {
        addSharedNodes(2,3,dt) ;
        break ;
    }
    case QUADRIC_TIME_LINEAR:
    {
        addSharedNodes(3,2,dt) ;
        break ;
    }
    case QUADRIC_TIME_QUADRATIC:
    {
        addSharedNodes(3,3,dt) ;
        break ;
    }
    case QUINTIC_TIME_LINEAR:
    {
        addSharedNodes(3,2,dt) ;
        break ;
    }
    case QUINTIC_TIME_QUADRATIC:
    {
        addSharedNodes(3,3,dt) ;
        break ;
    }
    default:
        break ;

    }
}

void DelaunayTree::refresh(TriElement *father)
{

    for(auto i = tree.begin() ; i != tree.end() ; ++i)
    {
        if((*i)->isAlive() && (*i)->isTriangle)
        {
            (static_cast<DelaunayTriangle *>(*i))->refresh(father) ;
        }
    }
}

size_t DelaunayTree::numPoints() const
{
    return this->global_counter ;
}

void DelaunayTreeItem::flatConflicts(std::valarray<bool> & visited,std::vector<DelaunayTreeItem *> &toTest, std::vector< DelaunayTriangle* > & ret, const Amie::Geometry* g)
{
    if(visited[index])
    {
        return ;
    }
    visited[index] = true ;

    if(!isConflicting(g) )
    {
        return ;
    }

    for (size_t i  = 0 ;  i < stepson.size() ; i++)
    {
        bool limit = false ;

        if( (!visited[stepson[i]] && getStepson(i)->isConflicting(g)) || limit)
        {
            toTest.push_back(getStepson(i)) ;
        }
    }

    for (size_t i  = 0 ;  i < son.size() ; i++)
    {
        bool limit = false ;

        if( (!visited[son[i]] && getSon(i)->isConflicting(g)) || limit)
        {
            toTest.push_back(getSon(i)) ;
        }
    }

    if(isAlive() && isTriangle && !isDeadTriangle && isConflicting(g))
    {
        ret.push_back(static_cast<DelaunayTriangle *>(this)) ;
    }

    for (size_t i  = 0 ;  i < neighbour.size() ; i++)
    {
        bool limit = false ;

        if( (!visited[neighbour[i]] && getNeighbour(i)->isConflicting(g)) || limit)
        {
            toTest.push_back(getNeighbour(i)) ;
        }
    }
}

// void DelaunayTreeItem::conflicts(std::valarray<bool> & visited, std::vector<DelaunayTriangle *> & ret, const Geometry *g)
// {
//     if(visited[index])
//     {
//         return ;
//     }
//     visited[index] = true ;
// 
//     if(!isConflicting(g) )
//     {
//         return ;
//     }
// 
//     for (size_t i  = 0 ;  i < stepson.size() ; i++)
//     {
//         if( (!visited[stepson[i]] && getStepson(i)->isConflicting(g)) )
//         {
//             getStepson(i)->conflicts(visited, ret, g) ;
//         }
//     }
// 
//     for (size_t i  = 0 ;  i < son.size() ; i++)
//     {
//         if( (!visited[son[i]] && getSon(i)->isConflicting(g)))
//         {
//             getSon(i)->conflicts(visited, ret,g) ;
//         }
//     }
// 
//     if(isAlive() && isTriangle && !isDeadTriangle && isConflicting(g))
//     {
//         ret.push_back(static_cast<DelaunayTriangle *>(this)) ;
//     }
// 
//     for (size_t i  = 0 ;  i < neighbour.size() ; i++)
//     {
// 
//         if( (!visited[neighbour[i]] && getNeighbour(i)->isConflicting(g)) )
//         {
//             getNeighbour(i)->conflicts(visited,ret, g) ;
//         }
//     }
// }
// 
// void DelaunayTreeItem::conflicts(std::valarray<bool> & visited,std::vector<DelaunayTreeItem *>& ret,const Point *p)
// {
// // 	print() ;
//     if(visited[index])
//         return  ;
//     visited[index] = true ;
// 
// 
//     for (size_t i  = 0 ;  i < stepson.size() ; i++)
//     {
//         if( !visited[stepson[i]] && getStepson(i)->inCircumCircle(*p) )
//         {
//             getStepson(i)->conflicts(visited,ret,p) ;
//         }
//     }
// 
// 
//     for (size_t i  = 0 ;  i < son.size() ; i++)
//     {
//         if( !visited[son[i]] && getSon(i)->inCircumCircle(*p) )
//         {
//             getSon(i)->conflicts(visited,ret,p) ;
//         }
// 
//     }
// 
//     if(!inCircumCircle(*p))
//         return  ;
// 
//     if(isAlive())
//     {
//         ret.push_back(this) ;
//     }
// 
//     for (size_t i  = 0 ;  i < neighbour.size() ; i++)
//     {
// 
//         if( !visited[neighbour[i]] && getNeighbour(i)->inCircumCircle(*p))
//         {
//             getNeighbour(i)->conflicts(visited,ret,p) ;
//         }
//     }
// }

void DelaunayTreeItem::flatConflicts(std::valarray<bool> & visited,std::vector<DelaunayTreeItem *> & toTest,  std::vector<DelaunayTreeItem *> & ret,const Point *p)
{
    if(visited[index])
        return  ;
    visited[index] = true ;

    if(!inCircumCircle(*p))
        return  ;

    for (size_t i  = 0 ;  i < stepson.size() ; i++)
    {
        bool limit = false ;
        if(!visited[stepson[i]])
        {
            if(getStepson(i)->isTriangle && !getStepson(i)->isDeadTriangle)
            {
                DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getStepson(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }
            if(getStepson(i)->isDeadTriangle)
            {
                DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getStepson(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }

            if( (getStepson(i)->inCircumCircle(*p)) || limit)
            {
                toTest.push_back(getStepson(i)) ;
            }
        }
    }


    for (size_t i  = 0 ;  i < son.size() ; i++)
    {
        bool limit = false ;
        if(!visited[son[i]])
        {
            if(getSon(i)->isTriangle && !getSon(i)->isDeadTriangle)
            {
                DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getSon(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }
            if(getSon(i)->isDeadTriangle)
            {
                DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getSon(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }

            if( (getSon(i)->inCircumCircle(*p)) || limit)
            {
                toTest.push_back(getSon(i)) ;
            }
        }
    }

    if(isAlive())
    {
        ret.push_back(this) ;
    }

    for (size_t i  = 0 ;  i < neighbour.size() ; i++)
    {

        bool limit = false ;
        if(!visited[neighbour[i]])
        {
            if(getNeighbour(i)->isTriangle && !getNeighbour(i)->isDeadTriangle)
            {
                DelaunayTriangle * t = static_cast<DelaunayTriangle *>(getNeighbour(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),*p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }
            if(getNeighbour(i)->isDeadTriangle)
            {
                DelaunayDeadTriangle * t = static_cast<DelaunayDeadTriangle *>(getNeighbour(i)) ;
                limit = std::abs(squareDist2D(t->getCircumCenter(),p)-t->getRadius()*t->getRadius())
                        < 20000.*POINT_TOLERANCE*POINT_TOLERANCE ;
            }
            // 		limit = true ;
            if( (getNeighbour(i)->inCircumCircle(*p)) || limit)
            {
                toTest.push_back(getNeighbour(i)) ;
            }
        }
    }
}

void DelaunayTreeItem::removeNeighbour(DelaunayTreeItem * t)
{
    unsigned int* e = std::find(&neighbour[0], &neighbour[neighbour.size()], t->index) ;

    if(e !=  &neighbour[neighbour.size()])
    {
        std::valarray<unsigned int>  newneighbours(neighbour.size()-1) ;
        int iterator = 0 ;
        for(size_t i = 0 ; i < neighbour.size() ; i++)
        {
            if(neighbour[i] == t->index)
                continue ;

            newneighbours[iterator] = neighbour[i] ;
            iterator++ ;
        }
        /*

        		std::copy(&neighbour[0], e, &newneighbours[0]) ;
        		if(e+1 != &neighbour[neighbour.size()])
        			std::copy(e+1, &neighbour[neighbour.size()], &newneighbours[e-&neighbour[0]]) ;*/
        neighbour.resize(neighbour.size()-1) ;
        neighbour = newneighbours ;
    }
}

void DelaunayTriangle::removeNeighbourhood(DelaunayTriangle * t)
{
    unsigned int* e = std::find(&neighbourhood[0], &neighbourhood[neighbourhood.size()], t->index) ;

    if(e !=  &neighbourhood[neighbourhood.size()])
    {
        std::valarray<unsigned int>  newneighbours(neighbourhood.size()-1) ;
        std::copy(&neighbourhood[0], e, &newneighbours[0]) ;
        std::copy(e+1, &neighbourhood[neighbourhood.size()], &newneighbours[e-&neighbourhood[0]]) ;
        neighbourhood.resize(neighbourhood.size()-1) ;
        neighbourhood = newneighbours ;
    }
}

void DelaunayTreeItem::addNeighbour(DelaunayTreeItem * t)
{
    assert(t != this) ;

    if(std::find(&neighbour[0], &neighbour[neighbour.size()],t->index) !=  &neighbour[neighbour.size()])
    {
        return ;
    }

    if(t->isAlive())
    {
        std::valarray<unsigned int>  newneighbours(neighbour.size()+1) ;
        std::copy(&neighbour[0], &neighbour[neighbour.size()], &newneighbours[0]) ;
        newneighbours[neighbour.size()] = t->index ;
        neighbour.resize(neighbour.size()+1) ;
        neighbour = newneighbours ;
    }
}

void DelaunayTriangle::addNeighbourhood(DelaunayTriangle * t)
{
    assert(t != this) ;

    if(std::find(&neighbourhood[0], &neighbourhood[neighbourhood.size()],t->index) !=  &neighbourhood[neighbourhood.size()])
    {
        return ;
    }

    if(t->isAlive())
    {
        std::valarray<unsigned int>  newneighbours(neighbourhood.size()+1) ;
        std::copy(&neighbourhood[0], &neighbourhood[neighbourhood.size()], &newneighbours[0]) ;
        newneighbours[neighbourhood.size()] = t->index ;
        neighbourhood.resize(neighbourhood.size()+1) ;
        neighbourhood = newneighbours ;
    }
}

DelaunayTriangle * DelaunayTriangle::getNeighbourhood(size_t i) const
{
    return static_cast<DelaunayTriangle *>(tree->getInTree(neighbourhood[i])) ;
}

DelaunayTreeItem * DelaunayTreeItem::getNeighbour(size_t i) const
{
    return tree->getInTree(neighbour[i]) ;
}

DelaunayTreeItem * DelaunayTreeItem::getSon(size_t i) const
{
    return tree->getInTree(son[i]) ;
}

DelaunayTreeItem * DelaunayTreeItem::getStepson(size_t i) const
{
    if(i < stepson.size())
        return tree->getInTree(stepson[i]) ;
    return nullptr ;
}

DelaunayTreeItem * DelaunayTreeItem::getFather() const
{
    return father ;
}
DelaunayTreeItem * DelaunayTreeItem::getStepfather() const
{
    return stepfather ;
}

void DelaunayTreeItem::addStepson(DelaunayTreeItem * s)
{
    if(s == this)
        return ;

    std::valarray<unsigned int>  newstepson(stepson.size()+1) ;
    std::copy(&stepson[0], &stepson[stepson.size()], &newstepson[0]) ;
    newstepson[stepson.size()] = s->index ;
    stepson.resize(stepson.size()+1) ;
    stepson = newstepson ;
    s->setStepfather(this) ;
    addNeighbour(s) ;
}

void DelaunayTreeItem::addSon(DelaunayTreeItem * s)
{
    std::valarray<unsigned int>  newson(son.size()+1) ;
    std::copy(&son[0], &son[son.size()], &newson[0]) ;
    newson[son.size()] = s->index ;
    son.resize(son.size()+1) ;
    son = newson ;
}

void DelaunayTreeItem::removeSon(DelaunayTreeItem * t)
{
    if(!son.size())
        return ;
    unsigned int* e = std::find(&son[0], &son[son.size()], t->index) ;

    if(e !=  &son[son.size()])
    {
        std::valarray<unsigned int>  newson(son.size()-1) ;
        std::copy(&son[0], e, &newson[0]) ;
        std::copy(e+1, &son[son.size()], &newson[e-&son[0]]) ;
        son.resize(son.size()-1) ;
        son = newson ;
    }
}

void DelaunayTreeItem::removeStepson(DelaunayTreeItem * t)
{
    if(stepson.size() == 0)
        return ;
    unsigned int* e = std::find(&stepson[0], &stepson[stepson.size()], t->index) ;
    if(e !=  &stepson[stepson.size()])
    {
        std::valarray<unsigned int>  newstepson(stepson.size()-1) ;
        std::copy(&stepson[0], e, &newstepson[0]) ;
        std::copy(e+1, &stepson[stepson.size()], &newstepson[e-&stepson[0]]) ;
        stepson.resize(stepson.size()-1) ;
        stepson = newstepson ;
        if(newstepson.size() == 0)
            t->stepfather = nullptr ;
    }
}


void DelaunayTreeItem::erase(const Point * p)
{
    dead = true ;
    m_k  = p ;
}

void DelaunayTreeItem::kill(const Point * p)
{
    dead = true ;
    m_k  = p ;

    for(size_t i = 0 ; i < this->neighbour.size() ; i++)
    {
        this->getNeighbour(i)->removeNeighbour(this) ;
    }
}

bool DelaunayTreeItem::isAlive() const
{
    return !this->dead ;
}



void DelaunayTreeItem::setStepfather(DelaunayTreeItem * s)
{
    stepfather = s ;
    addNeighbour(s) ;
}


DelaunayTriangle::DelaunayTriangle(Mesh<DelaunayTriangle, DelaunayTreeItem> *t, DelaunayTreeItem * father,  Point *p0,  Point *p1, Point *p2,  Point * c) : TriElement(p0, p1, p2),DelaunayTreeItem(t, father, c),  neighbourhood(0)
{
    first = &getBoundingPoint(0) ;
    second = &getBoundingPoint(1) ;
    third = &getBoundingPoint(2) ;

    assert(in(this->getCenter())) ;

    isPlane = false ;
    isTriangle = true ;
    assert(first->getId() > -1) ;
    assert(second->getId() > -1) ;
    assert(third->getId() > -1) ;
}

DelaunayTriangle::DelaunayTriangle() :  TriElement(),DelaunayTreeItem(nullptr, nullptr, nullptr), neighbourhood(0)
{
    first = &getBoundingPoint(0) ;
    second = &getBoundingPoint(1) ;
    third = &getBoundingPoint(2) ;

    isPlane = false ;
    isTriangle = true ;
}

DelaunayDeadTriangle::DelaunayDeadTriangle( DelaunayTriangle * parent) : DelaunayTreeItem(*parent),
    center(parent->getCircumCenter()),
    radius(parent->getRadius()), sqradius(radius*radius)
{
    index = parent->index ;
    tree = parent->tree ;

    stepson.resize(parent->stepson.size()) ;
    stepson = parent->stepson ;

    neighbour.resize(parent->neighbour.size()) ;
    neighbour = parent->neighbour ;

    son.resize(parent->son.size()) ;
    son = parent->son ;

    father = parent->father ;
    stepfather = parent->stepfather ;

    isTriangle = true ;
    isPlane = false ;
    isDeadTriangle = true ;

    first = parent->first ;
    second = parent->second ;
    third = parent->third ;
    dead = true ;
}


double DelaunayDeadTriangle::getRadius() const
{
    return radius ;
}

const Point * DelaunayDeadTriangle::getCircumCenter() const
{
    return &center ;
}

bool DelaunayDeadTriangle::onCircumCircle(const Point & p) const
{
    if(p.getX() > center.getX()+1.01*radius)
        return false ;
    if(p.getX() < center.getX()-1.01*radius)
        return false ;
    if(p.getY() > center.getY()+1.01*radius)
        return false ;
    if(p.getY() < center.getY()-1.01*radius)
        return false ;

    double d = squareDist2D(center, p) ;
    return  std::abs(d/(radius*radius)-1) < POINT_TOLERANCE/radius ;
}

bool DelaunayDeadTriangle::inCircumCircle(const Point & p) const
{


    if(p.getX() > center.getX()+1.01*radius)
        return false ;
    if(p.getX() < center.getX()-1.01*radius)
        return false ;
    if(p.getY() > center.getY()+1.01*radius)
        return false ;
    if(p.getY() < center.getY()-1.01*radius)
        return false ;


    double d = squareDist2D(center, p) ;
    return  d/(radius*radius)-1 < POINT_TOLERANCE/radius ;

}

bool DelaunayDeadTriangle::isNeighbour( const DelaunayTreeItem * t) const
{
    size_t cv = this->numberOfCommonVertices(t) ;
    return (cv == 2);
}

bool DelaunayDeadTriangle::isVertex(const Point *p) const
{
    return (*p == *first) || (*p == *second) || (*p == *third)  ;
}

double DelaunayDeadTriangle::distanceToVertex(const Point *p) const
{
    if(this->isVertex(p))
        return 0 ;
    return sqrt(std::min( squareDist2D(*p, *first), std::min( squareDist2D(*p, *second),  squareDist2D(*p, *third) ) ) ) ;
}

bool DelaunayDeadTriangle::isVertexByID(const Point *p) const
{
    return p == first || p == second || p == third  ;
}

bool DelaunayDeadTriangle::in( const Point & p) const
{
    return Triangle(third, first, second).in(p) ;
}

void DelaunayDeadTriangle::print() const
{

    std::cout << "(" << first->getX() << ", " << first->getY() <<  ") " ;
    std::cout << "(" << second->getX() << ", " << second->getY() <<  ") " ;
    std::cout << "(" << third->getX() << ", " << third->getY() <<  ") " ;
    std::cout <<  ":: "<< isAlive() << std::endl ;
}

inline std::pair<Point*, Point*> DelaunayDeadTriangle::commonEdge(const DelaunayTreeItem * t) const
{
    if(t == this)
    {
        return std::pair<Point*, Point*>(nullptr, nullptr) ;
    }

    if(t->isTriangle)
    {
        if(this->isVertex(t->first))
        {
            if(this->isVertex(t->second))
                return std::pair< Point*,  Point*>(t->first , t->second) ;
            if(this->isVertex(t->third))
                return std::pair< Point*,  Point*>(t->first , t->third) ;

            return std::pair<Point*, Point*>(nullptr, nullptr) ;
        }

        if(this->isVertex(t->second))
        {
            if(this->isVertex(t->third))
                return std::pair< Point*,  Point*>(t->third , t->second) ;

            return std::pair<Point*, Point*>(nullptr, nullptr) ;
        }

    }
    else
    {
        if(isAligned(t->first, first, second) && isAligned(t->second, first, second))
            return std::pair< Point*,  Point*>(first, second ) ;
        if(isAligned(t->first, third, second) && isAligned(t->second, third, second))
            return std::pair< Point*,  Point*>(third, second ) ;
        if(isAligned(t->first, third, first) && isAligned(t->second, third, first))
            return std::pair< Point*,  Point*>(first, third ) ;
    }

    return std::pair< Point*,  Point*>(nullptr, nullptr) ;
}


DelaunayTriangle::~DelaunayTriangle()
{

}

DelaunayDemiPlane::~DelaunayDemiPlane()
{
// 	first = nullptr ;
// 	second = nullptr ;
// 	third = nullptr ;
}

inline bool DelaunayTriangle::isVertex(const Point * p) const
{
    return (*p == *first || *p == *second || *p == *third) ;
}

double DelaunayTriangle::distanceToVertex(const Point *p) const
{
    if(this->isVertex(p))
        return 0 ;
    return sqrt(std::min( squareDist2D(*p, *first), std::min( squareDist2D(*p, *second),  squareDist2D(*p, *third) ) ) ) ;
}

bool DelaunayTriangle::isVertexByID(const Point * p) const
{
    return p == first || p == second || p == third ;
}


inline bool DelaunayTriangle::hasVertexByID(const std::valarray<Point *> * p) const
{
    for(size_t i = 0 ; i < p->size() ; i++)
    {
        if((*p)[i]->getId() == first->getId() || (*p)[i]->getId() == second->getId() || (*p)[i]->getId() == third->getId())
            return true ;
    }
    return false ;
}

std::pair<Point*, Point*> DelaunayTriangle::commonEdge(const DelaunayTreeItem * t) const
{
    if(t == this)
    {
        return std::pair<Point*, Point*>(nullptr, nullptr) ;
    }

    if(t->isTriangle || t->isDeadTriangle)
    {
        if(this->isVertex(t->first))
        {
            if(this->isVertex(t->second))
                return std::pair< Point*,  Point*>(t->first , t->second) ;
            if(this->isVertex(t->third))
                return std::pair< Point*,  Point*>(t->first , t->third) ;

            return std::pair<Point*, Point*>(nullptr, nullptr) ;
        }

        if(this->isVertex(t->second))
        {
            if(this->isVertex(t->third))
                return std::pair< Point*,  Point*>(t->third , t->second) ;

            return std::pair<Point*, Point*>(nullptr, nullptr) ;
        }

    }
    else if (numberOfCommonVertices(t) == 2)
    {
        if(!isAligned(first, t->first, t->second))
            return std::pair< Point*,  Point*>(third, second ) ;
        if(!isAligned(second, t->first, t->second))
            return std::pair< Point*,  Point*>(third, first ) ;
        if(!isAligned(t->third, t->first, t->second))
            return std::pair< Point*,  Point*>(first, second ) ;

        //now, this means they are all "aligned"
        Segment seg(*t->first,*t->second) ;
        double d0 = squareDist2D(*first, seg.project(*first)) ;
        double d1 = squareDist2D(*second, seg.project(*second)) ;
        double d2 = squareDist2D(*third, seg.project(*third)) ;
        if(d0 >= d1 && d0 >=d2)
            return std::pair< Point*,  Point*>(third, second ) ;
        else if(d1 >= d0 && d1 >=d2)
            return std::pair< Point*,  Point*>(third, first ) ;

        return std::pair< Point*,  Point*>(first, second ) ;
    }

    return std::pair< Point*,  Point*>(nullptr, nullptr) ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::commonEdge(const DelaunayTreeItem * t) const
{
    if(t->isTriangle)
    {
        Point test = (*first)*0.5 + (*second)*0.5 ;

        if(isAligned(test, (*t->first), (*t->second)))
        {
            return std::pair< Point*,  Point*>(t->first, t->second ) ;
        }
        if(isAligned(test, (*t->third), (*t->second)))
        {
            return std::pair< Point*,  Point*>(t->third, t->second ) ;
        }
        if(isAligned(test, (*t->first), (*t->third)))
        {
            return std::pair< Point*,  Point*>( t->first, t->third ) ;
        }
    }
    else
    {
        if((t->first == first && t->second == second) ||
                (t->first == second &&  t->second == first))
            return std::pair< Point*,  Point*>(t->first , t->second) ;
    }

    return std::pair< Point*,  Point*>(nullptr, nullptr) ;
}


void DelaunayDemiPlane::merge(DelaunayDemiPlane * p)
{
    if(isAlive() && p != this && p->isAlive())
    {

        if( (first == p->first && second == p->second ) || (first == p->second && second == p->first) )
        {
            for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
            {
                for(size_t j = 0 ; j <  neighbour.size() ; j++)
                {
                    if(p->getNeighbour(i)->isNeighbour(getNeighbour(j)))
                        makeNeighbours(getNeighbour(j), p->getNeighbour(i)) ;
                }
            }
            p->kill(first) ;
            kill(first) ;
            return ;
        }

        if(isAligned(first, second, p->first) && isAligned(first, second, p->second))
        {

            for(size_t i = 0 ; i <  p->neighbour.size() ; i++)
            {
                makeNeighbours(this, p->getNeighbour(i)) ;
            }
            for(size_t i = 0 ; i <  p->son.size() ; i++)
            {
                p->getSon(i)->father = this ;
                addSon(p->getSon(i)) ;
            }
            for(size_t i = 0 ; i <  p->stepson.size() ; i++)
            {
                p->getStepson(i)->stepfather = this ;
                addStepson(p->getStepson(i)) ;
                makeNeighbours(this, p->getStepson(i)) ;
            }
            p->father->addSon(this) ;
            p->addSon(this) ;
            p->kill(first) ;
        }
    }
}


std::pair< Point*,  Point*> DelaunayDeadTriangle::nearestEdge(const Point & p) const
{
    std::map<double, Point> cen ;
    Point c0(((*first) + (*second))/2.) ;
    cen[squareDist2D(c0, p)] = c0 ;
    Point c1(((*third) + (*second))/2.) ;
    cen[squareDist2D(c1, p)] = c1 ;
    Point c2(((*third) + (*first))/2.) ;
    cen[squareDist2D(c2, p)] = c2 ;

    if(cen.begin()->second == c0)
        return std::pair< Point*,  Point*>(first, second) ;
    if(cen.begin()->second == c1)
        return std::pair< Point*,  Point*>(third, second) ;
    if(cen.begin()->second == c2)
        return std::pair< Point*,  Point*>(third, first) ;

    return std::pair< Point*,  Point*>(first, second) ;
}

std::pair< Point*,  Point*> DelaunayTriangle::nearestEdge(const Point & p) const
{
    std::map<double, Point> cen ;
    Point c0(((*first) + (*second))/2.) ;
    cen[squareDist2D(c0, p)] = c0 ;
    Point c1(((*third) + (*second))/2.) ;
    cen[squareDist2D(c1, p)] = c1 ;
    Point c2(((*third) + (*first))/2.) ;
    cen[squareDist2D(c2, p)] = c2 ;

    if(cen.begin()->second == c0)
        return std::pair< Point*,  Point*>(first, second) ;
    if(cen.begin()->second == c1)
        return std::pair< Point*,  Point*>(third, second) ;
    if(cen.begin()->second == c2)
        return std::pair< Point*,  Point*>(third, first) ;

    return std::pair< Point*,  Point*>(first, second) ;
}

bool DelaunayTriangle::inCircumCircle(const Point &p) const
{
    if(p.getX() > circumCenter.getX()+1.0001*radius)
        return false ;
    if(p.getX() < circumCenter.getX()-1.0001*radius)
        return false ;
    if(p.getY() > circumCenter.getY()+1.0001*radius)
        return false ;
    if(p.getY() < circumCenter.getY()-1.0001*radius)
        return false ;

    double d = dist(circumCenter, p) ;
    return  d/radius-1 < POINT_TOLERANCE ;
}

bool DelaunayTriangle::onCircumCircle(const Point &p) const
{
    if(p.getX() > circumCenter.getX()+1.0001*radius)
        return false ;
    if(p.getX() < circumCenter.getX()-1.0001*radius)
        return false ;
    if(p.getY() > circumCenter.getY()+1.0001*radius)
        return false ;
    if(p.getY() < circumCenter.getY()-1.0001*radius)
        return false ;

    double d = dist(getCircumCenter(), p) ;
    return  std::abs(d/radius-1) < POINT_TOLERANCE ;
}


void DelaunayTriangle::insert(std::valarray<bool> & visited,std::vector<DelaunayTreeItem *> &ret, Point *p,  Star* s)
{
    if (visited[index])
        return ;

    visited[index] = true ;
    for (size_t i = 0 ; i < neighbour.size() ; i++)
    {
        if (!visited[neighbour[i]] && (!getNeighbour(i)->inCircumCircle(*p) || getNeighbour(i)->onCircumCircle(*p)))
        {
            std::pair< Point*,  Point*> pp = commonEdge(getNeighbour(i)) ;
            if(!isAligned(p,pp.first,pp.second))
            {
                DelaunayTriangle *ss = new DelaunayTriangle(tree, this, p, pp.first, pp.second, p) ;
                addSon(ss) ;

                getNeighbour(i)->addStepson(ss) ;
                ret.push_back(ss) ;
            }
        }
    }

}


bool DelaunayTriangle::isNeighbour(const DelaunayTreeItem * t)  const
{
// 	if(t->isPlane)
// 		return ((int)isAligned(t->first, t->second, first) + (int)isAligned(t->first, t->second, second) + (int)isAligned(t->first, t->second, third)) == 2;

    return numberOfCommonVertices(t) == 2 ;
}

bool DelaunayTriangle::isInNeighbourhood(const DelaunayTriangle * t) const
{
    return this->numberOfCommonVertices(t) > 0 ;
}

void DelaunayTriangle::print() const
{
    for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
    {
        std::cout << "(" << getBoundingPoint(i).getX() <<  ";" << getBoundingPoint(i).getY() <<") " ;
    }
    std::cout <<  ":: "<< isAlive() << std::endl ;
}

DelaunayDemiPlane::DelaunayDemiPlane(Mesh<DelaunayTriangle, DelaunayTreeItem> *t, DelaunayTreeItem * father,  Point  * _begin,  Point  * _end,  Point  * p,  Point * c) : DelaunayTreeItem(t, father, c)
{
    second  = _end ;
    first = _begin ;
    third = p ;
    dead = false ;
    vector =(*second)- (*first) ;
    Point pseudonormal = (*third) - (*first);
    direction = (vector.getX()*pseudonormal.getY() - vector.getY()*pseudonormal.getX()) ;
    isPlane =true ;
    isTriangle = false ;
}

std::pair< Point*,  Point*> DelaunayDemiPlane::nearestEdge(const Point & p) const
{
    return std::pair< Point*,  Point*>(first, second) ;
}

bool DelaunayDemiPlane::inCircumCircle(const Point &p) const
{
    return fma(vector.getX(),(p.getY() - first->getY()), - vector.getY()*(p.getX() - first->getX())) * direction < -POINT_TOLERANCE ;
}

bool DelaunayDemiPlane::onCircumCircle(const Point &p) const
{
    return std::abs(fma(vector.getX(),(p.getY() - first->getY()), - vector.getY()*(p.getX() - first->getX())) * direction) < POINT_TOLERANCE || isAligned(p, *first, *second);
}

bool DelaunayDemiPlane::isVertex(const Point *p) const
{
    return ( (*p) == (*first) || (*p) == (*second) ) || isAligned(p, first, second) ;
}

double DelaunayDemiPlane::distanceToVertex(const Point *p) const
{
    if(this->isVertex(p))
        return 0 ;
    Line l(*first, *second) ;
    Point p_ = l.projection(*p) ;
    return dist(*p, p_) ;
}


bool  DelaunayDemiPlane::isNeighbour(const DelaunayTreeItem * t) const
{

    if(t->isTriangle)
    {
        return t->isNeighbour(this) ;
    }
    else
    {
        return (first ==  t->first || first == t->first ||
                second == t->first || second == t->second) ;
    }
    return false ;
}

void DelaunayDemiPlane::insert(std::valarray<bool> & visited,std::vector<DelaunayTreeItem *> & ret, Point *p, Star *s)
{
    if (visited[index])
        return ;
    visited[index] = true ;

    std::vector<Point*> lims ;

    for(size_t i = 0 ; i < neighbour.size() ; i++)
    {
        if(numberOfCommonVertices(getNeighbour(i)) == 2)
        {
            std::pair< Point*,  Point*> pp = getNeighbour(i)->commonEdge(this) ;
            if(std::find(lims.begin(), lims.end(), pp.first) != lims.end())
                lims.erase(std::find(lims.begin(), lims.end(), pp.first)) ;
            else
                lims.push_back(pp.first) ;

            if(std::find(lims.begin(), lims.end(), pp.second) != lims.end())
                lims.erase(std::find(lims.begin(), lims.end(), pp.second)) ;
            else
                lims.push_back(pp.second) ;

            if (!visited[neighbour[i]] && (!getNeighbour(i)->inCircumCircle(*p) || getNeighbour(i)->onCircumCircle(*p)))
            {

                assert(getNeighbour(i)->isNeighbour(this) );

                if(sqrt(Triangle(*p, *pp.first, *pp.second).area()) > static_cast<DelaunayTriangle *>(getNeighbour(i))->getRadius()*std::numeric_limits<double>::epsilon())
                {
                    DelaunayTriangle *ss = new DelaunayTriangle(this->tree, this, p, pp.first, pp.second , p) ;

                    addSon(ss) ;

                    getNeighbour(i)->addStepson(ss) ;

                    ret.push_back(ss) ;
                }
            }
        }
    }


    if(lims.empty())
    {
        lims.push_back(first) ;
        lims.push_back(second) ;
    }

    if(!isAligned(lims[0], p, lims[1]))
    {
        this->kill(p) ;
        DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this->tree,this, lims[0], p, lims[1], p) ;

        DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this->tree,this, lims[1], p, lims[0], p) ;

        addSon(p0) ;
        addSon(p1) ;

        ret.push_back(p0) ;
        ret.push_back(p1) ;

        for(size_t i = 0 ; i < son.size()-2 ; i++)
        {
            if(getSon(i)->isNeighbour(p0))
                getSon(i)->addStepson(p0) ;

            if(getSon(i)->isNeighbour(p1))
                getSon(i)->addStepson(p1) ;

            for(size_t j = i ; j < son.size()-2 ; j++)
                if(getSon(i)->isNeighbour(getSon(j)))
                    makeNeighbours(getSon(i), getSon(j)) ;
        }
    }
    else
    {
// 		this->kill(p) ;
        DelaunayDemiPlane *p0 = new DelaunayDemiPlane(this->tree,this, lims[0], third, lims[1], p) ;
        DelaunayDemiPlane *p1 = new DelaunayDemiPlane(this->tree,this, lims[1], third, lims[0], p) ;

        addSon(p0) ;
        addSon(p1) ;

        ret.push_back(p0) ;
        ret.push_back(p1) ;

        for(size_t i = 0 ; i < son.size()-2 ; i++)
        {
            if(getSon(i)->isNeighbour(p0))
                getSon(i)->addStepson(p0) ;

            if(getSon(i)->isNeighbour(p1))
                getSon(i)->addStepson(p1) ;

            for(size_t j = i ; j < son.size()-2 ; j++)
                if(getSon(i)->isNeighbour(getSon(j)))
                    makeNeighbours(getSon(i), getSon(j)) ;
        }
    }

// 	s->updateNeighbourhood() ;

    return  ;
}

void DelaunayDemiPlane::print() const
{
    std::cerr << "###############(" << first->getX() << ", " << first->getY() << ") (" <<
              second->getX() << ", " << second->getY() << ")" << "X (" << third->getX() << ", " << third->getY() << ") :: " << isAlive()  << std::endl ;
}

void makeNeighbours(DelaunayTreeItem *t0, DelaunayTreeItem *t1 )
{
    if(t0 == t1 || !t0->isAlive() || !t1->isAlive())
        return ;
    if(t0->isPlane && t1->isPlane)
        return ;
    if(!t0->isNeighbour(t1))
        return ;
    t0->addNeighbour(t1) ;
    t1->addNeighbour(t0) ;
}

void updateNeighbours(std::vector<DelaunayTreeItem *> * t)
{
    for(size_t i = 0 ; i < t->size() ; i++)
    {
        for(size_t j = i ; j < t->size() ; j++)
        {
            if((*t)[i]->isNeighbour((*t)[j]))
            {
                makeNeighbours( (*t)[i],(*t)[j] ) ;
            }
        }
    }
}

DelaunayRoot::DelaunayRoot(Mesh<DelaunayTriangle, DelaunayTreeItem> *tt, Point * p0, Point * p1, Point * p2) : DelaunayTreeItem(tt, nullptr, nullptr)
{
    isPlane = false ;
    isTriangle =false ;
    isDeadTriangle =false ;
    this->father = nullptr ;
    DelaunayTriangle *t = new DelaunayTriangle(tt, this, p0, p1, p2, nullptr) ;
    DelaunayDemiPlane * pl0 = new DelaunayDemiPlane(tt,this, p0, p1, p2, nullptr);
    DelaunayDemiPlane * pl1 = new DelaunayDemiPlane(tt,this, p0, p2, p1, nullptr);
    DelaunayDemiPlane * pl2 = new DelaunayDemiPlane(tt,this, p1, p2, p0, nullptr);

    makeNeighbours(t,pl0 ) ;
    makeNeighbours(t,pl1) ;
    makeNeighbours(t,pl2) ;

    addSon(t) ;
    addSon(pl0) ;
    addSon(pl1) ;
    addSon(pl2) ;
    kill(p0) ;
}

void DelaunayRoot::print() const
{
    std::cerr << "I am root !" << std::endl ;
}

bool DelaunayRoot::isVertex(const Point *p) const
{
    return false ;
}

bool DelaunayRoot::inCircumCircle(const Point & p) const
{
    return true ;
}

bool DelaunayRoot::onCircumCircle(const Point & p) const
{
    return true ;
}


std::pair< Point*,  Point*> DelaunayRoot::nearestEdge(const Point & p) const
{
    assert(false) ;
    return std::pair< Point*,  Point*>(nullptr, nullptr) ;
}

void DelaunayRoot::insert(std::valarray<bool> & visited, std::vector<DelaunayTreeItem *> & ret,Point *p, Star *s)
{

    for (size_t i  = 0 ;  i < son.size() ; i++)
    {

        std::vector<DelaunayTreeItem *> temp ;
        getSon(i)->insert(visited, temp, p, s) ;

        for(size_t j = 0 ; j< temp.size() ; j++)
        {
            ret.push_back(temp[j]) ;
        }

    }

// 	updateNeighbours(&ret) ;
    return  ;
}



void DelaunayRoot::conflicts(std::valarray<bool> & visited, std::vector<DelaunayTriangle *> & ret, const Geometry *g)
{

    visited[index] = true ;

    for (size_t i  = 0 ;  i < son.size() ; i++)
    {
        std::vector<DelaunayTreeItem *> toTest ;
        getSon(i)->flatConflicts(visited,toTest,ret,g) ;
        while(!toTest.empty())
        {
            std::vector<DelaunayTreeItem *> tempToTest ;
            for(size_t j  = 0 ;  j < toTest.size() ; j++)
            {
                toTest[j]->flatConflicts(visited,tempToTest,ret,g) ;
            }

            toTest = tempToTest ;
        }

    }

    return ;
}


void DelaunayRoot::conflicts(std::valarray<bool> & visited, std::vector<DelaunayTreeItem *> & ret,const Point *p)
{

    visited[index] = true ;

// #pragma omp parallel for
    for (size_t i  = 0 ;  i < son.size() ; i++)
    {
// 		std::valarray<bool> localVisitedItems = visitedItems ;
// 		std::vector<DelaunayTreeItem *> localret ;
        if(!visited[son[i]])
        {
            std::vector<DelaunayTreeItem *> toTest ;
            getSon(i)->flatConflicts(visited,toTest,ret,p) ;
            while(!toTest.empty())
            {
                std::vector<DelaunayTreeItem *> tempToTest ;
                for(size_t j  = 0 ;  j < toTest.size() ; j++)
                {
                    toTest[j]->flatConflicts(visited,tempToTest,ret,p) ;
                }
                toTest = tempToTest ;
            }
        }
// 		#pragma omp critical
// 		{
// 			ret.insert(ret.end(), localret.begin(), localret.end());
// 		}
    }
}

Star::Star(std::vector<DelaunayTreeItem *> *t, const Point *p) : creator(p)
{
    treeitem.insert(treeitem.end(), t->begin(), t->end()) ;
    for(size_t i = 0 ; i < t->size() ; i++)
    {
        if((*t)[i]->isTriangle)
        {
            this->edge.push_back((*t)[i]->first) ;
            this->edge.push_back((*t)[i]->second) ;
            this->edge.push_back((*t)[i]->third) ;
        }
        else if((*t)[i]->isPlane)
        {
            this->edge.push_back((*t)[i]->first) ;
            this->edge.push_back((*t)[i]->second) ;
        }
    }
    std::sort(edge.begin(), edge.end()) ;
    edge.erase(unique(edge.begin(), edge.end()), edge.end()) ;
}

bool Star::contains(DelaunayTreeItem * t) const 
{
   return (std::find(treeitem.begin(), treeitem.end(), t) != treeitem.end() ) ;
}

size_t Star::size()
{
    return edge.size() ;
}

const Point * Star::getEdge(size_t i) const
{
    assert(i < edge.size()) ;
    return this->edge[i] ;
}


DelaunayTree::DelaunayTree(Point * p0, Point *p1, Point *p2): Mesh< Amie::DelaunayTriangle, Amie::DelaunayTreeItem >(SPACE_TWO_DIMENSIONAL)
{
    neighbourhood = false ;
    this->global_counter = 3;
    p0->setId(0) ;
    p1->setId(1)  ;
    p2->setId(2);
    DelaunayRoot *root = new DelaunayRoot( this, p0, p1, p2) ;
    plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(1))) ;
    plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(2))) ;
    plane.push_back(static_cast<DelaunayDemiPlane *>(root->getSon(3))) ;
}

DelaunayTree::~DelaunayTree()
{
    PointArray nularray ;
    for(size_t i = 0 ;  i < this->tree.size() ; i++)
    {
        if(this->tree[i] && this->tree[i]->isTriangle && !this->tree[i]->isDeadTriangle )
        {
            DelaunayTriangle * t = static_cast<DelaunayTriangle *>(tree[i]) ;

            t->setBoundingPoints(nularray) ;
        }

    }

    for(size_t i = 0 ;  i < this->tree.size() ; i++)
    {
        delete this->tree[i] ;
    }

    for(size_t i = 0 ; i < additionalPoints.size() ; i++)
        delete additionalPoints[i] ;

} 

void DelaunayTree::insertIf( Point *p, std::vector<SamplingCriterion *> v, double minScore )
{
    std::vector<DelaunayTreeItem *> cons = conflicts(p) ;
    std::valarray<bool> visited(false, size()) ;

    //We store all the pointers to the affected elements of the tree, and make copies of those
    std::vector<DelaunayTreeItem *> backup = tree ;


    Star * s = new Star(&cons, p) ;

    std::vector<DelaunayTreeItem *> ret ;

    for(size_t i = 0 ; i < cons.size() ; i++)
    {
        std::vector<DelaunayTreeItem *> temp ;
        cons[i]->insert(visited,temp,p, s) ;
        ret.insert(ret.end(), temp.begin(), temp.end()) ;
    }

    for(size_t i =  0 ; i < ret.size() ; i++)
    {
        if(ret[i]->isTriangle)
        {
            double score = 0;
            for(size_t j = 0 ; j < v.size() ; j++)
                if (v[j]->meetsCriterion((DelaunayTriangle *)(ret[i])))
                    score += 1/v.size() ;

            if (score < minScore)
            {
                tree = backup ;
                delete s ;
                for(size_t j =  0 ; j < ret.size() ; j++)
                {
                    delete ret[j] ;
                }

                return ;
            }
            else
            {
                p->setId(this->global_counter++) ;
                neighbourhood = false ;
            }
        }
    }

    bool weGotPlanes = false ;

    for(size_t j = 0 ; j< ret.size() ; j++)
    {
        if(ret[j]->isAlive() && ret[j]->isPlane )
        {
            weGotPlanes = true ;

            if(ret[j]->neighbour.size() > 1 || ret[j]->neighbour.size() == 0)
            {
                ret[j]->kill(p) ;
            }
            else
            {
                plane.push_back((DelaunayDemiPlane*)(ret[j])) ;
            }
        }
    }

    if(weGotPlanes)
    {
        for(int k = 0 ; k< (int)plane.size()-1 ; k++)
        {
            for(int j = k ; j< (int)plane.size() ; j++)
            {
                plane[j]->merge(plane[k]) ;
            }
        }
    }

    for(size_t i = 0 ; i < ret.size() ; i++)
    {

        assert(ret[i]->isPlane  || ret[i]->neighbour.size() == 3) ;

    }

    std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
    plane.clear() ;
    plane.insert(plane.end(), hull->begin(), hull->end()) ;
    delete hull ;
    delete s ;
}

std::vector<DelaunayTreeItem *> DelaunayTree::addElements(std::vector<DelaunayTreeItem *> & cons, Point * p)
{
    p->setId(global_counter++) ;
    std::valarray<bool> visited(false, size()) ;

    Star * s = new Star(&cons, p) ;

    std::vector<DelaunayTreeItem *> ret ;

    std::valarray<bool> in(cons.size()) ;
    std::valarray<bool> out(cons.size()) ;
    in = false ;
    out = false ;

    for(size_t i = 0 ; i < cons.size() ; i++)
        in[i] = !cons[i]->onCircumCircle(*p) ;

    if(cons.size() > 3)
    {
        for(size_t i = 0 ; i < cons.size() ; i++)
        {
            if(!in[i] && cons[i]->isTriangle)
            {
                for(size_t j = 0 ; j < cons.size() ; j++)
                {
                    if( i != j && cons[j]->isTriangle && cons[i]->isNeighbour(cons[j]))
                    {
                        bool alone = true ;
                        for(size_t k = 0 ; k < cons.size() && alone ; k++)
                        {
                            if( i != j && i != k && cons[k]->isTriangle && cons[j]->isNeighbour(cons[k]) )
                                alone = false ;
                        }
                        if(alone)
                            out[j] = true ;
                    }
                }
            }
        }
    }
    if(cons.size() == 1)
       in[0] = true ;

    size_t realInCount = 0 ;
    for(size_t i = 0 ; i < in.size() ; i++)
        realInCount += (in[i] && !out[i]) ;

    if(realInCount > 0)
    {
        for(size_t i = 0 ; i < in.size() ; i++)
        {
            if(in[i] && out[i] && cons[i]->isTriangle && cons.size() > 1)
                in[i] = false ;
        }
    }

    for(size_t i = 0 ; i < cons.size() ; i++)
    {

        if(in[i])
        {
            std::vector<DelaunayTreeItem *> temp ;
            cons[i]->insert(visited,temp,p, s) ;
            ret.insert(ret.end(), temp.begin(), temp.end()) ;
        }
    }

    s->updateNeighbourhood() ;

    bool weGotPlanes = false ;

    for(size_t j = 0 ; j< ret.size() ; j++)
    {
        if(!ret[j]->isPlane && ret[j]->isAlive())
        {
            if(ret[j]->neighbour.size() != 3 )
            {
                falseTopology = true ;
            }
        }
        if(ret[j]->isAlive() && ret[j]->isPlane )
        {
            weGotPlanes = true ;
            plane.push_back((DelaunayDemiPlane*)(ret[j])) ;
        }
    }

    if(weGotPlanes)
    {
        for(int k = 0 ; k< (int)plane.size()-1 ; k++)
        {
            for(int j = k ; j< (int)plane.size() ; j++)
            {
                plane[j]->merge(plane[k]) ;
            }
        }
    }

    for(size_t i = 0 ; i < cons.size() ; i++)
    {
        if(in[i])
            cons[i]->kill(p) ;
    }
    std::vector<DelaunayDemiPlane *> * hull = this->getConvexHull() ;
    plane.clear() ;
    plane.insert(plane.end(), hull->begin(), hull->end()) ;

    for(size_t i = 0 ; i < cons.size() ; i++)
    {
        if(!cons[i]->isAlive() && cons[i]->isTriangle && !cons[i]->isDeadTriangle)
        {
            DelaunayDeadTriangle* dt = new DelaunayDeadTriangle(static_cast<DelaunayTriangle *>(cons[i])) ;
            std::valarray<Point *> nularray(0) ;
            static_cast<DelaunayTriangle *>(cons[i])->setBoundingPoints(nularray) ;
            tree[cons[i]->index] = dt ;
            delete cons[i] ;
        }
    }
    delete hull ;
    delete s ;
    return ret ;
}

void DelaunayTree::insert(Point *p, double minDist)
{
    if(falseTopology)
    {
        std::cout << "Warning: false mesh topology" << std::endl ;
        exit(0) ;
        return ;
    }
    neighbourhood = false ;
    std::vector<DelaunayTreeItem *> cons = this->conflicts(p) ;

    if(cons.empty())
    {
        std::cout << "Failed insertion : in nothing !" << std::endl ;
        return ;
    }
    for(size_t i = 0 ; i < cons.size() ; i++)
    {
        if(cons[i]->isVertex(p))
        {
            return ;
        }

        if(cons[i]->isTriangle && cons[i]->neighbour.size() != 3)
        {
            falseTopology = true ;
            return ;
        }

        else if(minDist > POINT_TOLERANCE)
        {
            if(cons[i]->distanceToVertex(p) < minDist)
            {
                return ;
            }
        }

    }

    addElements(cons, p) ;

}



std::vector<DelaunayTreeItem *> DelaunayTree::conflicts( const Point *p)
{
    if(visited.size() != size())
      visited.resize(size(), false);
    else
      visited = false ;
    
    std::vector<DelaunayTreeItem *> cons ;

    static_cast<DelaunayRoot *>(tree[0])->conflicts(visited,cons,p) ;

    for(size_t i = 0 ; i < plane.size() ; i++)
    {

        if(!visited[plane[i]->index])
        {

	    std::vector<DelaunayTreeItem *> toTest ;
            plane[i]->flatConflicts(visited,toTest,cons,p) ;
            while(!toTest.empty())
            {
                std::vector<DelaunayTreeItem *> tempToTest ;
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


std::vector<DelaunayTriangle *> DelaunayTree::conflicts(const Geometry *g)
{
    buildNeighbourhoods() ;
    std::vector<DelaunayTreeItem *> cons = conflicts(&g->getCenter()) ;
    DelaunayTriangle * origin = nullptr ;
    for(size_t i = 0 ; i < cons.size() ; i++)
    {
        if(cons[i]->isTriangle && (static_cast<DelaunayTriangle *>(cons[i])->in(g->getCenter()) || static_cast<DelaunayTriangle *>(cons[i])->intersects(g)))
        {
            origin = static_cast<DelaunayTriangle *>(cons[i]) ;
            break ;
        }
    }
    if(!origin)
    {
        std::vector<DelaunayTriangle *> cons = getTriangles(true) ;

        for(size_t i = 0 ; i < cons.size() ; i++)
        {
            if(cons[i]->in(g->getCenter()) || cons[i]->intersects(g))
            {
                origin = cons[i] ;
                break ;
            }
        }
    }

    return getNeighbouringElementsInGeometry(origin, g) ;

}

std::vector<DelaunayDemiPlane *> * DelaunayTree::getConvexHull()
{
    std::vector<DelaunayDemiPlane *> * ret = new std::vector<DelaunayDemiPlane *> ;

    for(size_t i = 0 ; i < plane.size() ; i++)
    {
        if(plane[i]->isAlive())
            ret->push_back(plane[i]) ;
    }

    return ret ;
}

void DelaunayTree::buildNeighbourhoods()
{
    if(visited.size() != size())
      visited.resize(size(), false);
    else
      visited = false ;

    if(!neighbourhood)
    {
        neighbourhood = true ;
//      std::cerr << "\r building neighbourhood... element 0/" << ret.size() << std::flush ;
        for( size_t i = 0 ; i < tree.size() ; i++)
        {
            if(tree[i]->isTriangle && tree[i]->isAlive())
                static_cast<DelaunayTriangle *>(tree[i])->neighbourhood.resize(0) ;
        }



        for( size_t i = 0 ; i < tree.size() ; i++)
        {
            if(!tree[i]->isTriangle || !tree[i]->isAlive())
                continue ;
//          if(i%100 == 0)
//              std::cerr << "\r building neighbourhood... element " << i <<"/" << ret.size() << std::flush ;

            std::vector<DelaunayTriangle *> tocheck ;
            std::vector<DelaunayTriangle *> toclean ;
            for(size_t j = 0 ; j < tree[i]->neighbour.size() ; j++)
            {
                if(tree[i]->getNeighbour(j)->isTriangle && !visited[tree[i]->neighbour[j]])
                {
                    tocheck.push_back(static_cast<DelaunayTriangle *>(tree[i]->getNeighbour(j)));
                    visited[tree[i]->neighbour[j]] = true ;
                    toclean.push_back(*tocheck.rbegin()) ;
                    static_cast<DelaunayTriangle *>(tree[i])->addNeighbourhood(*tocheck.rbegin()) ;
                }
            }

            while(!tocheck.empty())
            {
                std::vector<DelaunayTriangle *> tocheck_temp ;
                for(size_t k = 0 ; k < tocheck.size() ; k++)
                {
                    for(size_t j = 0 ; j < tocheck[k]->neighbour.size() ; j++)
                    {
                        if(
                            tocheck[k]->getNeighbour(j)->isTriangle
                            && !visited[tocheck[k]->neighbour[j]]
                            && tocheck[k]->getNeighbour(j) != tree[i]
                            && static_cast<DelaunayTriangle *>(tocheck[k]->getNeighbour(j))->isInNeighbourhood(static_cast<DelaunayTriangle *>(tree[i]))
                        )
                        {
                            tocheck_temp.push_back(static_cast<DelaunayTriangle *>(tocheck[k]->getNeighbour(j)));
                            visited[tocheck[k]->neighbour[j]] = true ;
                            toclean.push_back(*tocheck_temp.rbegin()) ;
                            static_cast<DelaunayTriangle *>(tree[i])->addNeighbourhood(*tocheck_temp.rbegin()) ;
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


std::vector<DelaunayTriangle *> DelaunayTree::getTriangles(bool buildNeighbourhood)
{
    std::vector<DelaunayTriangle *> ret ;
    for(size_t i = 0 ; i < tree.size() ; i++)
    {
        if(tree[i]->isTriangle && !tree[i]->isDeadTriangle)
        {
            ret.push_back((DelaunayTriangle *)(tree[i])) ;
        }
    }
    buildNeighbourhoods() ;
    

    return ret ;
}

void DelaunayTree::print() const
{
    size_t alive = 0 ;
    std::cerr << "we have a total of " << tree.size() << " elements" << std::endl ;
    for(size_t i = 0 ; i < tree.size() ; i++)
    {
        if (tree[i]->isAlive())
        {
            alive++ ;
        }
    }
    std::cerr << "of which " << alive << " are alive" << std::endl ;
// 	#ifdef DEBUG
    for(size_t i = 0 ; i < tree.size() ; i++)
    {
        if (tree[i]->isAlive())
        {
            tree[i]->print() ;
            for(size_t j = 0 ; j< tree[i]->neighbour.size() ; j++)
            {
                std::cerr << "\t --> " ;
                tree[i]->getNeighbour(j)->print() ;
            }
        }
    }
// 	#endif
}

void DelaunayTriangle::refresh(const TriElement * father)
{
    TriElement::refresh(father) ;

// 	if(!cachedGPs)
// 		cachedGPs = new GaussPointArray(getSubTriangulatedGaussPoints()) ;

// 	this->computeCenter() ;
}

std::valarray<std::valarray<Matrix> > & DelaunayTriangle::getElementaryMatrix(VirtualMachine * vm )
{
    size_t dofCount = getShapeFunctions().size()+getEnrichmentFunctions().size() ;

    if(!behaviourUpdated && !enrichmentUpdated && cachedElementaryMatrix.size() != 0 && cachedElementaryMatrix[0].size() == dofCount)
    {
        return cachedElementaryMatrix ;
    }
    needAssembly = true ;
    clearBoundaryConditions() ;
    size_t ndofs = getBehaviour()->getNumberOfDegreesOfFreedom() ;

    if(        cachedElementaryMatrix.size() == 0
            || cachedElementaryMatrix.size() != dofCount
            ||  (cachedElementaryMatrix.size() && cachedElementaryMatrix[0].size() != dofCount))
    {

        std::valarray< Matrix > v_j(Matrix(ndofs, ndofs), dofCount) ;
        cachedElementaryMatrix.resize(dofCount,v_j) ;
        getSubTriangulatedGaussPoints() ;
    }
    bool cleanup = false ;
    if(!vm)
    {
        cleanup = true ;
        vm = new VirtualMachine() ;
    }

    if(!getState().JinvCache)
        getState().updateInverseJacobianCache(Point( 1./3.,1./3.)) ;   
//     std::cout <<"\n -> "<< getGaussPoints().gaussPoints.size() << std::endl ;
    
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

    if(behaviour->isSymmetric())
    {
        Matrix elem(getShapeFunctions().size()*2., getShapeFunctions().size()*2.);
        for(size_t i = start ; i < getShapeFunctions().size() ; i++)
        {
            behaviour->apply(getShapeFunction(i), getShapeFunction(i),getGaussPoints(), Jinv, cachedElementaryMatrix[i][i], vm) ;

            for(size_t j = start ; j < i ; j++)
            {
                behaviour->apply(getShapeFunction(i), getShapeFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i][j], vm) ;
                cachedElementaryMatrix[i][j].transpose(cachedElementaryMatrix[j][i]) ;
                
            }
            for(size_t j = startEnriched ; j < getEnrichmentFunctions().size() ; j++)
            {  
                behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                cachedElementaryMatrix[i][j+getShapeFunctions().size()].transpose(cachedElementaryMatrix[j+getShapeFunctions().size()][i]) ;
            }
        }

        for(size_t i = startEnriched ; i < getEnrichmentFunctions().size() ; i++)
        {
            behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = startEnriched ; j < i ; j++)
            {
                behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()].transpose(cachedElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()]) ;
            }
        }
   
    }
    else
    {
    
        for(size_t i = start ; i < getShapeFunctions().size() ; i++)
        {

            behaviour->apply(getShapeFunction(i), getShapeFunction(i),getGaussPoints(), Jinv, cachedElementaryMatrix[i][i], vm) ;

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

    if(behaviour->isViscous())
        getViscousElementaryMatrix(vm) ;

    enrichmentUpdated = false ;
    behaviourUpdated = false ;
    
    if(behaviour->hasInducedForces())
        cachedForces.resize(0) ;

    if(cleanup)
        delete vm ;

    return cachedElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > & DelaunayTriangle::getViscousElementaryMatrix(VirtualMachine * vm)
{
    size_t dofCount = getShapeFunctions().size()+getEnrichmentFunctions().size() ;

    if( !behaviourViscoUpdated && !enrichmentUpdated && cachedViscousElementaryMatrix.size() && cachedViscousElementaryMatrix[0].size() == dofCount)
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

    getState().updateInverseJacobianCache(Point(1./3, 1./3)) ;
    
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
// 		#pragma omp parallel for
        for(size_t i = start ; i < getShapeFunctions().size() ; i++)
        {

            
            behaviour->applyViscous(getShapeFunction(i), getShapeFunction(i),getGaussPoints(), Jinv, cachedViscousElementaryMatrix[i][i], vm) ;

            for(size_t j = 0 ; j < i ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getShapeFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j], vm) ;
                cachedViscousElementaryMatrix[i][j].transpose(cachedViscousElementaryMatrix[j][i]) ;
            }
            for(size_t j = startEnriched ; j < getEnrichmentFunctions().size() ; j++)
            {
                behaviour->applyViscous(getShapeFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i][j+getShapeFunctions().size()], vm) ;
                cachedViscousElementaryMatrix[i][j+getShapeFunctions().size()].transpose(cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i]) ;
            }
        }
// 		#pragma omp parallel for
        for(size_t i = startEnriched ; i < getEnrichmentFunctions().size() ; i++)
        {

            behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(i),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], vm) ;

            for(size_t j = startEnriched ; j <  i ; j++)
            {
                behaviour->applyViscous(getEnrichmentFunction(i), getEnrichmentFunction(j),getGaussPoints(), Jinv,cachedViscousElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], vm) ;
                cachedViscousElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()].transpose(cachedViscousElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()]) ;
            }
        }
    }
    else
    {
// 		#pragma omp parallel for
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

// 		#pragma omp parallel for
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
    behaviourViscoUpdated = false ;

    if(behaviour->hasInducedForces())
        cachedForces.resize(0) ;

    if(cleanup)
        delete vm ;
    
    return cachedViscousElementaryMatrix ;
}

void DelaunayTriangle::scaleCachedViscousElementaryMatrix(double s)
{
    for(size_t i = 0 ; i < cachedViscousElementaryMatrix.size() ; i++)
    {
        for(size_t j = 0 ; j < cachedViscousElementaryMatrix[i].size() ; j++)
        {
            cachedViscousElementaryMatrix[i][j] *= s ;
        }
    }
}

void DelaunayTriangle::adjustElementaryMatrix(double previousTimeStep, double nextTimeStep)
{
    if(std::abs(previousTimeStep- nextTimeStep) < POINT_TOLERANCE )
        return ;

    if( getBehaviour() && !getBehaviour()->isViscous() )
        return ;

    getState().updateInverseJacobianCache(Point( 1./3.,1./3.)) ;

//    getState().

    if( ! this->getBehaviour()->timeDependent() && ! this->getBehaviour()->spaceDependent() && getEnrichmentFunctions().size() == 0)
    {
//        scaleCachedElementaryMatrix( previousTimeStep / nextTimeStep) ;
        scaleCachedViscousElementaryMatrix( /*(previousTimeStep / nextTimeStep) */ (previousTimeStep / nextTimeStep)) ;
    }
    else
    {
        clearElementaryMatrix() ;
// 		getElementaryMatrix() ;
// 		getViscousElementaryMatrix() ;
    }
}

void DelaunayTriangle::scaleCachedElementaryMatrix(double s)
{
    for(size_t i = 0 ; i < cachedElementaryMatrix.size() ; i++)
    {
        for(size_t j = 0 ; j < cachedElementaryMatrix[i].size() ; j++)
        {
            cachedElementaryMatrix[i][j] *= s ;
        }
    }
}

std::vector<Point *> DelaunayTriangle::getIntegrationHints() const
{
    std::vector<Point *> to_add ;
    to_add.push_back(new Point(0,1)) ;
    to_add.push_back(new Point(0,0)) ;
    to_add.push_back(new Point(1,0)) ;

    for(size_t i = 0 ; i <  getEnrichmentFunctions().size() ; i++)
    {
        for(size_t j = 0 ; j < getEnrichmentFunction(i).getIntegrationHint().size() ; j++)
        {
            to_add.push_back(new Point(getEnrichmentFunction(i).getIntegrationHint(j))) ;
        }
    }

    return to_add ;
}

std::vector<std::pair<Point, double> > monteCarloGaussPoints(int nPoints, DelaunayTriangle *t, double fact  =1)
{
        std::vector<std::pair<Point, double> > gp_alternative ;
        TriElement father(LINEAR) ;

        size_t target = nPoints ;

        while(gp_alternative.size() < target)
        {

            Point test = Point((double)rand()/RAND_MAX,(double)rand()/RAND_MAX);

            if( father.in( test ) )
            {
                gp_alternative.push_back( std::make_pair( test, 0.5 ) ) ;
            }
        }

        for( size_t i = 0 ; i < gp_alternative.size() ; i++ )
            gp_alternative[i].second = fact*t->area()/gp_alternative.size() ;
        
        return gp_alternative ;
}

std::vector<std::pair<Point, double> > monteCarloGaussPoints(size_t nPoints, DelaunayTriangle *t, const Geometry * enr, double fact  =1)
{
        std::vector<std::pair<Point, double> > gp_alternative_in ;
        std::vector<std::pair<Point, double> > gp_alternative_out ;
        TriElement father(LINEAR) ;

        Function xtrans = t->getXTransform() ;
        Function ytrans = t->getYTransform() ;
        VirtualMachine vm ;
        double infrac = 0 ;
        double count = 0 ;
        while(count < 32000 || (gp_alternative_out.size() < nPoints && gp_alternative_in.size() < nPoints) )
        {
            Point test = Point((double)std::rand()/RAND_MAX,(double)std::rand()/RAND_MAX);

            if( father.in( test ) )
            {
                count++ ;
 
                Point gtest = Point(vm.eval(xtrans, test.getX(), test.getY()), vm.eval(ytrans, test.getX(), test.getY())) ;
                if(enr->in(gtest))
                {
                    infrac++ ;
                    if(gp_alternative_in.size() < nPoints)
                    {
                        gp_alternative_in.push_back( std::make_pair( test, 0.5 ) ) ;
                    }
                }
                else if(gp_alternative_out.size() < nPoints)
                    gp_alternative_out.push_back( std::make_pair( test, 0.5 ) ) ;
            }
        }
        
        double f = infrac/count ;

        for( size_t i = 0 ; i < gp_alternative_in.size() ; i++ )
            gp_alternative_in[i].second = f*fact/( gp_alternative_in.size()) ;
        for( size_t i = 0 ; i < gp_alternative_out.size() ; i++ )
            gp_alternative_out[i].second = (1.-f)*fact/(gp_alternative_out.size()) ;
        
        gp_alternative_in.insert(gp_alternative_in.end(), gp_alternative_out.begin(),gp_alternative_out.end());
        double sum = 0 ;
        for( size_t i = 0 ; i < gp_alternative_in.size() ; i++ )
            sum += gp_alternative_in[i].second ;
        double j = t->area() ;
        for( size_t i = 0 ; i < gp_alternative_in.size() ; i++ )
            gp_alternative_in[i].second *= j/sum ;
        
        return gp_alternative_in ;
}

const GaussPointArray & DelaunayTriangle::getSubTriangulatedGaussPoints()
{
    if(!enrichmentUpdated && getCachedGaussPoints())
    {
        return *getCachedGaussPoints() ;
    }
    size_t numberOfRefinements = 0 ;
    
//     if(getEnrichmentFunctions().size() > getBoundingPoints().size())
//         numberOfRefinements = 4 ;
    VirtualMachine vm ;
    if(getEnrichmentFunctions().size() > 0)
    {
        std::vector<std::pair<Point, double> > gp_alternative ;
        if( order >= CONSTANT_TIME_LINEAR )
        {

            if((getCachedGaussPoints() != nullptr) && (getCachedGaussPoints()->regularGrid ))
                return *getCachedGaussPoints() ;


            Point A(0,1) ;
            Point B(0,0) ;
            Point C(1,0) ;
            TriElement father (LINEAR) ;
            TriangularInclusion trg(A,B,C) ;

            trg.sample(1./16., 1.) ;
            DelaunayTree * dt = new DelaunayTree(&A,&B,&C) ;

            for(size_t i = 0 ; i < trg.getBoundingPoints().size() ; i++)
            {
                dt->insert(&trg.getBoundingPoint(i), 0);
            }
            for(size_t i = 0 ; i < trg.getInPoints().size() ; i++)
            {
                dt->insert(&trg.getInPoint(i), 0);
            }

            std::vector<DelaunayTriangle *> tris = dt->getElements() ;
            double fsum = 0 ;
            for(size_t i = 0 ; i < tris.size() ; i++)
            {
                auto gpl = tris[i]->getGaussPoints() ;
                tris[i]->refresh(&father) ;
                Function xtr = tris[i]->getXTransform() ;
                Function ytr = tris[i]->getYTransform() ;
                for(size_t j = 0 ; j < gpl.gaussPoints.size() ; j++)
                {
                    Point c(vm.eval(xtr,gpl.gaussPoints[j].first.getX(), gpl.gaussPoints[j].first.getY()), vm.eval(ytr,gpl.gaussPoints[j].first.getX(), gpl.gaussPoints[j].first.getY())) ;
                    double ar = 2.*tris[i]->area() ;
                    gp_alternative.push_back(std::make_pair(Point(c.getX(),c.getY(),0,-1), gpl.gaussPoints[j].second*ar*1./3.));
                    gp_alternative.push_back(std::make_pair(Point(c.getX(),c.getY(),0, 0), gpl.gaussPoints[j].second*ar*4./3.));
                    gp_alternative.push_back(std::make_pair(Point(c.getX(),c.getY(),0, 1), gpl.gaussPoints[j].second*ar*1./3.));

                }

            }

//            double jac = 2.*area() ;

//             for(size_t i = 0 ; i < gp_alternative.size() ; i++)
//             {
//                 gp_alternative[i].second *= jac ;
//                 fsum += gp_alternative[i].second ;
// 
//             }

            double originalSum = 0 ;
            for(auto & i : getGaussPoints().gaussPoints)
                originalSum += i.second ;

            for(size_t i = 0 ; i < gp_alternative.size() ; i++)
            {
                gp_alternative[i].second *= originalSum/fsum ;
            }

            delete dt ;
            if(cachedGps)
                *cachedGps = gp_alternative ;
            else
                cachedGps = new GaussPointArray(gp_alternative) ;
            cachedGps->regularGrid = true ;      
            return *getCachedGaussPoints();

        }
        else
        {

            if( false )
            {
                if(getCachedGaussPoints() && getCachedGaussPoints()->regularGrid)
                    return *getCachedGaussPoints() ;

                delete cachedGps ;
                cachedGps = new GaussPointArray(monteCarloGaussPoints(256, this)) ;
                cachedGps->regularGrid = true ;
                return *getCachedGaussPoints() ;
            }

            VirtualMachine vm ;
            std::vector<Point *> to_add = getIntegrationHints();
            std::vector<Point *> pointsToCleanup = to_add;
            std::vector<DelaunayTriangle *> triangleToCleanup;
            std::vector<DelaunayTriangle *> tri ;

            DelaunayTree * dt = new DelaunayTree(to_add[0], to_add[1], to_add[2]) ;
            
            TriElement f(QUADRATIC) ;
            if(to_add.size() > 4)
                std::random_shuffle(to_add.begin()+3, to_add.end()) ;

            for(size_t i = 3 ; i < to_add.size() ; i++)
            {
                dt->insert(to_add[i], .05) ;
            }

            dt->setElementOrder(QUADRATIC); 
            dt->refresh(&f);
            tri = dt->getTriangles() ;

            Function gx = getXTransform() ;
            Function gy = getYTransform() ;
            double parentArea = 0 ;
            for(size_t i = 0 ; i < tri.size() ; i++)
            {
                
                Function x = tri[i]->getXTransform() ;
                Function y = tri[i]->getYTransform() ;
                               
                GaussPointArray gp_temp = tri[i]->getGaussPoints() ;
//                  GaussPointArray gp_temp = monteCarloGaussPoints(512, tri[i]) ;

                for(size_t j = 0 ; j < gp_temp.gaussPoints.size() ; j++)
                {
                    gp_temp.gaussPoints[j].first.set(vm.eval(x, gp_temp.gaussPoints[j].first), vm.eval(y, gp_temp.gaussPoints[j].first)) ;
                    gp_alternative.push_back(gp_temp.gaussPoints[j]) ;
                    parentArea += gp_temp.gaussPoints[j].second ;
                }
            }
            if(tri.size() < 3 || std::abs(parentArea - 0.5) > 1e-6)
            {
                if(getCachedGaussPoints() && getCachedGaussPoints()->regularGrid)
                    return *getCachedGaussPoints() ;
                delete cachedGps ;
                cachedGps = new GaussPointArray(monteCarloGaussPoints(512, this)) ;
                cachedGps->regularGrid = true ;
                return *getCachedGaussPoints() ;
            }

            double originalSum = area() ;

            for(size_t i = 0 ; i < gp_alternative.size() ; i++)
                gp_alternative[i].second *= originalSum / parentArea ;

            delete dt ;

            if(numberOfRefinements)
            {

                std::valarray<Point *> nularray(0) ;

                for(size_t i = 0 ; i < triangleToCleanup.size() ; i++)
                {
                    triangleToCleanup[i]->setBoundingPoints(nularray) ;
                    delete triangleToCleanup[i];
                }

            }
            for(size_t i = 0 ; i < pointsToCleanup.size() ; i++)
                delete pointsToCleanup[i] ;
        }

        if(cachedGps)
            *cachedGps = gp_alternative ;
        else
            cachedGps = new GaussPointArray(gp_alternative) ;
        return *getCachedGaussPoints();

    }


    setCachedGaussPoints(new GaussPointArray(getGaussPoints()));
    return *getCachedGaussPoints();
}


DelaunayDeadTriangle::~DelaunayDeadTriangle() { }

std::pair<std::vector<DelaunayTriangle *>, std::vector<Point *> > quad(const DelaunayTriangle * t)
{
    std::vector<DelaunayTriangle* > tris ;
    std::vector<Point *> points ;

    points.push_back(new Point(*t->first + (*t->second - *t->first)*.5) ) ;
    points.push_back(new Point(*t->first + (*t->third - *t->first)*.5 )) ;
    points.push_back(new Point(*t->second + (*t->third - *t->second)*.5 )) ;

    tris.push_back(new DelaunayTriangle(nullptr,nullptr, points[0], points[1], t->first, nullptr)) ;
    tris.push_back(new DelaunayTriangle(nullptr,nullptr,points[0], points[2], t->second, nullptr)) ;
    tris.push_back(new DelaunayTriangle(nullptr,nullptr,points[1], points[2], t->third, nullptr)) ;
    tris.push_back(new DelaunayTriangle(nullptr,nullptr,points[0], points[1], points[2], nullptr)) ;
    return std::make_pair(tris, points) ;
}

}
