
//
// C++ Implementation: countact
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "contactboundarycondition.h"

namespace Amie
{
    
ContactBoundaryCondition::ContactBoundaryCondition(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh, Geometry * geo) : baseGeometry(geo)
{
    for(auto element = mesh->begin() ; element != mesh->end() ; element++)
    {
        DelaunayTreeItem * VoidItem = nullptr;
        bool border = false ;

        for ( size_t j = 0 ; j < element->neighbour.size() ; j++ )
        {
            bool voidNeighbour = ( element->getNeighbour ( j )->isTriangle
                                   && dynamic_cast<DelaunayTriangle *> ( element->getNeighbour ( j ) )->getBehaviour()->type == VOID_BEHAVIOUR ) ;
            border = border || element->getNeighbour ( j )->isPlane
                     || voidNeighbour ;

            if ( voidNeighbour )
            {
                VoidItem = element->getNeighbour ( j ) ;
            }

            if ( element->getNeighbour ( j )->isPlane )
            {
                VoidItem = element->getNeighbour ( j ) ;
            }
        }

        if ( element->getBehaviour()->type == VOID_BEHAVIOUR )
        {
            border = false ;
        }

        if ( border )
        {
            edgeElements.push_back(element) ;
            std::pair<Point *, Point*> commonSurface = element->commonEdge ( VoidItem ) ;
            double d = dist(commonSurface.first, commonSurface.second) ;
            auto pt0 = contactPointsAndTributary.find(commonSurface.first) ;
            auto pt1 = contactPointsAndTributary.find(commonSurface.second) ;
            if(pt0 == contactPointsAndTributary.end())
            {
                contactPointsAndTributary[commonSurface.first] = d*.5 ;
                affectedElements[commonSurface.first] = element ;
            }
            else
            {
               pt0->second += d*.5 ; 
            }
            
            if(pt1 == contactPointsAndTributary.end())
            {
                contactPointsAndTributary[commonSurface.second] = d*.5 ;
                affectedElements[commonSurface.second] = element ;
            }
            else
            {
               pt1->second += d*.5 ; 
            }

        }
    }
    
    //now we have all the relevant points. Onto the normals
    double mindist = 1e12 ;
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        Point proj(*pt->first) ;
        baseGeometry->project(&proj);
        double inSign = (baseGeometry->in(*pt->first))?-1.:1 ;
        double pdist = inSign*dist(&proj, pt->first) ;
        mindist = std::min(mindist, pdist) ;
        Point nbase = proj - *pt->first ;
        double nbasenorm = nbase.norm() ;
        if(nbasenorm  < POINT_TOLERANCE)
        {
           nbase =  proj - (*pt->first+baseGeometry->getCenter())*.5 ;
           nbasenorm = nbase.norm() ;
           inSign = -1 ;
        }
        normalVectors[pt->first] = nbase/nbasenorm*inSign ;
        
    }
    baselength -= mindist ;
    
    //now we can define the referencePoints
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        referencePoints[pt->first] = *pt->first+normalVectors[pt->first]*(distances[pt->first]+baselength) ;
        distances[pt->first] += baselength ;
        stiffnesses[pt->first] = pt->second*1e9 ;
    }
}

void ContactBoundaryCondition::update(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh)
{
    Vector disps(2) ;
    VirtualMachine vm ;
    for(auto element: edgeElements)
    {
        for(size_t i = 0 ; i < element->getBoundingPoints().size() ; i++)
        {
            auto pt = contactPointsAndTributary.find( &element->getBoundingPoint(i)) ;
            if(pt != contactPointsAndTributary.end())
            {
                element->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
                Point movedPoint = *pt->first+disps;
                stiffnesses[pt->first] = stiffnesses[pt->first]*(dist(referencePoints[pt->first], movedPoint))/distances[pt->first] ;
            }
        }
    }
}

std::vector<BoundaryCondition *> ContactBoundaryCondition::getBoundaryconditions(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) 
{
    std::vector<BoundaryCondition *> ret ;
    for(auto pt : contactPointsAndTributary)
    {
        double d = dist(referencePoints[pt.first], *pt.first) ;
        double force  = d*stiffnesses[pt.first] ;
        
        std::valarray<Matrix> Jinv(*affectedElements[pt.first]->getState().JinvCache, affectedElements[pt.first]->getGaussPoints().gaussPoints.size());
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA,(ElementarySurface*)affectedElements[pt.first],affectedElements[pt.first]->getGaussPoints(), Jinv, pt.first->getId(), force*normalVectors[pt.first].getY())) ;
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, (ElementarySurface*)affectedElements[pt.first],affectedElements[pt.first]->getGaussPoints(), Jinv, pt.first->getId(), force*normalVectors[pt.first].getX()));                                                    
    }
    
    return ret ;
}

}
