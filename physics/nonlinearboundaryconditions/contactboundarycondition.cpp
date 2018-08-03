
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
    
ContactBoundaryCondition::ContactBoundaryCondition(Geometry * geo) : baseGeometry(geo), conv(false)
{

}

void ContactBoundaryCondition::setScale(double s)
{
    scale = s;
}

void ContactBoundaryCondition::initialise(Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh)
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
            if(pt0 == contactPointsAndTributary.end() && baseGeometry->in(*commonSurface.first))
            {
                contactPointsAndTributary[commonSurface.first] = d*.5 ;
                affectedElements[commonSurface.first] = element ;
            }
            else if(baseGeometry->in(*commonSurface.first))
            {
               pt0->second += d*.5 ; 
            }
            
            if(pt1 == contactPointsAndTributary.end()&& baseGeometry->in(*commonSurface.second))
            {
                contactPointsAndTributary[commonSurface.second] = d*.5 ;
                affectedElements[commonSurface.second] = element ;
            }
            else if(baseGeometry->in(*commonSurface.second))
            {
               pt1->second += d*.5 ; 
            }

        }
    }
    
    //now we have all the relevant points. Onto the normals
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        Point proj(*pt->first) ;
        baseGeometry->project(&proj);
        double inSign = (baseGeometry->in(*pt->first))?-1.:1 ;
        double pdist = inSign*dist(&proj, pt->first) ;
        Point nbase = proj - *pt->first ;
        double nbasenorm = nbase.norm() ;
        if(nbasenorm  < POINT_TOLERANCE)
        {
           nbase =  proj - (*pt->first+baseGeometry->getCenter())*.5 ;
           nbasenorm = nbase.norm() ;
           inSign = -1 ;
        }
        normalVectors[pt->first] = nbase/nbasenorm*inSign ;        
        referencePoints[pt->first] = *pt->first+normalVectors[pt->first]*(baselength) ;
        stiffnesses[pt->first] = pt->second ;
        
    }
    
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<<stiffnesses[pt->first] << std::endl;
    }
    
}


void ContactBoundaryCondition::reInitialise()
{
    contactPointsAndTributary.clear() ;
    affectedElements.clear() ;
    stiffnesses.clear() ;
    normalVectors.clear() ;
    referencePoints.clear() ;
    
    for(auto element : edgeElements)
    {
        DelaunayTreeItem * VoidItem = nullptr;
        for ( size_t j = 0 ; j < element->neighbour.size() ; j++ )
        {
            bool voidNeighbour = ( element->getNeighbour ( j )->isTriangle
                                   && dynamic_cast<DelaunayTriangle *> ( element->getNeighbour ( j ) )->getBehaviour()->type == VOID_BEHAVIOUR ) ;

            if ( voidNeighbour )
            {
                VoidItem = element->getNeighbour ( j ) ;
            }

            if ( element->getNeighbour ( j )->isPlane )
            {
                VoidItem = element->getNeighbour ( j ) ;
            }
        }
        
        std::pair<Point *, Point*> commonSurface = element->commonEdge ( VoidItem ) ;
        double d = dist(commonSurface.first, commonSurface.second) ;
        auto pt0 = contactPointsAndTributary.find(commonSurface.first) ;
        auto pt1 = contactPointsAndTributary.find(commonSurface.second) ;
        if(pt0 == contactPointsAndTributary.end() && baseGeometry->in(*commonSurface.first))
        {
            contactPointsAndTributary[commonSurface.first] = d*.5 ;
            affectedElements[commonSurface.first] = element ;
        }
        else if(baseGeometry->in(*commonSurface.first))
        {
            pt0->second += d*.5 ; 
        }
        
        if(pt1 == contactPointsAndTributary.end()&& baseGeometry->in(*commonSurface.second))
        {
            contactPointsAndTributary[commonSurface.second] = d*.5 ;
            affectedElements[commonSurface.second] = element ;
        }
        else if(baseGeometry->in(*commonSurface.second))
        {
            pt1->second += d*.5 ; 
        }
    }
    
    Vector disps(2) ;
    VirtualMachine vm ;
    //now we have all the relevant points. Onto the normals
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        affectedElements[pt->first]->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
        Point movedPoint = *pt->first+disps;
        Point proj(movedPoint) ;
        baseGeometry->project(&proj);
        double inSign = (baseGeometry->in(movedPoint))?-1.:1 ;
        double pdist = inSign*dist(&proj, pt->first) ;
        Point nbase = proj - *pt->first ;
        double nbasenorm = nbase.norm() ;
        if(nbasenorm  < POINT_TOLERANCE)
        {
           nbase =  proj - (movedPoint+baseGeometry->getCenter())*.5 ;
           nbasenorm = nbase.norm() ;
           inSign = -1 ;
        }
        normalVectors[pt->first] = nbase/nbasenorm*inSign ;
        referencePoints[pt->first] = movedPoint+normalVectors[pt->first]*baselength ;
        stiffnesses[pt->first] += pt->second ;

    }
    
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<<stiffnesses[pt->first] << std::endl;
    }

}

bool ContactBoundaryCondition::converged() const
{
    return conv ;
}

void ContactBoundaryCondition::update()
{
    Vector disps(2) ;
    VirtualMachine vm ;
    
    conv = true ;
    double threshold = 1e-12 ;
    
    for(auto element: edgeElements)
    {
        for(size_t i = 0 ; i < element->getBoundingPoints().size() ; i++)
        {
            auto pt = contactPointsAndTributary.find( &element->getBoundingPoint(i)) ;
            if(pt != contactPointsAndTributary.end())
            {
                element->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
                Point movedPoint = *pt->first+disps;
                Point proj(movedPoint) ;
                baseGeometry->project(&proj); ;
                double initstf = stiffnesses[pt->first] ;
                double fac = dist(referencePoints[pt->first], movedPoint)/dist(proj, referencePoints[pt->first]) ;
                double nxtstf = initstf*fac ;
                if(nxtstf <  threshold)
                {
                    stiffnesses[pt->first] = threshold ; 
                }
                else if(std::abs((initstf-nxtstf)) / nxtstf < threshold )
                {
                    stiffnesses[pt->first] = nxtstf ;
                }
                else
                {
                    stiffnesses[pt->first] = nxtstf ;
                    conv = false ;
                }
            }
        }
    }
    
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<< stiffnesses[pt->first] << std::endl;
    }
}

void ContactBoundaryCondition::applyBoundaryConditions( Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) 
{
    std::vector<BoundaryCondition *> ret ;
    for(auto pt : contactPointsAndTributary)
    {
        double d = dist(referencePoints[pt.first], *pt.first) ;
        double force  = -d*stiffnesses[pt.first]*scale ;
        
        int JinvSize = 3 ;
        if ( affectedElements[pt.first]->spaceDimensions() == SPACE_THREE_DIMENSIONAL && affectedElements[pt.first]->timePlanes() > 1 )
            JinvSize = 4 ;
        if ( affectedElements[pt.first]->spaceDimensions() == SPACE_TWO_DIMENSIONAL && affectedElements[pt.first]->timePlanes() == 1 )
            JinvSize = 2 ;
        
        std::valarray<Matrix> Jinv ( (bool) affectedElements[pt.first]->getState().JinvCache ? (*affectedElements[pt.first]->getState().JinvCache) : Matrix( JinvSize, JinvSize ),  affectedElements[pt.first]->getGaussPoints().gaussPoints.size()) ;

        if( ! affectedElements[pt.first]->getState().JinvCache )
        {
            affectedElements[pt.first]->getState().updateInverseJacobianCache(affectedElements[pt.first]->spaceDimensions() == SPACE_THREE_DIMENSIONAL? Point(1./3.) : Point(.25, .25, .25)) ;
        }
        
        DofDefinedBoundaryCondition ycond(SET_FORCE_ETA,(ElementarySurface*)affectedElements[pt.first],affectedElements[pt.first]->getGaussPoints(), Jinv, pt.first->getId(), force*normalVectors[pt.first].getY()) ;
        
        ycond.apply ( a, mesh ) ;
        DofDefinedBoundaryCondition xcond(SET_FORCE_XI, (ElementarySurface*)affectedElements[pt.first],affectedElements[pt.first]->getGaussPoints(), Jinv, pt.first->getId(), force*normalVectors[pt.first].getX());    
        xcond.apply ( a, mesh ) ;
    }

}

}