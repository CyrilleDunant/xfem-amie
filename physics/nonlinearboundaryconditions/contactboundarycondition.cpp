
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
            
            if(VoidItem)
                break ;
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
        Point nbase = proj - *pt->first ;
        double nbasenorm = nbase.norm() ;
        if(nbasenorm  < POINT_TOLERANCE)
        {
           nbase =  proj - (*pt->first*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm = nbase.norm() ;
           inSign = -1 ;
        }
        normalVectors[pt->first] = nbase/nbasenorm*inSign ;        
        referencePoints[pt->first] = proj-normalVectors[pt->first] ;
        stiffnesses[pt->first] = 0 ;
        distances[pt->first] = dist(*pt->first, referencePoints[pt->first]) ;
    }
    
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<<stiffnesses[pt->first]*pt->second << std::endl;
    }
    
}


void ContactBoundaryCondition::reInitialise()
{
//     contactPointsAndTributary.clear() ;
//     affectedElements.clear() ;
//     stiffnesses.clear() ;
//     normalVectors.clear() ;
//     referencePoints.clear() ;
    
    Vector disps(2) ;
    VirtualMachine vm ;
    
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
            
            if(VoidItem)
                break ;
        }
        
        std::pair<Point *, Point*> commonSurface = element->commonEdge ( VoidItem ) ;
        
        auto pt0 = contactPointsAndTributary.find(commonSurface.first) ;
        auto pt1 = contactPointsAndTributary.find(commonSurface.second) ;
        
        element->getState().getField(DISPLACEMENT_FIELD, *commonSurface.first ,disps, false, &vm, 0) ;
        Point movedPoint0 = *commonSurface.first+disps;
        element->getState().getField(DISPLACEMENT_FIELD, *commonSurface.second ,disps, false, &vm, 0) ;
        Point movedPoint1 = *commonSurface.second+disps;
        double d = dist(movedPoint0, movedPoint1) ;
        
        if(pt0 == contactPointsAndTributary.end() && baseGeometry->in(movedPoint0))
        {
            contactPointsAndTributary[commonSurface.first] = d*.5 ;
            affectedElements[commonSurface.first] = element ;
        }
        else if(baseGeometry->in(movedPoint0))
        {
            pt0->second += d*.5 ; 
        }
        
        if(pt1 == contactPointsAndTributary.end() && baseGeometry->in(movedPoint1))
        {
            contactPointsAndTributary[commonSurface.second] = d*.5 ;
            affectedElements[commonSurface.second] = element ;
        }
        else if(baseGeometry->in(movedPoint1))
        {
            pt1->second += d*.5 ; 
        }
    }
    
//     for(auto pt = stiffnesses.begin() ; pt != stiffnesses.end() ; pt++)
//     {
//         auto search  = contactPointsAndTributary.find(pt->first) ;
//         if(search == contactPointsAndTributary.end())
//         {
//             stiffnesses.erase(pt) ;
//             pt = stiffnesses.begin() ;
//         }
//     }
    
    //now we have all the relevant points. Onto the normals
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        affectedElements[pt->first]->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
        Point movedPoint = *pt->first+disps;
        Point proj(movedPoint) ;
        baseGeometry->project(&proj);
        double inSign = (baseGeometry->in(movedPoint))?-1.:1 ;
        Point nbase = proj - movedPoint ;
        double nbasenorm = nbase.norm() ;
        if(nbasenorm  < POINT_TOLERANCE)
        {
           nbase =  proj - (movedPoint*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm = nbase.norm() ;
           inSign = -1 ;
        }
        
        normalVectors[pt->first] = nbase/nbasenorm*inSign ;
        
        referencePoints[pt->first] = proj-normalVectors[pt->first] ;
        distances[pt->first] = dist(*pt->first, referencePoints[pt->first]) ;
        
        if(stiffnesses.find(pt->first) == stiffnesses.end())
            stiffnesses[pt->first] = 0 ;
        else
            stiffnesses[pt->first] *= 1.+ dist(movedPoint,movedPoint)/sqrt(disps[0]*disps[0]+disps[1]*disps[1]);
  
    }
    
//     update() ;
    currentError = 0 ;
    
    errors.clear() ;
    positions.clear() ;
    previousStiffnesses.clear() ;
    previousErrors.clear() ;
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
    double threshold = 1e-6 ;
    std::set<const Point * > done ;
    currentError = 0 ;
// std::cout << std::endl;
    for(auto element: affectedElements)
    {
        for(size_t i = 0 ; i < element.second->getBoundingPoints().size() ; i++)
        {
            auto pt = contactPointsAndTributary.find( &element.second->getBoundingPoint(i)) ;
            if(pt != contactPointsAndTributary.end() && done.find(element.first) == done.end())
            {                
                element.second->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
                Point movedPoint = *pt->first+disps;
                Point proj(movedPoint) ;
                baseGeometry->project(&proj);

                double nerr = dist(proj, movedPoint) ;
                double pstiff = stiffnesses[pt->first] ;               
                double perr = errors[pt->first] ;  
                double dsderr = 1.;
                
                if(pstiff < perr || std::abs(pstiff) < 1e-7)
                {
                    pstiff = 1. ;
                }
                else if(previousStiffnesses.find(pt->first) != previousStiffnesses.end() )
                {
                    if(std::abs(perr-previousErrors[pt->first]) > 1e-7)
                        dsderr = (pstiff-previousStiffnesses[pt->first])/(perr-previousErrors[pt->first]) ;  
                }                
                previousStiffnesses[pt->first] = pstiff ;
                previousErrors[pt->first] = perr ;

                bool isIn = baseGeometry->in(movedPoint) ;
                
                //Newton
                if(isIn)
                    stiffnesses[pt->first] -= nerr/dsderr ;
                else
                    stiffnesses[pt->first] += nerr/dsderr ; 
                
                if(stiffnesses[pt->first] > 0)
                    stiffnesses[pt->first] = 0 ;

                errors[pt->first] = nerr ;
                currentError  += nerr ;
                done.insert(pt->first) ;

            }
        }
    }
    
    conv = currentError < threshold ;
  
}


void ContactBoundaryCondition::print()
{
    
    Vector disps(2) ;
    VirtualMachine vm ;
    std::cout << std::endl;
    for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
    {
        affectedElements[pt->first]->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
        Point movedPoint = *pt->first+disps;
        double d = dist(referencePoints[pt->first], movedPoint) ;
        double force  = d*stiffnesses[pt->first]*scale*pt->second ;
        
        
        std::cout << movedPoint.getX() << "  " << movedPoint.getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<< force << "  "<< errors[pt->first] << std::endl;
    }
}

double ContactBoundaryCondition::error() const
{
    return currentError ;
}

void ContactBoundaryCondition::applyBoundaryConditions( Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) 
{
    std::vector<BoundaryCondition *> ret ;
    Vector disps(2) ;
    VirtualMachine vm ;
    for(auto pt : contactPointsAndTributary)
    {
        affectedElements[pt.first]->getState().getField(DISPLACEMENT_FIELD, *pt.first,disps, false, &vm, 0) ;
        double d = dist(referencePoints[pt.first], *pt.first+disps) ;
        double force  = d*stiffnesses[pt.first]*scale/**pt.second*/ ;
        
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
