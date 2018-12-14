
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
            affectedElements[std::make_pair(commonSurface.first,commonSurface.second)] = element ;
            double d = dist(commonSurface.first, commonSurface.second) ;
            auto pt0 = contactPointsAndTributary.find(commonSurface.first) ;
            auto pt1 = contactPointsAndTributary.find(commonSurface.second) ;
            if(pt0 == contactPointsAndTributary.end() && baseGeometry->in(*commonSurface.first))
            {
                contactPointsAndTributary[commonSurface.first] = d*.5 ;
            }
            else if(baseGeometry->in(*commonSurface.first))
            {
               pt0->second += d*.5 ; 
            }
            
            if(pt1 == contactPointsAndTributary.end()&& baseGeometry->in(*commonSurface.second))
            {
                contactPointsAndTributary[commonSurface.second] = d*.5 ;
            }
            else if(baseGeometry->in(*commonSurface.second))
            {
               pt1->second += d*.5 ; 
            }

        }
    }
    
    //now we have all the relevant points. Onto the normals
    for(auto pt = affectedElements.begin() ; pt != affectedElements.end() ; pt++)
    {
        Point proj0(*pt->first.first) ;
        baseGeometry->project(&proj0);
        Point proj1(*pt->first.second) ;
        baseGeometry->project(&proj1);
        double inSign0 = (baseGeometry->in(*pt->first.first) )?-.1:.1 ;
        double inSign1 = (baseGeometry->in(*pt->first.second))?-.1:.1 ;

        Point cbase = (*pt->first.first+*pt->first.second)*.5-pt->second->getCircumCenter() ;        
        Point nbase0 = (proj0 - *pt->first.first)*.5 - cbase*.5 ;
        Point nbase1 = (proj1 - *pt->first.second)*.5 - cbase*.5 ;
        double nbasenorm0 = nbase0.norm() ;
        double nbasenorm1 = nbase1.norm() ;
        if(nbasenorm0  < POINT_TOLERANCE)
        {
           nbase0 =  proj0 - (*pt->first.first*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm0 = nbase0.norm() ;
           inSign0 = -.1 ;
        }
        
        if(nbasenorm1  < POINT_TOLERANCE)
        {
           nbase1 =  proj1 - (*pt->first.second*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm1 = nbase1.norm() ;
           inSign1 = -.1 ;
        }
        
        if(inSign0 > 0 && inSign1 > 0)
            continue ;
                
        if(nbasenorm0 < POINT_TOLERANCE)
            nbasenorm0 = 1. ;
        
        if(nbasenorm1 < POINT_TOLERANCE)
            nbasenorm1 = 1. ;
        
//         if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0 ;    
            double dx = proj0.getX()-normalVectors[pt->first.first].getX() ;
            referenceX[pt->first.first] = dx ;
            stiffnessesX[pt->first.first] = 0 ;
            double dy = proj0.getY()-normalVectors[pt->first.first].getY() ;
            referenceY[pt->first.first] = dy ;
            stiffnessesY[pt->first.first] = 0 ;
            
            normalVectors[pt->first.second] = nbase1/nbasenorm1*inSign1 ; 
            dx = proj1.getX()-normalVectors[pt->first.second].getX() ;
            referenceX[pt->first.second] = dx ;
            stiffnessesX[pt->first.second] = 0 ;
            dy = proj1.getY()-normalVectors[pt->first.second].getY() ;
            referenceY[pt->first.second] = dy ;
            stiffnessesY[pt->first.second] = 0 ;
        }
        
//         if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) == contactPointsAndTributary.end())
//         {
//             normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;   
//             double dx = proj0.getX()-normalVectors[pt->first.first].getX() ;
//             referenceX[pt->first.first] = dx ;
//             stiffnessesX[pt->first.first] = 0 ;
//             double dy = proj0.getY()-normalVectors[pt->first.first].getY() ;
//             referenceY[pt->first.first] = dy ;
//             stiffnessesY[pt->first.first] = 0 ;
//         }
//         
//         if(contactPointsAndTributary.find(pt->first.first) == contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
//         {
//             normalVectors[pt->first.second] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;    
//             double dx = proj1.getX()-normalVectors[pt->first.second].getX() ;
//             referenceX[pt->first.second] = dx ;
//             stiffnessesX[pt->first.second] = 0 ;
//             double dy = proj1.getY()-normalVectors[pt->first.second].getY() ;
//             referenceY[pt->first.second] = dy ;
//             stiffnessesY[pt->first.second] = 0 ;
//         }
        
    }
//     std::cout << "padoum" << std::endl ;
//     for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
//     {
//         std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<<stiffnesses[pt->first]*pt->second << std::endl;
//     }
    
}


void ContactBoundaryCondition::reInitialise()
{
//     if(affectedElements.size() > 3)
//         return ;
    contactPointsAndTributary.clear() ;
    affectedElements.clear() ;
    auto stiffCopyX = stiffnessesX ;
    stiffnessesX.clear() ;
    previousStiffnessesX.clear() ;
    auto stiffCopyY = stiffnessesY ;
    stiffnessesY.clear() ;
    previousStiffnessesY.clear() ;
    normalVectors.clear() ;
    referenceX.clear() ;
    referenceY.clear() ;
    
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
        affectedElements[std::make_pair(commonSurface.first, commonSurface.second)] = element ;
        
        if(pt0 == contactPointsAndTributary.end() && baseGeometry->in(movedPoint0))
        {
            contactPointsAndTributary[commonSurface.first] = d*.5 ;
        }
        else if(baseGeometry->in(movedPoint0))
        {
            pt0->second += d*.5 ; 
        }
        
        if(pt1 == contactPointsAndTributary.end() && baseGeometry->in(movedPoint1))
        {
            contactPointsAndTributary[commonSurface.second] = d*.5 ;
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
    for(auto pt = affectedElements.begin() ; pt != affectedElements.end() ; pt++)
    {
        affectedElements[pt->first]->getState().getField(DISPLACEMENT_FIELD, *pt->first.first,disps, false, &vm, 0) ;
        Point movedPoint0= *pt->first.first+disps;
        Point proj0(movedPoint0) ;
        baseGeometry->project(&proj0);
        if(stiffCopyX.find(pt->first.first) == stiffCopyX.end() && contactPointsAndTributary.find(pt->first.first) !=  contactPointsAndTributary.end())
            stiffnessesX[pt->first.first] = 0 ;
        else
            stiffnessesX[pt->first.first] = stiffCopyX[pt->first.first];

        if(stiffCopyY.find(pt->first.first) == stiffCopyY.end() && contactPointsAndTributary.find(pt->first.first) !=  contactPointsAndTributary.end())
            stiffnessesY[pt->first.first] = 0 ;
        else
            stiffnessesY[pt->first.first] = stiffCopyY[pt->first.first];
//         else if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end())
//         {
// //             stiffnesses[pt->first.first] *= 1.+ dist(movedPoint0,proj0)/sqrt(disps[0]*disps[0]+disps[1]*disps[1]);
// //             if(stiffnesses[pt->first.first] > 0)
//                 stiffnesses[pt->first.first] = 0 ;
//         }
        
        affectedElements[pt->first]->getState().getField(DISPLACEMENT_FIELD, *pt->first.second,disps, false, &vm, 0) ;
        Point movedPoint1 = *pt->first.second+disps;
        Point proj1(movedPoint1) ;
        baseGeometry->project(&proj1);
        if(stiffnessesX.find(pt->first.second) == stiffnessesX.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
            stiffnessesX[pt->first.second] = 0 ;
        else
            stiffnessesX[pt->first.second] = stiffCopyX[pt->first.second];
        
        if(stiffnessesY.find(pt->first.second) == stiffnessesY.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
            stiffnessesY[pt->first.second] = 0 ;
        else
            stiffnessesY[pt->first.second] = stiffCopyY[pt->first.second];
//         else  if(contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
//         {
// //             stiffnesses[pt->first.second] *= 1.+ dist(movedPoint1,proj1)/sqrt(disps[0]*disps[0]+disps[1]*disps[1]);
// //             if(stiffnesses[pt->first.second] > 0)
//                 stiffnesses[pt->first.second] = 0 ;
//         }
        
        double inSign0 = (baseGeometry->in(movedPoint0) )?-.1:.1 ;
        double inSign1 = (baseGeometry->in(movedPoint1))?-.1:.1 ;
        Point cbase = (*pt->first.first+*pt->first.second)*.5-pt->second->getCircumCenter() ;        
        Point nbase0 = (proj0 - movedPoint0)*.5 - cbase*.5;
        Point nbase1 = (proj1 - movedPoint1)*.5 - cbase*.5;
        double nbasenorm0 = nbase0.norm() ;
        double nbasenorm1 = nbase1.norm() ;
        if(nbasenorm0  < POINT_TOLERANCE)
        {
           nbase0 =  proj0 - (movedPoint0*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm0 = nbase0.norm() ;
           inSign0 = -1 ;
        }
        
        if(nbasenorm1  < POINT_TOLERANCE)
        {
           nbase1 =  proj1 - (movedPoint1*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm1 = nbase1.norm() ;
           inSign1 = -1 ;
        }
        
        if(nbasenorm0 < POINT_TOLERANCE)
            nbasenorm0 = 1. ;
        
        if(nbasenorm1 < POINT_TOLERANCE)
            nbasenorm1 = 1. ;
        
        if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0 ;  
            referenceX[pt->first.first] = proj0.getX()-normalVectors[pt->first.first].getX() ;
            referenceY[pt->first.first] = proj0.getY()-normalVectors[pt->first.first].getY() ;
            
            normalVectors[pt->first.second] = nbase1/nbasenorm1*inSign1 ;     
            referenceX[pt->first.second] = proj1.getX()-normalVectors[pt->first.second].getX() ;
            referenceY[pt->first.second] = proj1.getY()-normalVectors[pt->first.second].getY() ;
        }
        else if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) == contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;   
            referenceX[pt->first.first] = proj0.getX()-normalVectors[pt->first.first].getX() ;
            referenceY[pt->first.first] = proj0.getY()-normalVectors[pt->first.first].getY() ;
        }
        else if(contactPointsAndTributary.find(pt->first.first) == contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.second] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;    
            referenceX[pt->first.second] = proj1.getX()-normalVectors[pt->first.second].getX() ;
            referenceY[pt->first.second] = proj1.getY()-normalVectors[pt->first.second].getY() ;
        }
        
        

        
    }
    
//     print() ;
//     update() ;
    currentError = 0 ;
    
//     errors.clear() ;
//     positions.clear() ;
//     previousStiffnesses.clear() ;
//     previousErrors.clear() ;
}

bool ContactBoundaryCondition::converged() const
{
    return conv ;
}

void ContactBoundaryCondition::update()
{
    conv = true ;
    if(!active)
        return ;
    Vector disps(2) ;
    VirtualMachine vm ;
    
    
    std::set<const Point * > done ;
    currentError = 0 ;
// std::cout << std::endl;

    std::map<const Point *, double> computedErrorX ;
    std::map<const Point *, double> computedErrorY ;
    std::map<const Point *, double> smoothErrorX ;
    std::map<const Point *, double> smoothErrorY ;
    std::map<const Point *, bool> isIn ;
    double trb = 0 ;
    int cnt = 0 ;
    for(auto element: affectedElements)
    {
        for(size_t i = 0 ; i < element.second->getBoundingPoints().size() ; i++)
        {
            auto pt = contactPointsAndTributary.find( &element.second->getBoundingPoint(i)) ;
            
            if(done.find(pt->first) == done.end())
                done.insert(pt->first) ;
            else
                continue ; 
            
            if(pt == contactPointsAndTributary.end() )
                continue ;
            
            trb += pt->second ;
            cnt ++ ;
            element.second->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
            Point movedPoint = *pt->first+disps;
            Point proj(movedPoint) ;
            baseGeometry->project(&proj);

            double nerr = std::abs(proj.getX()- movedPoint.getX()) ;
            isIn[pt->first] = baseGeometry->in(movedPoint) ;
            if(isIn[pt->first])
                computedErrorX[pt->first] = -nerr*(normalVectors.find(&element.second->getBoundingPoint(i)))->second.getX() ;
            else
               computedErrorX[pt->first] = nerr*(normalVectors.find(&element.second->getBoundingPoint(i)))->second.getX() ; 
            
            nerr = std::abs(proj.getY()- movedPoint.getY()) ;
            if(isIn[pt->first])
                computedErrorY[pt->first] = -nerr*(normalVectors.find(&element.second->getBoundingPoint(i)))->second.getY() ;
            else
               computedErrorY[pt->first] = nerr*(normalVectors.find(&element.second->getBoundingPoint(i)))->second.getY() ; 
        }
    }
    if(cnt)
        trb /= cnt ;
    
    //smooth the errors
    double maxerr = 0 ;
    for(auto er : computedErrorX)
    {
        double sum = 1 ;
        double nerr = computedErrorX[er.first] ;
        smoothErrorX[er.first] =  nerr/sum ;
       maxerr = std::max(maxerr, nerr/sum) ;
    }
    for(auto er : computedErrorY)
    {
        double sum = 1 ;
        double nerr = computedErrorY[er.first] ;
        smoothErrorY[er.first] =  nerr/sum ;
       maxerr = std::max(maxerr, nerr/sum) ;
    }
    done.clear() ;
    
//     std::cout << maxerr << std::endl ;
    
    double fac = 1 ;
    double acc = 25. ;
//     if(maxerr*acc > 0.002)
//     {
//         fac *= 0.002/(maxerr*acc) ;
//     }
    
    
    for(auto element: affectedElements)
    {
        for(size_t i = 0 ; i < element.second->getBoundingPoints().size() ; i++)
        {
 
            auto pt = contactPointsAndTributary.find( &element.second->getBoundingPoint(i)) ;
            
            if(done.find(pt->first) == done.end())
                done.insert(pt->first) ;
            else
                continue ;
            
            if(pt == contactPointsAndTributary.end() )
                continue ;
            
            if(errorsX.find(pt->first) == errorsX.end() )
            {
                double nerr = smoothErrorX[pt->first] ;
                errorsX[pt->first] = nerr ;
                currentError += nerr*nerr ;
                stiffnessesY[pt->first] += acc*fac*nerr ;
                
                nerr = smoothErrorY[pt->first] ;
                errorsY[pt->first] = nerr ;
                currentError += nerr*nerr ;
                stiffnessesY[pt->first] += acc*fac*nerr ;

                continue ;
            }
            
            if(previousErrorsX.find(pt->first) == previousErrorsX.end())
            {
                double nerr = smoothErrorX[pt->first] ;
                double perr = errorsX[pt->first] ;  
                previousErrorsX[pt->first] = perr;
                errorsX[pt->first] = nerr ;
                currentError += nerr*nerr ;
                double pstiff = stiffnessesX[pt->first] ;                
                previousStiffnessesX[pt->first] = pstiff ;
                //Newton
                stiffnessesX[pt->first] += acc*fac*nerr ;
                
                nerr = smoothErrorY[pt->first] ;
                perr = errorsY[pt->first] ;  
                previousErrorsY[pt->first] = perr;
                errorsY[pt->first] = nerr ;
                currentError += nerr*nerr ;
                pstiff = stiffnessesY[pt->first] ;                
                previousStiffnessesY[pt->first] = pstiff ;
                //Newton
                stiffnessesY[pt->first] += acc*fac*nerr ;
                
//                 stiffnesses[pt->first] = stiffnesses[pt->first]*.75 + pstiff*.25 ;
                
                continue ;
            }   
              
            double nerr = smoothErrorX[pt->first] ;
            double pstiff = stiffnessesX[pt->first] ;               
            double perr = errorsX[pt->first] ;  
         
            previousStiffnessesX[pt->first] = pstiff ;
           
            previousErrorsX[pt->first] = perr ;

            stiffnessesX[pt->first] += acc*fac*nerr-(nerr-perr)*0.1*acc*fac*std::abs(nerr);
//             stiffnesses[pt->first] = stiffnesses[pt->first]*.75 + pstiff*.25 ;

//             std::cout << stiffnesses[pt->first] << std::endl ;
            errorsX[pt->first] = nerr ;
            currentError += nerr*nerr ;
            
            nerr = smoothErrorY[pt->first] ;
            pstiff = stiffnessesY[pt->first] ;               
            perr = errorsY[pt->first] ;  
         
            previousStiffnessesY[pt->first] = pstiff ;
           
            previousErrorsY[pt->first] = perr ;

            stiffnessesY[pt->first] += acc*fac*nerr-(nerr-perr)*0.1*acc*fac*std::abs(nerr);
//             stiffnesses[pt->first] = stiffnesses[pt->first]*.75 + pstiff*.25 ;

//             std::cout << stiffnesses[pt->first] << std::endl ;
            errorsY[pt->first] = nerr ;
            currentError += nerr*nerr ;

        }
    }
    
    
    for(auto pt : affectedElements)
    {
        if(stiffnessesX.find(pt.first.first) != stiffnessesX.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.first,disps, false, &vm, 0) ;
            
            double dx = std::abs(referenceX[pt.first.first]- pt.first.first->getX()-disps[0]) ;
            double force  = dx*stiffnessesX[pt.first.first]*scale/**pt.second*/ ;
            if(force > 0)
            {
                stiffnessesX[pt.first.first] *= .1 ;
            }
            double dy = std::abs(referenceY[pt.first.first]- pt.first.first->getY()-disps[1]) ;
            force  = dy*stiffnessesY[pt.first.first]*scale/**pt.second*/ ;
            if(force > 0)
            {
                stiffnessesY[pt.first.first] *= .1 ;
            }
        }
    }
    

    currentError = sqrt(currentError) ;

//     std::cout << currentError << "  " << std::flush ;
    conv = currentError < threshold ;
  
}

bool ContactBoundaryCondition::verifyConvergence() const
{
    if(!active)
        return true ;
    Vector disps(2) ;
    VirtualMachine vm ;
    
    std::set<const Point * > done ;
    double lcurrentError = 0 ;
// std::cout << std::endl;
    for(auto element: affectedElements)
    {
        for(size_t i = 0 ; i < element.second->getBoundingPoints().size() ; i++)
        {
            auto pt = contactPointsAndTributary.find( &element.second->getBoundingPoint(i)) ;
            if(pt != contactPointsAndTributary.end() && done.find(pt->first) == done.end())
            {                
                element.second->getState().getField(DISPLACEMENT_FIELD, *pt->first,disps, false, &vm, 0) ;
                Point movedPoint = *pt->first+disps;
                Point proj(movedPoint) ;
                baseGeometry->project(&proj);
                done.insert(pt->first) ;

                double nerr = std::abs(proj.getX()- movedPoint.getX()) ;
                Point p = (normalVectors.find(&element.second->getBoundingPoint(i)))->second ;
                double d = std::abs(referenceX.find(pt->first)->second- pt->first->getX()-disps[0]) ;
                double force  = d*stiffnessesX.find(pt->first)->second*scale/**pt.second*/ ;

                if(force < 0)
                    lcurrentError += nerr*nerr*p.x ;
                
                nerr = std::abs(proj.getY()- movedPoint.getY()) ;
                d = std::abs(referenceY.find(pt->first)->second- pt->first->getY()-disps[1]) ;
                force  = d*stiffnessesY.find(pt->first)->second*scale/**pt.second*/ ;

                if(force < 0)
                    lcurrentError += nerr*nerr*p.y ;

            }
        }
    }
    
    lcurrentError = sqrt(lcurrentError) ;
//     std::cout << lcurrentError << std::endl ;
//     print() ;

    bool pos = true ;

    
    return lcurrentError < threshold && pos;
  
}


void ContactBoundaryCondition::setActive(bool act)
{
    active = act ;
}

void ContactBoundaryCondition::print() const 
{
    
    Vector disps(2) ;
    VirtualMachine vm ;
    std::cout << std::endl;
    for(const auto pt : affectedElements)
    {
        if(contactPointsAndTributary.find(pt.first.first) != contactPointsAndTributary.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.first,disps, false, &vm, 0) ;
            Point movedPoint = *pt.first.first+disps;
            double dx = std::abs(referenceX.find(pt.first.first)->second- movedPoint.getX()) ;
            double forceX  = dx*stiffnessesX.find(pt.first.first)->second*scale*contactPointsAndTributary.find(pt.first.first)->second ;
            double dy = std::abs(referenceY.find(pt.first.first)->second- movedPoint.getY()) ;
            double forceY  = dy*stiffnessesY.find(pt.first.first)->second*scale*contactPointsAndTributary.find(pt.first.first)->second ;
            
            std::cout << movedPoint.getX() << "  " << movedPoint.getY() << "  "<< normalVectors.find(pt.first.first)->second.getX() << "  "<< normalVectors.find(pt.first.first)->second.getY() << "  "<< std::min(forceX, 0.) << "  "<< std::min(forceY, 0.) << "  "<< errorsX.find(pt.first.first)->second << "  "<< errorsY.find(pt.first.first)->second <<std::endl;
        }
        
        if(contactPointsAndTributary.find(pt.first.second) != contactPointsAndTributary.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.second,disps, false, &vm, 0) ;
            Point movedPoint = *pt.first.second+disps;
            double dx = std::abs(referenceX.find(pt.first.second)->second- movedPoint.getX()) ;
            double forceX  = dx*stiffnessesX.find(pt.first.second)->second*scale*contactPointsAndTributary.find(pt.first.second)->second ;
            double dy = std::abs(referenceY.find(pt.first.second)->second- movedPoint.getY()) ;
            double forceY  = dy*stiffnessesY.find(pt.first.second)->second*scale*contactPointsAndTributary.find(pt.first.second)->second ;
            
            std::cout << movedPoint.getX() << "  " << movedPoint.getY() << "  "<< normalVectors.find(pt.first.second)->second.getX() << "  "<< normalVectors.find(pt.first.second)->second.getY() << "  "<< std::min(forceX, 0.) << "  "<< std::min(forceY, 0.) << "  "<< errorsX.find(pt.first.second)->second << "  "<< errorsY.find(pt.first.second)->second <<std::endl;
        }
    }
}

double ContactBoundaryCondition::error() const
{
    return currentError ;
}

void ContactBoundaryCondition::postProcess()
{
//    if(conv)
//    {
//        double max = 0 ;
//         for(auto pt : affectedElements)
//         {
//             if(stiffnessesX.find(pt.first.first) != stiffnessesX.end())
//             {
//                 if(stiffnessesX[pt.first.first] < 0)
//                     max = std::max(max, std::abs(stiffnessesX[pt.first.first])) ;
//                 if(stiffnessesY[pt.first.first] < 0)
//                     max = std::max(max, std::abs(stiffnessesY[pt.first.first])) ;
//             }
//         }
//         for(auto pt : affectedElements)
//         {
//             if(stiffnessesX.find(pt.first.first) != stiffnessesX.end())
//             {
//                 stiffnessesX[pt.first.first] += 1e-2*max/stiffnessesX.size() ;
//                 stiffnessesY[pt.first.first] += 1e-2*max/stiffnessesX.size() ;
//             }
//         }
//    }
}

void ContactBoundaryCondition::applyBoundaryConditions( Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) 
{
    if(!active)
        return ;
//     print() ;
    std::vector<BoundaryCondition *> ret ;
    Vector disps(2) ;
    VirtualMachine vm ;
    
    //    if(conv)
//    {
       double max = 0 ;
        for(auto pt : affectedElements)
        {
            if(stiffnessesX.find(pt.first.first) != stiffnessesX.end())
            {
                if(stiffnessesX[pt.first.first] < 0)
                    max = std::max(max, std::abs(stiffnessesX[pt.first.first])) ;
                if(stiffnessesY[pt.first.first] < 0)
                    max = std::max(max, std::abs(stiffnessesY[pt.first.first])) ;
            }
        }

//    }
    
    for(auto pt : affectedElements)
    {
        if(stiffnessesX.find(pt.first.first) != stiffnessesX.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.first,disps, false, &vm, 0) ;
//             if(!baseGeometry->in(*pt.first.first+disps))
//             {
//                 continue ;
//             }
            double dx = std::abs(referenceX[pt.first.first]- pt.first.first->getX()-disps[0]) ;
            double forcex  = /*0.1**/std::pow(dx,0.1)*stiffnessesX[pt.first.first]*scale  + 1e-1*max;
            if(forcex > 0)
                forcex = 0 ;
            double dy = std::abs(referenceY[pt.first.first]- pt.first.first->getY()-disps[1]) ;
            double forcey  = /*0.1**/std::pow(dy,0.1)*stiffnessesY[pt.first.first]*scale + 1e-1*max;
            if(forcey > 0)
                forcey = 0 ;    
//             std::cout << (*pt.first.first).getY() << "  "<< force << "  "<< errors[pt.first.first]<< "  "<< currentError<< std::endl ;

            int JinvSize = 3 ;
            if ( pt.second->spaceDimensions() == SPACE_THREE_DIMENSIONAL && pt.second->timePlanes() > 1 )
                JinvSize = 4 ;
            if ( pt.second->spaceDimensions() == SPACE_TWO_DIMENSIONAL && pt.second->timePlanes() == 1 )
                JinvSize = 2 ;
            
            std::valarray<Matrix> Jinv ( (bool) pt.second->getState().JinvCache ? (*pt.second->getState().JinvCache) : Matrix( JinvSize, JinvSize ),  pt.second->getGaussPoints().gaussPoints.size()) ;

            if( ! pt.second->getState().JinvCache )
            {
                pt.second->getState().updateInverseJacobianCache(pt.second->spaceDimensions() == SPACE_THREE_DIMENSIONAL? Point(1./3.) : Point(.25, .25, .25)) ;
            }
            
            DofDefinedBoundaryCondition ycond(SET_FORCE_ETA,(ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.first->getId(), forcey*normalVectors[pt.first.first].getY()) ;
            ycond.apply ( a, mesh, &vm ) ;
            DofDefinedBoundaryCondition xcond(SET_FORCE_XI, (ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.first->getId(), forcex*normalVectors[pt.first.first].getX());    
            xcond.apply ( a, mesh, &vm ) ;
//             std::cout << force << "  "<< force*normalVectors[pt.first.first].getX() << "  " << force*normalVectors[pt.first.first].getY() << std::endl ;
        }
        
        if(stiffnessesX.find(pt.first.second) != stiffnessesX.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.second,disps, false, &vm, 0) ;
//             if(!baseGeometry->in(*pt.first.first+disps))
//                 continue ;
            double dx = std::abs(referenceX[pt.first.second]- pt.first.second->getX()-disps[0]) ;
            double forcex  = /*0.1**/std::pow(dx,0.1)*stiffnessesX[pt.first.second]*scale + 1e-1*max;
            if(forcex > 0)
                forcex = 0 ;
            double dy = std::abs(referenceY[pt.first.second]- pt.first.second->getY()-disps[1]) ;
            double forcey  = /*0.1**/std::pow(dy,0.1)*stiffnessesY[pt.first.second]*scale+ 1e-1*max;
            if(forcey > 0)
                forcey = 0 ;  

            int JinvSize = 3 ;
            if ( pt.second->spaceDimensions() == SPACE_THREE_DIMENSIONAL && pt.second->timePlanes() > 1 )
                JinvSize = 4 ;
            if ( pt.second->spaceDimensions() == SPACE_TWO_DIMENSIONAL && pt.second->timePlanes() == 1 )
                JinvSize = 2 ;
            
            std::valarray<Matrix> Jinv ( (bool) pt.second->getState().JinvCache ? (*pt.second->getState().JinvCache) : Matrix( JinvSize, JinvSize ),  pt.second->getGaussPoints().gaussPoints.size()) ;

            if( ! pt.second->getState().JinvCache )
            {
                pt.second->getState().updateInverseJacobianCache(pt.second->spaceDimensions() == SPACE_THREE_DIMENSIONAL? Point(1./3.) : Point(.25, .25, .25)) ;
            }
            
            DofDefinedBoundaryCondition ycond(SET_FORCE_ETA,(ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.second->getId(), forcey*normalVectors[pt.first.second].getY()) ;
            ycond.apply ( a, mesh, &vm ) ;
            DofDefinedBoundaryCondition xcond(SET_FORCE_XI, (ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.second->getId(), forcex*normalVectors[pt.first.second].getX());    
            xcond.apply ( a, mesh, &vm ) ;
//             std::cout << force << "  "<< force*normalVectors[pt.first.second].getX() << "  " << force*normalVectors[pt.first.second].getY() << std::endl ;
        }
    }

}

}
