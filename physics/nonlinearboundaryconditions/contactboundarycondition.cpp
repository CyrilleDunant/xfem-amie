
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
        double inSign0 = (baseGeometry->in(*pt->first.first) )?-1.:1 ;
        double inSign1 = (baseGeometry->in(*pt->first.second))?-1.:1 ;
        Point nbase0 = proj0 - *pt->first.first ;
        Point nbase1 = proj1 - *pt->first.second ;
        double nbasenorm0 = nbase0.norm() ;
        double nbasenorm1 = nbase1.norm() ;
        if(nbasenorm0  < POINT_TOLERANCE)
        {
           nbase0 =  proj0 - (*pt->first.first*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm0 = nbase0.norm() ;
           inSign0 = -1 ;
        }
        
        if(nbasenorm1  < POINT_TOLERANCE)
        {
           nbase1 =  proj1 - (*pt->first.second*.95-baseGeometry->getCenter()*.05) ;
           nbasenorm1 = nbase1.norm() ;
           inSign1 = -1 ;
        }
        
        if(inSign0 > 0 && inSign1 > 0)
            continue ;
        
        
        if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0 ;        
            referencePoints[pt->first.first] = proj0-normalVectors[pt->first.first] ;
            stiffnesses[pt->first.first] = 0 ;
            distances[pt->first.first] = dist(*pt->first.first, referencePoints[pt->first.first]) ;
            
            normalVectors[pt->first.second] = nbase1/nbasenorm1*inSign1 ;        
            referencePoints[pt->first.second] = proj1-normalVectors[pt->first.second] ;
            stiffnesses[pt->first.second] = 0 ;
            distances[pt->first.second] = dist(*pt->first.second, referencePoints[pt->first.second]) ;
        }
        
        if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) == contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;        
            referencePoints[pt->first.first] = proj0-normalVectors[pt->first.first] ;
            stiffnesses[pt->first.first] = 0 ;
            distances[pt->first.first] = dist(*pt->first.first, referencePoints[pt->first.first]) ;
        }
        
        if(contactPointsAndTributary.find(pt->first.first) == contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.second] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;        
            referencePoints[pt->first.second] = proj1-normalVectors[pt->first.second] ;
            stiffnesses[pt->first.second] = 0 ;
            distances[pt->first.second] = dist(*pt->first.second, referencePoints[pt->first.second]) ;
        }
        
    }
//     std::cout << "padoum" << std::endl ;
//     for(auto pt = contactPointsAndTributary.begin() ; pt != contactPointsAndTributary.end() ; pt++)
//     {
//         std::cout << (pt->first)->getX() << "  " << (pt->first)->getY() << "  "<< normalVectors[pt->first].getX() << "  "<< normalVectors[pt->first].getY() << "  "<<stiffnesses[pt->first]*pt->second << std::endl;
//     }
    
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
        if(stiffnesses.find(pt->first.first) == stiffnesses.end() && contactPointsAndTributary.find(pt->first.first) !=  contactPointsAndTributary.end())
            stiffnesses[pt->first.first] = 0 ;
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
        if(stiffnesses.find(pt->first.second) == stiffnesses.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
            stiffnesses[pt->first.second] = 0 ;
//         else  if(contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
//         {
// //             stiffnesses[pt->first.second] *= 1.+ dist(movedPoint1,proj1)/sqrt(disps[0]*disps[0]+disps[1]*disps[1]);
// //             if(stiffnesses[pt->first.second] > 0)
//                 stiffnesses[pt->first.second] = 0 ;
//         }
        
        double inSign0 = (baseGeometry->in(*pt->first.first) )?-1.:1 ;
        double inSign1 = (baseGeometry->in(*pt->first.second))?-1.:1 ;
        Point nbase0 = proj0 - movedPoint0 ;
        Point nbase1 = proj1 - movedPoint1 ;
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
        
        
        if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0 ;        
            referencePoints[pt->first.first] = proj0-normalVectors[pt->first.first] ;
            stiffnesses[pt->first.first] = 0 ;
            distances[pt->first.first] = dist(*pt->first.first, referencePoints[pt->first.first]) ;
            
            normalVectors[pt->first.second] = nbase1/nbasenorm1*inSign1 ;        
            referencePoints[pt->first.second] = proj1-normalVectors[pt->first.second] ;
            stiffnesses[pt->first.second] = 0 ;
            distances[pt->first.second] = dist(*pt->first.second, referencePoints[pt->first.second]) ;
        }
        else if(contactPointsAndTributary.find(pt->first.first) != contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) == contactPointsAndTributary.end())
        {
            normalVectors[pt->first.first] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;        
            referencePoints[pt->first.first] = proj0-normalVectors[pt->first.first] ;
            stiffnesses[pt->first.first] = 0 ;
            distances[pt->first.first] = dist(*pt->first.first, referencePoints[pt->first.first]) ;
        }
        else if(contactPointsAndTributary.find(pt->first.first) == contactPointsAndTributary.end() && contactPointsAndTributary.find(pt->first.second) != contactPointsAndTributary.end())
        {
            normalVectors[pt->first.second] = nbase0/nbasenorm0*inSign0*.5 + nbase1/nbasenorm1*inSign1*.5;        
            referencePoints[pt->first.second] = proj1-normalVectors[pt->first.second] ;
            stiffnesses[pt->first.second] = 0 ;
            distances[pt->first.second] = dist(*pt->first.second, referencePoints[pt->first.second]) ;
        }
        
        

        
    }
    
//     print() ;
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
    conv = true ;
    if(!active)
        return ;
    Vector disps(2) ;
    VirtualMachine vm ;
    
    
    std::set<const Point * > done ;
    currentError = 0 ;
// std::cout << std::endl;

    std::map<const Point *, double> computedError ;
    std::map<const Point *, double> smoothError ;
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

            double nerr = dist(proj, movedPoint) ;

            computedError[pt->first] = nerr*(normalVectors.find(&element.second->getBoundingPoint(i)))->second.norm() ;
            isIn[pt->first] = baseGeometry->in(movedPoint) ;
        }
    }
    trb /= cnt ;
    
    //smooth the errors
    double maxerr = 0 ;
    for(auto er : computedError)
    {
        double sum = 1 ;
        double nerr = computedError[er.first] ;
//         for( auto alter : computedError)
//         {
//             
//             double x = dist(er.first, alter.first) ;
//             double fac = exp(-x*x/(3.*trb*trb)) ;
//             sum += fac ;
//             nerr += fac*alter.second ;
//         }
        
        smoothError[er.first] =  nerr/sum ;
       maxerr = std::max(maxerr, nerr/sum) ;
    }
    done.clear() ;
    
//     std::cout << maxerr << std::endl ;
    
    double fac = 1 ;
    if(maxerr*5 > 0.004)
    {
        fac *= 0.004/(maxerr*5) ;
    }
    
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
            
            if(errors.find(pt->first) == errors.end() )
            {
                double nerr = smoothError[pt->first] ;
                errors[pt->first] = nerr ;
                currentError += nerr*nerr ;

                //Newton
                if(isIn[pt->first])
                    stiffnesses[pt->first] -= 5.*fac*nerr ;
                else
                    stiffnesses[pt->first] += 10.*fac*nerr ; 
                continue ;
            }
            
            if(previousErrors.find(pt->first) == previousErrors.end())
            {
                double nerr = smoothError[pt->first] ;
                previousErrors[pt->first] = errors[pt->first];
                errors[pt->first] = nerr ;
                currentError += nerr*nerr ;
                                
                previousStiffnesses[pt->first] = stiffnesses[pt->first] ;
                //Newton
                if(isIn[pt->first])
                    stiffnesses[pt->first] -= 5.*fac*nerr ;
                else
                    stiffnesses[pt->first] += 10.*fac*nerr ; 
                
                continue ;
            }   
              
            double nerr = smoothError[pt->first] ;
            double pstiff = stiffnesses[pt->first] ;               
            double perr = errors[pt->first] ;  
         
            previousStiffnesses[pt->first] = pstiff ;
           
            previousErrors[pt->first] = perr ;

            if(isIn[pt->first])
                stiffnesses[pt->first] -= 5.*fac*nerr ;
            else
                stiffnesses[pt->first] += 10.*fac*nerr ;  

            errors[pt->first] = nerr ;
            currentError += nerr*nerr ;

        }
    }
    
    for (auto s : stiffnesses)
    {
        if(s.second > 0)
            s.second = 0 ;
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

                double nerr = dist(proj, movedPoint) ;
                Point p = (normalVectors.find(&element.second->getBoundingPoint(i)))->second ;
                lcurrentError += nerr*nerr*sqrt(p.y*p.y+p.x*p.x) ;
            }
        }
    }
    
    lcurrentError = sqrt(lcurrentError) ;
//     std::cout << lcurrentError << std::endl ;
//     print() ;

    return lcurrentError < threshold ;
  
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
            double d = dist(referencePoints.find(pt.first.first)->second, movedPoint) ;
            double force  = d*stiffnesses.find(pt.first.first)->second*scale*contactPointsAndTributary.find(pt.first.first)->second ;
//             if(!baseGeometry->in(movedPoint))
//                 force *= .01 ;
            
            std::cout << movedPoint.getX() << "  " << movedPoint.getY() << "  "<< normalVectors.find(pt.first.first)->second.getX() << "  "<< normalVectors.find(pt.first.first)->second.getY() << "  "<< force << "  "<< errors.find(pt.first.first)->second << std::endl;
        }
        
        if(contactPointsAndTributary.find(pt.first.second) != contactPointsAndTributary.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.second,disps, false, &vm, 0) ;
            Point movedPoint = *pt.first.second+disps;
            double d = dist(referencePoints.find(pt.first.second)->second, movedPoint) ;
            double force  = d*stiffnesses.find(pt.first.second)->second*scale*contactPointsAndTributary.find(pt.first.second)->second ;
            
//             if(!baseGeometry->in(movedPoint))
//                 force *= .01 ;
            
            std::cout << movedPoint.getX() << "  " << movedPoint.getY() << "  "<< normalVectors.find(pt.first.second)->second.getX() << "  "<< normalVectors.find(pt.first.second)->second.getY() << "  "<< force << "  "<< errors.find(pt.first.second)->second << std::endl;
        }
    }
}

double ContactBoundaryCondition::error() const
{
    return currentError ;
}

void ContactBoundaryCondition::applyBoundaryConditions( Assembly * a, Mesh<DelaunayTriangle,DelaunayTreeItem> * mesh) 
{
    if(!active)
        return ;
//     print() ;
    std::vector<BoundaryCondition *> ret ;
    Vector disps(2) ;
    VirtualMachine vm ;
    for(auto pt : affectedElements)
    {
        if(stiffnesses.find(pt.first.first) != stiffnesses.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.first,disps, false, &vm, 0) ;
            double d = dist(referencePoints[pt.first.first], *pt.first.first+disps) ;
            double force  = d*stiffnesses[pt.first.first]*scale/**pt.second*/ ;
//             if(!baseGeometry->in(*pt.first.first+disps))
//                 force *= .01 ;
            
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
            
            DofDefinedBoundaryCondition ycond(SET_FORCE_ETA,(ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.first->getId(), force*normalVectors[pt.first.first].getY()) ;
            ycond.apply ( a, mesh ) ;
            DofDefinedBoundaryCondition xcond(SET_FORCE_XI, (ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.first->getId(), force*normalVectors[pt.first.first].getX());    
            xcond.apply ( a, mesh ) ;
//             std::cout << force << "  "<< force*normalVectors[pt.first.first].getX() << "  " << force*normalVectors[pt.first.first].getY() << std::endl ;
        }
        
        if(stiffnesses.find(pt.first.second) != stiffnesses.end())
        {
            pt.second->getState().getField(DISPLACEMENT_FIELD, *pt.first.second,disps, false, &vm, 0) ;
            double d = dist(referencePoints[pt.first.second], *pt.first.second+disps) ;
            double force  = d*stiffnesses[pt.first.second]*scale/**pt.second*/ ;
//             if(!baseGeometry->in(*pt.first.second+disps))
//                 force *= .01 ;
            
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
            
            DofDefinedBoundaryCondition ycond(SET_FORCE_ETA,(ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.second->getId(), force*normalVectors[pt.first.second].getY()) ;
            ycond.apply ( a, mesh ) ;
            DofDefinedBoundaryCondition xcond(SET_FORCE_XI, (ElementarySurface*)pt.second,pt.second->getGaussPoints(), Jinv, pt.first.second->getId(), force*normalVectors[pt.first.second].getX());    
            xcond.apply ( a, mesh ) ;
//             std::cout << force << "  "<< force*normalVectors[pt.first.second].getX() << "  " << force*normalVectors[pt.first.second].getY() << std::endl ;
        }
    }

}

}
