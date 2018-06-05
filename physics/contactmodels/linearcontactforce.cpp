//
// C++ Implementation: linearcontact force
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "linearcontactforce.h"
#include "../collisiondetectors/collisiondetector.h"

namespace Amie {

LinearContactForce::LinearContactForce(Geometry  *geo, double stiffness, double tangentStiffness) : geo(geo), stiffness(stiffness), tangentStiffness(tangentStiffness)
{
	getState(true).resize(1, 0.) ;
	isNull = false ;
	state = 0 ;
    change = false ;
    setConvergenceType(CONSERVATIVE);
}

std::pair< Vector, Vector > LinearContactForce::computeDamageIncrement(ElementState &s)
{
    if(!es)
        es = &s ;
//     change = false ;
    int dim = s.getParent()->spaceDimensions() ;
    
    if(deltaForce.size() != s.getParent()->getBoundingPoints().size()*dim)
        deltaForce.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(forces.size() != deltaForce.size())
        forces.resize(deltaForce.size(), 0.) ;
    if(tangentDeltaForce.size() != s.getParent()->getBoundingPoints().size()*dim)
        tangentDeltaForce.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(tangentForces.size() != tangentDeltaForce.size())
        tangentForces.resize(deltaForce.size(), 0.) ;   
    
    
    if(s.getDeltaTime() > POINT_TOLERANCE)
    {
        tangentForces = 0 ;
    }
    
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint() )
    {

//         forces *= .95 ;
//         tangentForces *= .95 ;
//         getState(true)[0] = .05 ;
//         deltaForce = forces ;
//         tangentDeltaForce = tangentForces ;
//         forces += deltaForce*(getState()[0]) ;
//         tangentForces += tangentDeltaForce*(getState()[0]) ;
//         getState(true)[0] = 0 ;
        
        VirtualMachine vm ;
        
        Vector disp(dim) ;

        int count = 0 ;
        
//         double mindist = -1 ;
//         for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
//         {
//             s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
//             Point test(s.getParent()->getBoundingPoint(i) + disp) ;
//             Point base(test) ;
//             geo->project(&test);                    
//             
//             double dx = test.x- base.x ;
//             double dy = test.y- base.y ;
//             double dz = test.z- base.z ;
//             
//             double num = sqrt(dx*dx+dy*dy+dz*dz) ;
//             
//             if(geo->in(base))
//             {
//                 if (num > mindist)
//                     mindist = num ;
//             }
//         }
//         std::cout << mindist << "  " << std::flush ;
        
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            disp = {s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1]};
            if(dim == 3)
                disp = {s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1], s.getDisplacements()[i*dim+2]};
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            Point base(test) ;
            geo->project(&test);                    
            
            double dx = test.x- base.x ;
            double dy = test.y- base.y ;
            double dz = test.z- base.z ;
            
            double num = sqrt(dx*dx+dy*dy+dz*dz) ;
//             dx /= num ;
//             dy /= num ;
//             dz /= num ;

            if(geo->in(base))
            {  
//                 std::cout << num << "  " << dx<< "  "<<std::flush ;
                count++ ;
                deltaForce[i*dim] = dx ;
                if(std::abs(dx) < 1e-4 && std::abs(forces[i*dim]) >=  1e-4)
                    deltaForce[i*dim] = forces[i*dim] ;
                deltaForce[i*dim+1] = dy ;
                if(std::abs(dy) < 1e-4 && std::abs(forces[i*dim+1]) >=  1e-4)
                    deltaForce[i*dim+1] = forces[i*dim+1] ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = dz ;
                    if(std::abs(dz) < 1e-4 && std::abs(forces[i*dim+2]) >=  1e-4)
                        deltaForce[i*dim+2] = forces[i*dim+2] ;
                }

                tangentDeltaForce[i*dim] = num*disp[0] ;
                tangentDeltaForce[i*dim+1] = num*disp[1] ;
                if(dim == 3)
                    tangentDeltaForce[i*dim+2] = num*disp[2] ;
                
            }
            else //if(std::abs(forces[i*dim])+std::abs(forces[i*dim+1])+((dim == 3)?std::abs(forces[i*dim+2]):0.) > POINT_TOLERANCE)
            { 
                deltaForce[i*dim] = -forces[i*dim] ;
                deltaForce[i*dim+1] = -forces[i*dim+1] ;
                if(dim == 3)
                    deltaForce[i*dim+2] = -forces[i*dim+2] ;
                
                tangentDeltaForce[i*dim] = 0 ;
                tangentDeltaForce[i*dim+1] = 0 ;
                if(dim == 3)
                    tangentDeltaForce[i*dim+2] = 0 ;
                
            }
        }
        
//         std::cout << std::endl ;
        
//         if(count == 1)
//         {
//             deltaForce = 0 ;
//             tangentDeltaForce = 0 ;
//         }

//         if(std::abs(tangentDeltaForce).max() > POINT_TOLERANCE || std::abs(deltaForce).max() > POINT_TOLERANCE)
//             change = true ;

    }
   
	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    if(converged && state[0] > POINT_TOLERANCE)
    {

//         Vector disp(es->getParent()->spaceDimensions()) ;
        for(size_t i = 0 ; i < es->getParent()->getBoundingPoints().size() ; i++)
        {
//             disp = {es->getDisplacements()[i*2], es->getDisplacements()[i*2+1]};
//             if(es->getParent()->spaceDimensions() == 3)
//                 disp = {es->getDisplacements()[i*3], es->getDisplacements()[i*3+1], es->getDisplacements()[i*3+2]};
//             Point test(es->getParent()->getBoundingPoint(i) + disp) ;
//             
//              if(geo->in(test))
//              {
                if(es->getParent()->spaceDimensions() == 2)
                {
                    forces[i*2] += deltaForce[i*2]*(getState()[0]) ;
                    tangentForces[i*2] += tangentDeltaForce[i*2]*(getState()[0]) ;
                    forces[i*2+1] += deltaForce[i*2+1]*(getState()[0]) ;
                    tangentForces[i*2+1] += tangentDeltaForce[i*2+1]*(getState()[0]) ;
                }
                else
                {
                    forces[i*3] += deltaForce[i*3]*(getState()[0]) ;
                    tangentForces[i*3] += tangentDeltaForce[i*3]*(getState()[0]) ;
                    forces[i*3+1] += deltaForce[i*3+1]*(getState()[0]) ;
                    tangentForces[i*3+1] += tangentDeltaForce[i*3+1]*(getState()[0]) ;
                    forces[i*3+2] += deltaForce[i*3+2]*(getState()[0]) ;
                    tangentForces[i*3+2] += tangentDeltaForce[i*3+2]*(getState()[0]) ;
                }
//              }
//              else 
//              {
//                 if(es->getParent()->spaceDimensions() == 2)
//                 {
//                     forces[i*2] += deltaForce[i*2]*(getState()[0]) ;
//                     tangentForces[i*2] += tangentDeltaForce[i*2]*(getState()[0]) ;
//                     forces[i*2+1] += deltaForce[i*2+1]*(getState()[0]) ;
//                     tangentForces[i*2+1] += tangentDeltaForce[i*2+1]*(getState()[0]) ;
//                     
// //                     forces[i*2] *=.9995 ;
// //                     tangentForces[i*2] *=.9995 ;
// //                     forces[i*2+1] *=.9995 ;
// //                     tangentForces[i*2+1] *=.9995 ;
//                 }
//                 else
//                 {
//                     forces[i*3] += deltaForce[i*3]*(getState()[0]) ;
//                     tangentForces[i*3] += tangentDeltaForce[i*3]*(getState()[0]) ;
//                     forces[i*3+1] += deltaForce[i*3+1]*(getState()[0]) ;
//                     tangentForces[i*3+1] += tangentDeltaForce[i*3+1]*(getState()[0]) ;
//                     forces[i*3+2] += deltaForce[i*3+2]*(getState()[0]) ;
//                     tangentForces[i*3+2] += tangentDeltaForce[i*3+2]*(getState()[0]) ;
//                     
// //                     forces[i*3] *=.9995 ;
// //                     tangentForces[i*3] *=.9995 ;
// //                     forces[i*3+1] *=.9995 ;
// //                     tangentForces[i*3+1] *=.9995 ;
// //                     forces[i*3+2] *=.9995 ;
// //                     tangentForces[i*3+2] *=.9995 ;
//                 }
//              }
        }
        

        getState(true)[0] = 0 ;
    }
//     else if(converged && es->getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet())
//     {
//         Vector disp(es->getParent()->spaceDimensions()) ;
//         for(size_t i = 0 ; i < es->getParent()->getBoundingPoints().size() ; i++)
//         {
//             disp = {es->getDisplacements()[i*2], es->getDisplacements()[i*2+1]};
//             if(es->getParent()->spaceDimensions() == 3)
//                 disp = {es->getDisplacements()[i*3], es->getDisplacements()[i*3+1], es->getDisplacements()[i*3+2]};
//             Point test(es->getParent()->getBoundingPoint(i) + disp) ;
//             
//              if(!geo->in(test))
//              {
//                 if(es->getParent()->spaceDimensions() == 2)
//                 {
//                     forces[i*2] *=.995 ;
//                     tangentForces[i*2] *=.995 ;
//                     forces[i*2+1] *=.995 ;
//                     tangentForces[i*2+1] *=.995 ;
//                 }
//                 else
//                 {
//                     forces[i*3] *=.995 ;
//                     tangentForces[i*3] *=.995 ;
//                     forces[i*3+1] *=.995 ;
//                     tangentForces[i*3+1] *=.995 ;
//                     forces[i*3+2] *=.995 ;
//                     tangentForces[i*3+2] *=.995 ;
//                 }
//              }
//         }
//     }


}

void LinearContactForce::computeDelta(ElementState & s)
{
	delta = 1 ;
}

Matrix LinearContactForce::apply(const Matrix & m, const Point & p , const IntegrableEntity * e , int g) const
{
    return m ;
}

std::vector<BoundaryCondition * > LinearContactForce::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{

    std::vector<BoundaryCondition * > ret ;    
    if(forces.size() == 0)
        return ret ;
    
    double factor = s.getParent()->getRadius();

    int dim = s.getParent()->spaceDimensions() ;

    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        if(id != s.getParent()->getBoundingPoint(i).getId())
            continue ;
        
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                        (   forces[dim*i]*stiffness+
                                                            tangentForces[dim*i]*tangentStiffness+
                                                            (deltaForce[dim*i]*stiffness+tangentDeltaForce[dim*i]*tangentStiffness)*getState()[0] 
                                                        )*factor
                                                        )
                        );
        
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                        (   forces[dim*i+1]*stiffness+
                                                            tangentForces[dim*i+1]*tangentStiffness+
                                                            (deltaForce[dim*i+1]*stiffness+tangentDeltaForce[dim*i+1]*tangentStiffness)*getState()[0] 
                                                        )*factor
                                                        )
                        );
        
        if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
        {
                        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                        (   forces[dim*i+2]*stiffness+
                                                            tangentForces[dim*i+2]*tangentStiffness+
                                                            (deltaForce[dim*i+2]*stiffness+tangentDeltaForce[dim*i+2]*tangentStiffness)*getState()[0] 
                                                        )*factor
                                                        )
                        );
        }
    }


    return ret ;
}

LinearContactForce::~LinearContactForce()
{
}

DamageModel * LinearContactForce::getCopy() const
{
    LinearContactForce * ret = new LinearContactForce(geo, stiffness, tangentStiffness) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}
