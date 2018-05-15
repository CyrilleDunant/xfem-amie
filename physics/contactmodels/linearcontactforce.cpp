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
    
    
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint())
    {
        VirtualMachine vm ;
        
        Vector disp(dim) ;

        int count = 0 ;
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            Point base(test) ;
            geo->project(&test);                    
            
            double dx = test.x- base.x ;
            double dy = test.y- base.y ;
            double dz = test.z- base.z ;
            
            double num = sqrt(dx*dx+dy*dy+dz*dz) ;

            if(geo->in(base))
            {  
                count++ ;
                deltaForce[i*dim] = dx ;
                deltaForce[i*dim+1] = dy ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = dz ;
                }
                
                tangentDeltaForce[i*dim] = -num*disp[0] ;
                tangentDeltaForce[i*dim+1] = -num*disp[1] ;
                if(dim == 3)
                    tangentDeltaForce[i*dim+2] = num*disp[2] ;
                
            }
            else
            { 

                deltaForce[i*dim] = 0 ;
                deltaForce[i*dim+1] = 0 ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = 0 ;
                }
                
                tangentDeltaForce[i*dim] =0 ;
                tangentDeltaForce[i*dim+1] = 0 ;
                if(dim == 3)
                    tangentDeltaForce[i*dim+2] = 0 ;

            }
        }
        
        if(count == 1)
        {
            deltaForce = 0 ;
            tangentDeltaForce = 0 ;
        }

//         if(std::abs(tangentDeltaForce).max() > POINT_TOLERANCE || std::abs(deltaForce).max() > POINT_TOLERANCE)
//             change = true ;

    }
   
	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    change = false ;
    if(converged && state[0] > 0)
    {
        forces += deltaForce*(getState()[0]) ;
        tangentForces += tangentDeltaForce*(getState()[0]) ;
        getState(true)[0] = 0 ;

    }


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
    
    double factor = 1.*s.getParent()->getRadius();
    if(!s.getParent()->getBehaviour()->getCollisionDetection()->met())
        return ret ;

    int dim = s.getParent()->spaceDimensions() ;
    if(getState()[0] > POINT_TOLERANCE)
//     std::cout << "conditions! " <<  getState()[0] << "  "<< (forces[0]*stiffness+
//                                                               tangentForces[0]*tangentStiffness+
//                                                               (deltaForce[0]*stiffness+tangentDeltaForce[0]*tangentStiffness)*getState()[0] 
//                                                           )*factor << std::endl ;
    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {

        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
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
        }
    }
    if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
    {
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {  
//             if( s.getParent()->getBoundingPoint(i).getId() != id)
//                 continue ;
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), (forces[dim*i]+deltaForce[dim*i]*getState()[0])*stiffness));
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(),(forces[dim*i+1]+deltaForce[dim*i+1]*getState()[0])*stiffness));
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(),(forces[dim*i+2]+deltaForce[dim*i+2]*getState()[0])*stiffness));
        }

    }
    return ret ;
}

LinearContactForce::~LinearContactForce()
{
}

DamageModel * LinearContactForce::getCopy() const
{
    LinearContactForce * ret = new LinearContactForce(geo, stiffness) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}
