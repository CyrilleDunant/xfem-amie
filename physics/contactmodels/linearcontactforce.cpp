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
    change = false ;
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
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getCenter(), disp,false,  &vm);
        Point test(s.getParent()->getCenter() + disp) ;
        Point base(test) ;
        geo->project(&test);
        double mdx = test.x- base.x ;
        double mdy = test.y- base.y ;
        double mdz = test.z- base.z ;
        int count = 0 ;
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
  
            s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            Point base(test) ;
            geo->project(&test);         
            
            if(geo->in(base))
            { 
                
                double dx = test.x- base.x ;
                double dy = test.y- base.y ;
                double dz = test.z- base.z ;
                
                double num = sqrt((dx*dx+dy*dy+dz*dz)*(mdx*mdx+mdy*mdy+mdz*mdz)) ;
//                 std::cout << num << std::endl ;
                if(num <  1e-12)
                    continue ;
                count++ ;
                double normalFactor = std::abs((dx*mdx+dy*mdy+dz*mdz)/num);
                double tangentFactor = 1.-normalFactor ;
                
                deltaForce[i*dim] = -dx*normalFactor ;
                deltaForce[i*dim+1] = -dy*normalFactor ;
                if(dim == 3)
                {
                    deltaForce[i*dim+2] = -dz*normalFactor ;
                }
                
                tangentDeltaForce[i*dim] = -disp[0]*tangentFactor ;
                tangentDeltaForce[i*dim+1] = -disp[1]*tangentFactor ;
                if(dim == 3)
                    tangentDeltaForce[i*dim+2] = -disp[2]*tangentFactor ;
                
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
            tangentDeltaForce = 0 ;
            deltaForce = 0 ;
        }

        if(std::abs(tangentDeltaForce).max() > POINT_TOLERANCE || std::abs(deltaForce).max() > POINT_TOLERANCE)
            change = true ;
        
//         std::cout << std::abs(deltaForce).max() << "  "<< std::abs(tangentDeltaForce).max() << "  " << count << std::endl ;
    }
   
	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    
    if(converged /*&& !es->getParent()->getBehaviour()->getCollisionDetection()->met(-1e-6)*/)
    {
//         if(es->getParent()->getBehaviour()->getCollisionDetection()->met())
//         {
            forces += deltaForce*(getState()[0]) ;
            tangentForces += tangentDeltaForce*(getState()[0]) ;
//         }
        getState(true)[0] = 0 ;
    }
    
//     if(converged && !es->getParent()->getBehaviour()->getCollisionDetection()->met(1e-6))
//     {
//         forces *= .5 ;
//         tangentForces *= .5 ;
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
    if(!s.getParent()->getBehaviour()->getCollisionDetection()->met())
        return ret ;

    int dim = s.getParent()->spaceDimensions() ;
//     std::cout << "conditions! " <<  getState()[0] << "  "<< (forces[0]+deltaForce[0]*getState()[0])*stiffness << std::endl ;
    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {

        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if( s.getParent()->getBoundingPoint(i).getId() != id)
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
                                                              
                                                          )));
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
