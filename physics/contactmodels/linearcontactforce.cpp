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
    if(displacements.size() != s.getParent()->getBoundingPoints().size()*dim)
        displacements.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(forces.size() != deltaForce.size())
        forces.resize(deltaForce.size(), 0.) ;
    if(tangentDeltaForce.size() != s.getParent()->getBoundingPoints().size()*dim)
        tangentDeltaForce.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(tangentForces.size() != tangentDeltaForce.size())
        tangentForces.resize(deltaForce.size(), 0.) ;   
    
    if(s.getDeltaTime() > POINT_TOLERANCE)
        displacements = 0 ;
    
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint() )
    {   
        if((s.getParent()->getBehaviour()->getFractureCriterion() && s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState() > s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) || !s.getParent()->getBehaviour()->getFractureCriterion())
        {
            Vector disp(dim) ;

            std::multimap<double, Point> mmap ;
            size_t npoints = s.getParent()->getBoundingPoints().size() ;
            active = false ;
            deltaForce = 0 ;
            tangentDeltaForce = 0 ;
            refscore = s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
            forces *= scorefactor ;
            scorefactor = 1. ;

            for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
            {

                Point test(.5*s.getParent()->getBoundingPoint(i) + .5*s.getParent()->getBoundingPoint((i+1)%npoints) + 
                Vector({.5*s.getDisplacements()[i*dim]+.5*s.getDisplacements()[((i+1)%npoints)*dim], .5*s.getDisplacements()[i*dim+1]+.5*s.getDisplacements()[((i+1)%npoints)*dim+1]})) ;
                Point base(test) ;
                geo->project(&test);             
                
                double dx = test.x- base.x ;
                double dy = test.y- base.y ;

                Point norm(-dy, dx, 0) ;
                double no = norm.norm() ;
                mmap.insert(std::make_pair(no, base)) ;
                norm /= no ;
                
                
                disp = {.5*s.getDisplacements()[i*dim]+.5*s.getDisplacements()[((i+1)%npoints)*dim], .5*s.getDisplacements()[i*dim+1]+.5*s.getDisplacements()[((i+1)%npoints)*dim+1]};
            
            
                if(geo->in(base))
                {  
                    double d = dist(s.getParent()->getBoundingPoint(i) ,s.getParent()->getBoundingPoint((i+1)%npoints))*.5 ;

                    deltaForce[i*dim] += dx*d ;
                    deltaForce[i*dim+1] += dy*d ;

                    tangentDeltaForce[i*dim] += -disp[0]*norm.getX()*.5 ;
                    tangentDeltaForce[i*dim+1] += -disp[1]*norm.getY()*.5 ;
                    
                    deltaForce[((i+1)%npoints)*dim] += dx*d ;
                    deltaForce[((i+1)%npoints)*dim+1] += dy*d ;

                    tangentDeltaForce[((i+1)%npoints)*dim] += -disp[0]*norm.getX()*.5 ;
                    tangentDeltaForce[((i+1)%npoints)*dim+1] += -disp[1]*norm.getY()*.5 ;


                }

            }
            active = true ;
        }
        else
        {
          scorefactor = (s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState()+1)/(refscore+1) ;  
          active = false ;  
        }

    }
    else if( s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint())
    {
        scorefactor = (s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState()+1)/(refscore+1) ;  
        active = false ;
    }

	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    if(converged && state[0] > POINT_TOLERANCE)
    {
        refscore = es->getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
        
        for(size_t i = 0 ; i < es->getParent()->getBoundingPoints().size() ; i++)
        {
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
        }

        getState(true)[0] = 0 ;
        displacements = es->getDisplacements() ;
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
    int dim = s.getParent()->spaceDimensions() ;
    if( active || true )
    {
        double factor = 1.;
        if(!active)
            factor *= scorefactor ;


        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if(id != s.getParent()->getBoundingPoint(i).getId())
                continue ;
            
    //         Vector disp(dim) ;
    //         disp = {s.getDisplacements()[i*dim]-displacements[i*dim], s.getDisplacements()[i*dim+1]-displacements[i*dim+1]};
    //         if(dim == 3)
    //             disp = {s.getDisplacements()[i*dim]-displacements[i*dim], s.getDisplacements()[i*dim+1]-displacements[i*dim+1], s.getDisplacements()[i*dim+2]-displacements[i*dim+2]};
            
    //         if(!geo->in(s.getParent()->getBoundingPoint(i)+disp))
    //             continue ;
                

            
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            (   forces[dim*i]*stiffness-
                                                                tangentForces[dim*i]*stiffness*tangentStiffness*forces[dim*i+1]+
                                                                (deltaForce[dim*i]*stiffness+deltaForce[dim*i]*tangentDeltaForce[dim*i]*tangentStiffness*stiffness)*getState()[0] 
                                                            )*factor
                                                            )
                            );
            
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            (   forces[dim*i+1]*stiffness+
                                                                forces[dim*i]*stiffness*tangentForces[dim*i+1]*tangentStiffness+
                                                                (deltaForce[dim*i+1]*stiffness+deltaForce[dim*i]*tangentDeltaForce[dim*i+1]*tangentStiffness*stiffness)*getState()[0] 
                                                            )*factor
                                                            )
                            );
            
//             if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
//             {
//                             ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
//                                                             (   forces[dim*i+2]*stiffness+
//                                                                 tangentForces[dim*i+2]*stiffness*tangentStiffness+
//                                                                 (deltaForce[dim*i+2]*stiffness+tangentDeltaForce[dim*i+2]*stiffness*tangentStiffness)*getState()[0] 
//                                                             )*factor
//                                                             )
//                             );
//             }
        }
    }
    else if(std::abs(forces).max() > POINT_TOLERANCE)
    {
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if(id != s.getParent()->getBoundingPoint(i).getId())
                continue ;
            
    //         Vector disp(dim) ;
    //         disp = {s.getDisplacements()[i*dim]-displacements[i*dim], s.getDisplacements()[i*dim+1]-displacements[i*dim+1]};
    //         if(dim == 3)
    //             disp = {s.getDisplacements()[i*dim]-displacements[i*dim], s.getDisplacements()[i*dim+1]-displacements[i*dim+1], s.getDisplacements()[i*dim+2]-displacements[i*dim+2]};
            
    //         if(!geo->in(s.getParent()->getBoundingPoint(i)+disp))
    //             continue ;
                
            double factor =0 ;

            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            (   forces[dim*i]*stiffness-
                                                                tangentForces[dim*i]*tangentStiffness*forces[dim*i+1]+
                                                                (deltaForce[dim*i]*stiffness+deltaForce[dim*i]*tangentDeltaForce[dim*i]*stiffness*tangentStiffness)*getState()[0] 
                                                            )*factor
                                                            )
                            );
            
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            (   forces[dim*i+1]*stiffness+
                                                                forces[dim*i]*tangentForces[dim*i+1]*tangentStiffness+
                                                                (deltaForce[dim*i+1]*stiffness+deltaForce[dim*i]*tangentDeltaForce[dim*i+1]*stiffness*tangentStiffness)*getState()[0] 
                                                            )*factor
                                                            )
                            );
            
            if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
            {
                            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            (   forces[dim*i+2]*stiffness+
                                                                tangentForces[dim*i+2]*stiffness*tangentStiffness+
                                                                (deltaForce[dim*i+2]*stiffness+tangentDeltaForce[dim*i+2]*stiffness*tangentStiffness)*getState()[0] 
                                                            )*factor
                                                            )
                            );
            }
            
            
//             ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), displacements[i*dim]
//                                                             )
//                             );
//             
//             ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), displacements[i*dim+1]
//                                                             )
//                             );
//             
//             if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
//             {
//                             ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_ZETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(),displacements[i*dim+2]
//                                                             )
//                             );
//             }
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
