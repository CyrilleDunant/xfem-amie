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
    active = false ;
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
    
    scorefactor = (s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState()+1)/(refscore+1) ;
    
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint())
    {   

        Vector disp(dim) ;

        std::multimap<double, Point> mmap ;
        size_t npoints = s.getParent()->getBoundingPoints().size() ;
        active = true ;
        deltaForce = 0 ;
        tangentDeltaForce = 0 ;

        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            disp = {.5*s.getDisplacements()[i*dim]+.5*s.getDisplacements()[((i+1)%npoints)*dim], .5*s.getDisplacements()[i*dim+1]+.5*s.getDisplacements()[((i+1)%npoints)*dim+1]};

            Point test(.5*s.getParent()->getBoundingPoint(i) + .5*s.getParent()->getBoundingPoint((i+1)%npoints) + disp) ;
            Point base(test) ;
            geo->project(&test);             
            
            double dx = test.x- base.x ;
            double dy = test.y- base.y ;

            Point norm(-dy, dx, 0) ;
            double no = norm.norm() ;
            norm /= no ;

            if(geo->in(base))
            {  
                double d = dist(s.getParent()->getBoundingPoint(i) ,s.getParent()->getBoundingPoint((i+1)%npoints))*.5 ;
                
                Point A(s.getParent()->getBoundingPoint(i)+Vector({s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1]})) ;
                Point AA(A) ;
                geo->project(&A);  
                
                dx = A.x- AA.x ;
                dy = A.y- AA.y ;

                norm.set(-dy, dx, 0) ;
                no = norm.norm() ;
                norm /= no ;
                
//                     if(geo->in(AA))
//                     {
                    deltaForce[i*dim] += dx*d ;
                    deltaForce[i*dim+1] += dy*d ;

                    tangentDeltaForce[i*dim] = disp[0]*norm.getX() ;
                    tangentDeltaForce[i*dim+1] = disp[1]*norm.getY() ;
//                     }
                
                Point B(s.getParent()->getBoundingPoint(((i+1)%npoints))+Vector({s.getDisplacements()[((i+1)%npoints)*dim], s.getDisplacements()[((i+1)%npoints)*dim+1]})) ;
                Point BB(B) ;
                geo->project(&B);
                
                dx = B.x- BB.x ;
                dy = B.y- BB.y ;

                norm.set(-dy, dx, 0) ;
                no = norm.norm() ;
                norm /= no ;
                
//                     if(geo->in(BB))
//                     {
                    deltaForce[((i+1)%npoints)*dim] += dx*d ;
                    deltaForce[((i+1)%npoints)*dim+1] += dy*d ;

                    tangentDeltaForce[((i+1)%npoints)*dim] = disp[0]*norm.getX() ;
                    tangentDeltaForce[((i+1)%npoints)*dim+1] = disp[1]*norm.getY() ;
//                     }

            }

        }

        deltaForce = 0.9 *deltaForce+0.1*forces ;
        tangentDeltaForce = 0.9*tangentDeltaForce+0.1*tangentForces ;
    }
//     else if (s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint())
//     {
// //         if(active)
// //         {
//             forces *=.0 ;
//             tangentForces *=.0 ;
//             deltaForce *=0 ;
//             tangentDeltaForce *=0;
// //         }
// //         active  = false ;
//     }


	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    if(converged && es && es->getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && state[0] > 1e-6)
    {
        std::cout << "more contact " << state[0] << std::endl ;
        refscore = es->getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState() ;
        forces += deltaForce*getState()[0] ;
        tangentForces += tangentDeltaForce*getState()[0] ;

        getState(true)[0] = 0 ;


        active = false ;
        displacements = es->getDisplacements() ;
    }
    else if(converged&& es->getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet())
    {
//         forces -= deltaForce*getState()[0] ;
//         tangentForces -= tangentDeltaForce*getState()[0] ;
        active = false ;
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
    double factor = 1.;
    if(!active)
    {
        factor = std::max((s.getParent()->getBehaviour()->getCollisionDetection()->getScoreAtState()+1)/(refscore+1), .1) ;
    }
//     if(!active)
//     {
//         return ret ;
//     }
//         factor = scorefactor ;


    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        if(id != s.getParent()->getBoundingPoint(i).getId())
            continue ;
//         
        
//         Vector disp = {s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1]};
// 
//         Point test(s.getParent()->getBoundingPoint(i) + disp) ;
//         if(!geo->in(test))
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
