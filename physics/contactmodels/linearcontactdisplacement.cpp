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
#include "linearcontactdisplacement.h"
#include "../collisiondetectors/collisiondetector.h"

namespace Amie {

LinearContactDisplacement::LinearContactDisplacement(Geometry  *geo,  double tangentStiffness) : geo(geo), tangentStiffness(tangentStiffness)
{
	getState(true).resize(1, 0.) ;
	isNull = false ;
	state = 0 ;
    change = false ;
}

std::pair< Vector, Vector > LinearContactDisplacement::computeDamageIncrement(ElementState &s)
{
    if(!es)
        es = &s ;
//     change = false ;
    int dim = s.getParent()->spaceDimensions() ;
    
    if(deltaPosition.size() != s.getParent()->getBoundingPoints().size()*dim)
        deltaPosition.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(initialPosition.size() != s.getParent()->getBoundingPoints().size()*dim)
        initialPosition.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(tangentForces.size() != s.getParent()->getBoundingPoints().size()*dim)
        tangentForces.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;   
    if(tangentDeltaForces.size() != s.getParent()->getBoundingPoints().size()*dim)
        tangentDeltaForces.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ; 
    
    
    if( s.getParent()->getBehaviour()->getCollisionDetection()->isInDamagingSet() && s.getParent()->getBehaviour()->getCollisionDetection()->isAtCheckpoint() )
    {

        VirtualMachine vm ;
        
        Vector disp(dim) ;

        int count = 0 ;
        
        double dx_ = 0;
        double dy_ = 0 ;
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            disp = {s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1]};
            if(dim == 3)
                disp = {s.getDisplacements()[i*dim], s.getDisplacements()[i*dim+1], s.getDisplacements()[i*dim+2]};
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            Point base(test) ;
            geo->project(&test);             
            if(geo->in(base))
            {  
                deltaPosition[i*dim] = test.getX()-base.getX() ;
                dx_ += deltaPosition[i*dim] ;
                deltaPosition[i*dim+1] = test.getY()-base.getY() ;
                dy_ += deltaPosition[i*dim+1] ;
                if(dim == 3)
                    deltaPosition[i*dim+2] = test.getZ()-base.getZ() ;
                    
                double dx = test.x- base.x ;
                double dy = test.y- base.y ;
                double dz = test.z- base.z ;
                
                double num = sqrt(dx*dx+dy*dy+dz*dz) ;

                count++ ;

                tangentDeltaForces[i*dim] = num*disp[0] ;
                tangentDeltaForces[i*dim+1] = num*disp[1] ;
                if(dim == 3)
                    tangentDeltaForces[i*dim+2] = num*disp[2] ;
                
            }

        }
        if(std::abs(dy_) > std::abs(dx_))
            vert = true ;
        else
            vert = false ;
    }
   
	return std::make_pair(Vector(0., 1), Vector(1., 1)) ;
}

void LinearContactDisplacement::postProcess()
{
    if(converged && state[0] > POINT_TOLERANCE)
    {

        Vector disp(es->getParent()->spaceDimensions()) ;
        for(size_t i = 0 ; i < es->getParent()->getBoundingPoints().size() ; i++)
        {
            disp = {es->getDisplacements()[i*2], es->getDisplacements()[i*2+1]};
            if(es->getParent()->spaceDimensions() == 3)
                disp = {es->getDisplacements()[i*3], es->getDisplacements()[i*3+1], es->getDisplacements()[i*3+2]};
            Point test(es->getParent()->getBoundingPoint(i) + disp) ;
            
             if(geo->in(test))
             {
                if(es->getParent()->spaceDimensions() == 2)
                {
                    tangentForces[i*2] += tangentDeltaForces[i*2]*(getState()[0]) ;
                    tangentForces[i*2+1] += tangentDeltaForces[i*2+1]*(getState()[0]) ;
                }
                else
                {
                    tangentForces[i*3] += tangentDeltaForces[i*3]*(getState()[0]) ;
                    tangentForces[i*3+1] += tangentDeltaForces[i*3+1]*(getState()[0]) ;
                    tangentForces[i*3+2] += tangentDeltaForces[i*3+2]*(getState()[0]) ;
                }
             }

        }
        
        initialPosition += deltaPosition*getState()[0] ;

        getState(true)[0] = 0 ;
    }


}

void LinearContactDisplacement::computeDelta(ElementState & s)
{
	delta = 1 ;
}

Matrix LinearContactDisplacement::apply(const Matrix & m, const Point & p , const IntegrableEntity * e , int g) const
{
    return m ;
}

std::vector<BoundaryCondition * > LinearContactDisplacement::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{

    std::vector<BoundaryCondition * > ret ;    
    if(tangentForces.size() == 0)
        return ret ;
    
    int dim = s.getParent()->spaceDimensions() ;
    Vector disp(dim) ;
    VirtualMachine vm ;

    if(!vert)
    {
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if(id != s.getParent()->getBoundingPoint(i).getId())
                continue ;
            
            
            es->getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            if(!geo->in(test))
            {
                ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            initialPosition[i*dim])) ;
                                                                        
                ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            tangentForces[dim*i+1]*tangentStiffness)
                            );
                continue ;
            }

            ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            +deltaPosition[i*dim]*getState()[0]+initialPosition[i*dim])) ;
                                                                        
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            tangentForces[dim*i+1]*tangentStiffness+tangentDeltaForces[dim*i+1]*getState()[0]*tangentStiffness)
                            );

            if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
            {

            }
        }
    }
    else
    {
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if(id != s.getParent()->getBoundingPoint(i).getId())
                continue ;
            
            
            es->getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
            Point test(s.getParent()->getBoundingPoint(i) + disp) ;
            if(!geo->in(test))
            {
                ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            initialPosition[i*dim+1])) ;
                                                            
                ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            tangentForces[dim*i]*tangentStiffness)
                            );
                continue ;
            }

            ret.push_back(new DofDefinedBoundaryCondition(SET_ALONG_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            +deltaPosition[i*dim+1]*getState()[0]+initialPosition[i*dim+1])) ;
                                                            
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), 
                                                            tangentForces[dim*i]*tangentStiffness+tangentDeltaForces[dim*i]*getState()[0]*tangentStiffness)
                            );

            
            if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
            {

            }
        }
    }


    return ret ;
}

LinearContactDisplacement::~LinearContactDisplacement()
{
}

DamageModel * LinearContactDisplacement::getCopy() const
{
    LinearContactDisplacement * ret = new LinearContactDisplacement(geo, tangentStiffness) ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}


}
