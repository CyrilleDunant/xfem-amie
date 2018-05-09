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
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

LinearContactForce::LinearContactForce(Geometry  *geo, double stiffness) : geo(geo), stiffness(stiffness)
{
	getState(true).resize(1, 0.) ;
	isNull = false ;
	state = 0 ;
    change = false ;
}

std::pair< Vector, Vector > LinearContactForce::computeDamageIncrement(ElementState &s)
{
    change = false ;
    int dim = s.getParent()->spaceDimensions() ;
    deltaForce.resize(s.getParent()->getBoundingPoints().size()*dim, 0.) ;
    if(forces.size() != deltaForce.size())
        forces.resize(deltaForce.size(), 0.) ;
    VirtualMachine vm ;
    for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
    {
        Vector disp(dim) ;
        s.getField(DISPLACEMENT_FIELD,s.getParent()->getBoundingPoint(i), disp,false,  &vm);
        Point test(s.getParent()->getBoundingPoint(i) + disp) ;
        Point base(test) ;
        geo->project(&test);
        double dx = test.x- base.x ;
        double dy = test.y- base.y ;
        double dz = test.z- base.z ;
        if(geo->in(base))
        {
            if(std::abs(dx) > POINT_TOLERANCE || std::abs(dy) >  POINT_TOLERANCE || std::abs(dz) > POINT_TOLERANCE)
                change = true ; 
            deltaForce[i*dim] = -dx ;
            deltaForce[i*dim+1] = -dy ;
            if(dim == 3)
                deltaForce[i*dim+2] = -dz ;
            
        }
        else
        { 
            deltaForce[i*dim] = 0 ;
            deltaForce[i*dim+1] = 0 ;
            if(dim == 3)
                deltaForce[i*dim+2] = 0 ;
        }
    }
        
	return std::make_pair(Vector(.0, 1), Vector(1., 1)) ;
}

void LinearContactForce::postProcess()
{
    if(converged && state[0] > 0)
    {
        forces += deltaForce*state[0] ;
        state[0] = 0 ;
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
//     std::cout << "conditions! " << (forces[0]+deltaForce[0]*getState()[0])*stiffness << std::endl ;
    if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        for(size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++)
        {
            if( s.getParent()->getBoundingPoint(i).getId() != id)
                continue ;
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(), (forces[dim*i]+deltaForce[dim*i]*getState()[0])*stiffness*s.getParent()->getRadius()));
            ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp, Jinv, s.getParent()->getBoundingPoint(i).getId(),(forces[dim*i+1]+deltaForce[dim*i+1]*getState()[0])*stiffness*s.getParent()->getRadius()));
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
