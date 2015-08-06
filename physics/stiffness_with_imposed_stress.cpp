//
// C++ Implementation: stiffness_with_imposed_deformation
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_imposed_stress.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

StiffnessWithImposedStress::StiffnessWithImposedStress(const Matrix & rig, const Vector & imposedDef) : LinearForm(rig, false, false, rig.numRows()/3+1) , imposed(imposedDef)
{
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);
    this->time_d = false ;
}

StiffnessWithImposedStress::StiffnessWithImposedStress(double E, double nu, double alpha, SpaceDimensionality dim, planeType pt) : LinearForm(Tensor::cauchyGreen(std::make_pair(E,nu), true,dim,pt), false, false, dim),v(2)
{
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);

    if(dim == SPACE_TWO_DIMENSIONAL)
    {
        imposed.resize(3, 0.) ;
        for(size_t i = 0 ; i < 2 ; i++)
            imposed[i] = alpha ;
    }
    else if(dim == SPACE_THREE_DIMENSIONAL)
    {
        imposed.resize(6, 0.) ;
        for(size_t i = 0 ; i < 3 ; i++)
            imposed[i] = alpha ;
    }

    this->time_d = false ;
}

StiffnessWithImposedStress::~StiffnessWithImposedStress() { }

void StiffnessWithImposedStress::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithImposedStress::fractured() const
{
    return false ;
}

Form * StiffnessWithImposedStress::getCopy() const
{
    StiffnessWithImposedStress * copy = new StiffnessWithImposedStress(param, imposed) ;

    return copy ;
}

void StiffnessWithImposedStress::step(double timestep, ElementState & currentState, double maxscore)
{
}

Vector StiffnessWithImposedStress::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return (param * imposed) ;
}

Vector StiffnessWithImposedStress::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    return imposed ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedStress::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{

    std::vector<BoundaryCondition * > ret ;
    if(v.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[1]));
    }
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[2]));
    }
    return ret ;
}


