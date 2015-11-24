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

#include "stiffness_with_variable_imposed_deformation_and_fracture.h"
#include <limits>
#include "../features/boundarycondition.h"

using namespace Amie ;

StiffnessWithVariableImposedDeformationAndFracture::StiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * c) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef),criterion(c)
{
    dfunc = new IsotropicLinearDamage() ;
    frac = false ;
    init = param[0][0] ;
    change  = false ;
    previousDamage = 0 ;
    damage = 0 ;

    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36)
        v.push_back(ZETA);
}

StiffnessWithVariableImposedDeformationAndFracture::~StiffnessWithVariableImposedDeformationAndFracture()
{
    delete dfunc ;
}

void StiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

void StiffnessWithVariableImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    currentState.getParent()->behaviourUpdated = false ;

    dfunc->step(currentState, maxscore) ;
    change = dfunc->changed() ;
    currentState.getParent()->behaviourUpdated = true ;
    frac = dfunc->fractured() ;

    previousDamage = damage ;

    damage = dfunc->getState()[0] ;

    if(!currentState.getParent()->behaviourUpdated && timestep > std::numeric_limits<double>::epsilon())
    {
        double randomVar = 1 ; //(double)rand()/(double)RAND_MAX ;
        imposed[0] = timestep*randomVar ;
        imposed[1] = timestep*randomVar ;
        currentState.getParent()->behaviourUpdated = true ;
    }

    if(frac)
    {
        imposed = 0 ;
    }

    change = currentState.getParent()->behaviourUpdated  ;
    currentState.getParent()->needAssembly = change  ;
}

bool StiffnessWithVariableImposedDeformationAndFracture::changed() const
{
    return change ;
}

bool StiffnessWithVariableImposedDeformationAndFracture::fractured() const
{
    return frac;
}

Vector StiffnessWithVariableImposedDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return (param * imposed) ;
}

Vector StiffnessWithVariableImposedDeformationAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    return imposed ;
}

Form * StiffnessWithVariableImposedDeformationAndFracture::getCopy() const
{
    StiffnessWithVariableImposedDeformationAndFracture * copy = new StiffnessWithVariableImposedDeformationAndFracture(param, imposed, criterion->getCopy()) ;
    copy->damage = damage ;

    return copy ;
}

std::vector<BoundaryCondition * > StiffnessWithVariableImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    Vector f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;

    std::vector<BoundaryCondition * > ret ;
    if(f.size() == 2)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[1]));
    }
    if(f.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[2]));
    }
    return ret ;
}
