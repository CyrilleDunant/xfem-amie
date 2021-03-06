//
// C++ Interface: stiffness_and_fracture
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"
#include "damagemodels/prandtlgrauertplasticstrain.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

using namespace Amie ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit, DamageModel * d) : LinearForm(rig, false, true, rig.numRows()/3+1)
{
    if(!d)
        dfunc = new FiberBasedIsotropicLinearDamage() ;
    else
        dfunc = d ;
    criterion = crit ;

    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36 )
    {
        v.push_back(ZETA);
    }
// 	v.push_back(TIME_VARIABLE);
}

StiffnessAndFracture::StiffnessAndFracture(double E, double nu, FractureCriterion * crit, DamageModel * d, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke) : LinearForm(Tensor::cauchyGreen(E, nu, dim, pt, hooke), false, true, dim)
{
    if(!d)
        dfunc = new FiberBasedIsotropicLinearDamage() ;
    else
        dfunc = d ;
    criterion = crit ;

    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() == 36 )
    {
        v.push_back(ZETA);
    }
// 	v.push_back(TIME_VARIABLE);
}

StiffnessAndFracture::~StiffnessAndFracture()
{
    delete criterion ;
    delete dfunc ;
}

FractureCriterion * StiffnessAndFracture::getFractureCriterion() const
{
    return criterion ;
}

DamageModel * StiffnessAndFracture::getDamageModel() const
{
    return dfunc ;
}

void StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

void StiffnessAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
//   if(timestep >= 0 && maxscore > 0)
//   {
    dfunc->step(currentState, maxscore) ;
    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
    currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;
//   }

}

bool StiffnessAndFracture::changed() const
{
    return dfunc->changed() ;
}

bool StiffnessAndFracture::fractured() const
{
    return dfunc->fractured() ;
}

Form * StiffnessAndFracture::getCopy() const
{
    StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy(), dfunc->getCopy()) ;
    copy->dfunc->getState(true).resize(dfunc->getState().size());
    copy->dfunc->getState(true) = dfunc->getState() ;
    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());

    return copy ;
}

Matrix StiffnessAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return dfunc->apply(param, p, e, g) ;
}

void StiffnessAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        delete criterion ;
        criterion = frac ;
    }

}


