//
// C++ Implementation: weibull_distributed_stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "weibull_distributed_stiffness.h"
#include "physics_base.h"
#include "stiffness.h"
#include "stiffness_and_fracture.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/confinedmohrcoulomb.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/confinedmohrcoulombwithstrain.h"
#include "fracturecriteria/maxstrain.h"
#include "fracturecriteria/ruptureenergy.h"

using namespace Amie ;

WeibullDistributedStiffness::WeibullDistributedStiffness(double E, double nu, SpaceDimensionality dim, double down, double up, planeType pt , double var, double r, IsotropicMaterialParameters hooke) : LinearForm(Tensor::cauchyGreen(E, nu,dim, pt, hooke), true, true, dim), variability(var), down(down), up(up), E(E), nu(nu), dim(dim)
{
    materialRadius = r;

    v.push_back(XI);
    v.push_back(ETA);
    if(dim == SPACE_THREE_DIMENSIONAL)
        v.push_back(ZETA);
    damageModel = nullptr ;
}

WeibullDistributedStiffness::~WeibullDistributedStiffness() { } 

void WeibullDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool WeibullDistributedStiffness::fractured() const
{
    return false ;
}

Form * WeibullDistributedStiffness::getCopy() const
{
    std::default_random_engine generator(std::rand());
    std::weibull_distribution< double > distribution(5, 1);
    double weib = distribution(generator) ;
    double factor = 1. - variability + variability*weib ;
    StiffnessAndFracture * copy = new StiffnessAndFracture(
        param*factor,
        new NonLocalMCFT(
            down*factor ,
            E,
            materialRadius, UPPER_BOUND)
    ) ;
    if(damageModel)
    {
        delete copy->dfunc ;
        copy->dfunc = damageModel->getCopy() ;
        copy->dfunc->setThresholdDamageDensity(damageModel->getThresholdDamageDensity());
        copy->dfunc->setSecondaryThresholdDamageDensity(damageModel->getSecondaryThresholdDamageDensity());
    }
    copy->criterion->setMaterialCharacteristicRadius(materialRadius);


    return copy ;
}

void WeibullDistributedStiffness::setDamageModel(DamageModel * newmodel)
{
    damageModel = newmodel ;
}

