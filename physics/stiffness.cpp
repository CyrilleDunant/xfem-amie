//
// C++ Implementation: stiffness
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/vonmises.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/nonlocalvonmises.h"
#include <valarray>


namespace Amie {

Stiffness::Stiffness(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() > 9)
    {
        v.push_back(ZETA);
    }
}

Stiffness::Stiffness(double E, double nu, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke) : LinearForm( Tensor::cauchyGreen(E, nu, dim, pt, hooke), false, false, dim)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(dim > 2)
        v.push_back(ZETA);
}

Stiffness::~Stiffness() { } 

void Stiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

bool Stiffness::fractured() const
{
    return false ;
}

Form * Stiffness::getCopy() const
{
    Stiffness* copy = new Stiffness( param) ;

    return copy ;
}

Form * WeibullDistributedElasticStiffness::getCopy() const
{
    std::default_random_engine generator(std::rand());
    std::weibull_distribution< double > distribution(5, 1);
    double weib = distribution(generator) ;
    double factor = 1. - variability + variability*weib ;
    return new Stiffness( param*factor ) ;
}

DerivedStiffness::DerivedStiffness(const Matrix & rig) : LinearForm(Matrix(), false, false, rig.numRows()/3+1), rig(rig)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(rig.size() > 9)
    {
        v.push_back(ZETA);
    }
} 

DerivedStiffness::~DerivedStiffness() { } 

void DerivedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

    vm->ieval(Gradient(p_i) * rig * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

bool DerivedStiffness::fractured() const
{
    return false ;
}

Form * DerivedStiffness::getCopy() const
{
    DerivedStiffness* copy = new DerivedStiffness( param) ;

    return copy ;
}


PseudoPlastic::PseudoPlastic(const Amie::Matrix& rig, double E, double limitStrain, double radius): LinearForm(rig, false, true, rig.numRows()/3+1), alpha(0), change(true), radius(radius), limitStrain(limitStrain)
{
    stiffness = E ;
    vm = new NonLocalDeviatoricVonMises(limitStrain, radius) ;
    vm->setMaterialCharacteristicRadius(radius);
    initialised = false ;
    lastDamage = alpha ;
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() > 9)
        v.push_back(ZETA);

    frac = false ;
    fixedfrac = false ;
}

bool PseudoPlastic::changed() const
{
    return change ;
}

void PseudoPlastic::fixLastDamage()
{
    lastDamage = alpha ;
    fixedfrac = frac ;
}

PseudoPlastic::~PseudoPlastic() {
    delete vm ;
} 


void PseudoPlastic::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * (param*(1.-alpha)) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void PseudoPlastic::step(double timestep, ElementState & currentState, double maxscore)
{
    if(timestep > POINT_TOLERANCE)
    {
        fixLastDamage() ;
    }

    frac = fixedfrac ;
    change = false ;
    double lastalpha = alpha ;

    Vector str = vm->getSmoothedField( PRINCIPAL_REAL_STRESS_FIELD, currentState) ;
    double maxStress = sqrt( ( ( str[0] - str[1] ) * ( str[0] - str[1] ) + str[0] * str[0] + str[1] * str[1] ) / 2. ) ;

    if(maxStress > POINT_TOLERANCE)
    {
        alpha = std::max(1.-(vm->threshold/maxStress)*(1.-alpha), lastDamage) ;
        change = std::abs(alpha-lastalpha) > 1e-6 ;
        currentState.getParent()->behaviourUpdated = change ;
        currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;
    }
}

Matrix PseudoPlastic::getTensor(const Point & p, IntegrableEntity *, int) const
{
    return (param*(1.-alpha)) ;
}


Matrix PseudoPlastic::getPreviousTensor(const Point & p) const
{
    return (param*lastDamage) ;
}

FractureCriterion * PseudoPlastic::getFractureCriterion() const
{
    return vm ;
}

bool PseudoPlastic::fractured() const
{
    return frac ;
}

Form * PseudoPlastic::getCopy() const
{
    return new PseudoPlastic(param, stiffness, limitStrain, radius) ;

}

}

