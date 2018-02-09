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

#include "stiffness_with_large_deformation.h"
#include "../features/boundarycondition.h"

namespace Amie {

StiffnessWithLargeDeformation::StiffnessWithLargeDeformation(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1) , largeDeformationTransformc(2,2)
{
    change = false ;
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);
    this->time_d = false ;
    
    largeDeformationTransformc[0][0] = 1 ;
    largeDeformationTransformc[1][1] = 1 ;
}

StiffnessWithLargeDeformation::StiffnessWithLargeDeformation(double E, double nu, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke) : LinearForm(Tensor::cauchyGreen(E, nu,dim, pt, hooke), false, false, dim),v(2), largeDeformationTransformc(2,2)
{
    change = false ;
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);
    
    largeDeformationTransformc[0][0] = 1 ;
    largeDeformationTransformc[1][1] = 1 ;

    this->time_d = false ;
}

StiffnessWithLargeDeformation::~StiffnessWithLargeDeformation() { }

void StiffnessWithLargeDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    std::valarray<Matrix> JJ(largeDeformationTransformc,Jinv.size()) ;
//     std::cout << 0.5/det(Jinv[0]) << "  " << gp.gaussPoints[0].second << std::endl ;
    GaussPointArray gpj = gp ;
    for(size_t i = 0 ;  i < gpj.gaussPoints.size() ; i++ )
        gpj.gaussPoints[i].second *= det(largeDeformationTransformc)/det(Jinv[i]) ;
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gpj, JJ,v,ret) ;
//     ret *= .5 ;
//     ret = largeDeformationTransformc*ret*largeDeformationTransformc.transpose()/**largeDeformationTransformc*/ ;
//     ret  = Tensor::rotate2ndOrderTensor2D(ret, globalTransformAngle) ;
}

bool StiffnessWithLargeDeformation::fractured() const
{
    return false ;
}

Form * StiffnessWithLargeDeformation::getCopy() const
{
    StiffnessWithLargeDeformation * copy = new StiffnessWithLargeDeformation(param) ;

    return copy ;
}

void StiffnessWithLargeDeformation::getInverseFiniteDeformationJacobianMatrix(const Point & p, Matrix & ret, ElementState & currentState)
{

    if(ret.isNull() || ret.size() != 4)
        ret.resize(2,2) ;
    
    VirtualMachine vm ;
    Vector deltaBaseDisplacements = lastConvergedDisplacements*0.5+(currentState.getDisplacements()-lastConvergedDisplacements)*0.5 ; 
        
    ret.array() = 0 ;
    for(size_t i = 0 ; i < currentState.getParent()->getBoundingPoints().size() ; i++)
    {
        double dxi = vm.deval(currentState.getParent()->getShapeFunction(i), XI, p) ;
        double deta = vm.deval(currentState.getParent()->getShapeFunction(i), ETA, p) ;
        
        ret[0][0] += dxi*(currentState.getParent()->getBoundingPoint(i).getX()+deltaBaseDisplacements[i*2]) ;
        ret[0][1] += dxi*(currentState.getParent()->getBoundingPoint(i).getY()+deltaBaseDisplacements[i*2+1]) ;

        ret[1][0] += deta*(currentState.getParent()->getBoundingPoint(i).getX()+deltaBaseDisplacements[i*2]) ;
        ret[1][1] += deta*(currentState.getParent()->getBoundingPoint(i).getY()+deltaBaseDisplacements[i*2+1]) ;
    }

    invert2x2Matrix(ret) ;
}


void StiffnessWithLargeDeformation::step(double timestep, ElementState & currentState, double maxscore)
{
    Matrix previousTransform = largeDeformationTransformc ;
    if(lastConvergedDisplacements.size() == 0)
    {
        currentState.getInverseJacobianMatrix(Point(1./3., 1./3.),previousTransform) ;
        norm = sqrt(std::inner_product(&previousTransform.array()[0], &previousTransform.array()[previousTransform.array().size()], &previousTransform.array()[0], 0.)) ;
        lastConvergedDisplacements.resize(currentState.getDisplacements().size(), 0.) ;
    }
    currentState.getInverseJacobianMatrix(Point(1./3., 1./3.), largeDeformationTransformc) ;
    
    change = false ;
    Vector delta = previousTransform.array()-largeDeformationTransformc.array() ;
    double err = std::inner_product(&delta[0], &delta[delta.size()], &delta[0], 0.) ;

    if(sqrt(err)/norm > 1e-8 || timestep > POINT_TOLERANCE)
        change = true ;
    
    currentState.getParent()->behaviourUpdated = change ;
    currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;
    
}

bool StiffnessWithLargeDeformation::changed() const
{
    return change ;
}

}
