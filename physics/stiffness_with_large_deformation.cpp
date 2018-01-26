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

StiffnessWithLargeDeformation::StiffnessWithLargeDeformation(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1) , largeDeformationTransformc(2,2), globalTransformAngle(0)
{
    poisson = 0.2 ;
    change = false ;
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);
    this->time_d = false ;
    
    largeDeformationTransformc[0][0] = 1 ;
    largeDeformationTransformc[1][1] = 1 ;
    
    if(v.size() == 2)
    {
        imposed.resize(3, 0.) ;
    }
    else 
    {
        imposed.resize(6, 0.) ;
    }
}

StiffnessWithLargeDeformation::StiffnessWithLargeDeformation(double E, double nu, SpaceDimensionality dim, planeType pt, IsotropicMaterialParameters hooke) : LinearForm(Tensor::cauchyGreen(E, nu,dim, pt, hooke), false, false, dim),v(2), largeDeformationTransformc(2,2), globalTransformAngle(0)
{
    poisson = nu ;
    change = false ;
    v.push_back(XI) ;
    v.push_back(ETA) ;
    if(param.size() == 36)
        v.push_back(ZETA);
    
    largeDeformationTransformc[0][0] = 1 ;
    largeDeformationTransformc[1][1] = 1 ;

    if(dim == SPACE_TWO_DIMENSIONAL)
    {
        imposed.resize(3, 0.) ;
    }
    else if(dim == SPACE_THREE_DIMENSIONAL)
    {
        imposed.resize(6, 0.) ;
    }

    this->time_d = false ;
}

StiffnessWithLargeDeformation::~StiffnessWithLargeDeformation() { }

void StiffnessWithLargeDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
    ret *= largeDeformationTransformc ;
    ret  = Tensor::rotate2ndOrderTensor2D(ret, globalTransformAngle) ;
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

void StiffnessWithLargeDeformation::step(double timestep, ElementState & currentState, double maxscore)
{
    double previousAngle = globalTransformAngle ;
    Matrix previousTransform = largeDeformationTransformc ;
    VirtualMachine vm ;
    Vector d(0.,2) ; 
    Point centre(1./3, 1./3.) ;
    Point gcentre = currentState.getParent()->getCenter() ;
    double bs = 0 ;
    double as = 0 ;
    Matrix Jinv(2,2) ;
    if(currentState.JinvCache)
        Jinv = *currentState.JinvCache ;
    else 
        currentState.getParent()->getInverseJacobianMatrix(gcentre, Jinv) ;
            
    for(size_t id = 0 ; id < currentState.getParent()->getBoundingPoints().size() ; id++)
    {
        double f_xi  = vm.deval ( currentState.getParent()->getShapeFunction( id ),  XI, centre ) ;
        double f_eta = vm.deval ( currentState.getParent()->getShapeFunction( id ), ETA, centre ) ;
        
        d[0] =  Jinv[1][0]*f_xi+ Jinv[1][1]*f_eta ;
        d[1] = - Jinv[0][0]*f_xi- Jinv[0][1]*f_eta ;
        
        Vector xy = {currentState.getDisplacements()[id*2+0] + currentState.getParent()->getBoundingPoint(id).getX() - gcentre.getX(), 
                        currentState.getDisplacements()[id*2+1] + currentState.getParent()->getBoundingPoint(id).getY() - gcentre.getY()} ;

        as += -d[1]*xy[0] +d[0]*xy[1] ;
        bs +=  d[0]*xy[0] +d[1]*xy[1] ;          
        
    }   
    globalTransformAngle = atan2(-bs, as) ;
       
    Matrix rot(2,2) ;
    rot[0][0] =  cos(globalTransformAngle) ; rot[0][1] = sin(globalTransformAngle) ; 
    rot[1][0] = -sin(globalTransformAngle) ; rot[1][1] = cos(globalTransformAngle) ; 
    
    Vector baseDisplacements = currentState.getDisplacements() ;
//     for(size_t i = 0 ; i < currentState.getParent()->getBoundingPoints().size() ; i++)
//     {
//         Vector ldisp = {baseDisplacements[i*2], baseDisplacements[i*2+1]} ;
//         ldisp = ldisp*rot ;
//         baseDisplacements[i*2]   = ldisp[0] ;
//         baseDisplacements[i*2+1] = ldisp[1] ;
//     }
    
    Matrix largeDeformationTransform(2,2) ;
//     largeDeformationTransform[0][0] = XdXTransform(baseDisplacements ,currentState.getParent()->getShapeFunctions(), XI , gcentre, 2) ;
//     largeDeformationTransform[0][1] = XdYTransform(baseDisplacements ,currentState.getParent()->getShapeFunctions(), XI , gcentre, 2) ;
//     largeDeformationTransform[1][0] = XdXTransform(baseDisplacements ,currentState.getParent()->getShapeFunctions(), ETA, gcentre, 2) ;
//     largeDeformationTransform[1][1] = XdYTransform(baseDisplacements ,currentState.getParent()->getShapeFunctions(), ETA, gcentre, 2) ;
    
    PointArray rotatedPoints(currentState.getParent()->getBoundingPoints().size()) ;
    for(size_t i = 0 ; i < rotatedPoints.size() ; i++)
    {
        rotatedPoints[i] = new Point(currentState.getParent()->getBoundingPoint(i)+Point(baseDisplacements[i*2],baseDisplacements[i*2+1] )*rot) ;
    }
    
    largeDeformationTransform[0][0] = dXTransform(rotatedPoints ,currentState.getParent()->getShapeFunctions(), XI , gcentre) ;
    largeDeformationTransform[0][1] = dYTransform(rotatedPoints ,currentState.getParent()->getShapeFunctions(), XI , gcentre) ;
    largeDeformationTransform[1][0] = dXTransform(rotatedPoints ,currentState.getParent()->getShapeFunctions(), ETA, gcentre) ;
    largeDeformationTransform[1][1] = dYTransform(rotatedPoints ,currentState.getParent()->getShapeFunctions(), ETA, gcentre) ;

    double j = 2.*currentState.getParent()->area() ;
     for(size_t i = 0 ; i < rotatedPoints.size() ; i++)
         delete rotatedPoints[i] ;
    
    invert2x2Matrix(Jinv) ; 
    largeDeformationTransformc[0][0] = 1. + largeDeformationTransform[0][0]* Jinv[0][0] + largeDeformationTransform[0][1]* Jinv[0][1] ;
    largeDeformationTransformc[0][1] =    + largeDeformationTransform[0][0]* Jinv[1][0] + largeDeformationTransform[0][1]* Jinv[1][1] ;
    largeDeformationTransformc[1][0] =    + largeDeformationTransform[1][0]* Jinv[0][0] + largeDeformationTransform[1][1]* Jinv[0][1] ;
    largeDeformationTransformc[1][1] = 1. + largeDeformationTransform[1][0]* Jinv[1][0] + largeDeformationTransform[1][1]* Jinv[1][1] ;
    double J = det(largeDeformationTransformc) ;
    invert2x2Matrix(largeDeformationTransformc) ;
    
    
//     Vector d0(2) ;
//     Vector d1(2) ;
//     Vector d2(2) ;
//     currentState.getField ( DISPLACEMENT_FIELD, currentState.getParent()->getBoundingPoint(0), d0, false, &vm, 0 ) ;
//     currentState.getField ( DISPLACEMENT_FIELD, currentState.getParent()->getBoundingPoint(1), d1, false, &vm, 0 ) ;
//     currentState.getField ( DISPLACEMENT_FIELD, currentState.getParent()->getBoundingPoint(2), d2, false, &vm, 0 ) ;
//     Triangle test(currentState.getParent()->getBoundingPoint(0)+Point(d0[0], d0[1]), 
//                   currentState.getParent()->getBoundingPoint(1)+Point(d1[0], d1[1]), 
//                   currentState.getParent()->getBoundingPoint(2)+Point(d2[0], d2[1])) ;
                
    largeDeformationTransformc /= J;
    
//     double dv = 1.-currentState.getParent()->area()/test.area() ;
//     imposed[0] = -largeDeformationTransform[0][0]* Jinv[0][0]*J - largeDeformationTransform[0][1]* Jinv[0][1]*J ;
//     imposed[1] = -largeDeformationTransform[1][0]* Jinv[1][0]*J - largeDeformationTransform[1][1]* Jinv[1][1]*J ;
//     imposed[2] = 0 ;
//     imposed *= .002 ;
    
    imposed = Tensor::rotate2ndOrderTensor2D( imposed, -globalTransformAngle ) ;

    change = false ;
    if(std::abs(previousAngle-globalTransformAngle) > 1e-6 || std::abs(previousTransform.array()-largeDeformationTransformc.array()).max() > 1e-8)
        change = true ;
    
}

Vector StiffnessWithLargeDeformation::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    return imposed*0. ;
}

bool StiffnessWithLargeDeformation::changed() const
{
    return change ;
}

Vector StiffnessWithLargeDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    return imposed ;
}

std::vector<BoundaryCondition * > StiffnessWithLargeDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    
//     if(v.size() == 2)
//     {
// //             Vector istress = VirtualMachine().ieval(Gradient(p_i)*(param * imposed),gp,Jinv, v)   ;
//         Vector istress = getTensor(Point(1./3., 1./3.)) * imposed   ;
//         ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
//         ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
//     }
    return ret ;
    
}

}
