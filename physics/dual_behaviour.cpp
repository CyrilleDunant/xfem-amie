// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "dual_behaviour.h"
#include "homogeneised_behaviour.h"
#include "viscoelasticity.h"
#include "fracturecriteria/fracturecriterion.h"
#include "../geometry/space_time_geometry_2D.cpp"
#include <typeinfo>

using namespace Amie ;

BimaterialInterface::BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour) : LinearForm(Matrix(outbehaviour->param), true, true, outbehaviour->getNumberOfDegreesOfFreedom()), inGeometry(in),inBehaviour(inbehaviour), outBehaviour(outbehaviour)  { }

BimaterialInterface::~BimaterialInterface() { }

void BimaterialInterface::transform(ElementarySurface * e)
{
    xtransform = e->getXTransform() ;
    ytransform = e->getYTransform() ;
    ztransform = Function("0") ;
    if(inGeometry->timePlanes() > 1)
    {
        ttransform = e->getTTransform() ;
    }
    else
        ttransform = Function("0") ;
}

void BimaterialInterface::transform(ElementaryVolume * e)
{
    xtransform = e->getXTransform() ;
    ytransform = e->getYTransform() ;
    ztransform = e->getZTransform() ;
    if(inGeometry->timePlanes() > 1)
        ttransform = e->getTTransform() ;
    else
        ttransform = Function("0") ;
}

Matrix BimaterialInterface::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

    if(inGeometry->in(test))
        return inBehaviour->getTensor(p, e, g) ;

    return outBehaviour->getTensor(p, e, g) ;
}

Matrix BimaterialInterface::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

    if(inGeometry->in(test))
        return inBehaviour->getViscousTensor(p, e, g) ;

    return outBehaviour->getViscousTensor(p, e, g) ;
}

Vector BimaterialInterface::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ytransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ztransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ttransform, p.getX(), p.getY(), p.getZ(), p.getT())) ;


    if(inGeometry->in(test))
    {
        return inBehaviour->getImposedStress(p,e,g) ;
    }
    
    return outBehaviour->getImposedStress(p,e,g) ;
}

Vector BimaterialInterface::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()), 
                       vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

    if(inGeometry->in(test))
        inBehaviour->getImposedStrain(p,e,g) ;
    
    return outBehaviour->getImposedStrain(p,e,g) ;
}

void BimaterialInterface::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    ret = 0 ;

    bool allin = true ;
    bool allout = true ;
    Vector x = vm->eval(xtransform,gp) ;
    Vector y = vm->eval(ytransform,gp) ;
    Vector z = vm->eval(ztransform,gp) ;
    Vector t = vm->eval(ttransform,gp) ;
    std::valarray<bool> inIn(false, x.size()) ;
    int inCount = 0;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
        {
            allout = false ;
            inIn[i] = true ;
            inCount++ ;
        }
        else
        {
            allin = false ;
        }

    }

    if(allin)
    {
        inBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
        return ;
    }
    else if(allout)
    {
        outBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
        return ;
    }

    std::valarray<std::pair<Point, double> > inArray(inCount) ;
    std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()), inCount) ;
    std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
    std::valarray<Matrix> outMatrixArray(Matrix(Jinv[0].numRows(), Jinv[0].numCols()) , gp.gaussPoints.size()-inCount) ;
    GaussPointArray gpIn(inArray, -1) ;
    GaussPointArray gpOut(outArray, -1) ;

    int outIterator = 0;
    int inIterator = 0 ;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inIn[i])
        {
            inMatrixArray[inIterator] = Jinv[i] ;
            gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
        }
        else
        {
            outMatrixArray[outIterator] = Jinv[i] ;
            gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
        }
    }

    Matrix retIn(ret) ;
    inBehaviour->apply(p_i, p_j, gpIn, inMatrixArray, retIn,vm) ;
    outBehaviour->apply(p_i, p_j, gpOut, outMatrixArray,ret,vm) ;
    ret += retIn ;
}

void BimaterialInterface::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    ret = 0 ;

    bool allin = true ;
    bool allout = true ;
    Vector x = vm->eval(xtransform,gp) ;
    Vector y = vm->eval(ytransform,gp) ;
    Vector z = vm->eval(ztransform,gp) ;
    Vector t = vm->eval(ttransform,gp) ;
//	std::cout << sqrt(x[0]*x[0] + y[0]*y[0]) << "," << t[0] << " || " ;
    std::valarray<bool> inIn(false, x.size()) ;
    int inCount = 0;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
        {
            allout = false ;
            inIn[i] = true ;
            inCount++ ;
        }
        else
        {
            allin = false ;
        }

    }
//		std::cout << t[0] << "\t" << inCount << "\t" << gp.gaussPoints.size() << std::endl ;

    if(allin)
    {
        inBehaviour->applyViscous(p_i, p_j, gp, Jinv,ret,vm) ;
        return ;
    }
    else if(allout)
    {
//		std::cout << "all aout !" << std::endl ;
        outBehaviour->applyViscous(p_i, p_j, gp, Jinv,ret,vm) ;
        return ;
    }

    std::valarray<std::pair<Point, double> > inArray(inCount) ;
    std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()), inCount) ;
    std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
    std::valarray<Matrix> outMatrixArray(Matrix(Jinv[0].numRows(), Jinv[0].numCols()) , gp.gaussPoints.size()-inCount) ;
    GaussPointArray gpIn(inArray, -1) ;
    GaussPointArray gpOut(outArray, -1) ;

    int outIterator = 0;
    int inIterator = 0 ;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inIn[i])
        {
            inMatrixArray[inIterator] = Jinv[i] ;
            gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
        }
        else
        {
            outMatrixArray[outIterator] = Jinv[i] ;
            gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
        }
    }
    Matrix retIn(ret) ;
    inBehaviour->applyViscous(p_i, p_j, gpIn, inMatrixArray, retIn,vm) ;
    outBehaviour->applyViscous(p_i, p_j, gpOut, outMatrixArray,ret,vm) ;
    ret += retIn ;

}

Form * BimaterialInterface::getCopy() const
{
    BimaterialInterface * copy = new BimaterialInterface(*this) ;

    return copy ;
}

bool BimaterialInterface::changed() const
{
    if(getDamageModel())
        return getDamageModel()->changed() ;
    return false ;
}

bool BimaterialInterface::fractured() const
{
//     if(getDamageModel())
//         return getDamageModel()->fractured() ;
    return false ;
}

std::vector<BoundaryCondition * > BimaterialInterface::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;

    if(Jinv.size() == 0)
        return ret ;
     
    VirtualMachine vm ;
    Vector x = vm.eval(xtransform,gp) ;
    Vector y = vm.eval(ytransform,gp) ;
    Vector z = vm.eval(ztransform,gp) ;
    Vector t = vm.eval(ttransform,gp) ;
    std::valarray<bool> inIn(false, x.size()) ;
    size_t inCount = 0;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
        {
            inIn[i] = true ;
            inCount++ ;
        }
    }

    if(inCount == 0)
        return outBehaviour->getBoundaryConditions(s,id, p_i, gp, Jinv) ;
    if(inCount == gp.gaussPoints.size())
        return inBehaviour->getBoundaryConditions(s,id, p_i, gp, Jinv) ;
    
    std::valarray<std::pair<Point, double> > inArray(inCount) ;
    std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),inCount) ;
    std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
    std::valarray<Matrix> outMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),gp.gaussPoints.size()-inCount) ;
    GaussPointArray gpIn(inArray, -1) ;
    GaussPointArray gpOut(outArray, -1) ;

    int outIterator = 0;
    int inIterator = 0 ;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inIn[i])
        {
            inMatrixArray[inIterator] = Jinv[i] ;
            gpIn.gaussPoints[inIterator] = gp.gaussPoints[i] ;
            inIterator++ ;
        }
        else
        {
            outMatrixArray[outIterator] = Jinv[i] ;
            gpOut.gaussPoints[outIterator] = gp.gaussPoints[i] ;
            outIterator++ ;
        }
    }

    std::vector<BoundaryCondition * > temp = inBehaviour->getBoundaryConditions(s,id, p_i, gpIn, inMatrixArray) ;

    ret.insert(ret.end(), temp.begin(), temp.end()) ;
    temp.clear() ;
    temp = outBehaviour->getBoundaryConditions(s,id, p_i, gpOut, outMatrixArray) ;

    ret.insert(ret.end(), temp.begin(), temp.end()) ;

    return ret ;
}


void BimaterialInterface::step(double timestep, ElementState & currentState, double maxScore)
{
    inBehaviour->step(timestep, currentState,maxScore) ;
    outBehaviour->step(timestep, currentState,maxScore) ;

}

DamageModel * BimaterialInterface::getDamageModel() const
{
    DamageModel * inDamage = inBehaviour->getDamageModel() ;
    DamageModel * outDamage = outBehaviour->getDamageModel() ;
    if(inDamage && !outDamage  /*&& !inBehaviour->fractured()*/)
        return inDamage ;
//     else if (inDamage && !outDamage  && inBehaviour->fractured())
//         return nullptr ;
    
    if(outDamage && !inDamage /*&& !outBehaviour->fractured()*/)
        return outDamage ;
//     else if (outDamage && !inDamage&& outBehaviour->fractured())
//         return nullptr ;
//     
    if(!inDamage && !outDamage)
        return nullptr ;
    
    if(inBehaviour->getFractureCriterion()->getScoreAtState() > outBehaviour->getFractureCriterion()->getScoreAtState()) 
        return inDamage ;
        
    return outDamage ;
  
}

FractureCriterion * BimaterialInterface::getFractureCriterion() const
{
    FractureCriterion * inCriterion = inBehaviour->getFractureCriterion() ;
    FractureCriterion * outCriterion = outBehaviour->getFractureCriterion() ;
    if(inCriterion && !outCriterion /*&& !inBehaviour->fractured()*/)
        return inCriterion ;
//     else if (inCriterion && !outCriterion  && inBehaviour->fractured())
//         return nullptr ;
    
    if(outCriterion && !inCriterion /*&& !outBehaviour->fractured()*/)
        return outCriterion ;
//     else if (outCriterion && !inCriterion&& outBehaviour->fractured())
//         return nullptr ;
    
    if(!inCriterion && !outCriterion)
        return nullptr ;
    
    if(inCriterion->getScoreAtState() > outCriterion->getScoreAtState()) 
        return inCriterion ;
        
    return outCriterion ;
   
}


Vector BimaterialInterface::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal)
{
    Vector ret(v.size(), 0.) ;

    if(Jinv.size() == 0)
        return ret ;
     
    VirtualMachine vm ;
    Vector x = vm.eval(xtransform,gp) ;
    Vector y = vm.eval(ytransform,gp) ;
    Vector z = vm.eval(ztransform,gp) ;
    Vector t = vm.eval(ttransform,gp) ;
    std::valarray<bool> inIn(false, x.size()) ;
    size_t inCount = 0;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
        {
            inIn[i] = true ;
            inCount++ ;
        }
    }

    if(inCount == 0)
        return outBehaviour->getForcesFromAppliedStress(data, shape, gp, Jinv, v, isVolumic, normal) ;
    if(inCount == gp.gaussPoints.size())
        return inBehaviour->getForcesFromAppliedStress(data, shape, gp, Jinv, v, isVolumic, normal) ;
    
    std::valarray<std::pair<Point, double> > inArray(inCount) ;
    std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),inCount) ;
    std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
    std::valarray<Matrix> outMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),gp.gaussPoints.size()-inCount) ;
    GaussPointArray gpIn(inArray, -1) ;
    GaussPointArray gpOut(outArray, -1) ;

    int outIterator = 0;
    int inIterator = 0 ;

    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inIn[i])
        {
            inMatrixArray[inIterator] = Jinv[i] ;
            gpIn.gaussPoints[inIterator] = gp.gaussPoints[i] ;
            inIterator++ ;
        }
        else
        {
            outMatrixArray[outIterator] = Jinv[i] ;
            gpOut.gaussPoints[outIterator] = gp.gaussPoints[i] ;
            outIterator++ ;
        }
    }
    
    return inBehaviour->getForcesFromAppliedStress(data, shape, gpIn, inMatrixArray, v, isVolumic, normal) + 
           outBehaviour->getForcesFromAppliedStress(data, shape, gpOut, outMatrixArray, v, isVolumic, normal);
    
//     VirtualMachine vm ;
//     double x = vm.eval(xtransform,gp.gaussPoints[0].first) ;
//     double y = vm.eval(ytransform,gp.gaussPoints[0].first) ;
//     double z = vm.eval(ztransform,gp.gaussPoints[0].first) ;
//     double t = vm.eval(ttransform,gp.gaussPoints[0].first) ;
// 
//     if(inGeometry->in( Point(x,y,z,t) ))
//         return inBehaviour->getForcesFromAppliedStress( data, shape, gp, Jinv, v, isVolumic, normal) ;
//     return outBehaviour->getForcesFromAppliedStress( data, shape, gp, Jinv, v, isVolumic, normal) ;
}

Vector BimaterialInterface::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic,  const Vector & normal)
{
    Vector ret(v.size(), 0.) ;

    if(Jinv.size() == 0)
        return ret ;
     
    VirtualMachine vm ;
    Vector x = vm.eval(xtransform,gp) ;
    Vector y = vm.eval(ytransform,gp) ;
    Vector z = vm.eval(ztransform,gp) ;
    Vector t = vm.eval(ttransform,gp) ;
    std::valarray<bool> inIn(false, x.size()) ;
    size_t inCount = 0;
    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
        {
            inIn[i] = true ;
            inCount++ ;
        }
    }

    if(inCount == 0)
        return outBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;
    if(inCount == gp.gaussPoints.size())
        return inBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;
    
    std::valarray<std::pair<Point, double> > inArray(inCount) ;
    std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),inCount) ;
    std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
    std::valarray<Matrix> outMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()),gp.gaussPoints.size()-inCount) ;
    GaussPointArray gpIn(inArray, -1) ;
    GaussPointArray gpOut(outArray, -1) ;

    int outIterator = 0;
    int inIterator = 0 ;

    for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
    {
        if(inIn[i])
        {
            inMatrixArray[inIterator] = Jinv[i] ;
            gpIn.gaussPoints[inIterator] = gp.gaussPoints[i] ;
            inIterator++ ;
        }
        else
        {
            outMatrixArray[outIterator] = Jinv[i] ;
            gpOut.gaussPoints[outIterator] = gp.gaussPoints[i] ;
            outIterator++ ;
        }
    }
    
    return inBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gpIn, inMatrixArray, v, isVolumic, normal) + 
           outBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e,  gpOut, outMatrixArray, v, isVolumic, normal);
    
/*    
    
    VirtualMachine vm ;
    double x = vm.eval(xtransform,gp.gaussPoints[0].first) ;
    double y = vm.eval(ytransform,gp.gaussPoints[0].first) ;
    double z = vm.eval(ztransform,gp.gaussPoints[0].first) ;
    double t = vm.eval(ttransform,gp.gaussPoints[0].first) ;
    if(inGeometry->in( Point(x,y,z,t) ))
        return inBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;
    return outBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;*/
}







