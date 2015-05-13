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
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

    if(inGeometry->in(test))
        return inBehaviour->getTensor(p, e, g) ;

    return outBehaviour->getTensor(p, e, g) ;
}

Matrix BimaterialInterface::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;
    Point test = Point(vm.eval(xtransform, p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ytransform,  p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ztransform,  p.getX(), p.getY(), p.getZ(), p.getT()), vm.eval(ttransform, p.getX(),p.getY(),p.getZ(),p.getT())) ;

    if(inGeometry->in(test))
        return inBehaviour->getViscousTensor(p, e, g) ;

    return outBehaviour->getViscousTensor(p, e, g) ;
}

Vector BimaterialInterface::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{

    return outBehaviour->getImposedStress(p,e,g) ;
}

Vector BimaterialInterface::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{

    return outBehaviour->getImposedStrain(p,e,g) ;
}

void BimaterialInterface::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
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
//		std::cout << t[0] << "\t" << inCount << "\t" << gp.gaussPoints.size() << std::endl ;

    if(allin)
    {
        inBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
        return ;
    }
    else if(allout)
    {
//		std::cout << "all aout !" << std::endl ;
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

bool BimaterialInterface::fractured() const
{
    return inBehaviour->fractured() && outBehaviour->fractured();
}

Form * BimaterialInterface::getCopy() const
{
    BimaterialInterface * copy = new BimaterialInterface(*this) ;

    return copy ;
}

std::vector<BoundaryCondition * > BimaterialInterface::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;

    if(Jinv.size() == 0)
        return ret ;

    Vector x = VirtualMachine().eval(xtransform,gp) ;
    Vector y = VirtualMachine().eval(ytransform,gp) ;
    Vector z = VirtualMachine().eval(ztransform,gp) ;
    Vector t = VirtualMachine().eval(ttransform,gp) ;
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
        return inBehaviour->getBoundaryConditions(s,id, p_i, gp, Jinv) ;
    if(inCount == gp.gaussPoints.size())
        return outBehaviour->getBoundaryConditions(s,id, p_i, gp, Jinv) ;

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
// 	for(size_t i = 0 ; i < temp.size() ; i++)
// 		temp[i]->setData(temp[i]->getData()*1.5) ;
    ret.insert(ret.end(), temp.begin(), temp.end()) ;
    temp.clear() ;
    temp = outBehaviour->getBoundaryConditions(s,id, p_i, gpOut, outMatrixArray) ;
//	std::cout << gpOut.gaussPoints.size() << "/" << gp.gaussPoints.size() << std::endl ;
    ret.insert(ret.end(), temp.begin(), temp.end()) ;

//
//
    return ret ;
}


void BimaterialInterface::step(double timestep, ElementState & currentState, double maxScore)
{
    inBehaviour->step(timestep, currentState,maxScore) ;
    outBehaviour->step(timestep, currentState,maxScore) ;

}

DamageModel * BimaterialInterface::getDamageModel() const
{
    return nullptr ;
    double max = -2 ;
    int ret = 0 ;

    double outScore = 0. ;
    DamageModel * inCriterion = inBehaviour->getDamageModel() ;
    DamageModel * outCriterion = outBehaviour->getDamageModel() ;
    if(inCriterion)
    {
        max = inCriterion->getState().max() ;
        ret = 1 ;
    }

    if(outCriterion)
    {
        outScore = outCriterion->getState().max() ;
        if(outScore > max || inCriterion == nullptr)
            ret = 2 ;
    }

    switch(ret)
    {
    case 0:
        return nullptr ;
    case 1:
        return inCriterion ;
    case 2:
        return outCriterion ;
    }
    return nullptr ;
}

FractureCriterion * BimaterialInterface::getFractureCriterion() const
{
    return nullptr ;
    double max = -2 ;
    int ret = 0 ;

    double outScore = 0. ;
    FractureCriterion * inCriterion = inBehaviour->getFractureCriterion() ;
    FractureCriterion * outCriterion = outBehaviour->getFractureCriterion() ;
    if(inCriterion)
    {
        max = inCriterion->getScoreAtState() ;
        ret = 1 ;
    }

    if(outCriterion)
    {
        outScore = outCriterion->getScoreAtState() ;
        if(outScore > max || inCriterion == nullptr)
            ret = 2 ;
    }

    switch(ret)
    {
    case 0:
        return nullptr ;
    case 1:
        return inCriterion ;
    case 2:
        return outCriterion ;
    }
    return nullptr ;
}


Vector BimaterialInterface::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal)
{
    Vector x = VirtualMachine().eval(xtransform,gp) ;
    Vector y = VirtualMachine().eval(ytransform,gp) ;
    Vector z = VirtualMachine().eval(ztransform,gp) ;
    Vector t = VirtualMachine().eval(ttransform,gp) ;
//  	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
//  	{
// 		 std::cout << std::sqrt(x[i]*x[i]+y[i]*y[i]) << "\t" ;
// 		 if(inGeometry->in(Point(x[i],y[i],z[i],t[i])))
// 			std::cout << "IN" << std::endl ;
// 		 else
// 			std::cout << "OUT" << std::endl ;
//  	}
    if(inGeometry->in( Point(x[0],y[0],z[0],t[0]) ))
        return inBehaviour->getForcesFromAppliedStress( data, shape, gp, Jinv, v, isVolumic, normal) ;
    return outBehaviour->getForcesFromAppliedStress( data, shape, gp, Jinv, v, isVolumic, normal) ;
}

Vector BimaterialInterface::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic,  const Vector & normal)
{
    Vector x = VirtualMachine().eval(xtransform,gp) ;
    Vector y = VirtualMachine().eval(ytransform,gp) ;
    Vector z = VirtualMachine().eval(ztransform,gp) ;
    Vector t = VirtualMachine().eval(ttransform,gp) ;
    if(inGeometry->in( Point(x[0],y[0],z[0],t[0]) ))
        return inBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;
    return outBehaviour->getForcesFromAppliedStress( data, index, externaldofs, shape, e, gp, Jinv, v, isVolumic, normal) ;
}







