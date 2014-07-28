#include "logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

LogarithmicCreep::LogarithmicCreep(const Matrix & rig) : Viscoelasticity(PURE_ELASTICITY, rig, 1), C(rig), E(rig*0), tau(0), isPurelyElastic(true), accumulatedStress(0), currentStress(0), updated(true)
{

}

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, const Matrix & vs, double e) : Viscoelasticity(MAXWELL, rig, vs), C(rig), E(vs), tau(e), isPurelyElastic(false), currentStress(0), accumulatedStress(0), updated(true)
{

}

void LogarithmicCreep::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{  
	Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
	Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	
	Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
    getBlockInMatrix(param, 0,0, buffer) ;

    vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
    vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
    a += b ;
			
    placeMatrixInBlock( a, 0,0, ret ) ;
    if(!isPurelyElastic)
    {
        b = a*(-1.) ;
        placeMatrixInBlock(b,1,0, ret);
        placeMatrixInBlock(b,0,1, ret);
        placeMatrixInBlock(a,1,1, ret);
    }
}

void LogarithmicCreep::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    if(!isPurelyElastic)
    {
        Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
        Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

        Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
        getBlockInMatrix(eta, 1,1, buffer) ;

        vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
        vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
        a += b ;
        placeMatrixInBlock( a, 1,1, ret ) ;
    }
}

Form * LogarithmicCreep::getCopy() const
{
    LogarithmicCreep * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreep( C ) ;
    else
        copy = new LogarithmicCreep( C, E, tau ) ;

    if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
	}
	return copy ; 
}

void LogarithmicCreep::print() const
{
    std::cout << "I am a logarithmic viscoelastic law" << std::endl ;
}

void LogarithmicCreep::preProcess(double timeStep, ElementState & currentState)
{
    updated = false ;
    if(E.array().max() < POINT_TOLERANCE_2D)
        return ;

    Vector stress(C.numRows()) ; stress = 0 ;
    currentState.getAverageField( REAL_STRESS_FIELD, stress, nullptr, -1, 1) ;
    Vector principalstress = toPrincipal(stress) ;

    currentStress = std::max(std::abs(principalstress.min()), principalstress.max()) ;
    accumulatedStress += currentStress*timeStep ;


    Matrix visc = E ;
    if(accumulatedStress > POINT_TOLERANCE_2D)
    {
        if(currentStress < POINT_TOLERANCE_2D)
            visc.array() = 0 ;
        else
            visc *= (1.+accumulatedStress / (tau*currentStress)) ;
    }

    placeMatrixInBlock(visc, 1,1, eta);

    isPurelyElastic = (visc.array().max() < POINT_TOLERANCE_2D) ;


    currentState.getParent()->behaviourUpdated = true ;

}

