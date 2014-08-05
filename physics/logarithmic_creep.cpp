#include "logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

LogarithmicCreep::LogarithmicCreep(const Matrix & rig) : Viscoelasticity(PURE_ELASTICITY, rig, 1), C(rig), E(rig*0), tau(0), isPurelyElastic(true), accumulatedStress(0), currentStress(0), updated(true), accumulator(LOGCREEP_FORWARD)
{

}

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, const Matrix & vs, double e) : Viscoelasticity(MAXWELL, rig, vs), C(rig), E(vs), tau(e), isPurelyElastic(false), currentStress(0), accumulatedStress(0), updated(true), accumulator(LOGCREEP_FORWARD)
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
        vm->ieval(GradientDotDot(p_i) * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
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

/*    if(getExtra2dMeshes())
	{
		for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
			copy->addMesh((*getExtra2dMeshes())[i]);
	}
	if(getExtra3dMeshes())
	{
		for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
			copy->addMesh((*getExtra3dMeshes())[i]);
    }*/
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

    Matrix visc = E ;

    switch(accumulator)
    {
        case LOGCREEP_CONSTANT:
            break ;
        case LOGCREEP_FORWARD:
        {
            Vector stress(C.numRows()) ; stress = 0 ;
            currentState.getAverageField( REAL_STRESS_FIELD, stress, nullptr, -1, 1) ;
            currentStress = std::max(std::abs(stress[0]), std::abs(stress[1])) ;
            accumulatedStress += currentStress*timeStep ;
            if(accumulatedStress > POINT_TOLERANCE_2D)
            {
                if(currentStress < POINT_TOLERANCE_2D)
                    visc.array() = 0 ;
                else
                    visc *= (1.+accumulatedStress / (tau*currentStress)) ;
            }
            break ;
        }
        case LOGCREEP_PREDICTED:
        {
            Vector buffer0(C.numRows()) ; buffer0 = 0 ;
            Vector buffer1(C.numRows()) ; buffer1 = 0 ;
            currentState.getAverageField( REAL_STRESS_FIELD, buffer0, nullptr, -1, -1) ;
            currentState.getAverageField( REAL_STRESS_FIELD, buffer1, nullptr, -1, 1) ;
            currentStress =  std::max(std::abs(buffer1[0]), std::abs(buffer1[1])) ;
            buffer1 = buffer1*2. - buffer0 ;
            double projStress = std::max(std::abs(buffer1[0]), std::abs(buffer1[1])) ;
            currentState.getAverageField( STRAIN_RATE_FIELD, buffer0, nullptr, -1, 1) ;
            buffer0 = E*buffer0 ;
            double projElasticStress = std::max(std::abs(buffer0[0]), std::abs(buffer0[1])) ;
            if(currentStress > POINT_TOLERANCE_2D)
            {
                if(projStress < POINT_TOLERANCE_2D || projElasticStress < POINT_TOLERANCE_2D)
                    visc.array() = 0 ;
                else
                    visc *= ( (1.+timeStep/(2.*tau)) + currentStress/projStress*(currentStress/projElasticStress - 1. + timeStep/(2.*tau))  ) ;
            }
            break ;
        }
    }


    placeMatrixInBlock(visc, 1,1, eta);
    isPurelyElastic = (visc.array().max() < POINT_TOLERANCE_2D) ;

    currentState.getParent()->behaviourUpdated = true ;

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp  ) : LogarithmicCreep(rig), imposed(imp)
{

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, double e, const Vector & imp  ) : LogarithmicCreep(rig, v, e), imposed(imp)
{

}

Form * LogarithmicCreepWithImposedDeformation::getCopy() const
{
    LogarithmicCreepWithImposedDeformation * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreepWithImposedDeformation( C, imposed ) ;
    else
        copy = new LogarithmicCreepWithImposedDeformation( C, E, tau, imposed ) ;

 /*   if(getExtra2dMeshes())
    {
        for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
            copy->addMesh((*getExtra2dMeshes())[i]);
    }
    if(getExtra3dMeshes())
    {
        for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
            copy->addMesh((*getExtra3dMeshes())[i]);
    }*/
    return copy ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
        return imposed ;
    return Vector(0., C.numCols()) ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
        return C*imposed ;
    return Vector(0., C.numCols()) ;
}

std::vector<BoundaryCondition * > LogarithmicCreepWithImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(imposed.size() == 0)
        return ret ;
    if(v.size() == 3)
    {
        Vector istress = C * imposed   ;
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
    }
    if(v.size() == 4)
    {
        Vector istress = C * imposed ;
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
    }
    return ret ;

}


LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp) : LogarithmicCreepWithImposedDeformation(rig, imp), dfunc(nullptr), criterion(nullptr), noFracture(false)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double t, const Vector & imp) : LogarithmicCreepWithImposedDeformation(rig, v, t, imp), dfunc(nullptr), criterion(nullptr), noFracture(false)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d) : LogarithmicCreepWithImposedDeformation(rig, imp), dfunc(d), criterion(c), noFracture(!d || !c)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double t, const Vector & imp, FractureCriterion * c , DamageModel * d) : LogarithmicCreepWithImposedDeformation(rig, v, t, imp), dfunc(d), criterion(c), noFracture(!d || !c)
{

}

Form * LogarithmicCreepWithImposedDeformationAndFracture::getCopy() const
{
    LogarithmicCreepWithImposedDeformationAndFracture * copy ;

    if(noFracture)
    {
	    if(isPurelyElastic)
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed ) ;
	    else
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, tau, imposed ) ;
    } else {

	    if(isPurelyElastic)
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed, criterion->getCopy(), dfunc->getCopy() ) ;
	    else
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, tau, imposed, criterion->getCopy(), dfunc->getCopy() ) ;
	    copy->dfunc->getState(true).resize(dfunc->getState().size());
	    copy->dfunc->getState(true) = dfunc->getState() ;
	    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
    }

/*    if(getExtra2dMeshes())
    {
        for(size_t i = 0 ; i < getExtra2dMeshes()->size() ; i++)
            copy->addMesh((*getExtra2dMeshes())[i]);
    }
    if(getExtra3dMeshes())
    {
        for(size_t i = 0 ; i < getExtra3dMeshes()->size() ; i++)
            copy->addMesh((*getExtra3dMeshes())[i]);
    }*/

    return copy ;
}

void LogarithmicCreepWithImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
    Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

    Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
    getBlockInMatrix(param, 0,0, buffer) ;
    if(!noFracture)
        buffer = dfunc->apply(buffer) ;

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

void LogarithmicCreepWithImposedDeformationAndFracture::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    if(!isPurelyElastic)
    {
        Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
        Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

        Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
        getBlockInMatrix(eta, 1,1, buffer) ;
        if(!noFracture)
    	    buffer = dfunc->applyViscous(buffer) ;

        vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
        vm->ieval(GradientDotDot(p_i)    * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
        a += b ;
        placeMatrixInBlock( a, 1,1, ret ) ;
    }
}

void LogarithmicCreepWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    if(!noFracture)
    {
	    dfunc->step(currentState, maxscore) ;
	    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
    }
}

void LogarithmicCreepWithImposedDeformationAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        criterion = frac ;
    }
}

Matrix LogarithmicCreepWithImposedDeformationAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(noFracture)
	return param ;
    return  dfunc->apply(param) ;
}

Matrix LogarithmicCreepWithImposedDeformationAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(noFracture)
	return eta ;
    return dfunc->applyViscous(eta) ;
}

