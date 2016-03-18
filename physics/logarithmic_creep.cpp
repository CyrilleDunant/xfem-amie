#include "logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, LogCreepAccumulator * acc) : Viscoelasticity(PURE_ELASTICITY, rig, 2), C(rig), E(rig*0), R(rig*0), tau(0), reducedTimeStep(-1.),  accumulator(acc), isPurelyElastic(true), updated(true), timeDependentIntegration(true),fixCreepVariable(false), prevParam(param), prevEta(eta)
{
		model = GENERALIZED_KELVIN_VOIGT ;
    makeBlockConnectivity() ;
}

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, const Matrix & vs, const Matrix & vr, double t, LogCreepAccumulator * acc) : Viscoelasticity(BURGER, rig, vr, vr*t, vs*t), C(rig), E(vs), R(vr), tau(t), reducedTimeStep(-1.), accumulator(acc), isPurelyElastic(false), updated(true),timeDependentIntegration(true), fixCreepVariable(false),  prevParam(param), prevEta(eta)
{
		model = GENERALIZED_KELVIN_VOIGT ;
    makeBlockConnectivity() ;
}

LogarithmicCreep::LogarithmicCreep(double young, double poisson, double E_creep, double nu_creep, double t, double E_recovery, double nu_recovery, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt) : Viscoelasticity(PURE_ELASTICITY, Tensor::cauchyGreen(young, poisson, dim, pt, YOUNG_POISSON), 2), tau(t), reducedTimeStep(-1.), accumulator(acc), updated(true), timeDependentIntegration(true),fixCreepVariable(false), prevParam(param), prevEta(eta)
{
    C = Tensor::cauchyGreen(young, poisson, dim, pt, YOUNG_POISSON) ;
    double nu = poisson ;
    if(nu_creep > 0 && E_creep > 0) { nu = nu_creep ; }
    E = Tensor::cauchyGreen(E_creep, nu, dim, pt, YOUNG_POISSON) ; if(E_creep < 0) { E = 0 ; }
    if(nu_recovery > 0 && tau > 0) { nu = nu_recovery ; }
    double y = E_creep ;
    if(E_recovery > 0) { y = E_recovery ; }
    R = Tensor::cauchyGreen(y, nu, dim, pt, YOUNG_POISSON) ; if(y < 0) { R = 0 ; }

    if(tau < 0 || E_creep < 0) 
    { 
        tau = 0 ;
        E = 0 ;
        R = 0 ;
        isPurelyElastic = true ;
    }
    else
    {
		double a = R[1][1]/E[1][1] ;
		tau /= exp(-a) ;
    }

        param = 0. ;
        eta = 0. ;
        placeMatrixInBlock(C,0,0,param) ;
        if(!isPurelyElastic)
        {
            Matrix S = -C ;
            placeMatrixInBlock(S,0,1,param) ;
            placeMatrixInBlock(S,0,2,param) ;
            placeMatrixInBlock(S,1,0,param) ;
            placeMatrixInBlock(S,2,0,param) ;
            placeMatrixInBlock(C,1,2,param) ;
            placeMatrixInBlock(C,2,1,param) ;
		placeMatrixInBlock( C+E*accumulator->getKelvinVoigtSpringReduction(), 1,1, param) ;
		placeMatrixInBlock( E*tau*accumulator->getKelvinVoigtDashpotReduction(), 1,1, eta) ;
		placeMatrixInBlock( C+R, 2,2, param) ;
		placeMatrixInBlock( R*tau, 2,2, eta) ;

	model = GENERALIZED_KELVIN_VOIGT ;
        }

    makeBlockConnectivity() ;

}

Form * LogarithmicCreep::getCopy() const
{
    LogarithmicCreep * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreep( C, accumulator->getCopy() ) ;
    else
        copy = new LogarithmicCreep( C, E, R, tau, accumulator->getCopy() ) ;
    copy->setBlocks(blocks) ;
    return copy ;
}

void LogarithmicCreep::print() const
{
    std::cout << "I am a logarithmic viscoelastic law" << std::endl ;
}

void LogarithmicCreep::preProcess(double timeStep, ElementState & currentState)
{
    if(reducedTimeStep < 0)
    {
        prevParam = param ;
        prevEta = eta ;

        param = 0. ;
        eta = 0. ;
        placeMatrixInBlock(C,0,0,param) ;
	model = PURE_ELASTICITY ;
        if(!isPurelyElastic)
        {
            accumulator->preProcess(timeStep, currentState) ;
            Matrix S = -C ;
            placeMatrixInBlock(S,0,1,param) ;
            placeMatrixInBlock(S,0,2,param) ;
            placeMatrixInBlock(S,1,0,param) ;
            placeMatrixInBlock(S,2,0,param) ;
            placeMatrixInBlock(C,1,2,param) ;
            placeMatrixInBlock(C,2,1,param) ;
		placeMatrixInBlock( C+E*accumulator->getKelvinVoigtSpringReduction(), 1,1, param) ;
		placeMatrixInBlock( E*tau*accumulator->getKelvinVoigtDashpotReduction(), 1,1, eta) ;
		placeMatrixInBlock( C+R, 2,2, param) ;
		placeMatrixInBlock( R*tau, 2,2, eta) ;

		model = GENERALIZED_KELVIN_VOIGT ;
        }
	currentState.getParent()->behaviourUpdated = true ;
    }

}

void LogarithmicCreep::step(double timestep, ElementState &s, double maxScore)
{

}

std::vector<BoundaryCondition * > LogarithmicCreep::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret ;
    if(fixCreepVariable)
    {
        size_t k = 0 ;
        ElementarySurface * tri = dynamic_cast<ElementarySurface *>(s.getParent()) ;
        for(size_t i = tri->getBoundingPoints().size()/tri->timePlanes() ; i < tri->getBoundingPoints().size() ; i++)
        {
            if((int)id == tri->getBoundingPoint(i).getId())
            {
                k = i-tri->getBoundingPoints().size()/tri->timePlanes() ;
                break ;
            }
        }
        if(k != 0)
        {
            ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*6+2], 2)) ;
            ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*6+3], 3)) ;
            ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*6+4], 4)) ;
            ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*6+5], 5)) ;
        }
    }
    return ret ;
}


LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc   ) : LogarithmicCreep(rig, acc), imposed(imp), prevImposed(imp.size())
{

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc   ) : LogarithmicCreep(rig, v, r, e, acc), imposed(imp), prevImposed(imp.size())
{

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation(double young, double poisson, double alpha, double E_creep, double nu_creep, double tau, double E_recovery, double nu_recovery, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt) : LogarithmicCreep( young, poisson, E_creep, nu_creep, tau, E_recovery, nu_recovery, acc, dim, pt), imposed(C.size()), prevImposed(C.size())
{
    for(size_t i = 0 ; i < dim ; i++)
        imposed[i] = alpha ;
}


Form * LogarithmicCreepWithImposedDeformation::getCopy() const
{
    LogarithmicCreepWithImposedDeformation * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreepWithImposedDeformation( C, imposed, accumulator->getCopy() ) ;
    else
        copy = new LogarithmicCreepWithImposedDeformation( C, E, R, tau, imposed, accumulator->getCopy() ) ;
    copy->setBlocks(blocks) ;
    return copy ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
    {
        /*		if(timeDependentIntegration)
        		{
        			Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
        			Function f2 = accumulator->getKelvinVoigtNextFunction() ;
        			Point p_ = p ;
        			if(e && g > -1)
        				p_ = e->getGaussPoints().gaussPoints[g].first ;
        			return prevImposed*VirtualMachine().eval(f1, p_) + imposed*VirtualMachine().eval(f2, p_) ;
        		}*/
//        std::cout << imposed[0]*(p.getT()+1)*0.5 + prevImposed[0]*(1-p.getT())*0.5 << std::endl ;
	return imposed*(p.getT()+1)*0.5 + prevImposed*(1-p.getT())*0.5 ;
    }
    return Vector(0., C.numCols()) ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
        return C*getImposedStrain(p,e,g) ;
    return Vector(0., C.numCols()) ;
}

std::vector<BoundaryCondition * > LogarithmicCreepWithImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = LogarithmicCreep::getBoundaryConditions(s, id, p_i, gp, Jinv ) ;
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

void LogarithmicCreepWithImposedDeformation::preProcess( double timeStep, ElementState & currentState )
{
    prevImposed = imposed ;
    LogarithmicCreep::preProcess(timeStep, currentState) ;
}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, imp, acc), dfunc(nullptr), criterion(nullptr), noFracture(false)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double t, const Vector & imp, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, v, r, t, imp, acc), dfunc(nullptr), criterion(nullptr), noFracture(false)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c)
{
   timeDependentIntegration = true ;
}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double t, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, v, r, t, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c)
{
   timeDependentIntegration = true ;
}

Form * LogarithmicCreepWithImposedDeformationAndFracture::getCopy() const
{
    LogarithmicCreepWithImposedDeformationAndFracture * copy ;

    if(noFracture)
    {
        if(isPurelyElastic)
            copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed, accumulator->getCopy() ) ;
        else
            copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, R, tau, imposed, accumulator->getCopy() ) ;
    } else {

        if(isPurelyElastic)
            copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
        else
            copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, R, tau, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
        copy->dfunc->getState(true).resize(dfunc->getState().size());
        copy->dfunc->getState(true) = dfunc->getState() ;
        copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
        copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
        copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
    }
    copy->setBlocks(blocks) ;
    return copy ;
}



void LogarithmicCreepWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    if(!noFracture)
    {
        dfunc->step(currentState, maxscore) ;
        currentState.getParent()->behaviourUpdated = dfunc->changed() ;
        currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;
    }
    if(maxscore > 0)
    {
        reducedTimeStep = currentState.getNodalDeltaTime() * (1. - maxscore) ;
        prevImposed += (imposed-prevImposed)*(1.-maxscore) ;
    }
    else
    {
        prevImposed = imposed ;
        reducedTimeStep = -1. ;
    }
}

void LogarithmicCreepWithImposedDeformationAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        criterion = frac ;
    }
}

Vector LogarithmicCreepWithImposedDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size() && !fractured())
    {
        if(dfunc)
            return (dfunc->apply(C))*getImposedStrain(p,e,g) ;
        return C*getImposedStrain(p,e,g) ;
    }
    return Vector(0., C.numCols()) ;
}

Matrix LogarithmicCreep::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
}

Matrix LogarithmicCreep::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
}

Matrix LogarithmicCreepWithImposedDeformation::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
}

Matrix LogarithmicCreepWithImposedDeformation::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
}

Matrix LogarithmicCreepWithImposedDeformationAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(noFracture)
        return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
    Matrix realParam = param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
    return  dfunc->apply(realParam, p,e,g) ;
}

Matrix LogarithmicCreepWithImposedDeformationAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(noFracture)
        return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
    Matrix realeta = eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
    return dfunc->applyViscous(realeta, p,e,g) ;
}

LogarithmicCreepWithImposedDeformationAndFracture::~LogarithmicCreepWithImposedDeformationAndFracture()
{
    if(dfunc)
    {
        delete dfunc ;
    }
    if(criterion)
    {
        delete criterion ;
    }
}

std::vector<BoundaryCondition * > LogarithmicCreepWithImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = LogarithmicCreep::getBoundaryConditions(s, id, p_i, gp, Jinv ) ;
    if(imposed.size() == 0 || fractured())
        return ret ;
    Vector istress = C*imposed ;
    if(dfunc)
        istress = dfunc->apply(C) * imposed   ;
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
    }
    if(v.size() == 4)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
    }
    return ret ;

}
