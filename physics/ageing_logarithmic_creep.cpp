#include "ageing_logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

AgeingLogarithmicCreep::AgeingLogarithmicCreep(const Matrix & rig, LogCreepAccumulator * acc) : Viscoelasticity(PURE_ELASTICITY, rig, 2), C(rig), E(rig*0), R(rig*0), tau(0), reducedTimeStep(-1.),  accumulator(acc), isPurelyElastic(true), updated(true), timeDependentIntegration(true),fixCreepVariable(false), prevParam(param), prevEta(eta)
{
		//model = GENERALIZED_KELVIN_VOIGT ;
    makeBlockConnectivity() ;
}

AgeingLogarithmicCreep::AgeingLogarithmicCreep(const Matrix & rig, const Matrix & vs, const Matrix & vr, double t, LogCreepAccumulator * acc) : Viscoelasticity(BURGER, rig, vr, vr*t, vs*t), C(rig), E(vs), R(vr), tau(t), reducedTimeStep(-1.), accumulator(acc), isPurelyElastic(false), updated(true),timeDependentIntegration(true), fixCreepVariable(false),  prevParam(param), prevEta(eta)
{
		//model = GENERALIZED_KELVIN_VOIGT ;
    makeBlockConnectivity() ;
}

AgeingLogarithmicCreep::AgeingLogarithmicCreep(double young, double poisson, double E_creep, double nu_creep, double t, double E_recovery, double nu_recovery, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt) : Viscoelasticity(PURE_ELASTICITY, Tensor::cauchyGreen(young, poisson, dim, pt, YOUNG_POISSON), 2), tau(t), reducedTimeStep(-1.), accumulator(acc), updated(true), timeDependentIntegration(true),fixCreepVariable(false), prevParam(param), prevEta(eta)
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

Form * AgeingLogarithmicCreep::getCopy() const
{
    AgeingLogarithmicCreep * copy ;
    if(isPurelyElastic)
        copy = new AgeingLogarithmicCreep( C, accumulator->getCopy() ) ;
    else
        copy = new AgeingLogarithmicCreep( C, E, R, tau, accumulator->getCopy() ) ;
    copy->setBlocks(blocks) ;
    return copy ;
}

void AgeingLogarithmicCreep::print() const
{
    std::cout << "I am a logarithmic viscoelastic law" << std::endl ;
}

void AgeingLogarithmicCreep::preProcess(double timeStep, ElementState & currentState)
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
            placeMatrixInBlock( C+E*agefactor*accumulator->getKelvinVoigtSpringReduction(), 1,1, param) ;
            placeMatrixInBlock( E*agefactor*tau*accumulator->getKelvinVoigtDashpotReduction(), 1,1, eta) ;
            placeMatrixInBlock( C+R, 2,2, param) ;
            placeMatrixInBlock( R*tau, 2,2, eta) ;

		model = GENERALIZED_KELVIN_VOIGT ;
        }
	currentState.getParent()->behaviourUpdated = true ;
    }

}

void AgeingLogarithmicCreep::step(double timestep, ElementState &s, double maxScore)
{
   age+=timestep;
   agefactor=1-exp(-age/5);
}

std::vector<BoundaryCondition * > AgeingLogarithmicCreep::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
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


AgeingLogarithmicCreepWithImposedDeformation::AgeingLogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc   ) : AgeingLogarithmicCreep(rig, acc), imposed(imp), prevImposed(imp.size())
{

}

AgeingLogarithmicCreepWithImposedDeformation::AgeingLogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, const Matrix & r, double e, const Vector & imp, LogCreepAccumulator * acc   ) : AgeingLogarithmicCreep(rig, v, r, e, acc), imposed(imp), prevImposed(imp.size())
{

}

AgeingLogarithmicCreepWithImposedDeformation::AgeingLogarithmicCreepWithImposedDeformation(double young, double poisson, double alpha, double E_creep, double nu_creep, double tau, double E_recovery, double nu_recovery, LogCreepAccumulator * acc, SpaceDimensionality dim, planeType pt) : AgeingLogarithmicCreep( young, poisson, E_creep, nu_creep, tau, E_recovery, nu_recovery, acc, dim, pt), imposed(C.size()), prevImposed(C.size())
{
    for(size_t i = 0 ; i < dim ; i++)
        imposed[i] = alpha ;
}


Form * AgeingLogarithmicCreepWithImposedDeformation::getCopy() const
{
    AgeingLogarithmicCreepWithImposedDeformation * copy ;
    if(isPurelyElastic)
        copy = new AgeingLogarithmicCreepWithImposedDeformation( C, imposed, accumulator->getCopy() ) ;
    else
        copy = new AgeingLogarithmicCreepWithImposedDeformation( C, E, R, tau, imposed, accumulator->getCopy() ) ;
    copy->setBlocks(blocks) ;
    return copy ;
}

Vector AgeingLogarithmicCreepWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
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

Vector AgeingLogarithmicCreepWithImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
        return C*getImposedStrain(p,e,g) ;
    return Vector(0., C.numCols()) ;
}

std::vector<BoundaryCondition * > AgeingLogarithmicCreepWithImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = AgeingLogarithmicCreep::getBoundaryConditions(s, id, p_i, gp, Jinv ) ;
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

void AgeingLogarithmicCreepWithImposedDeformation::preProcess( double timeStep, ElementState & currentState )
{
    prevImposed = imposed ;
    AgeingLogarithmicCreep::preProcess(timeStep, currentState) ;
}

AgeingLogarithmicCreepWithImposedDeformationAndFracture::AgeingLogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc) : AgeingLogarithmicCreepWithImposedDeformation(rig, imp, acc), dfunc(nullptr), criterion(nullptr), noFracture(false), etaSet(false), paramSet(false)
{
    currentEta = eta ;
    currentParam = param ;
}

AgeingLogarithmicCreepWithImposedDeformationAndFracture::AgeingLogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double t, const Vector & imp, LogCreepAccumulator * acc) : AgeingLogarithmicCreepWithImposedDeformation(rig, v, r, t, imp, acc), dfunc(nullptr), criterion(nullptr), noFracture(false), etaSet(false), paramSet(false)
{
    currentEta = eta ;
    currentParam = param ;

}

AgeingLogarithmicCreepWithImposedDeformationAndFracture::AgeingLogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : AgeingLogarithmicCreepWithImposedDeformation(rig, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c), etaSet(false), paramSet(false)
{
    currentEta = eta ;
    currentParam = param ;
   timeDependentIntegration = true ;
}

AgeingLogarithmicCreepWithImposedDeformationAndFracture::AgeingLogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, const Matrix & r, double t, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : AgeingLogarithmicCreepWithImposedDeformation(rig, v, r, t, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c), etaSet(false), paramSet(false)
{
    currentEta = eta ;
    currentParam = param ;
   timeDependentIntegration = true ;
}

Form * AgeingLogarithmicCreepWithImposedDeformationAndFracture::getCopy() const
{
    AgeingLogarithmicCreepWithImposedDeformationAndFracture * copy ;

    if(noFracture)
    {
        if(isPurelyElastic)
            copy = new AgeingLogarithmicCreepWithImposedDeformationAndFracture( C, imposed, accumulator->getCopy() ) ;
        else
            copy = new AgeingLogarithmicCreepWithImposedDeformationAndFracture( C, E, R, tau, imposed, accumulator->getCopy() ) ;
    } else {

        if(isPurelyElastic)
            copy = new AgeingLogarithmicCreepWithImposedDeformationAndFracture( C, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
        else
            copy = new AgeingLogarithmicCreepWithImposedDeformationAndFracture( C, E, R, tau, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
        copy->dfunc->getState(true).resize(dfunc->getState().size());
        copy->dfunc->getState(true) = dfunc->getState() ;
        copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
        copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
        copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
    }
    copy->setBlocks(blocks) ;
    return copy ;
}



void AgeingLogarithmicCreepWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    if(!noFracture)
    {
        dfunc->step(currentState, maxscore) ;
        currentState.getParent()->behaviourUpdated = dfunc->changed() ;
        currentState.getParent()->needAssembly = currentState.getParent()->behaviourUpdated ;

            etaSet = false ;
            paramSet = false ;
            Point p( 0,0,0, 1.-2.*maxscore ) ;
            currentParam = getTensor(p) ;
            currentEta = getViscousTensor(p) ;
            etaSet = true ;
            paramSet = true ;

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

void AgeingLogarithmicCreepWithImposedDeformationAndFracture::setFractureCriterion(FractureCriterion * frac)
{
    if(frac)
    {
        criterion = frac ;
    }
}

Vector AgeingLogarithmicCreepWithImposedDeformationAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    Vector epsilon = AgeingLogarithmicCreepWithImposedDeformation::getImposedStrain(p,e,g) ;
    if( dfunc && dfunc->hasInducedForces() ) { epsilon += dfunc->getImposedStrain(p) ; }
    return epsilon ;
}

Vector AgeingLogarithmicCreepWithImposedDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
    Vector ret = AgeingLogarithmicCreepWithImposedDeformation::getImposedStrain(p,e,g) ;
    if(dfunc)
    {
        ret = dfunc->apply(C)*ret ;
        if(dfunc->hasInducedForces())
            ret += dfunc->getImposedStress(p) ;
    }
    else
        ret = C*ret ;

    return ret ;
}

Matrix AgeingLogarithmicCreep::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
}

Matrix AgeingLogarithmicCreep::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
}

Matrix AgeingLogarithmicCreepWithImposedDeformation::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
}

Matrix AgeingLogarithmicCreepWithImposedDeformation::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
}

Matrix AgeingLogarithmicCreepWithImposedDeformationAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(paramSet) { return currentParam ; }

    if(noFracture)
        return param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;
    Matrix realParam = param*(p.getT()+1.)*0.5 + prevParam*(1.-p.getT())*0.5 ;

    if(dfunc->getState().max() < POINT_TOLERANCE)
        return realParam ;

    if(dfunc->getState().size() == 1)
        return  dfunc->apply(realParam, p,e,g) ;

    Matrix Ctemp( C.numRows(), C.numCols() ) ;
    Matrix Etemp( C.numRows(), C.numCols() ) ;
    Matrix Rtemp( C.numRows(), C.numCols() ) ;
    getBlockInMatrix( realParam, 0,0, Ctemp ) ;
    getBlockInMatrix( realParam, 1,1, Etemp ) ; 
    getBlockInMatrix( realParam, 2,2, Rtemp ) ; 
    Ctemp = dfunc->apply( Ctemp, p,e,g ) ;
    Etemp = dfunc->apply( Etemp, p,e,g ) ;
    Rtemp = dfunc->apply( Rtemp, p,e,g ) ;

    realParam = 0 ;
    placeMatrixInBlock(Ctemp,0,0,realParam) ;
    if(!isPurelyElastic)
    {
        placeMatrixInBlock(-Ctemp,0,1,realParam) ;
        placeMatrixInBlock(-Ctemp,0,2,realParam) ;
        placeMatrixInBlock(-Ctemp,1,0,realParam) ;
        placeMatrixInBlock(-Ctemp,2,0,realParam) ;
        placeMatrixInBlock(Ctemp,1,2,realParam) ;
        placeMatrixInBlock(Ctemp,2,1,realParam) ;
        placeMatrixInBlock(Etemp, 1,1, realParam) ;
        placeMatrixInBlock(Rtemp, 2,2, realParam) ;
    }

    return realParam ;
}

Matrix AgeingLogarithmicCreepWithImposedDeformationAndFracture::getViscousTensor(const Point & p, IntegrableEntity * e, int g) const
{
    if(isPurelyElastic)
        return eta ;

    if(noFracture)
        return eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;
    Matrix realeta = eta*(p.getT()+1.)*0.5 + prevEta*(1.-p.getT())*0.5 ;

    if(dfunc->getState().max() < POINT_TOLERANCE)
        return realeta ;

    if(dfunc->getState().size() == 1)
        return dfunc->applyViscous(realeta, p,e,g) ;

    Matrix Etemp( C.numRows(), C.numCols() ) ;
    Matrix Rtemp( C.numRows(), C.numCols() ) ;
    getBlockInMatrix( eta, 1,1, Etemp ) ; 
    getBlockInMatrix( eta, 2,2, Rtemp ) ; 
    Etemp = dfunc->apply( Etemp, p,e,g ) ;
    Rtemp = dfunc->apply( Rtemp, p,e,g ) ;

    realeta = 0 ;
    placeMatrixInBlock(Etemp, 1,1, realeta) ;
    placeMatrixInBlock(Rtemp, 2,2, realeta) ;

    return realeta ;

}

AgeingLogarithmicCreepWithImposedDeformationAndFracture::~AgeingLogarithmicCreepWithImposedDeformationAndFracture()
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

std::vector<BoundaryCondition * > AgeingLogarithmicCreepWithImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    std::vector<BoundaryCondition * > ret = AgeingLogarithmicCreep::getBoundaryConditions(s, id, p_i, gp, Jinv ) ;
    if((imposed.size() == 0 && !(dfunc && dfunc->hasInducedForces())) || fractured())
        return ret ;
    Vector istress = C*imposed ;
    if(dfunc)
    {
        istress = dfunc->apply(C) * imposed   ;
        if(dfunc->hasInducedForces())
            istress += dfunc->getImposedStrain( gp.gaussPoints[0].first ) ;
    }
    if(v.size() == 3)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
    }
    if(v.size() == 4)
    {
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
        ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
    }
    return ret ;

}

void AgeingLogarithmicCreepWithImposedDeformationAndFracture::preProcess( double timeStep, ElementState & currentState )
{
    AgeingLogarithmicCreepWithImposedDeformation::preProcess(timeStep, currentState) ;
    if(dfunc) { dfunc->prepare() ; }
}

bool AgeingLogarithmicCreepWithImposedDeformationAndFracture::fractured() const
{
    if(noFracture) { return false ; }

    for( size_t i = 0 ; i < dfunc->getState().size() ; i++)
    {
        if(! dfunc->fractured(i)) { return false ; }
    }
    return dfunc->fractured() ;
}
