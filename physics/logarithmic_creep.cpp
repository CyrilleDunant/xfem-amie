#include "logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

LogarithmicCreep::LogarithmicCreep(const Matrix & rig) : Viscoelasticity(PURE_ELASTICITY, rig, 1), C(rig), E(rig*0), Eprev(rig*0), Edot(rig*0), tau(0), isPurelyElastic(true), accumulatedStress(0), currentStress(0), updated(true), accumulator(LOGCREEP_FORWARD), previousTimeStep(0.), reducedTimeStep(-1.), sign(true), fixCreepVariable(false), currentStressTensor(2,2), accumulatedStressTensor(2,2)
{

}

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, const Matrix & vs, double e) : Viscoelasticity(MAXWELL, rig, vs), C(rig), E(vs), Eprev(vs), Edot(vs*0), tau(e), isPurelyElastic(false), currentStress(0), accumulatedStress(0), updated(true), accumulator(LOGCREEP_FORWARD), previousTimeStep(0.), reducedTimeStep(-1.), sign(true), fixCreepVariable(false), currentStressTensor(2,2), accumulatedStressTensor(2,2)
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

	if(Edot.array().max() > POINT_TOLERANCE_2D)
	{
		b *= 0. ;
		vm->ieval(Gradient(p_i) * Edot * GradientDot(p_j, true), gp, Jinv, v, b) ;
		addMatrixInBlock( b, 1,1, ret ) ;
	}

    }
}

Form * LogarithmicCreep::getCopy() const
{
    LogarithmicCreep * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreep( C ) ;
    else
        copy = new LogarithmicCreep( C, E, tau ) ;
    copy->accumulator = accumulator ;
	return copy ; 
}

void LogarithmicCreep::print() const
{
    std::cout << "I am a logarithmic viscoelastic law" << std::endl ;
}

void LogarithmicCreep::setLogCreepAccumulator( LogCreepStressAccumulator acc)
{
    accumulator = acc ;
}

void LogarithmicCreep::accumulateStress(double timeStep, ElementState & currentState)
{

    if(previousTimeStep < POINT_TOLERANCE_2D)
        previousTimeStep = timeStep ;

    switch(accumulator)
    {
        case LOGCREEP_CONSTANT:
        {
            break ;
        }
	case LOGCREEP_EXPSTRAIN:
	{
		accumulatedStress = 1. ;
		Vector strain(C.numRows()) ; strain = 0 ;
                Vector stress(C.numRows()) ; stress = 0 ;
                currentState.getAverageField( STRAIN_FIELD, strain, nullptr, -1, 1) ;
		accumulatedStress = std::pow(1 + std::max( std::abs(strain[0]), std::abs(strain[1]) ), 4. ) ;
                currentState.getAverageField( REAL_STRESS_FIELD, stress, nullptr, -1, 1) ;
                currentStress = std::max(std::abs(stress[0]), std::abs(stress[1])) ;
		break ;
	}
        case LOGCREEP_FORWARD:
        {
            Vector stress(2+(C.numRows()==6)) ; stress = 0 ;
            if(reducedTimeStep < 0)
            {
                if(!fractured())
                {
                    dynamic_cast<GeneralizedSpaceTimeViscoElasticElementState &>(currentState).getAverageField( PRINCIPAL_REAL_STRESS_FIELD, stress, nullptr, -1, 1.) ;
//		    Vector stress = C*strain ;
//		    std::cout << strain[0] << " " << strain[1] << "     " ;
		    currentStressTensor[0][0] = stress[0] ;
		    currentStressTensor[0][1] = 0. ;//stress[2]/2 ;
		    currentStressTensor[1][0] = 0. ;//stress[2]/2 ;
		    currentStressTensor[1][1] = stress[1] ;
		    double dt = (timeStep+previousTimeStep)/2. ;
		    accumulatedStressTensor[0][0] += currentStressTensor[0][0]*dt ;
		    accumulatedStressTensor[1][0] += currentStressTensor[1][0]*dt ;
		    accumulatedStressTensor[0][1] += currentStressTensor[0][1]*dt ;
		    accumulatedStressTensor[1][1] += currentStressTensor[1][1]*dt ;
                }
                previousTimeStep = timeStep ;
            }
            else
            {
//                accumulatedStress += currentStress*reducedTimeStep ;
            }
            break ;
        }
        case LOGCREEP_PREDICTED:
        {
            Vector buffer0(C.numRows()) ; buffer0 = 0 ;
            Vector buffer1(C.numRows()) ; buffer1 = 0 ;
            currentState.getAverageField( REAL_STRESS_FIELD, buffer0, nullptr, -1, -1) ;
            currentState.getAverageField( REAL_STRESS_FIELD, buffer1, nullptr, -1, 1) ;
            Vector buffer2 = buffer1-buffer0 ;
            buffer2 /= 2.*currentState.getNodalDeltaTime() ;
            buffer1 += timeStep*buffer2 ;
            currentStress =  std::max(std::abs(buffer1[0]), std::abs(buffer1[1])) ;
            std::cout << currentStress << " ";
            if(reducedTimeStep < 0)
            {
                accumulatedStress += currentStress*(timeStep+previousTimeStep)/2 ;
                previousTimeStep = timeStep ;
            }
            else
            {
//                accumulatedStress += currentStress*(reducedTimeStep+previousTimeStep)/2 ;
//                previousTimeStep = reducedTimeStep ;
            }
            break ;
        }
        case LOGCREEP_AGEING:
        {
            break ;
        }
    }

}

void LogarithmicCreep::makeEquivalentViscosity(double timeStep, ElementState & currentState)
{
    updated = false ;
    if(E.array().max() < POINT_TOLERANCE_2D)
        return ;

    getBlockInMatrix(Eprev, 1,1, eta) ;
    Matrix visc = E ;

    switch(accumulator)
    {
        case LOGCREEP_CONSTANT:
        {
            break ;
        }
        case LOGCREEP_EXPSTRAIN:
        {
		if(currentStress < POINT_TOLERANCE_2D)
		    visc.array() = 0 ;
		else
		    visc *= ((accumulatedStress*tau) / (currentStress)) ;
            break ;
        }
        case LOGCREEP_FORWARD:
        {
            if(accumulatedStressTensor.array().max() > POINT_TOLERANCE_2D)
            {
		Vector str = currentStressTensor.array() ;
                if((abs(str).max()) < 1.)
		{
		    fixCreepVariable = true ;
                    visc.array() = 0 ;
		}
                else
		{
		    fixCreepVariable = false ;
		    double n = 0. ;
		    Matrix S(2,2) ;
		    S = inverse2x2Matrix(currentStressTensor) ;
/*		    accumulatedStressTensor.print() ;
		    currentStressTensor.print() ;
		    S.print() ;*/
		    double smax = std::abs(str).max() ;
		    for(size_t i = 0 ; i < 2 ; i++)
		    {
			for(size_t j = 0 ; j < 2 ; j++)
			{
			    if(std::abs(currentStressTensor[i][j]) > smax*1e-4)
				n += accumulatedStressTensor[i][j]*S[j][i] ;
			}
		    }
//		    std::cout << n << std::endl ;
                    visc *= 1+std::abs(n)/tau ;
		}
            }
            break ;
        }
        case LOGCREEP_PREDICTED:
        {
            if(accumulatedStress > POINT_TOLERANCE_2D)
            {
                if(currentStress < POINT_TOLERANCE_2D)
                    visc.array() = 0 ;
                else
                    visc *= (1.+accumulatedStress / (tau*currentStress)) ;
            }
            break ;
        }
        case LOGCREEP_AGEING:
        {
            double t = currentState.getNodalCentralTime() ;
            visc *= (1.+t/tau) ;
            break ;
        }
    }


    placeMatrixInBlock(visc, 1,1, eta);
    Edot = visc-Eprev ;
	if(timeStep > POINT_TOLERANCE_2D)
	    Edot /= (1.*timeStep) ;
	else
	    Edot *= 0. ;
//    isPurelyElastic = (visc.array().max() < POINT_TOLERANCE_2D) ;

    currentState.getParent()->behaviourUpdated = true ;

}

void LogarithmicCreep::preProcess(double timeStep, ElementState & currentState)
{
    if(reducedTimeStep < 0)
    {
        accumulateStress(timeStep, currentState);
        makeEquivalentViscosity(timeStep, currentState);
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
			if(id == tri->getBoundingPoint(i).getId())
			{
				k = i-tri->getBoundingPoints().size()/tri->timePlanes() ;
				break ;
			}
		}
		if(k != 0)
		{
			ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*4+2], 2)) ;
			ret.push_back(new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, s.getDisplacements()[k*4+3], 3)) ;
		}
	}
	return ret ;
}


LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp  ) : LogarithmicCreep(rig), imposed(imp), prevImposed(imp.size())
{

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, double e, const Vector & imp  ) : LogarithmicCreep(rig, v, e), imposed(imp), prevImposed(imp.size())
{

}

Form * LogarithmicCreepWithImposedDeformation::getCopy() const
{
    LogarithmicCreepWithImposedDeformation * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreepWithImposedDeformation( C, imposed ) ;
    else
        copy = new LogarithmicCreepWithImposedDeformation( C, E, tau, imposed ) ;
    copy->accumulator = accumulator ;
    return copy ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
    if(imposed.size())
        return imposed ; //prevImposed+(imposed-prevImposed)*(p.getT()+1)*0.5 ;
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

void LogarithmicCreepWithImposedDeformation::step(double timestep, ElementState & currentState, double maxscore)
{
    prevImposed = imposed ;
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
    copy->accumulator = accumulator ;
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

/*	if(Edot.array().max() > POINT_TOLERANCE_2D)
	{
		b *= 0. ;
		buffer = Edot ;
		if(!noFracture)
		    buffer = dfunc->applyViscous(Edot) ;
		vm->ieval(Gradient(p_i) * buffer * GradientDot(p_j, true), gp, Jinv, v, b) ;
		addMatrixInBlock( b, 1,1, ret ) ;
	}*/

    }
}

void LogarithmicCreepWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
    if(!noFracture)
    {
	    dfunc->step(currentState, maxscore) ;
	    currentState.getParent()->behaviourUpdated = dfunc->changed() ;
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
