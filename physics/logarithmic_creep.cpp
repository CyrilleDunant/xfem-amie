#include "logarithmic_creep.h"
#include "viscoelasticity.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, LogCreepAccumulator * acc) : Viscoelasticity(PURE_ELASTICITY, rig, 1), C(rig), E(rig*0), tau(0), reducedTimeStep(-1.), isPurelyElastic(true), updated(true), accumulator(acc), fixCreepVariable(false), prevParam(param), prevEta(eta), timeDependentIntegration(false)
{

}

LogarithmicCreep::LogarithmicCreep(const Matrix & rig, const Matrix & vs, double t, LogCreepAccumulator * acc) : Viscoelasticity(GENERALIZED_KELVIN_VOIGT, rig, vs, vs*t), C(rig), E(vs), tau(t), reducedTimeStep(-1.), isPurelyElastic(false), updated(true), accumulator(acc), fixCreepVariable(false), prevParam(param), prevEta(eta), timeDependentIntegration(false)
{

}

void LogarithmicCreep::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	if(!timeDependentIntegration)
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

			a = 0. ;
			b = 0. ;
			getBlockInMatrix(param, 1,1, buffer) ;

			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
		
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
	}
	else
	{
		std::vector<Matrix> mat(Jinv.size()) ;
		Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
		Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	
		Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
		getBlockInMatrix(param, 0,0, buffer) ;
		Matrix prevBuffer(prevParam.numRows()/blocks, prevParam.numCols()/blocks) ;
		getBlockInMatrix(prevParam, 0,0, prevBuffer) ;

		Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
		Function f2 = accumulator->getKelvinVoigtPreviousFunction() ;

		for(size_t i = 0; i < Jinv.size() ; i++)
			mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f1,gp.gaussPoints[i].first)) ;

		vm->ieval(GradientDot(p_i) * mat * Gradient(p_j, true),    gp, Jinv,v, a) ;
		vm->ieval(Gradient(p_i)    * mat * GradientDot(p_j, true), gp, Jinv,v, b) ;
		a += b ;
		
		placeMatrixInBlock( a, 0,0, ret ) ;
		if(!isPurelyElastic)
		{
			b = a*(-1.) ;
			placeMatrixInBlock(b,1,0, ret);
			placeMatrixInBlock(b,0,1, ret);

			a = 0. ;
			b = 0. ;
			getBlockInMatrix(param, 1,1, buffer) ;
			getBlockInMatrix(prevParam, 1,1, prevBuffer) ;

			for(size_t i = 0; i < Jinv.size() ; i++)
				mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f2,gp.gaussPoints[i].first)) ;

			vm->ieval(GradientDot(p_i) * mat * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * mat * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
		
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
	}
}

void LogarithmicCreep::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	if(!isPurelyElastic)
	{
		if(!timeDependentIntegration)
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
		else
		{
			std::vector<Matrix> mat(Jinv.size()) ;
			Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
			Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

			Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
			getBlockInMatrix(eta, 1,1, buffer) ;
			Matrix prevBuffer(prevParam.numRows()/blocks, prevParam.numCols()/blocks) ;
			getBlockInMatrix(prevEta, 0,0, prevBuffer) ;

			Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
			Function f2 = accumulator->getKelvinVoigtPreviousFunction() ;

			for(size_t i = 0; i < Jinv.size() ; i++)
				mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f2,gp.gaussPoints[i].first)) ;

			vm->ieval(GradientDot(p_i) * mat * GradientDot(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(GradientDotDot(p_i) * mat * Gradient(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
		
	}
}

Form * LogarithmicCreep::getCopy() const
{
	LogarithmicCreep * copy ;
	if(isPurelyElastic)
		copy = new LogarithmicCreep( C, accumulator->getCopy() ) ;
	else
		copy = new LogarithmicCreep( C, E, tau, accumulator->getCopy() ) ;
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
		if(!isPurelyElastic)
		{
			accumulator->preProcess(timeStep, currentState) ;
			Matrix R = C*-1. ;
			placeMatrixInBlock(R,0,1,param) ;
			placeMatrixInBlock(R,1,0,param) ;
			placeMatrixInBlock(C,1,1,param) ;
			R = E*(accumulator->getKelvinVoigtSpringReduction()) ;
			addMatrixInBlock(R,1,1,param) ;
			R = E*accumulator->getKelvinVoigtDashpotReduction() ;
			placeMatrixInBlock(R,1,1,eta) ;
		}
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


LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Vector & imp, LogCreepAccumulator * acc   ) : LogarithmicCreep(rig, acc), imposed(imp), prevImposed(imp.size())
{

}

LogarithmicCreepWithImposedDeformation::LogarithmicCreepWithImposedDeformation( const Matrix & rig, const Matrix & v, double e, const Vector & imp, LogCreepAccumulator * acc   ) : LogarithmicCreep(rig, v, e, acc), imposed(imp), prevImposed(imp.size())
{

}

Form * LogarithmicCreepWithImposedDeformation::getCopy() const
{
    LogarithmicCreepWithImposedDeformation * copy ;
    if(isPurelyElastic)
        copy = new LogarithmicCreepWithImposedDeformation( C, imposed, accumulator->getCopy() ) ;
    else
        copy = new LogarithmicCreepWithImposedDeformation( C, E, tau, imposed, accumulator->getCopy() ) ;
    return copy ;
}

Vector LogarithmicCreepWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	if(imposed.size())
	{
		if(timeDependentIntegration)
		{
			Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
			Function f2 = accumulator->getKelvinVoigtNextFunction() ;
			Point p_ = p ;
			if(e && g > -1)
				p_ = e->getGaussPoints().gaussPoints[g].first ;
			return prevImposed*VirtualMachine().eval(f1, p_) + imposed*VirtualMachine().eval(f2, p_) ;
		}
        	return imposed ;
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

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double t, const Vector & imp, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, v, t, imp, acc), dfunc(nullptr), criterion(nullptr), noFracture(false)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c)
{

}

LogarithmicCreepWithImposedDeformationAndFracture::LogarithmicCreepWithImposedDeformationAndFracture( const Matrix & rig, const Matrix & v, double t, const Vector & imp, FractureCriterion * c , DamageModel * d, LogCreepAccumulator * acc) : LogarithmicCreepWithImposedDeformation(rig, v, t, imp, acc), dfunc(d), criterion(c), noFracture(!d || !c)
{

}

Form * LogarithmicCreepWithImposedDeformationAndFracture::getCopy() const
{
    LogarithmicCreepWithImposedDeformationAndFracture * copy ;

    if(noFracture)
    {
	    if(isPurelyElastic)
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed, accumulator->getCopy() ) ;
	    else
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, tau, imposed, accumulator->getCopy() ) ;
    } else {

	    if(isPurelyElastic)
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
	    else
		copy = new LogarithmicCreepWithImposedDeformationAndFracture( C, E, tau, imposed, criterion->getCopy(), dfunc->getCopy(), accumulator->getCopy() ) ;
	    copy->dfunc->getState(true).resize(dfunc->getState().size());
	    copy->dfunc->getState(true) = dfunc->getState() ;
	    copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	    copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	    copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
    }
    return copy ;
}

void LogarithmicCreepWithImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	if(!timeDependentIntegration)
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

			a = 0. ;
			b = 0. ;
			getBlockInMatrix(param, 1,1, buffer) ;

			if(!noFracture)
				buffer = dfunc->apply(buffer) ;

			vm->ieval(GradientDot(p_i) * buffer * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * buffer * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
		
			placeMatrixInBlock( a, 1,1, ret ) ;
		}


	}
	else
	{
		std::vector<Matrix> mat(Jinv.size()) ;
		Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
		Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;
	
		Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
		getBlockInMatrix(param, 0,0, buffer) ;
		Matrix prevBuffer(prevParam.numRows()/blocks, prevParam.numCols()/blocks) ;
		getBlockInMatrix(prevParam, 0,0, prevBuffer) ;

		Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
		Function f2 = accumulator->getKelvinVoigtPreviousFunction() ;

		for(size_t i = 0; i < Jinv.size() ; i++)
			mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f1,gp.gaussPoints[i].first)) ;

		if(!noFracture)
		{
			for(size_t i = 0; i < Jinv.size() ; i++)
				mat[i] = dfunc->apply(mat[i]) ;
		}

		vm->ieval(GradientDot(p_i) * mat * Gradient(p_j, true),    gp, Jinv,v, a) ;
		vm->ieval(Gradient(p_i)    * mat * GradientDot(p_j, true), gp, Jinv,v, b) ;
		a += b ;
		
		placeMatrixInBlock( a, 0,0, ret ) ;
		if(!isPurelyElastic)
		{
			b = a*(-1.) ;
			placeMatrixInBlock(b,1,0, ret);
			placeMatrixInBlock(b,0,1, ret);

			a = 0. ;
			b = 0. ;
			getBlockInMatrix(param, 1,1, buffer) ;
			getBlockInMatrix(prevParam, 1,1, prevBuffer) ;

			for(size_t i = 0; i < Jinv.size() ; i++)
				mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f2,gp.gaussPoints[i].first)) ;

			if(!noFracture)
			{
				for(size_t i = 0; i < Jinv.size() ; i++)
					mat[i] = dfunc->apply(mat[i]) ;
			}

			vm->ieval(GradientDot(p_i) * mat * Gradient(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(Gradient(p_i)    * mat * GradientDot(p_j, true), gp, Jinv,v, b) ;
			a += b ;
		
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
	}
}

void LogarithmicCreepWithImposedDeformationAndFracture::applyViscous(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	if(!isPurelyElastic)
	{
		if(!timeDependentIntegration)
		{
			Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
			Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

			Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
			getBlockInMatrix(eta, 1,1, buffer) ;

			if(!noFracture)
				buffer = dfunc->applyViscous(buffer) ;

			vm->ieval(GradientDot(p_i) * buffer * GradientDot(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(GradientDotDot(p_i) * buffer * Gradient(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
		else
		{
			std::vector<Matrix> mat(Jinv.size()) ;
			Matrix a(ret.numRows()/blocks, ret.numCols()/blocks) ;
			Matrix b(ret.numRows()/blocks, ret.numCols()/blocks) ;

			Matrix buffer(param.numRows()/blocks, param.numCols()/blocks) ;
			getBlockInMatrix(eta, 1,1, buffer) ;
			Matrix prevBuffer(prevParam.numRows()/blocks, prevParam.numCols()/blocks) ;
			getBlockInMatrix(prevEta, 0,0, prevBuffer) ;

			Function f1 = accumulator->getKelvinVoigtPreviousFunction() ;
			Function f2 = accumulator->getKelvinVoigtPreviousFunction() ;

			for(size_t i = 0; i < Jinv.size() ; i++)
				mat[i] = prevBuffer*(vm->eval(f1,gp.gaussPoints[i].first)) + buffer*(vm->eval(f2,gp.gaussPoints[i].first)) ;

			if(!noFracture)
			{
				for(size_t i = 0; i < Jinv.size() ; i++)
					mat[i] = dfunc->applyViscous(mat[i]) ;
			}

			vm->ieval(GradientDot(p_i) * mat * GradientDot(p_j, true),    gp, Jinv,v, a) ;
			vm->ieval(GradientDotDot(p_i) * mat * Gradient(p_j, true), gp, Jinv,v, b) ;
			a += b ;
			placeMatrixInBlock( a, 1,1, ret ) ;
		}
		
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
