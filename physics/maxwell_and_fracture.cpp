#include "maxwell_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "../elements/integrable_entity.h"
#include "homogenization/composite.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

using namespace Mu ;




GeneralizedIterativeMaxwellAndFracture::GeneralizedIterativeMaxwellAndFracture(const Matrix & r, const std::vector<Matrix> & rig, const std::vector<double> & chartime, FractureCriterion * crit, DamageModel * d) : LinearForm(r, false, false, r.numRows()/3+1), r0(r)
{
	if(!d)
		dfunc = new /*NonLocal*/FiberBasedIsotropicLinearDamage() ;
	else
		dfunc = d ;
	criterion = crit ;  
  
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	imposedStressAtGaussPoints.resize(1, Vector(0.,r0.numRows())) ;
	  
	for(size_t i = 0 ; i < rig.size() ; i++)
		branches.push_back(new IterativeMaxwell(rig[i], chartime[i])) ;
	
	
}


GeneralizedIterativeMaxwellAndFracture::~GeneralizedIterativeMaxwellAndFracture() 
{ 
	for(size_t b = 0 ; b < branches.size() ; b++)
		delete branches[b] ;
	delete criterion ;
	delete dfunc ;  
} ;

void GeneralizedIterativeMaxwellAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * GeneralizedIterativeMaxwellAndFracture::getCopy() const 
{
	std::vector<Matrix> rig ;
	std::vector<double> eta ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		rig.push_back(branches[i]->param) ;
		eta.push_back(branches[i]->chartime) ;
	}
	GeneralizedIterativeMaxwellAndFracture * copy = new GeneralizedIterativeMaxwellAndFracture(r0, rig, eta, criterion->getCopy(), dfunc->getCopy()) ;
	copy->dfunc->getState(true).resize(dfunc->getState().size());
	copy->dfunc->getState(true) = dfunc->getState() ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	return copy ;
}

Vector GeneralizedIterativeMaxwellAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	if(g != -1)
		return imposedStressAtGaussPoints[g] ;
	if(e)
	{
		if(imposedStressAtGaussPoints.size() > 1)
		{
			int ga = Mu::isGaussPoint(p, e) ;
			if( ga != -1)
				return imposedStressAtGaussPoints[ga] ;
		  
			VirtualMachine vm ;
			Vector ret ; ret.resize(imposedStressAtGaussPoints[0].size()) ; ret = 0 ;
			for(size_t i = 0 ; i < e->getShapeFunctions().size() ; i++)
			{
				Function ffi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( imposedStressAtGaussPoints, e) ;
				ret += mi*vm.eval( ffi, p.x, p.y, p.z, p.t ) ;
			}
			for(size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++)
			{
				Function ffi = e->getEnrichmentFunction(i) ;
				Vector mi = vm.ieval( imposedStressAtGaussPoints, e) ;
				ret += mi*vm.eval( ffi, p.x, p.y, p.z, p.t ) ;
			}
			return ret ;
		}
	}
	return imposedStressAtGaussPoints[0] ;  
}

Vector GeneralizedIterativeMaxwellAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = dfunc->apply(param) ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > GeneralizedIterativeMaxwellAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::vector<BoundaryCondition * > ret ;
	Vector istress ;
	for ( size_t i = 0 ; i < s.getParent()->getBoundingPoints().size() ; i++ )
	{
		if ( id == s.getParent()->getBoundingPoint( i ).id )
		{
			  Point p = s.getParent()->inLocalCoordinates(s.getParent()->getBoundingPoint(i)) ;
			  Vector tmp = getImposedStress(p, s.getParent()) ;
			  istress.resize(tmp.size()) ;
			  istress = tmp ;
			  break ;
		}
	}
	
 	if(istress.size() == 3)
 	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
//		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[2]));
	}
	else if(istress.size() == 6)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[3]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[4]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[5]));	  
	}
	return ret ;
 
}

FractureCriterion * GeneralizedIterativeMaxwellAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * GeneralizedIterativeMaxwellAndFracture::getDamageModel() const
{
	return dfunc ;
}

void GeneralizedIterativeMaxwellAndFracture::setFractureCriterion(FractureCriterion * frac) 
{
	if(frac)
	{
		criterion = frac ;
	}	
}

void GeneralizedIterativeMaxwellAndFracture::step(double timestep, ElementState & currentState, double maxscore) 
{
	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
}

bool GeneralizedIterativeMaxwellAndFracture::changed() const
{
	return dfunc->changed() ;
} 

bool GeneralizedIterativeMaxwellAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

Matrix GeneralizedIterativeMaxwellAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param, p, e, g) ;
}


ElementState * GeneralizedIterativeMaxwellAndFracture::createElementState( IntegrableEntity * e) 
{
	setNumberOfGaussPoints( e->getGaussPoints().gaussPoints.size() ) ;
	return new ElementStateWithInternalVariables(e, 1+branches.size(), param.numRows() ) ;
}

void GeneralizedIterativeMaxwellAndFracture::updateElementState(double timestep, ElementState & currentState) const 
{
// 	branches[2]->param = r0 ;
	LinearForm::updateElementState(timestep, currentState) ;
	if(timestep < POINT_TOLERANCE_2D)
	{
		return ;
	}
	
	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector strain_next( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
	for(size_t g = 0 ; g < imposedStressAtGaussPoints.size() ; g++)
	{
		strain_next = 0 ; strain_prev = 0 ; alpha_next = 0 ; alpha_prev = 0 ;
		currentState.getFieldAtGaussPoint( STRAIN_FIELD, g, strain_next) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, strain_prev, 0) ;
		for(size_t i = 0 ; i < branches.size() ; i++)
		{
			alpha_prev = 0 ; 
			currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, alpha_prev, i+1) ;
			
			alpha_next = (strain_next*branches[i]->coeff_unext) ;
			alpha_next += (strain_prev*branches[i]->coeff_uprev) ;
			alpha_next += (alpha_prev*branches[i]->coeff_aprev) ;
			
			
			dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, g, i+1) ;
		} 
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain_next, g, 0) ;
//		std::cout << strain_prev[0] << "\t" << strain_next[0] << "\t" << alpha_prev[0] << "\t" << alpha_next[0] << "\n" ;
	}
}

void GeneralizedIterativeMaxwellAndFracture::preProcess( double timeStep, ElementState & currentState ) 
{
	if(currentState.getParent()->getGaussPoints().gaussPoints.size() != imposedStressAtGaussPoints.size())
	      this->syncNumberOfGaussPoints(currentState) ;
	if(timeStep < POINT_TOLERANCE_2D)
	{
//		std::cout << param[0][0] << " " ;
//		this->getInstantaneousCoefficients() ;
		return ;
	}
	else
	{
		this->getCoefficients(timeStep) ;
	}
	
	for(size_t j = 0 ; j < currentState.getParent()->getGaussPoints().gaussPoints.size() ; j++)
	      this->preProcessAtGaussPoint(timeStep, currentState, j) ;
	
	param = r0 ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		param += (Matrix) (branches[i]->param*(1.-branches[i]->coeff_unext)) ;
		
	currentState.getParent()->behaviourUpdated = true ;	
//	std::cout << param[0][0] << " " ;
}

void GeneralizedIterativeMaxwellAndFracture::syncNumberOfGaussPoints(ElementState & state)
{
	size_t g = state.getParent()->getGaussPoints().gaussPoints.size() ;
	Vector imposed = imposedStressAtGaussPoints[0] ;
	imposedStressAtGaussPoints.resize(g) ;
	for(size_t i = 0 ; i < g ; i++)
	{
		imposedStressAtGaussPoints[i].resize(imposed.size()) ;
		imposedStressAtGaussPoints[i] = imposed ;
	}
	dynamic_cast<ElementStateWithInternalVariables &>(state).setNumberOfGaussPoints(g) ;
	
}

void GeneralizedIterativeMaxwellAndFracture::preProcessAtGaussPoint(double timestep, ElementState & currentState, int j) 
{  
	imposedStressAtGaussPoints[j] = 0. ;

	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		alpha_prev = 0 ; strain_prev = 0 ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, strain_prev, 0) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, alpha_prev, i+1) ;
		strain_prev *= branches[i]->coeff_uprev ;
		alpha_prev *= branches[i]->coeff_aprev ;
		
		imposedStressAtGaussPoints[j] += (Vector)(dfunc->apply(branches[i]->param)*strain_prev) ;
		imposedStressAtGaussPoints[j] += (Vector)(dfunc->apply(branches[i]->param)*alpha_prev) ;
	}
//	std::cout << imposedStressAtGaussPoints[j][0] ;
}

void GeneralizedIterativeMaxwellAndFracture::setNumberOfGaussPoints(size_t n) 
{
	if(n == 1)
		return ;
	imposedStressAtGaussPoints.resize(n) ;
	for(size_t i = 0 ; i < n ; i++)
	{
		imposedStressAtGaussPoints[i].resize(param.numRows()) ;
	}
	for(size_t b = 0 ; b < branches.size() ; b++)
	{
		branches[b]->setNumberOfGaussPoints(n) ;
	}
}

void GeneralizedIterativeMaxwellAndFracture::getCoefficients(double timestep) 
{
	for(size_t b = 0 ; b < branches.size() ; b++)
	{
		branches[b]->getCoefficients(timestep) ;
// 		std::cout << branches[b]->coeff_unext << "\t" << branches[b]->coeff_uprev << "\t" <<  branches[b]->coeff_aprev << "\t" <<   std::endl ;
	}
}

void GeneralizedIterativeMaxwellAndFracture::getInstantaneousCoefficients() 
{
	for(size_t b = 0 ; b < branches.size() ; b++)
	{
		branches[b]->getInstantaneousCoefficients() ;
	}
}










