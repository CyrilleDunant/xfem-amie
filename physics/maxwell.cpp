#include "maxwell.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "../elements/integrable_entity.h"
#include "homogenization/composite.h"

using namespace Mu ;


IterativeMaxwell::IterativeMaxwell(const Matrix & rig, double e) : LinearForm(rig, false, false, rig.numRows()/3+1)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	chartime = e ;
	
	imposedStressAtGaussPoints.resize(1, Vector(0.,rig.numRows())) ;
	
	coeff_unext = 0. ;
	coeff_uprev = 0. ;
	coeff_aprev = 0. ;
}


IterativeMaxwell::~IterativeMaxwell() { } ;

void IterativeMaxwell::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	vm->ieval(Gradient(p_i) * (param*(1.-coeff_unext)) * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * IterativeMaxwell::getCopy() const 
{
	return new IterativeMaxwell(param, chartime) ;
}

Vector IterativeMaxwell::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	if(g != -1)
		return imposedStressAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			int ga = Mu::isGaussPoint(p, e) ;
			if( ga != -1)
				return imposedStressAtGaussPoints[ga] ;
		  
			VirtualMachine vm ;
			Vector ret(0., imposedStressAtGaussPoints[0].size()) ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function ffi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( Gradient(ffi) * imposedStressAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(ffi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			return ret ;
		}
	}
	return imposedStressAtGaussPoints[0] ;  
}

Vector IterativeMaxwell::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = param*(1-coeff_unext) ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > IterativeMaxwell::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedStressAtGaussPoints, gp, Jinv,v) ;

	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[2]));
	}
	return ret ;
 
}

void IterativeMaxwell::step(double timestep, ElementState & currentState) 
{
	currentState.getParent()->behaviourUpdated = true ;  
}

ElementState * IterativeMaxwell::createElementState( IntegrableEntity * e) 
{
	setNumberOfGaussPoints( e->getGaussPoints().gaussPoints.size() ) ;
	return new ElementStateWithInternalVariables(e, 2, param.numRows() ) ;
}

void IterativeMaxwell::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	
	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector strain_next( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
	for(size_t g = 0 ; g < imposedStressAtGaussPoints.size() ; g++)
	{
		currentState.getFieldAtGaussPoint( STRAIN_FIELD, g, strain_next) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, strain_prev, 0) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, alpha_prev, 1) ;
		
		alpha_next = (strain_next*coeff_unext) ;
		alpha_next += (strain_prev*coeff_uprev) ;
		alpha_next += (alpha_prev*coeff_aprev) ;
		
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain_next, g, 0) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, g, 1) ;
	} 
}

void IterativeMaxwell::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
	{
		this->getInstantaneousCoefficients() ;
	}
	else
	{
		this->getCoefficients(timeStep) ;
	}
	for(size_t j = 0 ; j < currentState.getParent()->getGaussPoints().gaussPoints.size() ; j++)
	      this->preProcessAtGaussPoint(timeStep, currentState, j) ;
	currentState.getParent()->behaviourUpdated = true ;  
}

void IterativeMaxwell::preProcessAtGaussPoint(double timestep, ElementState & currentState, int j) 
{  
	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, strain_prev, 0) ;
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, alpha_prev, 1) ;
	strain_prev *= coeff_uprev ;
	alpha_prev *= coeff_aprev ;
	
	imposedStressAtGaussPoints[j] = 0. ;
	imposedStressAtGaussPoints[j] += (Vector)(param*strain_prev) ;
	imposedStressAtGaussPoints[j] += (Vector)(param*alpha_prev) ;	
}

void IterativeMaxwell::setNumberOfGaussPoints(size_t n) 
{
	if(n == 1)
		return ;
	imposedStressAtGaussPoints.resize(n) ;
	for(size_t i = 0 ; i < n ; i++)
	{
		imposedStressAtGaussPoints[i].resize(param.numRows()) ;
	}
}

void IterativeMaxwell::getCoefficients(double timestep) 
{
	double alpha = exp(-timestep/chartime) ;
	double gamma = 1.-alpha ;
	double lambda = gamma*(chartime/timestep) ;
	coeff_unext = 1.-lambda ;
	coeff_uprev = gamma - coeff_unext ;
	coeff_aprev = alpha ;
}

void IterativeMaxwell::getInstantaneousCoefficients() 
{
	coeff_unext = 0. ;
	coeff_uprev = 0. ;
	coeff_aprev = 1. ;
}











GeneralizedIterativeMaxwell::GeneralizedIterativeMaxwell(const Matrix & r, const std::vector<Matrix> & rig, const std::vector<double> & chartime) : LinearForm(r, false, false, r.numRows()/3+1), r0(r)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	
	imposedStressAtGaussPoints.resize(1, Vector(0.,r0.numRows())) ;
	  
	for(size_t i = 0 ; i < rig.size() ; i++)
		branches.push_back(new IterativeMaxwell(rig[i], chartime[i])) ;
	
	
}


GeneralizedIterativeMaxwell::~GeneralizedIterativeMaxwell() 
{ 
	for(size_t b = 0 ; b < branches.size() ; b++)
		delete branches[b] ;
  
} ;

void GeneralizedIterativeMaxwell::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * GeneralizedIterativeMaxwell::getCopy() const 
{
	std::vector<Matrix> rig ;
	std::vector<double> eta ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		rig.push_back(branches[i]->param) ;
		eta.push_back(branches[i]->chartime) ;
	}
	return new GeneralizedIterativeMaxwell(r0, rig, eta) ;
}

Vector GeneralizedIterativeMaxwell::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	if(g != -1)
		return imposedStressAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			int ga = Mu::isGaussPoint(p, e) ;
			if( ga != -1)
				return imposedStressAtGaussPoints[ga] ;
		  
			VirtualMachine vm ;
			Vector ret(0., imposedStressAtGaussPoints[0].size()) ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function ffi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( Gradient(ffi) * imposedStressAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(ffi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			return ret ;
		}
	}
	return imposedStressAtGaussPoints[0] ;  
}

Vector GeneralizedIterativeMaxwell::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = param ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > GeneralizedIterativeMaxwell::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedStressAtGaussPoints, gp, Jinv,v) ;

	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, f[2]));
	}
	return ret ;
 
}

void GeneralizedIterativeMaxwell::step(double timestep, ElementState & currentState) 
{
	currentState.getParent()->behaviourUpdated = true ;  
}

ElementState * GeneralizedIterativeMaxwell::createElementState( IntegrableEntity * e) 
{
	setNumberOfGaussPoints( e->getGaussPoints().gaussPoints.size() ) ;
	return new ElementStateWithInternalVariables(e, 1+branches.size(), param.numRows() ) ;
}

void GeneralizedIterativeMaxwell::updateElementState(double timestep, ElementState & currentState) const 
{
// 	branches[2]->param = r0 ;
	LinearForm::updateElementState(timestep, currentState) ;
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	
	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector strain_next( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_next( 0., 3+3*(num_dof == 3)) ;
	for(size_t g = 0 ; g < imposedStressAtGaussPoints.size() ; g++)
	{
		currentState.getFieldAtGaussPoint( STRAIN_FIELD, g, strain_next) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, strain_prev, 0) ;
		for(size_t i = 0 ; i < branches.size() ; i++)
		{
			currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, g, alpha_prev, i+1) ;
			
			alpha_next = (strain_next*branches[i]->coeff_unext) ;
			alpha_next += (strain_prev*branches[i]->coeff_uprev) ;
			alpha_next += (alpha_prev*branches[i]->coeff_aprev) ;
			
// 			std::cout << strain_prev[0] << "\t" << strain_next[0] << "\t" << alpha_prev[0] << "\t" << alpha_next[0] << "\n" ;
			
			dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha_next, g, i+1) ;
		} 
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain_next, g, 0) ;
	}
}

void GeneralizedIterativeMaxwell::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
	{
		this->getInstantaneousCoefficients() ;
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
}

void GeneralizedIterativeMaxwell::preProcessAtGaussPoint(double timestep, ElementState & currentState, int j) 
{  
	imposedStressAtGaussPoints[j] = 0. ;

	Vector strain_prev( 0., 3+3*(num_dof == 3)) ;
	Vector alpha_prev( 0., 3+3*(num_dof == 3)) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, strain_prev, 0) ;
		currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, alpha_prev, i+1) ;
		strain_prev *= branches[i]->coeff_uprev ;
		alpha_prev *= branches[i]->coeff_aprev ;
		
		imposedStressAtGaussPoints[j] += (Vector)(branches[i]->param*strain_prev) ;
		imposedStressAtGaussPoints[j] += (Vector)(branches[i]->param*alpha_prev) ;
	}
}

void GeneralizedIterativeMaxwell::setNumberOfGaussPoints(size_t n) 
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

void GeneralizedIterativeMaxwell::getCoefficients(double timestep) 
{
	for(size_t b = 0 ; b < branches.size() ; b++)
	{
		branches[b]->getCoefficients(timestep) ;
// 		std::cout << branches[b]->coeff_unext << "\t" << branches[b]->coeff_uprev << "\t" <<  branches[b]->coeff_aprev << "\t" <<   std::endl ;
	}
}

void GeneralizedIterativeMaxwell::getInstantaneousCoefficients() 
{
	for(size_t b = 0 ; b < branches.size() ; b++)
	{
		branches[b]->getInstantaneousCoefficients() ;
	}
}














Maxwell::Maxwell(const Matrix & rig, const Matrix & e ) : LinearForm(e, false, false, e.numRows()/3+1), decay(e)
{
	Matrix s(e.numRows(), e.numRows()) ;
	s = rig ;
	Composite::invertTensor(s) ;
	decay = decay*s ;
  
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
} ;

Maxwell::~Maxwell() { } ;

ElementState * Maxwell::createElementState( IntegrableEntity * e) 
{
	return new KelvinVoightSpaceTimeElementState(e) ;  
}

void Maxwell::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix temp(ret) ;
	Matrix temp2(ret) ;
	
	vm->ieval(GradientDot(p_i) *  param * GradientDot(p_j, true), gp, Jinv,v,temp);
	vm->ieval(GradientDotDot(p_i) * param  * Gradient(p_j, true), gp, Jinv,v,temp2);
	
	ret = temp+temp2 ;
	
}

bool Maxwell::fractured() const
{
	return false ;
}

bool Maxwell::changed() const
{
	return false ;
} 

Form * Maxwell::getCopy() const 
{
	return new Maxwell(*this) ;
}

Vector Maxwell::getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	Vector imposed = data ;
	return VirtualMachine().ieval(GradientDot( shape ) * ( imposed ), gp, Jinv, v) ;
}

Vector Maxwell::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	VirtualMachine vm ;
	
	size_t n = e->getBoundingPoints().size() ;
	Vector field(0., n*externaldofs) ;
	Vector test(0., externaldofs) ;
	std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;

	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
//	for(size_t i = 0 ; i < g.size() ; i++)
//		g[i] = rig*g[i] ;
	
	Vector f = vm.ieval( GradientDot( shape ) * g, e, v) ;

	field = 0. ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
// 	for(size_t i = 0 ; i < g.size() ; i++)
// 		g[i] = rig*g[i] ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	field = 0. ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.ddeval( data, TIME_VARIABLE, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
 	for(size_t i = 0 ; i < g.size() ; i++)
 		g[i] = decay*g[i] ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	field = 0. ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
 	for(size_t i = 0 ; i < g.size() ; i++)
 		g[i] = decay*g[i] ;
	
	f += vm.ieval( GradientDot( shape ) * g, e, v) ;

	return f ;
}

StandardLinearSolid::StandardLinearSolid(const Matrix & riga, const Matrix & rigb, const Matrix & e ) : LinearForm(riga, false, false, riga.numRows()/3+1), rig0(riga), rig1(rigb), eta1(e), decay(e), etaeq(e)
{
	Matrix s(e.numRows(), e.numRows()) ;
	s = rig1 ;
	Composite::invertTensor(s) ;
	decay = decay*s ;
  
	etaeq += decay * rig0 ;
  
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
	
} ;

StandardLinearSolid::~StandardLinearSolid() { } ;

ElementState * StandardLinearSolid::createElementState( IntegrableEntity * e) 
{
	return new KelvinVoightSpaceTimeElementState(e) ;  
}

void StandardLinearSolid::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix temp(ret) ;
	Matrix temp0(ret) ;
	Matrix temp1(ret) ;
	Matrix temp2(ret) ;
	
	vm->ieval(GradientDot(p_i) * etaeq   * GradientDot(p_j, true), gp, Jinv,v,temp);
	vm->ieval(GradientDotDot(p_i) * etaeq   * Gradient(p_j, true), gp, Jinv,v,temp2);
	vm->ieval(GradientDot(p_i) * param * Gradient(p_j, true),    gp, Jinv,v,temp0) ;
	vm->ieval(Gradient(p_i)    * param * GradientDot(p_j, true), gp, Jinv,v,temp1);
	
	ret = temp+temp2+temp0+temp1 ;
}

bool StandardLinearSolid::fractured() const
{
	return false ;
}

bool StandardLinearSolid::changed() const
{
	return false ;
} 

Form * StandardLinearSolid::getCopy() const 
{
	return new StandardLinearSolid(*this) ;
}

Vector StandardLinearSolid::getForcesFromAppliedStress( Vector & data, Function & shape, const GaussPointArray & gp, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	Vector imposed = data ;
	return VirtualMachine().ieval(GradientDot( shape ) * ( imposed ), gp, Jinv, v) ;
}

Vector StandardLinearSolid::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e, std::valarray<Matrix> & Jinv, std::vector<Variable> & v) 
{
	VirtualMachine vm ;
	
	size_t n = e->getBoundingPoints().size() ;
	Vector field(0., n*externaldofs) ;
	std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;

	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
//	for(size_t i = 0 ; i < g.size() ; i++)
//		g[i] = rig1*g[i] ;
	
	Vector f = vm.ieval( GradientDot( shape ) * g, e, v) ;
	
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
//	for(size_t i = 0 ; i < g.size() ; i++)
//		g[i] = rig1*g[i] ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
	for(size_t i = 0 ; i < g.size() ; i++)
		g[i] = decay*g[i] ;
	
	f += vm.ieval( GradientDot( shape ) * g, e, v) ;

	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.ddeval( data, TIME_VARIABLE, TIME_VARIABLE, e->getBoundingPoint(i) ) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
	for(size_t i = 0 ; i < g.size() ; i++)
		g[i] = decay*g[i] ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	return f ;
}

