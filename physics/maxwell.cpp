#include "maxwell.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "../elements/integrable_entity.h"
#include "homogenization/composite.h"

using namespace Mu ;


NewmarkNumeroffMaxwell::NewmarkNumeroffMaxwell(const Matrix & rig, Vector d, int p, double g) : LinearForm(rig, false, false, rig.numRows()/3+1), stiffness(rig), decay(d), p(p), gamma(g)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	reducedStiffnessAtGaussPoints.resize(1, rig) ;
	imposedStressAtGaussPoints.resize(1, Vector(0.,rig.numRows())) ;
	fi.resize(1, Vector(0.,rig.numRows())) ;
	gi.resize(1, Vector(0.,rig.numRows())) ;
	li.resize(1, Vector(0.,rig.numRows())) ;
	pi.resize(1, Vector(0.,rig.numRows())) ;
	
	Function k("x") ;
	affine = k/p;
	constant = 1.-k/p ;
	
}


NewmarkNumeroffMaxwell::NewmarkNumeroffMaxwell(const Matrix & rig, Vector d, Function a, Function c, int p, double g) : LinearForm(rig, false, false, rig.numRows()/3+1), stiffness(rig), decay(d), p(p), gamma(g), affine(a), constant(c)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	reducedStiffnessAtGaussPoints.resize(1, rig) ;
	imposedStressAtGaussPoints.resize(1, Vector(0.,rig.numRows())) ;
	fi.resize(1, Vector(0.,rig.numRows())) ;
	gi.resize(1, Vector(0.,rig.numRows())) ;
	li.resize(1, Vector(0.,rig.numRows())) ;
	pi.resize(1, Vector(0.,rig.numRows())) ;
}

NewmarkNumeroffMaxwell::~NewmarkNumeroffMaxwell() { } ;

void NewmarkNumeroffMaxwell::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const 
{
	vm->ieval(Gradient(p_i) * reducedStiffnessAtGaussPoints * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * NewmarkNumeroffMaxwell::getCopy() const 
{
	return new NewmarkNumeroffMaxwell(stiffness, decay, affine, constant, p, gamma) ;
}

Vector NewmarkNumeroffMaxwell::getImposedStress(const Point & p, IntegrableEntity * e, int g) const 
{
	if(g > -1)
		return imposedStressAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			int ga = Mu::isGaussPoint(p, e) ;
			if( ga > -1)
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

Vector NewmarkNumeroffMaxwell::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = this->getTensor(p,e,g) ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > NewmarkNumeroffMaxwell::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::cout << "before " << imposedStressAtGaussPoints[0][0] << std::endl ;
	
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedStressAtGaussPoints, gp, Jinv,v) ;
	
	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
 
}

void NewmarkNumeroffMaxwell::step(double timestep, ElementState & currentState) 
{
	currentState.getParent()->behaviourUpdated = true ;  
}

ElementState * NewmarkNumeroffMaxwell::createElementState( IntegrableEntity * e) 
{
	setNumberOfGaussPoints( e->getGaussPoints().gaussPoints.size() ) ;
	return new ElementStateWithInternalVariables(e, 3, param.numRows() ) ;
}

void NewmarkNumeroffMaxwell::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	for(size_t g = 0 ; g < imposedStressAtGaussPoints.size() ; g++)
	{
		Vector strain( 0., 3+3*(num_dof == 3)) ;
		currentState.getFieldAtGaussPoint( STRAIN_FIELD, g, strain) ;
		std::cout << "after " << strain[0] << std::endl ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain, g, 0) ;
		Vector alpha = this->updateInternalStrain(g, strain) ;
		Vector alphadot = this->updateInternalStrainRate(g, strain) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha, g, 1) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot, g, 2) ;
	} 
}

Matrix NewmarkNumeroffMaxwell::getTensor(const Point & p, IntegrableEntity * e, int g) const 
{
	if(g > -1)
		return reducedStiffnessAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			int g = isGaussPoint(p, e) ;
			if( g > -1)
				return reducedStiffnessAtGaussPoints[g] ;
			
			VirtualMachine vm ;
			Matrix ret = reducedStiffnessAtGaussPoints[0] ;
			ret.array() = 0. ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function ffi = e->getShapeFunction(i) ;
				Matrix mi = vm.ieval( Gradient(ffi) * reducedStiffnessAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(ffi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			return ret ;
		}
	}
	return reducedStiffnessAtGaussPoints[0] ;
}

void NewmarkNumeroffMaxwell::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
		return ;
	for(size_t j = 0 ; j < currentState.getParent()->getGaussPoints().gaussPoints.size() ; j++)
	      this->preProcessAtGaussPoint(timeStep, currentState, j) ;
}

void NewmarkNumeroffMaxwell::preProcessAtGaussPoint(double timestep, ElementState & currentState, int j) 
{  
	reducedStiffnessAtGaussPoints[j] = stiffness ;
	Matrix reduction(stiffness.numRows(), stiffness.numRows()) ;
	
	Vector fp(0., decay.size()) ;
	Vector fn(0., decay.size()) ;
	Vector gp(0., decay.size()) ;
	Vector gn(0., decay.size()) ;
	Vector lp(0., decay.size()) ;
	Vector ln(0., decay.size()) ;
	Vector pp(0., decay.size()) ;
	Vector pn(0., decay.size()) ;
	
	Vector delta0(0., decay.size()) ;
	Vector delta1(0., decay.size()) ;
	Vector delta2(0., decay.size()) ;
	Vector delta3(0., decay.size()) ;
	
	Vector lambda0(0., decay.size()) ;
	Vector lambda1(0., decay.size()) ;
	Vector lambda2(0., decay.size()) ;
	Vector lambda3(0., decay.size()) ;
	
	Vector mu0(0., decay.size()) ;
	Vector mu1(0., decay.size()) ;
	Vector mu2(0., decay.size()) ;
	Vector mu3(0., decay.size()) ;
	
	Vector prev(0., decay.size()) ;
	
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, prev, 0) ; 
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, gp, 1) ; 
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, pp, 2) ; 
	
	for(size_t k = 1 ; k < p+1 ; k++)
	{	  
		this->nextDelta(k, timestep/p, delta0, delta1, delta2, delta3, prev) ;
		
		lambda0 = delta0 ;
		lambda1 = 1. + delta1 ;
		lambda2 = timestep/p + delta2 ;
		lambda3 = delta3 ;

		delta0 *= (double) p/(gamma*timestep) ;
		delta1 *= (double) p/(gamma*timestep) ;
		delta2 *= (double) p/(gamma*timestep) ;
		delta3 *= (double) p/(gamma*timestep) ;
		
		mu0 = delta0 ;
		mu1 = delta1 ;
		mu2 = 1. + delta2 ;
		mu3 = delta3 ;
		
		fn = lambda0 + lambda1*fp + lambda2*lp ;
		ln = mu0 + mu1*fp + mu2*lp ;

		gn = lambda3 + lambda1*gp + lambda2*pp ;
		pn = mu3 + mu1*gp + mu2*pp ;
				
		fp = fn ;
		gp = gn ;
		lp = ln ;
		pp = pn ;
		
	}
	

	fi[j] = fn ;
	gi[j] = gn ;
	li[j] = ln ;
	pi[j] = pn ;
	
	for(size_t n = 0 ; n < decay.size() ; n++)
		reduction[n][n] = 1. - fn[n] ;
	
	reducedStiffnessAtGaussPoints[j] = stiffness * reduction ;
	imposedStressAtGaussPoints[j] = (Vector) (stiffness * gn) ;
	
}

void NewmarkNumeroffMaxwell::nextDelta(int k, double tau, Vector & delta0, Vector & delta1, Vector & delta2, Vector & delta3, Vector & prev) 
{
	VirtualMachine vm ;
	Point pk((double) k, 0.) ;
	for(size_t n = 0 ; n < delta0.size() ; n++)
	{
		double z = decay[n] / ( 1./(gamma*tau) + decay[n] ) ;
		delta0[n] = z * vm.eval(affine, pk) ;
		delta1[n] = 0. - z ;
		delta2[n] = 0. - (1. + decay[n]*tau)/( 1./(gamma*tau) + decay[n] ) ;
		delta3[n] = z * vm.eval(constant, pk) * prev[n] ;
	}
}

Vector NewmarkNumeroffMaxwell::updateInternalStrain( size_t g, const Vector & eps) const
{
	return fi[g] * eps + gi[g] ;
}  

Vector NewmarkNumeroffMaxwell::updateInternalStrainRate( size_t g, const Vector & eps) const
{
	return li[g] * eps + pi[g] ;
}  

void NewmarkNumeroffMaxwell::setNumberOfGaussPoints(size_t n) 
{
	if(n == 1)
		return ;
	reducedStiffnessAtGaussPoints.resize(n) ;
	imposedStressAtGaussPoints.resize(n) ;
	fi.resize(n) ;
	gi.resize(n) ;
	pi.resize(n) ;
	li.resize(n) ;
	for(size_t i = 0 ; i < n ; i++)
	{
		reducedStiffnessAtGaussPoints[i].resize(decay.size(), decay.size()) ;
		reducedStiffnessAtGaussPoints[i] = stiffness ;
		imposedStressAtGaussPoints[i].resize(decay.size()) ;
		fi[i].resize(decay.size()) ;
		gi[i].resize(decay.size()) ;
		pi[i].resize(decay.size()) ;
		li[i].resize(decay.size()) ;
	}
}



