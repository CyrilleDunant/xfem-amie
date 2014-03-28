//
// C++ Implementation: kelvinvoight
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
#include "kelvinvoight.h"
#include "homogenization/composite.h"
#include "../utilities/matrixops.h"
#include "../features/boundarycondition.h"
#include "../elements/generalized_spacetime_viscoelastic_element_state.h"

using namespace Mu ;

KelvinVoight::KelvinVoight(const Matrix & rig, const Matrix & e, double characteristicTime ) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e), characteristicTime(characteristicTime)
{
//	rig.print() ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
} ;

KelvinVoight::~KelvinVoight() { } ;

ElementState * KelvinVoight::createElementState( IntegrableEntity * e) 
{
	return new KelvinVoightSpaceTimeElementState(e) ;  
}

void KelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix temp(ret) ;
	Matrix temp0(ret) ;
	Matrix temp1(ret) ;
	Matrix temp2(ret) ;
	
	vm->ieval(GradientDot(p_i) * eta   * GradientDot(p_j, true), gp, Jinv,v,temp);
	vm->ieval(GradientDotDot(p_i) * eta   * Gradient(p_j, true), gp, Jinv,v,temp2);
	vm->ieval(GradientDot(p_i) * param * Gradient(p_j, true),    gp, Jinv,v,temp0) ;
	vm->ieval(Gradient(p_i)    * param * GradientDot(p_j, true), gp, Jinv,v,temp1);
	ret =(temp0+temp1) + temp+temp2 ;
}

bool KelvinVoight::fractured() const
{
	return false ;
}

bool KelvinVoight::changed() const
{
	return false ;
} 

Form * KelvinVoight::getCopy() const 
{
	KelvinVoight *copy =  new KelvinVoight(*this) ;
	
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

Vector KelvinVoight::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal)
{
	return VirtualMachine().ieval(GradientDot( shape ) * ( data ), gp, Jinv, v) ;
}

Vector KelvinVoight::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e, std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal )
{
	VirtualMachine vm ;
	
	size_t n = e->getBoundingPoints().size() ;
	Vector field(0., n*externaldofs) ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;
	
	std::vector<Vector> g(e->getGaussPoints().gaussPoints.size(), Vector(0., externaldofs)) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;
	
	Vector f = vm.ieval( GradientDot( shape ) * g, e, v) ;
	
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.deval( data, TIME_VARIABLE, e->getBoundingPoint(i) ) ;

	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

	f += vm.ieval( Gradient( shape ) * g, e, v) ;
	
	return f ;
}


/*IncrementalKelvinVoight::IncrementalKelvinVoight(const Matrix &rig, const Matrix &e, double dt) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e), stiff(rig), tau(dt)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() > 9)
	    v.push_back(ZETA);

    Matrix I = Composite::I4(rig) ;
    Matrix K = eta ;
    if(K.size()==36)
	invert6x6Matrix(K) ;
    else
	invert3x3Matrix(K) ;
    K *= stiff ;
    K *= -dt ; // K = - dt*C/E

    Matrix L = exp(K) ;
    if(L.size()==36)
	invert6x6Matrix(L) ;
    else
	invert3x3Matrix(L) ;

    Matrix J = I - L ; // J = (1 - exp( -dt*C/E ) )
    N.resize(J.numRows(), J.numCols()) ;
    N = I - L ;// in fact, N = J

    if(K.size()==36)
	invert6x6Matrix(K) ;
    else
	invert3x3Matrix(K) ;
    // K = -E/(dt*C)
    J *= K ; // J = -E/(dt*C) * (1 - exp( -dt*C/E ) )
    J += I ; // J = (1 -E/(dt*C) * (1 - exp( -dt*C/E ) ) )

    if(J.size()==36)
	invert6x6Matrix(J) ;
    else
	invert3x3Matrix(J) ;

    param *= J ;

    phi.resize(param.numCols()) ; // initialize phi = 0 ;
    up.resize(param.numCols());
}

IncrementalKelvinVoight::~IncrementalKelvinVoight()
{

}

void IncrementalKelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}


Form * IncrementalKelvinVoight::getCopy() const
{
    return new IncrementalKelvinVoight(*this) ;
}

Vector IncrementalKelvinVoight::getImposedStress(const Point &p , IntegrableEntity * e) const

{
    Matrix S = param * N ;
    return S * phi ;
}

std::vector<BoundaryCondition * > IncrementalKelvinVoight::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    Vector phi_n(v.size()) ;
    int n = (id - id%v.size())/v.size() ; // n is point id
    for(int i = 0 ; i < param.numCols() ; i++)
	phi_n[i] = phi[n*param.numCols()+i] ;

    Matrix S = param * N ;

    Vector f = VirtualMachine().ieval(Gradient(p_i) * (S * phi_n), gp, Jinv,v) ;

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

void IncrementalKelvinVoight::step(double timestep, ElementState &s)
{
    PointArray bc = s.getParent()->getBoundingPoints() ;
    Vector sigma(0., bc.size() * (3+3*(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL))) ;
    Vector epsilon(0., sigma.size()) ;
    s.getField( REAL_STRESS_FIELD, STRAIN_FIELD, bc, sigma, epsilon, false) ;

    Matrix compliance = stiff ;
    if(compliance.size()==36)
	invert6x6Matrix(compliance) ;
    else
	invert3x3Matrix(compliance) ;

    for(size_t i = 0 ; i < bc.size() ; i++)
    {
	size_t id = bc[i]->id ;
	if(!up[id])
	{
	    Vector stress(param.numCols()) ;
	    Vector strain(param.numCols()) ;
	    Vector phi_n(param.numCols()) ;

	    for(size_t j = 0 ; j < param.numCols() ; j++)
	    {
		stress[j] = sigma[i*param.numCols() + j] ;
		strain[j] = epsilon[i*param.numCols() + j] ;
		phi_n[j] = phi[id*param.numCols() + j] ;
	    }

	    stress = compliance * stress ;
	    phi_n = phi_n - stress + strain ;

	    up[i] = true ;
	}
    }

    bool done = true ;
    int n = 0 ;
    while(done && n < up.size())
    {
	done = up[n] ;
	n++ ;
    }
    if(done && n==up.size())
	up.resize(up.size());

}

void IncrementalKelvinVoight::resize(size_t num_points)
{
    std::cerr << "----" << std::endl ;
    std::cerr << num_points << std::endl ;
    phi.resize(num_points*param.numCols()) ;
    up.resize(num_points) ;
}*/


NewmarkNumeroffKelvinVoigt::NewmarkNumeroffKelvinVoigt(const Matrix & rig, const Vector & d, const double a) : LinearForm(rig, false, false, rig.numRows()/3+1), stiffness(rig), decay(d), viscosity(rig), alpha(a)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
  
  
	for(size_t i = 0 ; i < d.size() ; i++)
	{
		for(size_t j = 0 ; j < d.size() ; j++)
			viscosity[i][j] *= decay[i] ;
	}
	
	imposedAtGaussPoints.resize(1) ;
	imposedAtGaussPoints[0].resize(d.size()) ;
	imposedAtGaussPoints[0] = 0. ;
	
}
	
NewmarkNumeroffKelvinVoigt::~NewmarkNumeroffKelvinVoigt() { } ;

void NewmarkNumeroffKelvinVoigt::apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const 
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;  
}

Form * NewmarkNumeroffKelvinVoigt::getCopy() const 
{
	return new NewmarkNumeroffKelvinVoigt(stiffness, decay, alpha) ;
}

Vector NewmarkNumeroffKelvinVoigt::getImposedStress( const Point &p , IntegrableEntity * e, int g) const 
{
	if(g > -1)
		return imposedAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			VirtualMachine vm ;
			Matrix ret(imposedAtGaussPoints[0].size(),1) ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function fi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( Gradient(fi) * imposedAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(fi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			Vector vec(0., imposedAtGaussPoints[0].size()) ;
			for(size_t i = 0 ; i < vec.size() ; i++)
				vec[i] = ret[i][0] ;
			return vec ;
		}
	}
	return imposedAtGaussPoints[0] ;
}

Vector NewmarkNumeroffKelvinVoigt::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = this->getTensor(p,e,g) ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > NewmarkNumeroffKelvinVoigt::getBoundaryConditions( const ElementState &s, size_t id, const Function &p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const 
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedAtGaussPoints, gp, Jinv,v) ;
	
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

void NewmarkNumeroffKelvinVoigt::step( double timestep, ElementState &s ) 
{
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	s.getParent()->behaviourUpdated = true ;  
}

ElementState * NewmarkNumeroffKelvinVoigt::createElementState( IntegrableEntity * e) 
{
	if(e->getGaussPoints().gaussPoints.size() > 1)
	{
		imposedAtGaussPoints.resize(e->getGaussPoints().gaussPoints.size()) ;
		for(size_t i = 0 ; i < imposedAtGaussPoints.size() ; i++)
		{
			imposedAtGaussPoints[i].resize(decay.size()) ;
			imposedAtGaussPoints[i] = 0. ;
		}
	}
	return new ElementStateWithInternalVariables(e, 2, decay.size()) ;
}

void NewmarkNumeroffKelvinVoigt::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	for(size_t g = 0 ; g < imposedAtGaussPoints.size() ; g++)
	{
		Vector strainp(0., decay.size()) ;
		Vector speed(0., decay.size()) ;
		Vector strain(0., decay.size()) ;
		currentState.getFieldAtGaussPoint( STRAIN_FIELD , g , strain ) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD , INTERNAL_VARIABLE_FIELD , g , speed , strainp , 0 , 1 ) ;
		
		speed = (strain - strainp) / (alpha*timestep) - speed * (1.-alpha)/alpha ;
		strainp = strain ;//+ (speed * alpha*timestep) ;
		
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(speed, g, 0) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strainp, g, 1) ;		
	}  
}

Matrix NewmarkNumeroffKelvinVoigt::getTensor(const Point & p, IntegrableEntity * e, int g) const 
{
	return param ;
}

void NewmarkNumeroffKelvinVoigt::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
		return ;
	
	param = stiffness + (viscosity / (alpha*timeStep)) ;
	
	for(size_t g = 0 ; g < imposedAtGaussPoints.size() ; g++)
	{
		Vector strainp(0., decay.size()) ;
		Vector speed(0., decay.size()) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD , INTERNAL_VARIABLE_FIELD , g , speed , strainp , 0 , 1 ) ;
		
		speed *= (1.-alpha)/alpha ;
		strainp /= (alpha * timeStep) ;
		
		imposedAtGaussPoints[g] = (Vector) (viscosity * (speed + strainp)) ;
	}
//	std::cout << "before step " << imposedAtGaussPoints[0][0] << std::endl ;
  
}

ExponentiallyPredictedKelvinVoigt::ExponentiallyPredictedKelvinVoigt(const Matrix & rig, const Vector & d) :  LinearForm(rig, false, false, rig.numRows()/3+1), stiffness(rig), decay(d), viscosity(rig), reduction(rig)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
  
  
	for(size_t i = 0 ; i < d.size() ; i++)
	{
		for(size_t j = 0 ; j < d.size() ; j++)
			viscosity[i][j] *= decay[i] ;
	}
	reduction.array() = 0. ;
	
	Vector expl = std::exp(-decay) ;
	for(size_t i = 0 ; i < decay.size() ; i++)
		reduction[i][i] = decay[i] * (1. - expl[i]) / ( 1. - expl[i] - decay[i] ) ;
		
	imposedAtGaussPoints.resize(1) ;
	imposedAtGaussPoints[0].resize(d.size()) ;
	imposedAtGaussPoints[0] = 0. ;
}

ExponentiallyPredictedKelvinVoigt::~ExponentiallyPredictedKelvinVoigt() { } ;

void ExponentiallyPredictedKelvinVoigt::apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const 
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;  
}


Form * ExponentiallyPredictedKelvinVoigt::getCopy() const 
{
	return new ExponentiallyPredictedKelvinVoigt( stiffness, decay) ;
}

Vector ExponentiallyPredictedKelvinVoigt::getImposedStress( const Point &p , IntegrableEntity * e , int g) const 
{
	if(g > -1)
		return imposedAtGaussPoints[g] ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			VirtualMachine vm ;
			Matrix ret(imposedAtGaussPoints[0].size(),1) ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function fi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( Gradient(fi) * imposedAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(fi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			Vector vec(0., imposedAtGaussPoints[0].size()) ;
			for(size_t i = 0 ; i < vec.size() ; i++)
				vec[i] = ret[i][0] ;
			return vec ;
		}
	}
	return imposedAtGaussPoints[0] ;  
}


Vector ExponentiallyPredictedKelvinVoigt::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const 
{
	Vector imposed = this->getImposedStress(p,e,g) ;
	Matrix m = this->getTensor(p,e,g) ;
	Composite::invertTensor(m) ;
	return (Vector) (m*imposed) ;
}

std::vector<BoundaryCondition * > ExponentiallyPredictedKelvinVoigt::getBoundaryConditions( const ElementState &s, size_t id, const Function &p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const 
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * imposedAtGaussPoints, gp, Jinv,v) ;
	
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

void ExponentiallyPredictedKelvinVoigt::step( double timestep, ElementState &s ) 
{
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	s.getParent()->behaviourUpdated = true ;    
}

ElementState * ExponentiallyPredictedKelvinVoigt::createElementState( IntegrableEntity * e) 
{
	if(e->getGaussPoints().gaussPoints.size() > 1)
	{
		imposedAtGaussPoints.resize(e->getGaussPoints().gaussPoints.size()) ;
		for(size_t i = 0 ; i < imposedAtGaussPoints.size() ; i++)
		{
			imposedAtGaussPoints[i].resize(decay.size()) ;
			imposedAtGaussPoints[i] = 0. ;
		}
	}
	return new ElementStateWithInternalVariables(e, 2, decay.size()) ;  
}

void ExponentiallyPredictedKelvinVoigt::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	if(timestep < POINT_TOLERANCE_2D)
		return ;
	for(size_t g = 0 ; g < imposedAtGaussPoints.size() ; g++)
	{
		Vector strain(0., decay.size()) ;
		Vector strainp(0., decay.size()) ;
		Vector straine(0., decay.size()) ;
		Vector strainv(0., decay.size()) ;
		Vector speed(0., decay.size()) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD , INTERNAL_VARIABLE_FIELD , g , strainp , speed , 0 , 1 ) ;
		currentState.getFieldAtGaussPoint( STRAIN_FIELD , g , strain ) ;
		
		strainv = strain - strainp - speed * timestep ;
		strainv /= ( Vector(1., decay.size()) - std::exp(decay) - decay ) ;
		
		straine = speed * timestep - decay * strainv ;
		
		speed = decay * strainv + straine ;
		speed /= timestep ;
		
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain, g, 0) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(speed, g, 1) ;
		
	}    
}

Matrix ExponentiallyPredictedKelvinVoigt::getTensor(const Point & p, IntegrableEntity * e , int g) const 
{
	return param ;
}

void ExponentiallyPredictedKelvinVoigt::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
		return ;

	Vector expl = std::exp(-decay) ;
	for(size_t i = 0 ; i < decay.size() ; i++)
		reduction[i][i] = decay[i] * (1. - expl[i]) / ( 1. - expl[i] - decay[i] ) ;

	Matrix m = viscosity * reduction ;
	m /= timeStep ;
	param = stiffness - m ;
	
	for(size_t g = 0 ; g < imposedAtGaussPoints.size() ; g++)
	{
		Vector strainp(0., decay.size()) ;
		Vector speedp(0., decay.size()) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD , INTERNAL_VARIABLE_FIELD , g , strainp , speedp , 0 , 1 ) ;
		
		Vector imposed = speedp + (Vector) (reduction*speedp) + (Vector) (reduction*strainp) / timeStep ;
		imposedAtGaussPoints[g] = - (Vector) (viscosity * imposed) ;		
	}
}


