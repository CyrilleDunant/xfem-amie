#include "generalized_fd_maxwell.h"
#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"

using namespace Mu ;

MaxwellBranch::MaxwellBranch(const Matrix & rig, Vector d, int i , Function a, Function c, int p, double g) : stiffness(rig), decay(d), p(p), gamma(g), index(i), affine(a), constant(c)
{
	reducedStiffnessAtGaussPoints.resize(1, rig) ;
	imposedStressAtGaussPoints.resize(1, Vector(0.,rig.numRows())) ;
	fi.resize(1, Vector(0.,rig.numRows())) ;
	gi.resize(1, Vector(0.,rig.numRows())) ;
	li.resize(1, Vector(0.,rig.numRows())) ;
	pi.resize(1, Vector(0.,rig.numRows())) ;
}
		
void MaxwellBranch::step(double timestep, ElementState & currentState) 
{
	for(size_t j = 0 ; j < currentState.getParent()->getGaussPoints().gaussPoints.size() ; j++)
	      this->stepAtGaussPoint(timestep, currentState, j) ;
  
}

void MaxwellBranch::stepAtGaussPoint(double timestep, ElementState & currentState, int j) 
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
	
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, gp, index*2+1) ; 
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, pp, index*2+2) ; 
	currentState.getFieldAtGaussPoint( INTERNAL_VARIABLE_FIELD, j, prev, 0) ; 
//  	std::cout << gp[0] << "\t" << pp[0] << "\t" << prev[0] << std::endl ;

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
				
//		std::cout << (gp + timestep*(1.-gamma)/p * pp + timestep*gamma/p * pn - gn)[0] << std::endl ;
		
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

void MaxwellBranch::nextDelta(int k, double tau, Vector & delta0, Vector & delta1, Vector & delta2, Vector & delta3, Vector & prev) 
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

Vector MaxwellBranch::updateInternalStrain( size_t g, const Vector & eps) const
{
	return fi[g] * eps + gi[g] ;
}  

Vector MaxwellBranch::updateInternalStrainRate( size_t g, const Vector & eps) const
{
	return li[g] * eps + pi[g] ;
}  

void MaxwellBranch::setNumberOfGaussPoints(size_t n) 
{
	reducedStiffnessAtGaussPoints.resize(n) ;
	imposedStressAtGaussPoints.resize(n) ;
	fi.resize(n) ;
	gi.resize(n) ;
	pi.resize(n) ;
	li.resize(n) ;
	for(size_t i = 0 ; i < n ; i++)
	{
		reducedStiffnessAtGaussPoints[i].resize(decay.size(), decay.size()) ;
		imposedStressAtGaussPoints[i].resize(decay.size()) ;
		fi[i].resize(decay.size()) ;
		gi[i].resize(decay.size()) ;
		pi[i].resize(decay.size()) ;
		li[i].resize(decay.size()) ;
	}
}



GeneralizedFDMaxwell::GeneralizedFDMaxwell(const Matrix & rig, const std::vector<std::pair<Matrix, Vector> > & br , int p, double g) : Stiffness(rig)
{
	Function k("x") ;
	Function affine = k/p;
	Function constant = 1.-k/p ;
	change = false ;

	for(size_t i = 0 ; i < br.size() ; i++)
	{
		branches.push_back(MaxwellBranch(br[i].first, br[i].second, i, affine, constant, p, g)) ;
	}
}

GeneralizedFDMaxwell::GeneralizedFDMaxwell(const Matrix & rig, const std::pair<Matrix, Vector> & br , int p, double g) : Stiffness(rig)
{
	Function k("x") ;
	Function affine = k/p;
	Function constant = 1.-k/p ;
	change = false ;

	branches.push_back(MaxwellBranch(br.first, br.second, 0, affine, constant, p, g)) ;
}

GeneralizedFDMaxwell::~GeneralizedFDMaxwell() { } ;

void GeneralizedFDMaxwell::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	std::vector<Matrix> stiffnessAtGaussPoints = this->makeStiffnessMatrixAtGaussPoints(gp.gaussPoints.size()) ;
	vm->ieval(Gradient(p_i) * stiffnessAtGaussPoints * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

std::vector<Matrix> GeneralizedFDMaxwell::makeStiffnessMatrixAtGaussPoints(size_t p) const
{
	std::vector<Matrix> ret(p, param) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		for(size_t g = 0 ; g < p ; g++)
		{
			ret[g] += branches[i].reducedStiffnessAtGaussPoints[g] ;
		}
	}
	return ret ;
}

std::vector<Vector> GeneralizedFDMaxwell::makeImposedStressAtGaussPoints(size_t p) const
{
	std::vector<Vector> ret(p, Vector(param.numRows())) ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		for(size_t g = 0 ; g < p ; g++)
		{
			ret[g] += branches[i].imposedStressAtGaussPoints[g] ;
		}
	}
	return ret ;
}

Form * GeneralizedFDMaxwell::getCopy() const 
{
	GeneralizedFDMaxwell * ret = new GeneralizedFDMaxwell(param, std::vector<std::pair<Matrix, Vector> >(), branches[0].p, branches[0].gamma) ;
	Function k("x") ;
	Function affine = k/branches[0].p;
	Function constant = 1.-k/branches[0].p ;
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		ret->branches.push_back(MaxwellBranch(branches[i].stiffness, branches[i].decay, i, affine, constant, branches[i].p, branches[i].gamma)) ;
	}
	return ret;// new GeneralizedFDMaxwell(param, branches, p, gamma) ;
}

Vector GeneralizedFDMaxwell::getImposedStress(const Point & p, IntegrableEntity * e) const
{
	std::vector<Vector> imposed = this->makeImposedStressAtGaussPoints(branches[0].imposedStressAtGaussPoints.size()) ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			VirtualMachine vm ;
			Vector ret(0., imposed[0].size()) ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function fi = e->getShapeFunction(i) ;
				Vector mi = vm.ieval( Gradient(fi) * imposed, e, v) ;
				ret += vm.geval( Gradient(fi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			return ret ;
		}
	}
	return imposed[0] ;
}

std::vector<BoundaryCondition * > GeneralizedFDMaxwell::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * this->makeImposedStressAtGaussPoints(gp.gaussPoints.size()), gp, Jinv,v) ;
	
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

void GeneralizedFDMaxwell::step(double timestep, ElementState & currentState)
{
//	std::cout << timestep << std::endl ;
	currentState.getParent()->behaviourUpdated = true ;
	change = true ;
}

Vector GeneralizedFDMaxwell::updateInternalStrain( size_t g, size_t i, Vector & eps) 
{
	return branches[i].updateInternalStrain(g, eps) ;
}  

Vector GeneralizedFDMaxwell::updateInternalStrainRate( size_t g, size_t i, Vector & eps) 
{
	return branches[i].updateInternalStrainRate(g, eps) ;
}  

ElementState * GeneralizedFDMaxwell::createElementState( IntegrableEntity * e) 
{
	for(size_t i = 0 ; i < branches.size() ; i++)
	{
		branches[i].setNumberOfGaussPoints(e->getGaussPoints().gaussPoints.size()) ;
	}
	return new ElementStateWithInternalVariables(e, 1+2*branches.size(), param.numRows()) ;
}

void GeneralizedFDMaxwell::updateElementState(double timestep, ElementState & currentState) const 
{
	LinearForm::updateElementState(timestep, currentState) ;
	size_t gmax = currentState.getParent()->getGaussPoints().gaussPoints.size() ;
	for(size_t g = 0 ; g < gmax ; g++)
	{
		Vector strain( 0., 3+3*(num_dof == 3)) ;
		currentState.getFieldAtGaussPoint( STRAIN_FIELD, g, strain) ;
		dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(strain, g, 0) ;
		for(size_t i = 0 ; i < branches.size() ; i++)
		{
			Vector alpha = branches[i].updateInternalStrain(g, strain) ;
			Vector alphadot = branches[i].updateInternalStrainRate(g, strain) ;
			dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alpha, g, i*2+1) ;
			dynamic_cast<ElementStateWithInternalVariables &>(currentState).setInternalVariableAtGaussPoint(alphadot, g, i*2+2) ;
		}
	}
}

Matrix GeneralizedFDMaxwell::getTensor(const Point & p, IntegrableEntity * e) const 
{
	std::vector<Matrix> stiffnessAtGaussPoints = this->makeStiffnessMatrixAtGaussPoints(branches[0].reducedStiffnessAtGaussPoints.size()) ;
	if(e)
	{
		if(e->getOrder() > LINEAR)
		{
			VirtualMachine vm ;
			Matrix ret = stiffnessAtGaussPoints[0] ;
			ret.array() = 0. ;
			for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
			{
				Function fi = e->getShapeFunction(i) ;
				Matrix mi = vm.ieval( Gradient(fi) * stiffnessAtGaussPoints, e, v) ;
				ret += vm.geval( Gradient(fi), e, v, p.x, p.y, p.z, p.t ) * mi ;
			}
			return ret ;
		}
	}
	return stiffnessAtGaussPoints[0] ;
}

void GeneralizedFDMaxwell::preProcess( double timeStep, ElementState & currentState ) 
{
	if(timeStep < POINT_TOLERANCE_2D)
		return ;
	for(size_t i = 0 ; i < branches.size() ; i++)
		branches[i].step(timeStep, currentState) ;  
}

