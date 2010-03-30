
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "stress_defined_stiffness.h"
#include "physics.h"
#include "../mesher/delaunay.h" 
#include "../polynomial/vm_base.h" 

using namespace Mu ;




StressDefinedStiffness::StressDefinedStiffness(Function f, double n, IntegrableEntity * parent) : NonLinearStiffness(f,n,parent)
{
	VirtualMachine vm ;

	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		FunctionMatrix m0(3,3) ;
	
		m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = Function() ;
		m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = Function() ; 
		m0[2][0] = Function() ; m0[2][1] = Function() ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
		param = vm.eval(m0, Point(0., 0.)) ;
		num_dof = 2 ;
	}

	if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix m1(6,6) ;
	
		double divNu = 1./((1.+nu)*(1.-2.*nu)) ;

		m1[0][0] = E*(1-nu)*divNu ; m1[0][1] = E*nu*divNu ; m1[0][2] = E*nu*divNu ; m1[0][3] = Function() ; m1[0][4] = Function() ; m1[0][5] = Function() ;
		m1[1][0] = E*nu*divNu ; m1[1][1] = E*(1-nu)*divNu ; m1[1][2] = E*nu*divNu ; m1[1][3] = Function() ; m1[1][4] = Function() ; m1[1][5] = Function() ;
		m1[2][0] = E*nu*divNu ; m1[2][1] = E*nu*divNu ; m1[2][2] = E*(1-nu)*divNu ; m1[2][3] = Function() ; m1[2][4] = Function() ; m1[2][5] = Function() ;
		m1[3][0] = Function() ; m1[3][1] = Function() ; m1[3][2] = Function() ; m1[3][3] = E*(0.5-nu)*0.99*divNu ; m1[3][4] = Function() ; m1[3][5] = Function() ;
		m1[4][0] = Function() ; m1[4][1] = Function() ; m1[4][2] = Function() ; m1[4][3] = Function() ; m1[4][4] = E*(0.5-nu)*0.99*divNu ; m1[4][5] = Function() ;
		m1[5][0] = Function() ; m1[5][1] = Function() ; m1[5][2] = Function() ; m1[5][3] = Function() ; m1[5][4] = Function() ; m1[5][5] = E*(0.5-nu)*0.99*divNu ;

		param = vm.eval(m1, Point(0., 0., 0.)) ;

		num_dof = 3 ;
	}

}


StressDefinedStiffness::StressDefinedStiffness(Function f, double n, bool twod) : NonLinearStiffness(f,n)
{
	VirtualMachine vm ;

	if(twod)
	{
		FunctionMatrix m0(3,3) ;
	
		m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = Function() ;
		m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = Function() ; 
		m0[2][0] = Function() ; m0[2][1] = Function() ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
		param = vm.eval(m0, Point(0., 0.)) ;
		num_dof = 2 ;
	}
	else
	{
		FunctionMatrix m1(6,6) ;
	
		double divNu = 1./((1.+nu)*(1.-2.*nu)) ;

		m1[0][0] = E*(1-nu)*divNu ; m1[0][1] = E*nu*divNu ; m1[0][2] = E*nu*divNu ; m1[0][3] = Function() ; m1[0][4] = Function() ; m1[0][5] = Function() ;
		m1[1][0] = E*nu*divNu ; m1[1][1] = E*(1-nu)*divNu ; m1[1][2] = E*nu*divNu ; m1[1][3] = Function() ; m1[1][4] = Function() ; m1[1][5] = Function() ;
		m1[2][0] = E*nu*divNu ; m1[2][1] = E*nu*divNu ; m1[2][2] = E*(1-nu)*divNu ; m1[2][3] = Function() ; m1[2][4] = Function() ; m1[2][5] = Function() ;
		m1[3][0] = Function() ; m1[3][1] = Function() ; m1[3][2] = Function() ; m1[3][3] = E*(0.5-nu)*0.99*divNu ; m1[3][4] = Function() ; m1[3][5] = Function() ;
		m1[4][0] = Function() ; m1[4][1] = Function() ; m1[4][2] = Function() ; m1[4][3] = Function() ; m1[4][4] = E*(0.5-nu)*0.99*divNu ; m1[4][5] = Function() ;
		m1[5][0] = Function() ; m1[5][1] = Function() ; m1[5][2] = Function() ; m1[5][3] = Function() ; m1[5][4] = Function() ; m1[5][5] = E*(0.5-nu)*0.99*divNu ;

		num_dof = 3 ;

		param = vm.eval(m1, Point(0., 0., 0.)) ;

	}

}


StressDefinedStiffness::~StressDefinedStiffness()
{

}

Matrix StressDefinedStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	GaussPointArray gp = e->getGaussPoints() ;
	double E_ = 0;

	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Matrix m0(3,3) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			Vector stress = e->getState().getStress(gp.gaussPoints[i].first,false) ;
			E_ += vm.eval(E, stress[0],stress[1])*gp.gaussPoints[i].second ;
			
		}
		
		m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
		m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 

		std::vector<Variable> v ;
		v.push_back(XI);
		v.push_back(ETA);
		
		return vm.ieval(Gradient(p_i) * m0 * Gradient(p_j, true), e,v) ;
	}
	Matrix m1(6,6) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Vector stress3d = e->getState().getStress(gp.gaussPoints[i].first) ;
		E_ += vm.eval(E, stress3d[0],stress3d[1],stress3d[2])*gp.gaussPoints[i].second ;
	}
	
	m1[0][0] = 1. - nu ; m1[0][1] = nu ; m1[0][2] = nu ;
	m1[1][0] = nu ; m1[1][1] = 1. - nu ; m1[1][2] = nu ;
	m1[2][0] = nu ; m1[2][1] = nu ; m1[2][2] = 1. - nu ;
	m1[3][3] = (0.5 - nu)*.99 ;
	m1[4][4] = (0.5 - nu)*.99 ;
	m1[5][5] = (0.5 - nu)*.99 ;
	m1 *= E_/((1.+nu)*(1.-2.*nu)) ;

	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA) ;	

	return vm.ieval(Gradient(p_i) * m1 * Gradient(p_j, true), e,v) ;

}

void StressDefinedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	
	double E_ = 0;

	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Matrix m0(3,3) ;

		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			Vector stress = parent->getState().getStress(gp.gaussPoints[i].first,false) ;
			E_ += vm->eval(E, stress[0],stress[1])*gp.gaussPoints[i].second ;
			
		}
		
		m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
		m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 

		std::vector<Variable> v ;
		v.push_back(XI);
		v.push_back(ETA);
		
		vm->ieval(Gradient(p_i) * m0 * Gradient(p_j, true), gp, Jinv,v, ret) ;
	} else {
		Matrix m1(6,6) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			Vector stress3d = parent->getState().getStress(gp.gaussPoints[i].first) ;
			E_ += vm->eval(E, stress3d[0],stress3d[1],stress3d[2])*gp.gaussPoints[i].second ;
		}
	
		m1[0][0] = 1. - nu ; m1[0][1] = nu ; m1[0][2] = nu ;
		m1[1][0] = nu ; m1[1][1] = 1. - nu ; m1[1][2] = nu ;
		m1[2][0] = nu ; m1[2][1] = nu ; m1[2][2] = 1. - nu ;
		m1[3][3] = (0.5 - nu)*.99 ;
		m1[4][4] = (0.5 - nu)*.99 ;
		m1[5][5] = (0.5 - nu)*.99 ;
		m1 *= E_/((1.+nu)*(1.-2.*nu)) ;

		std::vector<Variable> v ;
		v.push_back(XI);
		v.push_back(ETA);
		v.push_back(ZETA) ;	

		vm->ieval(Gradient(p_i) * m1 * Gradient(p_j, true), gp, Jinv,v, ret) ;
	}
}

void StressDefinedStiffness::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	Vector stress = s.getStress(gp.gaussPoints) ; 

	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	f = VirtualMachine().ieval(Gradient(p_i, true)*stress, gp, Jinv,v) ;
}

bool StressDefinedStiffness::isActive() const 
{

	GaussPointArray gp = parent->getGaussPoints() ;
	
	double E_ = 0;
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Vector stress = parent->getState().getStress(gp.gaussPoints[i].first) ;
		if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			E_ += vm.eval(E, stress[i*2],stress[i*2+1])*gp.gaussPoints[i].second ;
		else
			E_ += vm.eval(E, stress[i*3],stress[i*3+1],stress[i*3+2])*gp.gaussPoints[i].second ;
	}
	
	return E_ > 1e-6 ;
}

Form * StressDefinedStiffness::getCopy() const 
{
	return new StressDefinedStiffness(*this) ;
}




