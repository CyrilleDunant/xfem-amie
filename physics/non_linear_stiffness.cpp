
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "non_linear_stiffness.h"
#include "../mesher/delaunay.h" 
#include "../polynomial/vm_base.h" 

using namespace Mu ;

NonLinearStiffness::NonLinearStiffness(Function f, double n, SpaceDimensionality dim) 
{
	E = f ;
	nu = n ;
	this->time_d = false ;
	this->type = NON_LINEAR ;
	VirtualMachine vm ;
	
	if(dim == SPACE_TWO_DIMENSIONAL)
	{
		FunctionMatrix m0(3,3) ;
		m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = Function() ;
		m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = Function() ; 
		m0[2][0] = Function() ; m0[2][1] = Function() ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
		param = vm.eval(m0, Point(0.333333333, 0.3333333333)) ;

		num_dof = 2 ;
	} else {
		FunctionMatrix m1(6,6) ;
	
		double divNu = 1./((1.+nu)*(1.-2.*nu)) ;

		m1[0][0] = E*(1-nu)*divNu ; m1[0][1] = E*nu*divNu ; m1[0][2] = E*nu*divNu ; m1[0][3] = Function() ; m1[0][4] = Function() ; m1[0][5] = Function() ;
		m1[1][0] = E*nu*divNu ; m1[1][1] = E*(1-nu)*divNu ; m1[1][2] = E*nu*divNu ; m1[1][3] = Function() ; m1[1][4] = Function() ; m1[1][5] = Function() ;
		m1[2][0] = E*nu*divNu ; m1[2][1] = E*nu*divNu ; m1[2][2] = E*(1-nu)*divNu ; m1[2][3] = Function() ; m1[2][4] = Function() ; m1[2][5] = Function() ;
		m1[3][0] = Function() ; m1[3][1] = Function() ; m1[3][2] = Function() ; m1[3][3] = E*(0.5-nu)*0.99*divNu ; m1[3][4] = Function() ; m1[3][5] = Function() ;
		m1[4][0] = Function() ; m1[4][1] = Function() ; m1[4][2] = Function() ; m1[4][3] = Function() ; m1[4][4] = E*(0.5-nu)*0.99*divNu ; m1[4][5] = Function() ;
		m1[5][0] = Function() ; m1[5][1] = Function() ; m1[5][2] = Function() ; m1[5][3] = Function() ; m1[5][4] = Function() ; m1[5][5] = E*(0.5-nu)*0.99*divNu ;

		param = vm.eval(m1, Point(0.3333333333, 0.3333333333, 0.3333333333)) ;

		num_dof = 3 ;
	}
}


NonLinearStiffness::NonLinearStiffness(Function f, double n, IntegrableEntity * p) : parent(p) 
{
	E = f ;
	nu = n ;
	this->time_d = false ;
	this->type = NON_LINEAR ;
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

NonLinearStiffness::~NonLinearStiffness() { } ;

void NonLinearStiffness::setParent(IntegrableEntity * p)
{
	parent = p;
}

Matrix NonLinearStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Point> pts(gp.gaussPoints.size()) ;
	for(size_t i = 0; i < gp.gaussPoints.size() ; i++)
		pts[i] = gp.gaussPoints[i].first ;
	double E_ = 0;

	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector displacements = e->getState().getDisplacements( pts) ;

		Matrix m0(3,3) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			E_ += vm.eval(E, displacements[i*2],displacements[i*2+1])*gp.gaussPoints[i].second ;
		
		m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
		m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 

		std::vector<Variable> v ;
		v.push_back(XI);
		v.push_back(ETA);
		
		return vm.ieval(Gradient(p_i) * m0 * Gradient(p_j, true), e,v) ;
	}

	Vector displacements3d = e->getState().getDisplacements( pts) ;


	Matrix m1(6,6) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		E_ += vm.eval(E, displacements3d[i*3],displacements3d[i*3+1],displacements3d[i*3+2])*gp.gaussPoints[i].second ;
	
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

bool NonLinearStiffness::hasInducedForces() const
{
	return true ;
} ;

bool NonLinearStiffness::hasInducedMatrix() const 
{
	return true ;
}

void NonLinearStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	std::valarray<Point> pts(gp.gaussPoints.size()) ;
	for(size_t i = 0; i < gp.gaussPoints.size() ; i++)
		pts[i] = gp.gaussPoints[i].first ;
	
	double E_ = 0;

	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector displacements = parent->getState().getDisplacements( pts) ;

		Matrix m0(3,3) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			E_ += vm->eval(E, displacements[i*2],displacements[i*2+1])*gp.gaussPoints[i].second ;
		
		m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
		m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 

		std::vector<Variable> v ;
		v.push_back(XI);
		v.push_back(ETA);
		
		vm->ieval(Gradient(p_i) * m0 * Gradient(p_j, true), gp, Jinv,v, ret) ;
	} else {
		Vector displacements3d = parent->getState().getDisplacements( pts) ;


		Matrix m1(6,6) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			E_ += vm->eval(E, displacements3d[i*3],displacements3d[i*3+1],displacements3d[i*3+2])*gp.gaussPoints[i].second ;
	
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

void NonLinearStiffness::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	Vector stress = s.getStress(gp.gaussPoints) ; 

	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(num_dof == 3)
		v.push_back(ZETA) ;
	
	f = VirtualMachine().ieval(Gradient(p_i, true)*stress, gp, Jinv,v) ;
}

bool NonLinearStiffness::isActive() const 
{

	GaussPointArray gp = parent->getGaussPoints() ;
	std::valarray<Point> pts(gp.gaussPoints.size()) ;
	for(size_t i = 0; i < gp.gaussPoints.size() ; i++)
		pts[i] = gp.gaussPoints[i].first ;
	
	Vector displacements = parent->getState().getDisplacements( pts) ;
	double E_ = 0;
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			E_ += vm.eval(E, displacements[i*2],displacements[i*2+1])*gp.gaussPoints[i].second ;
		else
			E_ += vm.eval(E, displacements[i*3],displacements[i*3+1],displacements[i*3+2])*gp.gaussPoints[i].second ;
	}
	
	return E_ > 1e-6 ;


	
}

Form * NonLinearStiffness::getCopy() const 
{
	return new NonLinearStiffness(*this) ;
}


