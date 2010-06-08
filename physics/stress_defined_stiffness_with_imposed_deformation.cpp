
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "stress_defined_stiffness_with_imposed_deformation.h"
#include "non_linear_stiffness.h"
#include "../mesher/delaunay.h" 
#include "../polynomial/vm_base.h" 
#include "../features/boundarycondition.h"

using namespace Mu ;


StressDefinedStiffnessWithImposedDeformation::StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defIso, IntegrableEntity * parent) : StressDefinedStiffness(f,n,parent), imposedX(defIso), imposedY(defIso), imposedZ(defIso)
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		imposedZ = Function("0") ;

}

StressDefinedStiffnessWithImposedDeformation::StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defIso, SpaceDimensionality dim) : StressDefinedStiffness(f,n,dim), imposedX(defIso), imposedY(defIso), imposedZ(defIso)
{
	if(dim == SPACE_TWO_DIMENSIONAL)
		imposedZ = Function("0") ;

}

StressDefinedStiffnessWithImposedDeformation::StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defX, Function defY, Function defZ, IntegrableEntity * parent) : StressDefinedStiffness(f,n,parent), imposedX(defX), imposedY(defY), imposedZ(defZ)
{
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		imposedZ = Function("0") ;

}


StressDefinedStiffnessWithImposedDeformation::StressDefinedStiffnessWithImposedDeformation(Function f, double n, Function defX, Function defY, Function defZ, SpaceDimensionality dim) : StressDefinedStiffness(f,n,dim), imposedX(defX), imposedY(defY), imposedZ(defZ)
{
	if(dim == SPACE_TWO_DIMENSIONAL)
		imposedZ = Function("0") ;

}


StressDefinedStiffnessWithImposedDeformation::~StressDefinedStiffnessWithImposedDeformation()
{

}


void StressDefinedStiffnessWithImposedDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
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

std::vector<BoundaryCondition * > StressDefinedStiffnessWithImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	#warning implementation wrong
	Vector stress = s.getStress(gp.gaussPoints) ; 
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(num_dof == 3)
		v.push_back(ZETA) ;
	Vector f =  VirtualMachine().ieval(Gradient(p_i, true)*stress, gp, Jinv,v) ;
	
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

bool StressDefinedStiffnessWithImposedDeformation::isActive() const 
{

	GaussPointArray gp = parent->getGaussPoints() ;
	
	double E_ = 0;
	Vector def(3) ;
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Vector stress = parent->getState().getStress(gp.gaussPoints[i].first) ;
		if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		{
			E_ += vm.eval(E, stress[0],stress[1])*gp.gaussPoints[i].second ;
			def[0] += vm.eval(imposedX,stress[0])*gp.gaussPoints[i].second ;
			def[1] += vm.eval(imposedY,stress[1])*gp.gaussPoints[i].second ;
		}
		else
		{
			E_ += vm.eval(E, stress[0],stress[1],stress[2])*gp.gaussPoints[i].second ;
			def[0] += vm.eval(imposedX,stress[0])*gp.gaussPoints[i].second ;
			def[1] += vm.eval(imposedY,stress[1])*gp.gaussPoints[i].second ;
			def[2] += vm.eval(imposedZ,stress[2])*gp.gaussPoints[i].second ;
		}
	}
	
	return (E_ > 1e-6 || (std::max(std::abs(def[0]),std::max(std::abs(def[1]),std::abs(def[2])))) > 1e-6) ;
}

Form * StressDefinedStiffnessWithImposedDeformation::getCopy() const 
{
	return new StressDefinedStiffnessWithImposedDeformation(*this) ;
}

Vector StressDefinedStiffnessWithImposedDeformation::getImposedStress(const Point & p) const
{
	GaussPointArray gp = parent->getGaussPoints() ;
	
	double E_ = 0;
	VirtualMachine vm ;
	
	if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector def(3) ;
		Matrix m0(3,3) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			Vector stress = parent->getState().getStress(gp.gaussPoints[i].first) ;
			E_ += vm.eval(E, stress[0],stress[1])*gp.gaussPoints[i].second ;
			def[0] += vm.eval(imposedX,stress[0])*gp.gaussPoints[i].second ;
			def[1] += vm.eval(imposedY,stress[1])*gp.gaussPoints[i].second ;
		}
		m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
		m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 

		return m0*def ;
	} else {
		Vector def3d(3) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			Vector stress3d = parent->getState().getStress(gp.gaussPoints[i].first) ;
			E_ += vm.eval(E, stress3d[0],stress3d[1],stress3d[2])*gp.gaussPoints[i].second ;
			def3d[0] += vm.eval(imposedX,stress3d[0])*gp.gaussPoints[i].second ;
			def3d[1] += vm.eval(imposedY,stress3d[1])*gp.gaussPoints[i].second ;
			def3d[2] += vm.eval(imposedZ,stress3d[2])*gp.gaussPoints[i].second ;
		}
		Matrix m1(6,6) ;
		m1[0][0] = 1. - nu ; m1[0][1] = nu ; m1[0][2] = nu ;
		m1[1][0] = nu ; m1[1][1] = 1. - nu ; m1[1][2] = nu ;
		m1[2][0] = nu ; m1[2][1] = nu ; m1[2][2] = 1. - nu ;
		m1[3][3] = (0.5 - nu)*.99 ;
		m1[4][4] = (0.5 - nu)*.99 ;
		m1[5][5] = (0.5 - nu)*.99 ;
		m1 *= E_/((1.+nu)*(1.-2.*nu)) ;
		return m1*def3d ;
	}
}




