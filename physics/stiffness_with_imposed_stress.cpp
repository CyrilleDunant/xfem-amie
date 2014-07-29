//
// C++ Implementation: stiffness_with_imposed_deformation
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_imposed_stress.h"
#include "../features/boundarycondition.h"
#include "homogenization/homogenization_base.h"
#include "homogenization/composite.h"

using namespace Amie ;

StiffnessWithImposedStress::StiffnessWithImposedStress(const Matrix & rig) : LinearForm(rig, false, false, rig.numRows()/3+1)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

	imposed.resize(v.size()) ;
	imposed = 0. ;
} ;

StiffnessWithImposedStress::StiffnessWithImposedStress(const Matrix & rig, Vector imposedStress) : LinearForm(rig, false, false, rig.numRows()/3+1) , imposed(imposedStress)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

StiffnessWithImposedStress::StiffnessWithImposedStress(double E, double nu, double beta, SpaceDimensionality dim) : LinearForm(Material::cauchyGreen(std::make_pair(E,nu), true,dim), false, false, dim)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
		
	if(dim == SPACE_TWO_DIMENSIONAL)
	{
		imposed.resize(3, 0.) ;
		for(size_t i = 0 ; i < 2 ; i++)
			imposed[i] = beta ;
	}
	else if(dim == SPACE_THREE_DIMENSIONAL)
	{
		imposed.resize(6, 0.) ;
		for(size_t i = 0 ; i < 3 ; i++)
			imposed[i] = beta ;
	}
} ;

StiffnessWithImposedStress::~StiffnessWithImposedStress() { } ;

void StiffnessWithImposedStress::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithImposedStress::fractured() const
{
	return false ;
}

Form * StiffnessWithImposedStress::getCopy() const 
{
	
	StiffnessWithImposedStress * copy = new StiffnessWithImposedStress(*this) ;
	
	return copy ; 
}

void StiffnessWithImposedStress::step(double timestep, ElementState & currentState, double maxscore)
{
}

Vector StiffnessWithImposedStress::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	return imposed ;
}

Vector StiffnessWithImposedStress::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	Matrix s = param ;
	Composite::invertTensor(s) ;
	return s*imposed ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedStress::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	
	std::vector<BoundaryCondition * > ret ;
	if(v.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[1]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, imposed[2]));
	}
	if(v.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[2]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[3]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[4]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, imposed[5]));
	}
	return ret ;
}

void StiffnessWithImposedStress::resetImposedStress() 
{
	imposed = 0. ;
}

void StiffnessWithImposedStress::setImposedStress(Vector beta) 
{
	imposed = beta ;
}

void StiffnessWithImposedStress::addImposedStress(Vector beta) 
{
	if(imposed.size() == beta.size())
		imposed += beta ;
}





