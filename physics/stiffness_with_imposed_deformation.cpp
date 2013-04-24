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

#include "stiffness_with_imposed_deformation.h"
#include "../features/boundarycondition.h"
#include "homogenization/homogenization_base.h"

using namespace Mu ;

StiffnessWithImposedDeformation::StiffnessWithImposedDeformation(const Matrix & rig, Vector imposedDef) : LinearForm(rig, false, false, rig.numRows()/3+1) , imposed(imposedDef)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	this->time_d = false ;
} ;

StiffnessWithImposedDeformation::StiffnessWithImposedDeformation(double E, double nu, double alpha, SpaceDimensionality dim) : LinearForm(Material::cauchyGreen(std::make_pair(E,nu), true,dim), false, false, dim)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
		
	if(dim == SPACE_TWO_DIMENSIONAL)
	{
	imposed.resize(3, 0.) ;
	for(size_t i = 0 ; i < 2 ; i++)
		imposed[i] = alpha ;
	}
	else if(dim == SPACE_THREE_DIMENSIONAL)
	{
		imposed.resize(6, 0.) ;
		for(size_t i = 0 ; i < 3 ; i++)
			imposed[i] = alpha ;
	}
	
	this->time_d = false ;
} ;

StiffnessWithImposedDeformation::~StiffnessWithImposedDeformation() { } ;

void StiffnessWithImposedDeformation::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithImposedDeformation::fractured() const
{
	return false ;
}

Form * StiffnessWithImposedDeformation::getCopy() const 
{
	StiffnessWithImposedDeformation * copy = new StiffnessWithImposedDeformation(param, imposed) ;
	
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

void StiffnessWithImposedDeformation::step(double timestep, ElementState & currentState, double maxscore)
{
}

Vector StiffnessWithImposedDeformation::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	return (param * imposed) ;
}

Vector StiffnessWithImposedDeformation::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	return imposed ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedDeformation::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	
	std::vector<BoundaryCondition * > ret ;
	if(v.size() == 2)
	{
// 		Vector istress = VirtualMachine().ieval(Gradient(p_i)*(param * imposed),gp,Jinv, v)   ;
		Vector istress = param * imposed   ;
		ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_VOLUMIC_STRESS_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[2]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[0]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[1]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementarySurface *>(s.getParent()),gp,Jinv, id, istress[2]));
	}
	if(v.size() == 3)
	{
// 		Vector istress = VirtualMachine().ieval(Gradient(p_i)*(param * imposed),gp,Jinv, v)   ;
		Vector istress = param * imposed ;
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[0]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[2]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[3]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_XI_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[4]));
// 		ret.push_back(new DofDefinedBoundaryCondition(SET_STRESS_ETA_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()),gp,Jinv, id, istress[5]));
	}
	return ret ;
}


