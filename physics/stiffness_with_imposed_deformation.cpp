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
	return new StiffnessWithImposedDeformation(param, imposed) ;
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
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
	
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


