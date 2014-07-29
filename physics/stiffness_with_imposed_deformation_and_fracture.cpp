//
// C++ Implementation: stiffness_with_imposed_deformation_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
// Author : Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_imposed_deformation_and_fracture.h"
#include "../features/boundarycondition.h"

using namespace Amie ;

StiffnessWithImposedDeformationAndFracture::StiffnessWithImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * crit) : LinearForm(rig, false, false, rig.numRows()/3+1), imposed(imposedDef), criterion(crit), eps(0.2)
{
	dfunc = new IsotropicLinearDamage() ;
	crit->setMaterialCharacteristicRadius(eps) ;
	change  = false ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

StiffnessWithImposedDeformationAndFracture::~StiffnessWithImposedDeformationAndFracture() 
{ 
	delete dfunc ;
	delete criterion ;
}

FractureCriterion * StiffnessWithImposedDeformationAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessWithImposedDeformationAndFracture::getDamageModel() const 
{
	return dfunc ;
}

void StiffnessWithImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void StiffnessWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore) 
{

	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
}

bool StiffnessWithImposedDeformationAndFracture::changed() const
{
	return dfunc->changed() ;
} 

bool StiffnessWithImposedDeformationAndFracture::fractured() const
{
	return dfunc->fractured() ;
}


Form * StiffnessWithImposedDeformationAndFracture::getCopy() const 
{
	StiffnessWithImposedDeformationAndFracture * copy = new StiffnessWithImposedDeformationAndFracture(param, imposed, criterion->getCopy()) ;
	copy->damage = damage ;
	
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


Matrix StiffnessWithImposedDeformationAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param) ;
}


Vector StiffnessWithImposedDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	return (param * imposed) ;
}

Vector StiffnessWithImposedDeformationAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	return imposed ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * (param * imposed), gp, Jinv,v) ;
	
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
