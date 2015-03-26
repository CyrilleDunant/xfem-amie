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

#include "stiffness_with_imposed_stress_and_fracture.h"
#include "../features/boundarycondition.h"
#include "homogenization/composite.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"
#include "damagemodels/plasticstrain.h"


using namespace Amie ;

StiffnessWithImposedStressAndFracture::StiffnessWithImposedStressAndFracture(const Matrix & rig, Vector imposedStress, FractureCriterion * c, DamageModel * d) : LinearForm(rig, false, false, rig.numRows()/3+1) , imposed(imposedStress)
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);

	if(!d) { dfunc = new /*NonLocal*/IsotropicLinearDamage() ; } else { dfunc = d ; }
	criterion = c ;

} 

StiffnessWithImposedStressAndFracture::~StiffnessWithImposedStressAndFracture() 
{
	delete criterion ;
	delete dfunc ;
} 

FractureCriterion * StiffnessWithImposedStressAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessWithImposedStressAndFracture::getDamageModel() const
{
	return dfunc ;
}


void StiffnessWithImposedStressAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool StiffnessWithImposedStressAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

bool StiffnessWithImposedStressAndFracture::changed() const
{
	return dfunc->changed() ;
}

Form * StiffnessWithImposedStressAndFracture::getCopy() const 
{
	StiffnessWithImposedStressAndFracture * copy = new StiffnessWithImposedStressAndFracture(param, imposed, criterion->getCopy(), dfunc->getCopy()) ;
	copy->dfunc->getState(true).resize(dfunc->getState().size());
	copy->dfunc->getState(true) = dfunc->getState() ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	
	return copy ;
}

void StiffnessWithImposedStressAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
}

Vector StiffnessWithImposedStressAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	return imposed ;
}

Vector StiffnessWithImposedStressAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	Matrix s = getTensor(p,e,g) ;
	Composite::invertTensor(s) ;
	return s*imposed ;
}

Matrix StiffnessWithImposedStressAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param, p, e, g) ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedStressAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
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

void StiffnessWithImposedStressAndFracture::resetImposedStress() 
{
	imposed = 0. ;
}

void StiffnessWithImposedStressAndFracture::setImposedStress(Vector beta) 
{
	imposed = beta ;
}

void StiffnessWithImposedStressAndFracture::addImposedStress(Vector beta) 
{
	if(imposed.size() == beta.size())
		imposed += beta ;
}


