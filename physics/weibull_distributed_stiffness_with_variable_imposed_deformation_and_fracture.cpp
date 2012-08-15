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

#include "weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.h"
#include "../features/boundarycondition.h"

using namespace Mu ;

WeibullStiffnessWithVariableImposedDeformationAndFracture::WeibullStiffnessWithVariableImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * c) : LinearForm(rig, true, false, rig.numRows()/3+1) , imposed(imposedDef), dfunc(),criterion(c)
{
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previousDamage = 0 ;
	damage = 0 ;
	variability = .2 ;	
	materialRadius = .001;
	neighbourhoodRadius = .004 ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

WeibullStiffnessWithVariableImposedDeformationAndFracture::~WeibullStiffnessWithVariableImposedDeformationAndFracture() { } ;

void WeibullStiffnessWithVariableImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
	
	
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

bool WeibullStiffnessWithVariableImposedDeformationAndFracture::hasInducedForces() const 
{
	return true ; 
} 

Vector WeibullStiffnessWithVariableImposedDeformationAndFracture::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	return (param * imposed) ;
}

Vector WeibullStiffnessWithVariableImposedDeformationAndFracture::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	return imposed ;
}

void WeibullStiffnessWithVariableImposedDeformationAndFracture::step(double timestep, ElementState & currentState, double maxscore)
{
	change = false ;

	if(!frac  )
	{
		dfunc.step(currentState, maxscore) ;
		previousDamage = damage ;
		
		Vector state = dfunc.getState() ;
		damage = 0 ;
		for(size_t i = 0 ; i < state.size() ; i++)
			damage += state[i] ;
		change = dfunc.changed() ;//std::abs(damage-previousDamage) > 1e-12 ;
		if(damage > .9)
		{
			frac = true ;
			damage = .9999 ;
// 			param[0][1] = 0 ;param[0][1] = 0 ;
// 			param[2][2] *= 0.0001 ;
// 			this->type = VOID_BEHAVIOUR ;
		}
	}

	if(!change && !frac)
	{
		double randomVar = (double)rand()/(double)RAND_MAX ;
		imposed[0] += timestep*randomVar ;
		imposed[1] += timestep*randomVar ;
	}
	
	if(frac)
	{
		imposed = 0 ;
	}
}


Vector WeibullStiffnessWithVariableImposedDeformationAndFracture::getPreviousDamage()
{
	Vector previous(1) ;
	previous[0] = previousDamage ;
	return previous ;
}



bool WeibullStiffnessWithVariableImposedDeformationAndFracture::changed() const
{
	return change ;
} 

bool WeibullStiffnessWithVariableImposedDeformationAndFracture::fractured() const
{
	return frac;
}

Form * WeibullStiffnessWithVariableImposedDeformationAndFracture::getCopy() const 
{
	double randomVar = (double)rand()/(double)RAND_MAX ;
	randomVar = 1.*pow(-log(randomVar),1./2.) ;
	Matrix newTensor(param*(1.-variability)+param*randomVar*variability) ;
	
	StiffnessWithVariableImposedDeformationAndFracture * copy = new StiffnessWithVariableImposedDeformationAndFracture(newTensor, imposed, criterion->getCopy()) ;
	copy->damage = damage ;
	copy->criterion->setMaterialCharacteristicRadius(materialRadius);
	return copy ;
}

std::vector<BoundaryCondition * > WeibullStiffnessWithVariableImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	Vector f = VirtualMachine().ieval(Gradient(p_i) * ( (param * damage) * imposed), gp, Jinv,v) ;
	
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

