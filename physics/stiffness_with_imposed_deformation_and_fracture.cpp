//
// C++ Implementation: stiffness_with_imposed_deformation_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007 ; Alain Giorla <alain.giorla@epfl.ch>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_with_imposed_deformation_and_fracture.h"
#include "../features/boundarycondition.h"

using namespace Mu ;

StiffnessWithImposedDeformationAndFracture::StiffnessWithImposedDeformationAndFracture(const Matrix & rig, Vector imposedDef, FractureCriterion * crit) : LinearForm(rig, false, false, rig.numRows()/3+1), dfunc(rig.numRows()-1, .01), imposed(imposedDef), criterion(crit), eps(0.2)
{
	crit->setNeighbourhoodRadius(eps) ;
	frac = false ;
	init = param[0][0] ;
	change  = false ;
	previouschange = false ;
	previousDamage.resize(dfunc.damageState().size()) ; previousDamage =0 ;
	intermediateDamage.resize(dfunc.damageState().size()) ;intermediateDamage = 0 ;
	count = 0 ;
	previousPreviousDamage.resize(dfunc.damageState().size()) ;previousPreviousDamage = 0 ;
	damage = 0 ;	
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36)
		v.push_back(ZETA);
} ;

StiffnessWithImposedDeformationAndFracture::~StiffnessWithImposedDeformationAndFracture() 
{ 
	delete criterion ;
}

FractureCriterion * StiffnessWithImposedDeformationAndFracture::getFractureCriterion() const
{
	return criterion ;
}

void StiffnessWithImposedDeformationAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc.apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}

void StiffnessWithImposedDeformationAndFracture::stepBack()
{
// 	if(change)
// 	{
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal /= .95 ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->downVal /= .95 ;
// 	}
	change = previouschange ;
	damage.resize(previousDamage.size()) ;
	damage = previousDamage ;
	dfunc.damageState() = damage ;
	frac = dfunc.fractured() ;

	previousDamage.resize(previousPreviousDamage.size()) ;
	previousDamage = previousPreviousDamage ;
}

void StiffnessWithImposedDeformationAndFracture::step(double timestep, ElementState & currentState) 
{
	previouschange = change ;
	change = false ;
	currentState.getParent()->behaviourUpdated = false ;
	if(!frac && criterion->met(currentState) )
	{
		
		dfunc.step(currentState) ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->upVal *= .95 ;
// 		dynamic_cast<MohrCoulomb *>(criterion)->downVal *= .95 ;
		change = true ;
		currentState.getParent()->behaviourUpdated = true ;
		frac = dfunc.fractured() ;
	}
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d = dfunc.damageState() ;
	damage.resize(d.size()) ;
	damage = d ;
}

void StiffnessWithImposedDeformationAndFracture::artificialDamageStep(double d)
{
	previouschange = change ;
	change = false ;

	dfunc.artificialDamageStep(d) ;
	change = true ;
	frac = dfunc.fractured() ;
	previousPreviousDamage.resize(previousDamage.size()) ;
	previousPreviousDamage = previousDamage ;
	previousDamage.resize(damage.size()) ;
	previousDamage = damage ;

	Vector d_ = dfunc.damageState() ;
	damage.resize(d_.size()) ;
	damage = d ;
}


void StiffnessWithImposedDeformationAndFracture::artificialPreviousDamage(Vector previous, Vector previousprevious)
{
	previousDamage.resize(damage.size()) ;
	if(previous.size() < previousDamage.size())
	{
		for(size_t i = 0 ; i < previous.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
		for(size_t j = previous.size() ; j < previousDamage.size() ; j++)
			previousDamage[j] = std::min(damage[j],previous[previous.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousDamage.size() ; i++)
			previousDamage[i] = std::min(damage[i],previous[i]) ;
	}
	previousPreviousDamage.resize(damage.size()) ;
	if(previousprevious.size() < previousPreviousDamage.size())
	{
		for(size_t i = 0 ; i < previousprevious.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
		for(size_t j = previous.size() ; j < previousPreviousDamage.size() ; j++)
			previousPreviousDamage[j] = std::min(previousDamage[j],previousprevious[previousprevious.size() - 1]) ;
	} else {
		for(size_t i = 0 ; i < previousPreviousDamage.size() ; i++)
			previousPreviousDamage[i] = std::min(previousDamage[i],previousprevious[i]) ;
	}
}


bool StiffnessWithImposedDeformationAndFracture::changed() const
{
	return change ;
} 

bool StiffnessWithImposedDeformationAndFracture::fractured() const
{
	return dfunc.fractured() ;
}


Form * StiffnessWithImposedDeformationAndFracture::getCopy() const 
{
	StiffnessWithImposedDeformationAndFracture * copy = new StiffnessWithImposedDeformationAndFracture(param, imposed, criterion->getCopy()) ;
	copy->damage = damage ;
	return copy ;
}


Matrix StiffnessWithImposedDeformationAndFracture::getTensor(const Point & p) const
{
	return dfunc.apply(param) ;
}


Vector StiffnessWithImposedDeformationAndFracture::getImposedStress(const Point & p) const
{
	return (param * imposed) ;
}

std::vector<BoundaryCondition * > StiffnessWithImposedDeformationAndFracture::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
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
