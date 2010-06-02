//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mohrcoulomb.h"

using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit, double eps) : LinearForm(rig, false, true, rig.numRows()/3+1),dfunc(rig.numRows()-1, .01)/*dfunc(rig.numRows()-1)*/, eps(eps)
{
	criterion = crit ;
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
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

void StiffnessAndFracture::setNeighbourhoodRadius(double d)
{
	criterion->setNeighbourhoodRadius(d);
	eps = d ;
}

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
} ;

FractureCriterion * StiffnessAndFracture::getFractureCriterion() const
{
	return criterion ;
}

Matrix StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{

	VirtualMachine vm ;
	return vm.ieval(Gradient(p_i) * dfunc.apply(param) * Gradient(p_j, true), e,v) ;
}

void StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc.apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}


void StiffnessAndFracture::stepBack()
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

void StiffnessAndFracture::step(double timestep, ElementState & currentState) 
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

void StiffnessAndFracture::artificialDamageStep(double d)
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

void StiffnessAndFracture::artificialPreviousDamage(Vector previous, Vector previousprevious)
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



bool StiffnessAndFracture::changed() const
{
	return change ;
} 

bool StiffnessAndFracture::fractured() const
{
	return dfunc.fractured() ;
}

Form * StiffnessAndFracture::getCopy() const 
{
	StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy(), eps) ;
	copy->damage = damage ;
	copy->dfunc.setCharacteristicRadius(dfunc.getCharacteristicRadius());
	copy->dfunc.setDamageDensityIncrement(dfunc.getDamageDensityIncrement());
	copy->dfunc.setThresholdDamageDensity(dfunc.getThresholdDamageDensity());
	return copy ;
}

void StiffnessAndFracture::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
}

Matrix StiffnessAndFracture::getTensor(const Point & p) const
{
	return dfunc.apply(param) ;
}

Material StiffnessAndFracture::toMaterial()
{
	Material mat(param) ;
	mat.merge(criterion->toMaterial()) ;
	return mat ;
}

