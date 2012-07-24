//
// C++ Interface: stiffness_and_fracture
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "fraction_stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"

using namespace Mu ;


FractionStiffnessAndFracture::FractionStiffnessAndFracture(const Matrix & rig, const Matrix & rig0, double phi, FractureCriterion * crit, double eps) : LinearForm(rig, false, true, rig.numRows()/3+1), /*dfunc(rig.numRows()-1)*/ eps(eps)
{
	dfunc = new FractionLinearDamage(rig0, phi) ;
	criterion = crit ;
	crit->setMaterialCharacteristicRadius(eps) ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

void FractionStiffnessAndFracture::setNeighbourhoodRadius(double d)
{
	criterion->setMaterialCharacteristicRadius(d);
	eps = d ;
}

FractionStiffnessAndFracture::~FractionStiffnessAndFracture() 
{ 
	delete criterion ;
	delete dfunc ;
} ;

FractureCriterion * FractionStiffnessAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * FractionStiffnessAndFracture::getDamageModel() const
{
	return dfunc ;
}

void FractionStiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
}


void FractionStiffnessAndFracture::step(double timestep, ElementState & currentState, double maxscore) 
{

	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;
	
}


bool FractionStiffnessAndFracture::changed() const
{
	 return dfunc->changed() ;
} 

bool FractionStiffnessAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

Form * FractionStiffnessAndFracture::getCopy() const 
{
	FractionStiffnessAndFracture * copy = new FractionStiffnessAndFracture(param, dfunc->remnant, dfunc->phi, criterion->getCopy(), criterion->getMaterialCharacteristicRadius()) ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	copy->dfunc->setSecondaryThresholdDamageDensity(dfunc->getSecondaryThresholdDamageDensity());
	return copy ;
}

Matrix FractionStiffnessAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param) ;
}

Material FractionStiffnessAndFracture::toMaterial()
{
	Material mat(getTensor(Point(0,0))) ;
	mat.setProperties(criterion->toMaterial()) ;
	return mat ;
}

