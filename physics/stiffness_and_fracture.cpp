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

#include "stiffness_and_fracture.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/ruptureenergy.h"
#include "fracturecriteria/mcft.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "damagemodels/nonlocalisotropiclineardamage.h"
#include "damagemodels/plasticstrain.h"
#include "damagemodels/fiberbasedisotropiclineardamage.h"

using namespace Mu ;


StiffnessAndFracture::StiffnessAndFracture(const Matrix & rig, FractureCriterion * crit, DamageModel * d) : LinearForm(rig, false, true, rig.numRows()/3+1)
{
	if(!d)
		dfunc = new /*NonLocal*/FiberBasedIsotropicLinearDamage() ;
	else
		dfunc = d ;
	criterion = crit ;

	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() == 36 )
	{
		v.push_back(ZETA);
	}
// 	v.push_back(TIME_VARIABLE);
} ;

StiffnessAndFracture::~StiffnessAndFracture() 
{ 
	delete criterion ;
	delete dfunc ;
} ;

FractureCriterion * StiffnessAndFracture::getFractureCriterion() const
{
	return criterion ;
}

DamageModel * StiffnessAndFracture::getDamageModel() const
{
	return dfunc ;
}

void StiffnessAndFracture::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
// 	if(gp.gaussPoints.size() == 1)
		vm->ieval(Gradient(p_i) * dfunc->apply(param) * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	else
// 	{
// 		std::valarray<Matrix> Jinv_i(Jinv[0], 1) ;
// 		GaussPointArray gp_i(gp.gaussPoints[0]) ;
// 		Matrix acc(ret) ;
// 		vm->ieval(Gradient(p_i) * dfunc->apply(param,gp.gaussPoints[0].first, nullptr, 0) * Gradient(p_j, true), gp_i, Jinv_i,v, acc) ;
// 		ret = acc ;
// 		for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
// 		{
// 			Jinv_i = Jinv[i] ;
// 			gp_i.gaussPoints[0] = gp.gaussPoints[i] ;
// 			vm->ieval(Gradient(p_i) * dfunc->apply(param,gp.gaussPoints[i].first, nullptr, i) * Gradient(p_j, true), gp_i, Jinv_i,v, acc) ;
// 			ret += acc ;
// 		}
// 	}
}

void StiffnessAndFracture::step(double timestep, ElementState & currentState, double maxscore) 
{
	dfunc->step(currentState, maxscore) ;
	currentState.getParent()->behaviourUpdated = dfunc->changed() ;

}

bool StiffnessAndFracture::changed() const
{
	return dfunc->changed() ;
} 

bool StiffnessAndFracture::fractured() const
{
	return dfunc->fractured() ;
}

Form * StiffnessAndFracture::getCopy() const 
{
	StiffnessAndFracture * copy = new StiffnessAndFracture(param, criterion->getCopy(), dfunc->getCopy()) ;
	copy->dfunc->getState(true).resize(dfunc->getState().size());
	copy->dfunc->getState(true) = dfunc->getState() ;
	copy->criterion->setMaterialCharacteristicRadius(criterion->getMaterialCharacteristicRadius()) ;
	copy->dfunc->setDamageDensityTolerance(dfunc->getDamageDensityTolerance());
	copy->dfunc->setThresholdDamageDensity(dfunc->getThresholdDamageDensity());
	return copy ;
}

Matrix StiffnessAndFracture::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	return dfunc->apply(param, p, e, g) ;
}

void StiffnessAndFracture::setFractureCriterion(FractureCriterion * frac) 
{
	if(frac)
	{
// //		delete criterion ;
		criterion = frac ;
	}
	
}


