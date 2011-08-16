//
// C++ Implementation: isotropiclineardamage
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "nonlocalisotropiclineardamage.h"
#include "../../elements/integrable_entity.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {


NonLocalIsotropicLinearDamage::NonLocalIsotropicLinearDamage() 
{
	nonLocalDamage = 0 ;
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	thresholdDamageDensity = 1 ;
	isNull = false ;
	es = NULL ;
}

std::pair< Vector, Vector > NonLocalIsotropicLinearDamage::computeDamageIncrement( Mu::ElementState &s)
{
	es = &s ;
	Vector ret(1) ;
	ret[0] = 1 ;
// 	ret[0] = std::min(thresholdDamageDensity/fraction+POINT_TOLERANCE-state[0], state[0]) ;
// 	ret[0] = std::min(.99999, state[0]) ;
// 	ret[0] = std::max(0., state[0]) ;
	return std::make_pair(state, ret) ;

}

void NonLocalIsotropicLinearDamage::artificialDamageStep(double d)
{
	getState(true)[0] = std::min(getState()[0]+d,thresholdDamageDensity/fraction+POINT_TOLERANCE_2D) ;
}


Matrix NonLocalIsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	return ret*std::min(1./*+getState()[0]-2.**/-nonLocalDamage, 1.) ;
}


Matrix NonLocalIsotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	return ret*std::min(1./*+getState()[0]-2.**/-nonLocalDamage, 1.) ;
}

bool NonLocalIsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState()[0] >= thresholdDamageDensity || -getState()[0]+2.*nonLocalDamage >= thresholdDamageDensity;
}

void NonLocalIsotropicLinearDamage::postProcess()
{
	if(es)
	{
		double newnldamage = smoothedState(*es, false)[0] ;
		change = std::abs(newnldamage-nonLocalDamage) > POINT_TOLERANCE_2D || change;
		nonLocalDamage = newnldamage ;
	}
}

}
