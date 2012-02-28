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
	auxiliarystate.resize(1, 0.);
	previousauxiliarystate.resize(1, 0.);
	previouspreviousauxiliarystate.resize(1, 0.);
	getState(true).resize(1, 0.);
	getPreviousState().resize(1, 0.);
	previouspreviousstate.resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > NonLocalIsotropicLinearDamage::computeDamageIncrement( Mu::ElementState &s)
{
	return std::make_pair(Vector(std::max(std::min(-state[0]+2.*auxiliarystate[0], 1.), state[0]), 1), Vector(1., 1) ) ;
}

void NonLocalIsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = (Vector(1., 1) - Vector(std::max(std::min(-state[0]+2.*auxiliarystate[0], 1.), state[0]), 1)).max() ;
}

Matrix NonLocalIsotropicLinearDamage::apply(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	
	double omega = 1 - exp(-std::max(std::min((1.-2.)*state[0]+2.*auxiliarystate[0], 1.), state[0])) ;
	return ret*(1.-omega) ;
}

Matrix NonLocalIsotropicLinearDamage::applyPrevious(const Matrix & m) const
{
	Matrix ret(m) ;

	if(fractured())
		return ret*0 ;
	
	double omega = 1 - exp(-std::max(std::min((1.-2.)*previousstate[0]+2.*previousauxiliarystate[0], 1.), previousstate[0])) ;
	return ret*(1.-omega)  ;
}

bool NonLocalIsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return std::max(std::min((1.-2)*state[0]+2.*auxiliarystate[0], 1.), state[0]) >= thresholdDamageDensity;
}

void NonLocalIsotropicLinearDamage::postProcess()
{
	if(elementState)
	{
		double newnldamage = smoothedState(*elementState, false)[0] ;
		change = std::abs(newnldamage-auxiliarystate[0]) > POINT_TOLERANCE_2D || change;
		if(change)
			auxiliarystate[0] = newnldamage ;
	}
}

}
