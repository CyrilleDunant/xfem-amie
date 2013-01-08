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
#include "fiberbasedisotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

FiberBasedIsotropicLinearDamage::FiberBasedIsotropicLinearDamage(double f, double c)  : fibreFraction(f)
{
	thresholdDamageDensity = c ;
	getState(true).resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > FiberBasedIsotropicLinearDamage::computeDamageIncrement( Mu::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void FiberBasedIsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix FiberBasedIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

	if(fractured())
		return m*0.001 ;
	
//	std::cout << getState()[0] << "_" ;
	
	return m*(1.-state[0]) ;
}



bool FiberBasedIsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState().max() >= thresholdDamageDensity ;
}


FiberBasedIsotropicLinearDamage::~FiberBasedIsotropicLinearDamage()
{
}

void FiberBasedIsotropicLinearDamage::step( ElementState &s , double maxscore)
{
	elementState = &s ;
	
	if( fraction < 0 )
	{
		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			fraction = s.getParent()->area() ;
		else
			fraction = s.getParent()->volume();

	}
	
	change = false ;
	if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
	{
		converged = true ;
		return ;
	}
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;//maxscore ;
	double maxScoreInNeighbourhood = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood() ;
	if(!fractured() && score > 0 && score == maxScoreInNeighbourhood)
	{
		state += fibreFraction ;
		for(size_t i = 0 ; i < state.size() ; i++)
		{
			if(state[i] > 1)
				state[i] = 1. ;
		}
		change = true ;
	}
	converged = true ;
	return ;
}


} ;

