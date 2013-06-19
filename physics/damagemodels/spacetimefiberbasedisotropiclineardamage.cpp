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
#include "spacetimefiberbasedisotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Mu {

SpaceTimeFiberBasedIsotropicLinearDamage::SpaceTimeFiberBasedIsotropicLinearDamage(double f, double t, double c)  : fibreFraction(f), timeTolerance(t)
{
	thresholdDamageDensity = c ;
	getState(true).resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > SpaceTimeFiberBasedIsotropicLinearDamage::computeDamageIncrement( Mu::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void SpaceTimeFiberBasedIsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix SpaceTimeFiberBasedIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

  
	if(fractured())
		return m*1e-9 ;
	
	if(state.size() == 1)
		return m*(1.-state[0]) ;
	
	double i = (p.t + 1.) * state.size() / 2 ;
	if(i >= state.size())
		i = state.size() - 1 ;
	
	return m*(1.-state[i]) ;
	
//	std::cout << getState()[0] << "_" ;
	
}



bool SpaceTimeFiberBasedIsotropicLinearDamage::fractured() const 
{
	if(fraction < 0)
		return false ;
	return getState().max() >= thresholdDamageDensity ;
}


SpaceTimeFiberBasedIsotropicLinearDamage::~SpaceTimeFiberBasedIsotropicLinearDamage()
{
}

void SpaceTimeFiberBasedIsotropicLinearDamage::step( ElementState &s , double maxscore)
{
	elementState = &s ;
	
	if( fraction < 0 )
	{
		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			fraction = s.getParent()->area() ;
		else
			fraction = s.getParent()->volume();

		if(state.size() +1 != s.getParent()->timePlanes())
		{
			state.resize(s.getParent()->timePlanes()-1) ;
			s.getParent()->getBehaviour()->setTimeDependent( s.getParent()->timePlanes() > 2) ;
		}
		
	}
	
	
	
	change = false ;
	if(!s.getParent()->getBehaviour()->getFractureCriterion() || maxscore < 0)
	{
		converged = true ;
		return ;
	}
	
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;//maxscore ;
	double maxScoreInNeighbourhood = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood() ;
	double dt = s.getNodalDeltaTime() ;
//	if(! fractured() && score > 0)
//		std::cout << (maxscore - score) << "\t" << timeTolerance*dt << std::endl ;


	if(!fractured() && score >= 1 && score == maxScoreInNeighbourhood)
	{
		state[state.size()-1] += fibreFraction ;
		change = true ;
		converged = true ;
		s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
//		std::cout << " before beginning << " << score ;
	}
	else if(!fractured() && score > 0 && (maxscore - score) < timeTolerance*dt)
	{

		state[state.size() -1] += fibreFraction*(1.-(maxscore-score)/maxscore) ;
		for(size_t i = 0 ; i < state.size() ; i++)
		{
			if(state[i] > 1)
				state[i] = 1. ;
		}
		change = true ;
// 		if((1.-score)*dt < timeTolerance)
// 			std::cout << " too close to beginning << " << dt << " / " << score ;
//		if((score)*dt < timeTolerance)
//			std::cout << " too close to end << " << dt << " / " << score ;
		converged = true ; /*((1.-score)*dt > timeTolerance) && score*dt > timeTolerance ;*/
		s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
	}
	else
		converged = true ;
	return ;
}

void SpaceTimeFiberBasedIsotropicLinearDamage::postProcess() 
{
	for(size_t i = 0 ; i < state.size()-1 ; i++)
	{
		state[i] = state[i+1] ;
	}
//	converged = true ;
}


} ;

