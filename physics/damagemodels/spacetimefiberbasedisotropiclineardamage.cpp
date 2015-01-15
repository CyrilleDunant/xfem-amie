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
#include "../../polynomial/vm_function_extra.h"

namespace Amie {

SpaceTimeFiberBasedIsotropicLinearDamage::SpaceTimeFiberBasedIsotropicLinearDamage(double f, double t, double density)  : fibreFraction(f), timeTolerance(t), visc("x")
{
	thresholdDamageDensity = density ;
	getState(true).resize(1, 0.);
	isNull = false ;
}

std::pair< Vector, Vector > SpaceTimeFiberBasedIsotropicLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
	return std::make_pair(state, Vector(1., 1)) ;
}

void SpaceTimeFiberBasedIsotropicLinearDamage::computeDelta(const ElementState & s)
{
	delta = 1.-getState()[0] ;
}

Matrix SpaceTimeFiberBasedIsotropicLinearDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{  
	if(fractured())
        return m*1e-6 ;
	
	if(state.size() == 1)
	{
		if(state[0] > 0)
			return m*(1.-VirtualMachine().eval(visc,state[0])) ;
		return m ;
	}
	
	double i = (p.getT() + 1.) * state.size() / 2 ;
	if(i >= state.size())
		i = state.size() - 1 ;
	
	if(state[i] > 0)
		return m*(1.-VirtualMachine().eval(visc,state[i])) ;
	return m ;
	
}

Matrix SpaceTimeFiberBasedIsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

  
	if(fractured())
        return m*1e-6 ;

	
	if(state.size() == 1)
		return m*(1.-state[0]) ;
	
	double i = (p.getT() + 1.) * state.size() / 2 ;
	if(i >= state.size())
		i = state.size() - 1 ;
	
	return m*(1.-state[i]) ;
	
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
	
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() ;//maxscore ;
	double dt = s.getNodalDeltaTime() ;
//	if(! fractured() && score > 0)
//		std::cout << (maxscore - score) << "\t" << timeTolerance*dt << std::endl ;


	if(!fractured() && score >= 1)/* && score == maxScoreInNeighbourhood*/
	{
		double maxScoreInNeighbourhood = s.getParent()->getBehaviour()->getFractureCriterion()->getMaxScoreInNeighbourhood(s) ;
		if(score == maxScoreInNeighbourhood)
		{
			state[state.size()-1] += fibreFraction ;
			change = true ;
			converged = true ;
			s.getParent()->getBehaviour()->getFractureCriterion()->inIteration = true ;
		}
//		std::cout << " before beginning << " << score ;
	}
	else if(!fractured() && score > 0 && (maxscore - score) < timeTolerance*dt)
	{

		state[state.size() -1] += fibreFraction ;//*(1.-(maxscore-score)/maxscore) ;
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

void SpaceTimeFiberBasedIsotropicLinearDamage::setLogitViscousDamageLaw(double a, double b, double c) 
{
	double f1 = log(b/(1.-b)) ;
	Function f="x 1 x - /";
	f /= Function("x") ;
	f = f_log(f) ;
	f = f-f1 ;
	f = f*f ;
	f /= (-c) ;
	f = f_exp(f) ;
	f*= a ;
	f /= Function("1 x -") ;
	f += Function("x") ;

	visc = f ;
}


void SpaceTimeFiberBasedIsotropicLinearDamage::postProcess() 
{
	for(size_t i = 0 ; i < state.size()-1 ; i++)
	{
		state[i] = state[i+1] ;
	}
//	converged = true ;
}

DamageModel * SpaceTimeFiberBasedIsotropicLinearDamage::getCopy() const
{ 
	SpaceTimeFiberBasedIsotropicLinearDamage * dam = new SpaceTimeFiberBasedIsotropicLinearDamage(fibreFraction, timeTolerance, thresholdDamageDensity) ;
	dam->visc = visc ;
	return dam ;
}

} ;

