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
#include "isotropiclineardamage.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

IsotropicLinearDamage::IsotropicLinearDamage()
{
    getState(true).resize(1, 0.);
    isNull = false ;
    alternate = false ;
    newtonIteration = true ;
    damage = 0 ;
    es = nullptr ;
}

std::pair< Vector, Vector > IsotropicLinearDamage::computeDamageIncrement( Amie::ElementState &s)
{
  if(!es)
    es = &s ;
  return std::make_pair(state, Vector(1., 1)) ;
}

void IsotropicLinearDamage::computeDelta(ElementState & s)
{
    delta = 1.-getState()[0]-damage ;
}

Matrix IsotropicLinearDamage::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

    if(fractured())
        return m*residualStiffnessFraction ;

    return m*(1.-getState()[0]-damage) ;
}

Matrix IsotropicLinearDamage::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

    if(fractured())
        return m*residualStiffnessFraction ;

    return m*(1.-getState()[0]-damage) ;
}

bool IsotropicLinearDamage::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return getState()[0]+damage >= thresholdDamageDensity ;
}

void IsotropicLinearDamage::step( ElementState &s , double maxscore)
{
 
  if(!newtonIteration)
    DamageModel::step(s, maxscore) ;
  else
  {
    converged = true ;
    std::pair<double, double> delta = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s , maxscore) ;
    double mindelta = std::max(delta.first,delta.second) ;
    
    computeDamageIncrement(s) ;
    
    if(std::abs(maxscore) < .5*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() && 
      ((s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() && 
        std::abs(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()) < .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance()) || !s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()))
    {
      change = false ;
      return ;
    }
    
    if(s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() > .05*s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() && s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
    {
      
	change = true ;
	damage = std::min(std::max(damage+std::min(mindelta, s.getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState()*(1.-getState()[0]-damage)), 0.), 1.) ;
    }
  }
}

IsotropicLinearDamage::~IsotropicLinearDamage()
{
}

IsotropicLinearDamageRate::IsotropicLinearDamageRate()
{
    getState(true).resize(1, 0.);
    isNull = false ;
}

std::pair< Vector, Vector > IsotropicLinearDamageRate::computeDamageIncrement( Amie::ElementState &s)
{
    return std::make_pair(getState(), Vector(1., 1)) ;
}

void IsotropicLinearDamageRate::computeDelta(ElementState & s)
{
    delta = 1.-getState()[0] ;
}

Matrix IsotropicLinearDamageRate::apply(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

    if(fraction < 0)
        return m ;
    if(fractured())
        return m*residualStiffnessFraction ;

    double ratio = (p.getT()+1)*.5 ;

    return m*(1.-std::min(initalState[0] + ratio*getState()[0], 1.)) ;
}

Matrix IsotropicLinearDamageRate::applyViscous(const Matrix & m, const Point & p,const IntegrableEntity * e, int g) const
{

    if(fractured())
        return m*residualStiffnessFraction ;
    if(fraction < 0)
        return m ;

    double ratio = (p.getT()+1)*.5 ;

    return m*(1.-std::min(initalState[0] + ratio*getState()[0], 1.)) ;
}

bool IsotropicLinearDamageRate::fractured(int direction) const
{
    if(fraction < 0)
        return false ;
    return getState()[0] >= thresholdDamageDensity ;
}


void IsotropicLinearDamage::postProcess()
{
  if(converged && state[0] > 0 ||
      newtonIteration /*&& 
      es->getParent()->getBehaviour()->getFractureCriterion()->getScoreAtState() < .05*es->getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() && 
      es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet()*/
    )
  {
    getState(true)[0] += damage ;
    getState(true)[0] = std::min(getState(true)[0], 1.) ;
    damage = 0 ;
  }
}

IsotropicLinearDamageRate::~IsotropicLinearDamageRate()
{
}

DamageModel * IsotropicLinearDamageRate::getCopy() const
{
    IsotropicLinearDamageRate * ret = new IsotropicLinearDamageRate() ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

DamageModel * IsotropicLinearDamage::getCopy() const
{
    IsotropicLinearDamage * ret = new IsotropicLinearDamage() ;
    ret->copyEssentialParameters( this ) ;
    return ret ;
}

}

