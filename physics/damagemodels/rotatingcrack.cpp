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
#include "rotatingcrack.h"
#include "damagemodel.h"
#include "../../elements/integrable_entity.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../orthotropicstiffness.h"
#include "../stiffness.h"

namespace Mu
{

RotatingCrack::RotatingCrack( double E, double nu ):  E( E ), nu( nu )
{
	getState( true ).resize( 2, 0. );
	getPreviousState().resize( 2, 0. );
	isNull = false ;
	currentAngle = 0 ;
	currentDamage = 0 ;
	inTension = true ;
	damaging = false ;
	tensionFailure = false ;
	compressionFailure = false ;
	factor = 1 ;
	es = NULL ;
}


double damageAtAngle( const std::vector<std::pair<double, double> > & increments , double angle )
{
	double ret = 0 ;
	double rettest = 0 ;

	for ( size_t i = 0 ; i <  increments.size() ; i++ )
	{
		rettest += increments[i].second ;
		if ( cos( increments[i].first - angle ) > POINT_TOLERANCE_2D )
		{
			double a = cos(  increments[i].first - angle ) * cos(  increments[i].first - angle ) ;
			ret +=  a * increments[i].second ;
		}
	}

	return ret ;
}

int RotatingCrack::getMode() const 
{ 
	if(es && es->getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() &&
		!inTension && es->getParent()->getBehaviour()->getFractureCriterion()->metInTension )
	{
		std::cout << es->getParent()->getBehaviour()->getFractureCriterion()->metInTension << inTension << es->getParent()->getBehaviour()->getFractureCriterion()->metInCompression << std::endl  ;
		return 1 ;
	}
	return -1 ;
}

double RotatingCrack::getAngleShift() const
{
	if(!es)
		return 0 ;
	return std::abs(currentAngle-es->getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle()) ;
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( ElementState &s )
{
	Vector range( 1., 2 ) ;
// 	std::cout << s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle() << std::endl ;
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && 
		   s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
	{
		es = &s ;
		currentAngle = s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle();
	
		if ( s.getParent()->getBehaviour()->getFractureCriterion()->metInTension && !tensionFailure)
		{
			inTension = true ;
			range[1] = getState()[1] ;
		}
		else if(!compressionFailure)
		{
			inTension = false ;
			range[0] = getState()[0] ;
		}
	}

	return std::make_pair( state,  range) ;
}

Matrix RotatingCrack::apply( const Matrix &m ) const
{

	if ( getState().max() < POINT_TOLERANCE_2D)
		return m ;


	return OrthothropicStiffness( factor * E * ( 1. - getState()[0] ), 
																factor * E * ( 1. - getState()[1] ), 
																factor * E * ( 1. -  getState()[1] ) * ( 1. - nu ) * .5, 
																nu * ( 1. -  getState()[0]), 
																currentAngle ).getTensor( Point() ) ;


}


void  RotatingCrack::computeDelta(const ElementState &s)
{
	Vector range( 1., 2 ) ;
		
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->metInTension && !tensionFailure)
	{
		range[1] = getState()[1] ;
	}
	else if(!compressionFailure)
	{
		range[0] = getState()[0] ;
	}
	delta = (range-state).max() ;
}

Matrix RotatingCrack::applyPrevious( const Matrix &m ) const
{

	if ( fractured() )
		return m * 0 ;

	return m * ( 1. - getPreviousState()[0] ) ;
}

bool RotatingCrack::fractured() const
{
//  return false ;
	if ( fraction < 0 )
		return false ;

 return getState()[0] >= thresholdDamageDensity && getState()[1] >= thresholdDamageDensity;
}

void addAndConsolidate( std::vector<std::pair<double, double> > & target, std::vector<double> & weights, double a, double v, double tol = 1e-2 )
{

	for ( size_t i = 0 ; i < target.size() ; i++ )
	{
		if ( std::abs( target[i].first - a ) < tol )
		{
			double newa = ( target[i].first * weights[i] + a ) / ( weights[i] + 1 ) ;
			double newv = target[i].second + v ;
			target[i] = std::make_pair( newa, newv ) ;
			weights[i]++ ;
			return ;
		}
		if(target[i].first > a)
		{
			target.insert(target.begin()+i, std::make_pair( a, v )) ;
			weights.insert(weights.begin()+i, 1 ) ;
			return  ;
		}
	}
	target.push_back( std::make_pair( a, v ) ) ;
	weights.push_back( 1 ) ;
}

void RotatingCrack::postProcess()
{
	if(converged && getState()[0]>= thresholdDamageDensity)
	{
		tensionFailure = true ;
		getState(true)[0] = 1 ;
	}
	if(converged && getState()[1]>= thresholdDamageDensity)
	{
		compressionFailure = true ;
		getState(true)[1] = 1 ;
	}

}

RotatingCrack::~RotatingCrack()
{
}



}
