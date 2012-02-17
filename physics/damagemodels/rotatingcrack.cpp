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

RotatingCrack::RotatingCrack( double E, double nu ): tdamage( 0 ), cdamage( 0 ), E( E ), nu( nu )
{
	getState( true ).resize( 1, 0. );
	getPreviousState().resize( 1, 0. );
	isNull = false ;
	currentAngle = 0 ;
	currentDamage = 0 ;
	inTension = false ;
	damaging = false ;
	factor = 1 ;
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

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( ElementState &s )
{
	
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && 
		   s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
	{
		currentAngle = M_PI*.5 ;s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle();
		damaging = true ;

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->metInTension )
		{
			inTension = true ;

			cdamage = damageAtAngle( compressionAngles, currentAngle ) ;
			tdamage = damageAtAngle( tensionAngles, currentAngle ) ;
			currentDamage = damageAtAngle( tensionAngles, currentAngle ) ;
		}
		else
		{
			inTension = false ;

			tdamage = damageAtAngle( tensionAngles, currentAngle ) ;
			cdamage = damageAtAngle( compressionAngles, currentAngle ) ;
			currentDamage = damageAtAngle( compressionAngles, currentAngle ) ;
		}
		state = currentDamage ;
	}
	return std::make_pair( state, Vector( 1., 1 ) ) ;
}

Matrix RotatingCrack::apply( const Matrix &m ) const
{
	if ( fractured() )
		return m * 0. ;

	if ( getState()[0] < POINT_TOLERANCE_2D
	        && tdamage < POINT_TOLERANCE_2D
	        && cdamage < POINT_TOLERANCE_2D )
		return m ;

// 	if ( damaging )
// 	{
		if ( inTension )
		{
			return OrthothropicStiffness( factor * E * ( 1. - getState()[0] ), 
																		factor * E * ( 1. - cdamage ), 
																		factor * E * ( 1. - 0.5 * ( cdamage + getState()[0] ) ) * ( 1. - nu ) * .5, 
																		nu * ( 1. - std::max( getState()[0], cdamage ) ), 
																		currentAngle ).getTensor( Point() ) ;
		}
		else
		{
			return OrthothropicStiffness( factor * E * ( 1. - tdamage ), 
																		factor * E * ( 1. - getState()[0] ), 
																		factor * E * ( 1. - 0.5 * ( getState()[0] + tdamage ) ) * ( 1. - nu ) * .5, 
																		nu * ( 1. - std::max( getState()[0], tdamage ) ), 
																    currentAngle ).getTensor( Point() ) ;
		}
// 	}
// 
// 	return OrthothropicStiffness( factor * E * ( 1. - tdamage ), factor * E * ( 1. - cdamage ), factor * E * ( 1. - 0.5 * ( cdamage + tdamage ) ) * ( 1. - nu ) * .5, nu*( 1. - std::max( tdamage, cdamage ) ), currentAngle ).getTensor( Point() ) ;
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

	return getState()[0] >= thresholdDamageDensity ;
}

void addAndConsolidate( std::vector<std::pair<double, double> > & target, std::vector<double> & weights, double a, double v, double tol = 1e-2 )
{

// 	for ( size_t i = 0 ; i < target.size() ; i++ )
// 	{
// 		if ( std::abs( target[i].first - a ) < tol )
// 		{
// 			double newa = ( target[i].first * weights[i] + a ) / ( weights[i] + 1 ) ;
// 			double newv = target[i].second + v ;
// 			target[i] = std::make_pair( newa, newv ) ;
// 			weights[i]++ ;
// 			return ;
// 		}
// 	}
	target.push_back( std::make_pair( a, v ) ) ;

	weights.push_back( 1 ) ;
}

void RotatingCrack::postProcess()
{
	if ( converged && getState()[0] > 0 && damaging )
	{
		if ( inTension )
		{
			addAndConsolidate( tensionAngles, tensionweights, currentAngle, getState()[0] - currentDamage ) ;
		}
		else
		{
			addAndConsolidate( compressionAngles, compressionweights, currentAngle, getState()[0] - currentDamage ) ;
		}
		damaging = false ;
	}

//  if(tensionAngles.size() > 400)
//  {
//    for(double i = 0 ; i < M_PI ; i+= 0.01)
//    {
//      std::cout << i << "   " << damageAtAngle(tensionAngles, i) << std::endl ;
//    }
//    exit(0) ;
//  }
}

RotatingCrack::~RotatingCrack()
{
}



}
