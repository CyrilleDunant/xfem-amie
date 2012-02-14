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

		if ( cos( std::abs( increments[i].first - angle ) ) > 0 ||
		        cos( std::abs( increments[i].first - angle + M_PI ) ) > 0 ||
		        cos( std::abs( increments[i].first - angle - M_PI ) ) > 0 )
		{
			double a = cos( std::abs( increments[i].first - angle ) ) * cos( std::abs( increments[i].first - angle ) ) ;
			double b = cos( std::abs( increments[i].first - angle + M_PI ) ) * cos( std::abs( increments[i].first - angle + M_PI ) ) ;
			double c = cos( std::abs( increments[i].first - angle - M_PI ) ) * cos( std::abs( increments[i].first - angle - M_PI ) ) ;
			ret += std::max( a, std::max( b, c ) ) * increments[i].second ;
		}
	}

	return ret ;
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( ElementState &s )
{
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
	{
		currentAngle = s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle();

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() )
		{
			damaging = false ;

			if ( s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() )
				damaging = true ;
		}

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->metInTension )
		{
			if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() )
				inTension = true ;

			cdamage = damageAtAngle( compressionAngles, currentAngle ) ;

			tdamage = damageAtAngle( tensionAngles, currentAngle ) ;

			if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint()  || converged )
				currentDamage = damageAtAngle( tensionAngles, currentAngle ) ;
		}
		else
		{
			if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() )
				inTension = false ;

			tdamage = damageAtAngle( tensionAngles, currentAngle ) ;

			cdamage = damageAtAngle( compressionAngles, currentAngle ) ;

			if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint()  || converged )
				currentDamage = damageAtAngle( compressionAngles, currentAngle ) ;
		}

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() || converged )
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

	if ( damaging )
	{
		if ( inTension )
		{
			return OrthothropicStiffness( factor *E * ( 1. - getState()[0] ), factor * E * ( 1. - cdamage ), factor * E * ( 1. - 0.5 * ( cdamage + getState()[0] ) ) * ( 1. - nu ) * .5, nu *( 1. - std::max( getState()[0], cdamage ) ), currentAngle ).getTensor( Point() ) ;
		}
		else
		{
			return OrthothropicStiffness( factor * E * ( 1. - tdamage ), factor * E * ( 1. - getState()[0] ), factor * E * ( 1. - 0.5 * ( getState()[0] + tdamage ) ) * ( 1. - nu ) * .5, nu, currentAngle ).getTensor( Point() ) ;
		}
	}

	return OrthothropicStiffness( factor * E * ( 1. - tdamage ), factor * E * ( 1. - cdamage ), factor * E * ( 1. - 0.5 * ( cdamage + tdamage ) ) * ( 1. - nu ) * .5, nu*( 1. - std::max( tdamage, cdamage ) ), currentAngle ).getTensor( Point() ) ;
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
	}

	target.push_back( std::make_pair( a, v ) ) ;

	weights.push_back( 1 ) ;
}

void RotatingCrack::postProcess()
{
	if ( converged && getState()[0] > POINT_TOLERANCE_2D && damaging )
	{
		if ( inTension )
		{
			addAndConsolidate( tensionAngles, tensionweights, currentAngle, getState()[0] - currentDamage ) ;
			addAndConsolidate( compressionAngles, compressionweights, currentAngle, 0.01*(getState()[0] - currentDamage) ) ;
		}
		else
		{
			addAndConsolidate( compressionAngles, compressionweights, currentAngle, getState()[0] - currentDamage ) ;
			addAndConsolidate( tensionAngles, tensionweights, currentAngle, 0.01*(getState()[0] - currentDamage) ) ;
		}
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
