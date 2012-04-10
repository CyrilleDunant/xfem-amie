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
	getState( true ).resize( 4, 0. );
// 	getState( true )[0] = 0.998 ;
	getPreviousState().resize( 4, 0. );
	isNull = false ;
	currentAngle = 0 ;
	factor = 1 ;
	es = NULL ;
	firstTension = false ;
	secondTension = false ;
	firstTensionFailure = false ;
	secondTensionFailure = false ;
	firstCompressionFailure = false ;
	secondCompressionFailure = false ;
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
		(!firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0) 
		|| firstTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0)
		|| !secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInTension(1) 
		|| secondTension && es->getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(1))
		) 
	{
		return 1 ;
	}
	return -1 ;
}

double RotatingCrack::getAngleShift() const
{
	if(!es)
		return 0 ;
	return std::min(std::abs(currentAngle-es->getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle()), std::abs(std::abs(currentAngle-es->getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle())-M_PI)) ;
}

std::pair< Vector, Vector > RotatingCrack::computeDamageIncrement( ElementState &s )
{
	Vector range( 1., 4 ) ;
// 	std::cout << s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle() << std::endl ;
	
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() && 
		s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet())
	{
		es = &s ;
// 		if(getState().max() < POINT_TOLERANCE_2D)
			currentAngle = s.getParent()->getBehaviour()->getFractureCriterion()->getCurrentAngle();
// 		if(!fractured())

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0) )
		{
			firstTension = true ;
		}
		else
		{
			firstTension = false ;
		}

		if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(1) )
		{
			secondTension = true ;
		}
		else
		{
			secondTension = false ;
		}
			
		if(s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(0))
		{
			if ( firstTension )
			{
				range[1] = getState()[1] ;
			}
			else
			{
				range[0] = getState()[0] ;
			}
		}
		
		if(s.getParent()->getBehaviour()->getFractureCriterion()->directionMet(1))
		{
			if ( secondTension )
			{
				range[3] = getState()[3] ;
			}
			else
			{
				range[2] = getState()[2] ;
			}
		}
		
// 		if(tensionFailure)
// 		{
// 			inTension = false ;
// 			range[0] = getState()[0] ;
// 		}
// 		if(compressionFailure)
// 			range[1] = getState()[1] ;
	}

	return std::make_pair( state,  range) ;
}

Matrix RotatingCrack::apply( const Matrix &m ) const
{

	if ( getState().max() < POINT_TOLERANCE_2D)
		return m ;
	
	if(fractured())
		return m *0 ;
	
	double E_0 = factor*E ;
	double E_1 = factor*E ;
	double fs = getState()[0] ;
	double ss = getState()[2] ;
	if(firstTension)
		E_0 *= ( 1. - fs ) ;
	else
	{
		fs = getState()[1] ;
		E_0 *= ( 1. - fs ) ;
	}
	
	if(secondTension)
		E_1 *= ( 1. - ss ) ;
	else
	{
		ss = getState()[3] ;
		E_1 *= ( 1. - ss ) ;
	}
	
	return OrthothropicStiffness( E_0, 
																E_1, 
																factor * E * (1.-std::max(fs, ss)) * ( 1. - nu ) * .5, 
																nu * ( 1. -  std::max(getState()[1], getState()[2])), 
																currentAngle ).getTensor( Point() ) ;


}


void  RotatingCrack::computeDelta(const ElementState &s)
{
	Vector range( 1., 2 ) ;
		
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(0))
	{
		firstTension = true ;
		range[1] = getState()[1] ;
	}
	
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(0))
	{
		firstTension = false ;
		range[0] = getState()[0] ;
	}
	
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInTension(1))
	{
		secondTension = true ;
		range[3] = getState()[3] ;
	}
	
	if ( s.getParent()->getBehaviour()->getFractureCriterion()->directionInCompression(1))
	{
		secondTension = false ;
		range[2] = getState()[2] ;
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

// 	if(firstTension && firstTensionFailure)
// 		return true ;
// 	if(!firstTension && firstCompressionFailure)
// 		return true ;
// 	if(secondTension && secondTensionFailure)
// 		return true ;
// 	if(!secondTension && secondCompressionFailure)
// 		return true ;
	
 return getState().min() >= thresholdDamageDensity ;
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
	if(converged && getState()[0] >= thresholdDamageDensity)
	{
		firstTensionFailure = true ;
		getState(true)[0] = 1 ;
	}
	if(converged && getState()[1] >= thresholdDamageDensity)
	{
		firstCompressionFailure = true ;
		getState(true)[1] = 1 ;
	}
	if(converged && getState()[2] >= thresholdDamageDensity)
	{
		secondTensionFailure = true ;
		getState(true)[2] = 1 ;
	}
	if(converged && getState()[3] >= thresholdDamageDensity)
	{
		secondCompressionFailure = true ;
		getState(true)[3] = 1 ;
	}
}

RotatingCrack::~RotatingCrack()
{
}



}
