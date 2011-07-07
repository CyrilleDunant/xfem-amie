//
// C++ Interface: damagemodel
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "damagemodel.h"
#include "../fracturecriteria/fracturecriterion.h"
#include "../../mesher/delaunay.h"

namespace Mu
{

/** \brief Damage model interface */


void DamageModel::step( ElementState &s )
{
	elementState = &s ;
	double phi = ( 1. + sqrt( 5. ) ) * .5 ;
	double resphi = 2. - phi ;   //goldensearch
// 		resphi = .5 ;              //bisection
// 		resphi = .1 ;                //down bias

	if( fraction < 0 )
	{
		upState.resize( state.size(), 0. );
		downState.resize( state.size(), 0. );

		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			fraction = s.getParent()->area() ;
		else
			fraction = s.getParent()->volume();

	}

	change = false ;

	if( wasBroken )
	{
		converged = true ;
		return ;
	}

	
	std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s ) ;
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;
	bool isInDamagingSet = s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() ;

	if( !isInDamagingSet )
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		converged = true ;

		if( fractured() )
			wasBroken = true ;

		return ;
	}
	std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ;
	bool checkpoint = s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ;

	if( checkpoint ) // initiate iteration
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		states.clear() ;

		getPreviousState() = getState() ;

		if( !fractured() )
		{
			converged = false ;
			change = true ;

			downState = damageIncrement.first;
			upState = damageIncrement.second;
			if(needRestart)
			{
				trialRatio = 0.  ;
				getState( true ) = downState ;
				return ;
			}
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), -setChange.first,0., score, setChange.second ) ) ;
			
			trialRatio = 1.  ;
			getState( true ) = downState + ( upState - downState ) * trialRatio ;
			while(fractured())
			{
				trialRatio -= damageDensityTolerance*.25 ;
				getState( true ) = downState + ( upState - downState ) * trialRatio ;
			}

			if( ( upState - downState ).min() < 0 )
			{
				while(!fractured())
				{
					trialRatio += damageDensityTolerance*.25 ;
					getState( true ) = downState + ( upState - downState ) * trialRatio ;
				}
				converged = true ;
				wasBroken = true ;
			}

			if( ( upState - downState ).max() < 2.*damageDensityTolerance )
			{
				while(!fractured())
				{
					trialRatio += damageDensityTolerance*.25 ;
					getState( true ) = downState + ( upState - downState ) * trialRatio ;
				}
				converged = true ;
				wasBroken = true ;
			}
			
		}
		else
		{
			wasBroken = true ;
			converged = true ;
		}

	}
	else if( !converged )
	{
		double scoreTolerance = s.getParent()->getBehaviour()->getFractureCriterion()->getScoreTolerance() ;
		change = true ;
		
		if(needRestart && states.empty()) 
		{
			getState( true ) = downState ;
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, 0, score, setChange.second ) ) ;
			return ;
		}

		if(states.size() == 1 && std::abs(trialRatio - states[0].fraction) < POINT_TOLERANCE_2D)
		{
			trialRatio = 1 ;
			getState( true ) = downState + ( upState - downState ) * trialRatio ;
			while(fractured())
			{
				trialRatio -= damageDensityTolerance*.25 ;
				getState( true ) = downState + ( upState - downState ) * trialRatio ;
			}
			
			return ;
		}
		
		states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second ) ) ;
		std::stable_sort( states.begin(), states.end() ) ;

		double minFraction = states[0].fraction ;
		double maxFraction = states[1].fraction ;
		double prevDelta = states[0].delta ;
		double prevScore = states[0].score ;
		double currentDelta = states[0].delta ;
		double currentScore = states[0].score ;
		bool deltaRoot = false ;
		bool scoreRoot = false ;
		
// 		std::cout << states[0].fraction << "  " << states[0].score << "  " << states[0].delta << std::endl ;
		for( int i = 1 ; i < states.size() ; i++ )
		{
// 			std::cout << states[i].fraction << "  " << states[i].score << "  " << states[i].delta << std::endl ;
			currentDelta = states[i].delta ;
			currentScore = states[i].score ;
			minFraction = states[i - 1].fraction ;
			maxFraction = states[i].fraction ;

			if( currentDelta * prevDelta < 0 || currentScore * prevScore < 0 )
			{
				deltaRoot = currentDelta * prevDelta < 0 ;
				scoreRoot = currentScore * prevScore < 0 ;
				break ;
			}
			else
			{
				prevDelta = states[i].delta ;
				prevScore = states[i].score ;
			}
		}
		trialRatio = ( minFraction + maxFraction ) * .5 ;

		getState( true ) = downState + ( upState - downState ) * trialRatio ;
		
		if( std::abs( minFraction - maxFraction ) * std::abs( upState - downState ).max()  < damageDensityTolerance )
		{
			getState( true ) = downState + ( upState - downState ) * trialRatio ;

			if( states.size() < 6 )
			{
				trialRatio = 1 ;
				getState( true ) = upState ;
				wasBroken = true ;
			}
			
			if(!deltaRoot && !scoreRoot)
			{
				trialRatio = 1 ;
				getState( true ) = downState + ( upState - downState ) * trialRatio ;
			}
// 			else if(deltaRoot)
// 			{
// 				if(prevScore < 0)
// 					trialRatio = minFraction ;
// 				else
// 					trialRatio = maxFraction ;
// 			}
// 			else
// 			{
// 				if(prevDelta > 0)
// 					trialRatio = minFraction ;
// 				else
// 					trialRatio = maxFraction ;
// 			}

			converged = true ;

			
		}
	}
}

DamageModel::DamageModel()
{
	elementState = NULL ;
	wasBroken = false ;
	change = false ;
	isNull = true ;
	needRestart = false ;
	thresholdDamageDensity = 1.-1.e-6 ;
	secondaryThresholdDamageDensity = 1.-1.e-6 ;

	fraction = -1 ;
	converged = true ;

	// The exploration increment is crucial for finding
	// the correct distribution of damage: the effect
	// of damage increment on the distribution of
	// fracture criterion scores is non-monotonic.
	damageDensityTolerance = 0.5e-3 ; 1. / pow( 2., 64 );
} ;

double DamageModel::getThresholdDamageDensity() const
{
	return thresholdDamageDensity ;
}

double DamageModel::getSecondaryThresholdDamageDensity() const
{
	return secondaryThresholdDamageDensity ;
}

Vector &DamageModel::getState( bool )
{
	return state ;
}


Vector DamageModel::smoothedState( const ElementState &s ) const
{
	if( !s.getParent()->getBehaviour()->getFractureCriterion() )
		return getState() ;

	Vector stra = getState() ;

	double physicalCharacteristicRadius = s.getParent()->getBehaviour()->getFractureCriterion()->getMaterialCharacteristicRadius() ;
	MirrorState mirroring = s.getParent()->getBehaviour()->getFractureCriterion()->mirroring ;
	double delta_x =  s.getParent()->getBehaviour()->getFractureCriterion()->delta_x ;
	double delta_y =  s.getParent()->getBehaviour()->getFractureCriterion()->delta_y ;
	double delta_z =  s.getParent()->getBehaviour()->getFractureCriterion()->delta_z ;
	std::vector <unsigned int> cache = s.getParent()->getBehaviour()->getFractureCriterion()->cache ;
	std::vector <DelaunayTreeItem * >* mesh2d = s.getParent()->getBehaviour()->getFractureCriterion()->mesh2d ;
	std::vector <DelaunayTreeItem3D * >* mesh3d = s.getParent()->getBehaviour()->getFractureCriterion()->mesh3d ;

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double area = s.getParent()->area() ;

		stra *= area ;

		double fact = area ;

		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;

			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()
			        || !ci->getBehaviour()->getFractureCriterion()
			        || ci->getBehaviour()->type == VOID_BEHAVIOUR
			        || ci->getBehaviour()->fractured()
			        || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource()
			        || dc > 2. * physicalCharacteristicRadius * physicalCharacteristicRadius )
			{
				continue ;
			}

			//this is to eliminate scaling effects ;
			double factor = 1 ;

			if( std::abs( s.getParent()->getBehaviour()->param[0][0] ) > POINT_TOLERANCE_3D )
				factor = std::abs( ci->getBehaviour()->param[0][0] / s.getParent()->getBehaviour()->param[0][0] ) ;

			double d = exp( -dc / ( 2.* physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;

			stra += ci->getBehaviour()->getDamageModel()->getState() * area * d ;
			fact += area * d ;

			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				stra += ci->getBehaviour()->getDamageModel()->getState() * area * d ;
				fact += area * d ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				stra += ci->getBehaviour()->getDamageModel()->getState() * area * d ;
				fact += area * d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				stra += ci->getBehaviour()->getDamageModel()->getState() * area * d ;
				fact += area * d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				stra += ci->getBehaviour()->getDamageModel()->getState() * area * d ;
				fact += area * d ;
			}
		}

		stra /= fact ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double volume = s.getParent()->volume() ;

		stra *= volume ;
		double fact = volume ;

		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			const DelaunayTetrahedron *ci = static_cast<const DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;

			if( dynamic_cast<const IntegrableEntity *>( ci ) == s.getParent()
			        || ci->getBehaviour()->getFractureCriterion()
			        || ci->getBehaviour()->type == VOID_BEHAVIOUR
			        || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource()
			        || dc > 2.* physicalCharacteristicRadius * physicalCharacteristicRadius
			  )
			{
				continue ;
			}

			double volume = ci->volume() ;
			double factor = 1 ;

			if( std::abs( s.getParent()->getBehaviour()->param[0][0] ) > POINT_TOLERANCE_3D )
				factor = std::abs( ci->getBehaviour()->param[0][0] / s.getParent()->getBehaviour()->param[0][0] ) ;

			double d = exp( -dc / ( 2.*physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;


			if( !ci->getBehaviour()->fractured() )
			{
				stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
				fact += volume * d ;

				if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}

				if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
				{
					stra += ci->getBehaviour()->getDamageModel()->getState() * volume * d ;
					fact += volume * d ;
				}
			}
		}

		stra /= fact ;
	}

	return stra ;
}

void DamageModel::setThresholdDamageDensity( double d )
{
	thresholdDamageDensity = d ;
}

void DamageModel::setSecondaryThresholdDamageDensity( double d )
{
	secondaryThresholdDamageDensity = d ;
}

void DamageModel::setDamageDensityTolerance( double d )
{
	damageDensityTolerance = d ;
}

bool DamageModel::changed() const
{
	return change ;
}
} ;
