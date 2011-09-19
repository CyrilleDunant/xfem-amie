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


void DamageModel::stepBack()
{
	change = previouschange ;
	previousstate = previouspreviousstate ;
	state = previousstate ;
	previousauxiliarystate = previouspreviousauxiliarystate ;
	auxiliarystate = previousauxiliarystate ;
}

void DamageModel::step( ElementState &s )
{
	previouschange = change ;
	previouspreviousstate = previousstate ;
	previousstate = state;
	previouspreviousauxiliarystate = previousauxiliarystate ;
	previousauxiliarystate = auxiliarystate ;
	elementState = &s ;
	
	
	double phi = ( 1. + sqrt( 5. ) ) * .5 ;
	double resphi = 2. - phi ;   //goldensearch
// 		resphi = .5 ;              //bisection
// 		resphi = .1 ;                //down bias

	if( fraction < 0 )
	{
		upState.resize( state.size(), 0. );
		downState.resize( state.size(), 0. );
		limitState.resize( state.size(), 0. );

		if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			fraction = s.getParent()->area() ;
		else
			fraction = s.getParent()->volume();

	}

	change = false ;
	
	std::pair<double, double> setChange = s.getParent()->getBehaviour()->getFractureCriterion()->setChange( s ) ;
	double score = s.getParent()->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() ;
	bool isInDamagingSet = s.getParent()->getBehaviour()->getFractureCriterion()->isInDamagingSet() ;

	if( !isInDamagingSet )
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		
		// this is necessary because we want to trigger side-effects
		//for example, plasticstrain gets a pointer to s
		computeDamageIncrement( s ) ;
		converged = true ;

		return ;
	}
	
	std::pair<Vector, Vector> damageIncrement = computeDamageIncrement( s ) ;
	if( s.getParent()->getBehaviour()->getFractureCriterion()->isAtCheckpoint() ) // initiate iteration
	{
		s.getParent()->getBehaviour()->getFractureCriterion()->setCheckpoint( false );
		states.clear() ;
		if(!fractured())
		{
			converged = false ;
			change = true ;
			
			downState = damageIncrement.first;
			upState = damageIncrement.second;
			
			if(!haslimit)
			{
				
				double fraction0 = 0 ;
				double fraction1 = 0.25 ;
				double fraction2 = 0.5 ;
				double fraction3 = 0.75 ;
				double fraction4 = 1 ;
				double score0 = score ;
				getState( true ) = downState + ( upState - downState ) * fraction1 ;
				double score1 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
				getState( true ) = downState + ( upState - downState ) * fraction2 ;
				double score2 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
				getState( true ) = downState + ( upState - downState ) * fraction3 ;
				double score3 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
				getState( true ) = downState + ( upState - downState ) * fraction4 ;
				double score4 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
				
				while(std::abs(fraction0-fraction4) > damageDensityTolerance*.125)
				{
					if(score1 <= score2 && score1 <= score3)
					{
						fraction4 = fraction2 ;
						score4 = score2 ;
						fraction2 = fraction1 ;
						score2 = score1 ;

						fraction1 = (fraction2+fraction0)*.5 ;
						fraction3 = (fraction4+fraction2)*.5 ;
						getState( true ) = downState + ( upState - downState ) * fraction1 ;
						score1 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
						getState( true ) = downState + ( upState - downState ) * fraction3 ;
						score3 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
					}
					else if(score2 <= score1 && score2 <= score3)
					{
						fraction0 = fraction1 ;
						score0 = score1 ;
						fraction4 = fraction3 ;
						score4 = score3 ;
						fraction1 = (fraction2+fraction0)*.5 ;
						fraction3 = (fraction4+fraction2)*.5 ;
						getState( true ) = downState + ( upState - downState ) * fraction1 ;
						score1 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
						getState( true ) = downState + ( upState - downState ) * fraction3 ;
						score3 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
					}
					else
					{
						fraction0 = fraction2 ;
						score0 = score2 ;
						fraction2 = fraction3 ;
						score2 = score3 ;
						fraction1 = (fraction2+fraction0)*.5 ;
						fraction3 = (fraction4+fraction2)*.5 ;
						getState( true ) = downState + ( upState - downState ) * fraction1 ;
						score1 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
						getState( true ) = downState + ( upState - downState ) * fraction3 ;
						score3 = s.getParent()->getBehaviour()->getFractureCriterion()->grade(s) ;
					}
				}
				limitState = downState + ( upState - downState ) * fraction2 ;
				haslimit = true ;
			}

			if(needRestart)
			{
				trialRatio = 0.  ;
				getState( true ) = downState ;
				return ;
			}
			states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first,0., score, setChange.second ) ) ;
			trialRatio = 1 ;
			getState( true ) = downState + ( upState - downState ) * trialRatio ;
			while(getState().max() > thresholdDamageDensity)
			{
				upState -= damageDensityTolerance*.25 ;
				for(int i = 0 ; i < upState.size() ; i++)
					upState[i]  = std::max(upState[i], 0.) ;
				getState( true ) = downState + ( upState - downState ) * trialRatio ;
			}
			if( ( upState - downState ).min() < 0 )
			{
				for(int i = 0 ; i < downState.size() ; i++)
					downState[i]  = std::min(downState[i], 1.) ;
				for(int i = 0 ; i < upState.size() ; i++)
					upState[i]  = std::min(upState[i], 1.) ;
				while(( upState - downState ).min() < 0)
				{
					upState += damageDensityTolerance*.25 ;
					getState( true ) = downState + ( upState - downState ) * trialRatio ;
				}
				converged = true ;
			}
			if( ( upState - downState ).max() < 2.*damageDensityTolerance )
			{
				getState( true ) += damageDensityTolerance ;
// 				
// 				while(getState().max() < thresholdDamageDensity)
// 				{
// 					getState( true ) = 1 ;
// 					for(int i = 0 ; i < upState.size() ; i++)
// 						upState[i]  = std::min(upState[i], 1.) ;
// 					getState( true ) = downState + ( upState - downState ) * trialRatio ;
// 				}
// 				for(int i = 0 ; i < getState( true ).size() ; i++)
// 				{
// 					getState( true )[i]  = std::min(getState( true )[i], 1.) ;
// 				}
				
				converged = true ;
			}
			
		}
		else
		{
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

		if(needRestart && states.size() == 1 && std::abs(trialRatio - states[0].fraction) < POINT_TOLERANCE_2D)
		{
			trialRatio = 1 ;
			getState( true ) = downState + ( upState - downState ) * trialRatio ;
			while(getState().max() > thresholdDamageDensity)
			{
				trialRatio -= damageDensityTolerance*.25 ;
				getState( true ) = downState + ( upState - downState ) * trialRatio ;
			}
			
			return ;
		}
		
		states.push_back( PointState( s.getParent()->getBehaviour()->getFractureCriterion()->met(), setChange.first, trialRatio, score, setChange.second ) ) ;
		std::stable_sort( states.begin(), states.end() ) ;

// 		if(states.size() < 5 && (limitState-downState).max() > damageDensityTolerance)
// 		{
// 			if(states.size() == 2)
// 			{
// 				trialRatio = (limitState-downState).max()/(upState-downState).max() ;
// 				getState( true ) = downState + ( upState - downState ) *trialRatio ;
// 				return ;
// 			}
// 			if(states.size() == 3)
// 			{
// 				std::cout << score << std::endl ;
// 				trialRatio *= 0.5*((limitState-downState).max()/(upState-downState).max()+0.) ;
// 				getState( true ) = downState + ( upState - downState ) *trialRatio ;
// 				return ;
// 			}
// 			if(states.size() == 4)
// 			{
// 				trialRatio = 0.5*((limitState-downState).max()/(upState-downState).max()+1) ;
// 				getState( true ) = downState + ( upState - downState ) *trialRatio ;
// 				return ;
// 			}
// 		}
		
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
		

		trialRatio = minFraction*(1.-0.5) + maxFraction*0.5  ;
		getState( true ) = downState + ( upState - downState ) *trialRatio ;
		
		if( std::abs( minFraction - maxFraction ) < damageDensityTolerance )
		{
// 			for(int i = 0 ; i < state.size() ; i++)
// 				std::cout << state[i] << std::endl ;
			
			getState( true ) = downState + ( upState - downState) * trialRatio;
			for(int i = 0 ; i < getState().size() ; i++)
					 getState(true)[i]  = std::min(getState(true)[i], 1.) ;
			//we have found no root, so now, we disturb the system
			if(!deltaRoot && !scoreRoot) 
			{
				trialRatio = 1 ;
				getState( true ) = upState ;
				converged = true ;
				return ;
			}
			
			if( states.size() < 6 )
			{
				trialRatio = 1 ;
				getState( true ) = upState ;
			}
			converged = true ;
		}
	}
}

void DamageModel::postProcess()
{
}

DamageModel::DamageModel(): state(0), previousstate(0), previouspreviousstate(0), auxiliarystate(0), previousauxiliarystate(0),previouspreviousauxiliarystate(0)
{
	elementState = NULL ;
	previouschange = false ;
	change = false ;
	isNull = true ;
	haslimit = false ;
	needRestart = false ;
	thresholdDamageDensity = 1./*-1.e-8*/ ;
	secondaryThresholdDamageDensity = 1./*-1.e-8*/ ;

	fraction = -1 ;
	converged = true ;

	// The exploration increment is crucial for finding
	// the correct distribution of damage: the effect
	// of damage increment on the distribution of
	// fracture criterion scores is non-monotonic.
	damageDensityTolerance =  1. / pow( 2., 14 );
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


Vector DamageModel::smoothedState( const ElementState &s , bool setUpdate) const
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
			        || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource()
			        || dc > 3. * physicalCharacteristicRadius * physicalCharacteristicRadius )
			{
				continue ;
			}

			if(setUpdate)
			{
				ci->behaviourUpdated = true ;
				ci->getBehaviour()->getDamageModel()->computeDamageIncrement(ci->getState()) ;
			}
			//this is to eliminate scaling effects ;
			double factor = 1 ;
			area = ci->area() ;
// 			if( std::abs( s.getParent()->getBehaviour()->param[0][0] ) > POINT_TOLERANCE_3D )
// 				factor = std::abs( ci->getBehaviour()->param[0][0] / s.getParent()->getBehaviour()->param[0][0] ) ;

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
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;

			if( dynamic_cast<const IntegrableEntity *>( ci ) == s.getParent()
			        || ci->getBehaviour()->getFractureCriterion()
			        || ci->getBehaviour()->type == VOID_BEHAVIOUR
			        || ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource()
			        || dc > 3.* physicalCharacteristicRadius * physicalCharacteristicRadius
			  )
			{
				continue ;
			}

			if(setUpdate)
			{
				ci->behaviourUpdated = true ;
				ci->getBehaviour()->getDamageModel()->computeDamageIncrement(ci->getState()) ;
			}
			double volume = ci->volume() ;
			double factor = 1 ;

// 			if( std::abs( s.getParent()->getBehaviour()->param[0][0] ) > POINT_TOLERANCE_3D )
// 				factor = std::abs( ci->getBehaviour()->param[0][0] / s.getParent()->getBehaviour()->param[0][0] ) ;

			double d = exp( -dc / ( 2.*physicalCharacteristicRadius * physicalCharacteristicRadius ) ) * factor;



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
