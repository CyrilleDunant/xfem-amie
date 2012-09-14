// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011

#include "boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"

using namespace Mu ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double d, int a) : BoundaryCondition( t, d, a ), pos( pos ) { } ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d, int a ) : BoundaryCondition( t, d, a ), pos( pos ) { } ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, double d, int a ) : BoundaryCondition( t, d, a ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( zm ), zmax( zp )
{

}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, const Function & d, int a ) : BoundaryCondition( t, d, a ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( zm ), zmax( zp )
{

}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, double d, int a  ) : BoundaryCondition( t, d, a ), pos( pos ), nearest( p ) {} ;

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, const Function & d, int a  ) : BoundaryCondition( t, d, a ), pos( pos ), nearest( p ) {} ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double d, int a  ): BoundaryCondition( t, d, a ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( 0 ), zmax( 0 )
{

}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, const Function & d, int a  ): BoundaryCondition( t, d, a ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( 0 ), zmax( 0 )
{

}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition( ElementarySurface * surface ) : BoundaryCondition( GENERAL, 0 ), surface( surface ), volume( nullptr )
{
}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition( ElementaryVolume * volume ) : BoundaryCondition( GENERAL, 0 ), surface( nullptr ), volume( volume )
{
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, double d, int a  ) : BoundaryCondition( t, d, a ), id( id ), surface( surface ), volume( nullptr ), gp(gp), Jinv(Jinv)
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a  ) : BoundaryCondition( t, d, a ), id( id ), surface( surface ), volume( nullptr ), gp(gp), Jinv(Jinv)
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementaryVolume * volume, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, double d, int a  ) : BoundaryCondition( t, d, a ), id( id ), surface( nullptr ), volume( volume ), gp(gp), Jinv(Jinv)
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementaryVolume * volume, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a  ) : BoundaryCondition( t, d, a ), id( id ), surface( nullptr ), volume( volume ), gp(gp), Jinv(Jinv)
{

}

void apply2DBC( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;

	double nTimePlanes = 1 ;
	if(e->getOrder() > CONSTANT_TIME_LINEAR)
	{
		nTimePlanes = e->timePlanes() ;
	}
	
	VirtualMachine vm ;
	
	for ( size_t idit = 0 ; idit < id.size() ; idit++ )
	{
		switch ( condition )
		{

			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;

			case FIX_ALONG_XI:
				a->setPointAlong( XI, 0, id[idit] ) ;
				break ;

			case SET_ALONG_XI:
				a->setPointAlong( XI, data, id[idit] ) ;
				break ;

			case FIX_ALONG_ETA:
				a->setPointAlong( ETA, 0, id[idit] ) ;
				break ;

			case SET_ALONG_ETA:
				a->setPointAlong( ETA, data, id[idit] ) ;
				break ;
				
			case SET_ALONG_INDEXED_AXIS:
				a->setPointAlongIndexedAxis( axis, data, id[idit] ) ;
				break ;

			case SET_FORCE_XI:
			{
				if ( !e->getBehaviour()->fractured() )
					a->setForceOn( XI, data/nTimePlanes, id[idit]) ;

				break ;
			}

			case SET_FORCE_ETA:
				if ( !e->getBehaviour()->fractured() )
					a->setForceOn( ETA, data/nTimePlanes, id[idit] ) ;

				break ;

			case SET_STRESS_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
					
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
				  v.push_back(TIME_VARIABLE) ;
				}

				Vector imposed( 3 ) ;
				imposed[0] = data ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;

				for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[j], gp, Jinv, v) ;
				  					
					a->addForceOn( XI, forces[0], id[idit] ) ;
					a->addForceOn( ETA, forces[1], id[idit] ) ;
				}

				break ;
			}

			case SET_STRESS_ETA:

			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
					
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
				  v.push_back(TIME_VARIABLE) ;
				}
				
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = data ;
				imposed[2] = 0 ;
				
				for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[j], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[idit] ) ;
					a->addForceOn( ETA, forces[1], id[idit] ) ;
				}
				
				return ;
			}

			case SET_STRESS_XI_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = data ;

				for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[j], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], idit ) ;
					a->addForceOn( ETA, forces[1], idit ) ;
				}

				return ;
			}

			case SET_FLUX_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 2 ) ;
				imposed[0] = data ;
				imposed[1] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, idit ) ;
				}

				return ;
			}

			case SET_FLUX_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 2 ) ;
				imposed[0] = 0 ;
				imposed[1] = data ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, idit ) ;
				}

				return ;
			}
						
			default:
				break;
		}
	}
}

void apply3DBC( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
// 	std::cout << "splash" << std::endl ;
// 	std::cout << (size_t)(e) << std::endl ;
// 	std::cout << (size_t)(e->getBehaviour()) << std::endl ;
// 	std::cout << (size_t)(e->getBehaviour()->type) << std::endl ;
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;
	
	VirtualMachine vm ;
	
	for ( size_t i = 0 ; i < id.size() ; i++ )
	{
		switch ( condition )
		{

			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;

			case FIX_ALONG_XI:
				a->setPointAlong( XI, 0., id[i] ) ;
				break ;

			case SET_ALONG_XI:
				a->setPointAlong( XI, data, id[i] ) ;
				break ;

			case FIX_ALONG_ETA:
				a->setPointAlong( ETA, 0., id[i] ) ;
				break ;

			case SET_ALONG_ETA:
				a->setPointAlong( ETA, data, id[i] ) ;
				break ;

			case FIX_ALONG_ZETA:
				a->setPointAlong( ZETA, 0., id[i] ) ;
				break ;

			case SET_ALONG_ZETA:
				a->setPointAlong( ZETA, data, id[i] ) ;
				break ;

			case SET_ALONG_INDEXED_AXIS:
				a->setPointAlongIndexedAxis( axis, data, id[i] ) ;
				break ;

			case SET_FORCE_XI:

				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( XI, data, id[i] ) ;

				break ;

			case SET_FORCE_ETA:
				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( ETA, data, id[i] ) ;

				break ;

			case SET_FORCE_ZETA:
				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( ZETA, data, id[i] ) ;

				break ;

			case SET_FLUX_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = data ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i] ) ;
				}

				return ;
			}

			case SET_FLUX_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = data ;
				imposed[2] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i] ) ;
				}

				return ;
			}
			
			case SET_FLUX_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = data ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i] ) ;
				}

				return ;
			}
			case SET_STRESS_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = data ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

				return ;
			}

			case SET_STRESS_ETA:

			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = 0 ;
				imposed[1] = data ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
//					std::cout << "constant\t" << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;

					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}
//				exit(0) ;
				
				return ;
			}

			case SET_STRESS_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = data ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

				return ;
			}

			case SET_STRESS_XI_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = data ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

				return ;
			}

			case SET_STRESS_XI_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = data ;
				imposed[5] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

				return ;
			}

			case SET_STRESS_ETA_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j] == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j] == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = data ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( imposed, shapeFunctions[i], gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

				return ;
			}

			default:
				break;
		}
	}
}

void apply2DBC( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
	std::vector<size_t> ids ;

	for ( size_t i = 0 ; i < id.size() ; i++ )
		ids.push_back( id[i].id );

	apply2DBC( e,gp,Jinv, ids, condition, data, a, axis ) ;
}

void apply3DBC( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
	std::vector<size_t> ids ;
	
	for ( size_t i = 0 ; i < id.size() ; i++ )
		ids.push_back( id[i].id );

	apply3DBC( e, gp, Jinv, ids, condition, data, a, axis ) ;
}

void apply2DBC( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a, int axis = 0 )
{
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;

	VirtualMachine vm ;

//	std::cerr << id.size() << std::endl ;
	for ( size_t i = 0 ; i < id.size() ; i++ )
	{
		switch ( condition )
		{

			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;

			case FIX_ALONG_XI:
				a->setPointAlong( XI, 0, id[i].id ) ;
				break ;

			case SET_ALONG_XI:
				a->setPointAlong( XI, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case FIX_ALONG_ETA:
				a->setPointAlong( ETA, 0, id[i].id ) ;
				break ;

			case SET_ALONG_ETA:
				a->setPointAlong( ETA,  vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case SET_ALONG_INDEXED_AXIS:
				a->setPointAlongIndexedAxis( axis, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case SET_FORCE_XI:

				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( XI,  vm.eval( data, id[i] ), id[i].id ) ;

				break ;

			case SET_FORCE_ETA:
				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( ETA,  vm.eval( data, id[i] ), id[i].id ) ;

				break ;

			case SET_STRESS_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;
				Vector imposed( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 0, 3, shapeFunctions[i], e, gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
				}

				return ;
			}

			case SET_STRESS_ETA:

			{
				if ( e->getBehaviour()->fractured() )
					return ;

//				std::cout << vm.eval(data, id[i]) << std::endl ;
//				id[i].print() ;
				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;

				VirtualMachine vm ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 1, 3, shapeFunctions[i], e,gp, Jinv, v) ;
				
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
				}

				return ;
			}

			case SET_STRESS_XI_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 2, 3, shapeFunctions[i], e,gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
				}

				return ;
			}

			case SET_FLUX_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 2 ) ;
				imposed[0] = vm.eval( data, id[i] ) ; ;
				imposed[1] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i].id ) ;
				}

				return ;
			}

			case SET_FLUX_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = vm.eval( data, id[i] ) ; ;
				imposed[2] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i].id ) ;
				}

				return ;
			}
			

			default:
				break;
		}
	}
}

void apply3DBC( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a, int axis = 0 )
{
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;

	VirtualMachine vm ;

	for ( size_t i = 0 ; i < id.size() ; i++ )
	{
		switch ( condition )
		{

			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;

			case FIX_ALONG_XI:
				a->setPointAlong( XI, 0, id[i].id ) ;
				break ;

			case SET_ALONG_XI:
				a->setPointAlong( XI, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case FIX_ALONG_ETA:
				a->setPointAlong( ETA, 0, id[i].id ) ;
				break ;

			case SET_ALONG_ETA:
				a->setPointAlong( ETA, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case FIX_ALONG_ZETA:
				a->setPointAlong( ZETA, 0, id[i].id ) ;
				break ;

			case SET_ALONG_ZETA:
				a->setPointAlong( ZETA, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case SET_ALONG_INDEXED_AXIS:
				a->setPointAlongIndexedAxis( axis, vm.eval( data, id[i] ), id[i].id ) ;
				break ;

			case SET_FORCE_XI:

				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( XI, vm.eval( data, id[i] ), id[i].id ) ;

				break ;

			case SET_FORCE_ETA:
				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( ETA, vm.eval( data, id[i] ), id[i].id ) ;

				break ;

			case SET_FORCE_ZETA:
				if ( e->getBehaviour()->fractured() )
					break ;

				a->setForceOn( ZETA, vm.eval( data, id[i] ), id[i].id ) ;

				break ;

			case SET_STRESS_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 0, 6, shapeFunctions[i], e,gp, Jinv, v) ;

//					std::cout << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
					
				}

				return ;
			}

			case SET_STRESS_ETA:

			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 1, 6, shapeFunctions[i], e,gp, Jinv, v) ;
					
//					std::cout << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;
					
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}
//				exit(0) ;

				return ;
			}

			case SET_STRESS_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;
				
				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 2, 6, shapeFunctions[i], e,gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

				return ;
			}

			case SET_STRESS_XI_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 3, 6, shapeFunctions[i], e,gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

				return ;
			}

			case SET_STRESS_XI_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 4, 6, shapeFunctions[i], e,gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

				return ;
			}

			case SET_STRESS_ETA_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID() )
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 6 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces = e->getBehaviour()->getForcesFromAppliedStress( data, 5, 6, shapeFunctions[i], e,gp, Jinv, v) ;
					
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

				return ;
			}

			case SET_FLUX_XI:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID())
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = vm.eval( data, id[i] ) ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;


				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i].id ) ;
				}

				return ;
			}

			case SET_FLUX_ETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID())
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = vm.eval( data, id[i] ) ;
				imposed[2] = 0 ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i].id ) ;
				}

				return ;
			}
			
			case SET_FLUX_ZETA:
			{
				if ( e->getBehaviour()->fractured() )
					return ;

				std::vector<Function> shapeFunctions ;

				for ( size_t j = 0 ; j < id.size() ; j++ )
				{
					for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
					{
						if ( id[j].id == e->getBoundingPoint( i ).id )
							shapeFunctions.push_back( e->getShapeFunction( i ) ) ;
					}
					for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
					{
						if ( id[j].id == e->getEnrichmentFunction( i ).getDofID())
							shapeFunctions.push_back( e->getEnrichmentFunction( i ) ) ;
					}
				}

				std::vector<Variable> v( 3 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				v[2] = ZETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = vm.eval( data, id[i] ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					double f = 0. ;
					f = vm.ieval( VectorGradient( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v) ;
					a->addForceOn( XI, f, id[i].id ) ;
				}

				return ;
			}
			default:
				break;
		}
	}
}


void applyVerticalPlaneSections(double topY, double bottomY, Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t)
{
	std::vector<Point *> points ;
	std::vector<DelaunayTriangle *> triangles = t->getElements() ;
	
	for(size_t i  = 0 ; i < triangles.size() ; i++)
	{
		for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
		{
			if(triangles[i]->getBoundingPoint(j).y <= topY && triangles[i]->getBoundingPoint(j).y >= bottomY)
			{
				points.push_back(&triangles[i]->getBoundingPoint(j));
			}
		}
	}
	std::sort(points.begin(), points.end()) ;
	auto e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	for(size_t i = 0 ; i <  points.size() ; i++)
	{
		Point topPoint(points[i]->x, topY) ;
		Point bottomPoint(points[i]->x, bottomY) ;
		DelaunayTriangle * topElement = t->getUniqueConflictingElement(&topPoint) ;
		DelaunayTriangle * bottomElement = t->getUniqueConflictingElement(&bottomPoint) ;
		if(topElement && bottomElement)
		{
			std::vector<int> idstop ;
			std::vector<double> coefficientstop ;
			idstop.push_back(topElement->first->id*2);
			coefficientstop.push_back(((topElement->second->y-topElement->third->y)*(topPoint.x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topPoint.y-topElement->third->y))
			/((topElement->second->y-topElement->third->y)*(topElement->first->x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topElement->first->y-topElement->third->y)));
			
			idstop.push_back(topElement->second->id*2);
			coefficientstop.push_back(((topElement->third->y-topElement->first->y)*(topPoint.y-topElement->third->x)
			+(topElement->first->x-topElement->third->x)*(topPoint.y-topElement->third->y))
			/((topElement->second->y-topElement->third->y)*(topElement->first->x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topElement->first->y-topElement->third->y)));
			
			idstop.push_back(topElement->third->id*2);
			coefficientstop.push_back(1.-coefficientstop[0]-coefficientstop[1]) ;
			
			std::vector<int> idsbot ;
			std::vector<double> coefficientsbot ;
			idsbot.push_back(bottomElement->first->id*2);
			coefficientsbot.push_back(((bottomElement->second->y-bottomElement->third->y)*(bottomPoint.x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomPoint.y-topElement->third->y))
			/((bottomElement->second->y-bottomElement->third->y)*(bottomElement->first->x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomElement->first->y-bottomElement->third->y)));
			
			idsbot.push_back(bottomElement->second->id*2);
			coefficientsbot.push_back(((bottomElement->third->y-bottomElement->first->y)*(bottomPoint.x-bottomElement->third->x)
			+(bottomElement->first->x-bottomElement->third->x)*(bottomPoint.y-bottomElement->third->y))
			/((bottomElement->second->y-bottomElement->third->y)*(bottomElement->first->x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomElement->first->y-bottomElement->third->y)));
			
			idsbot.push_back(bottomElement->third->id*2);
			coefficientsbot.push_back(1.-coefficientsbot[0]-coefficientsbot[1]) ;
			
			std::valarray<unsigned int> allids(7) ;
			Vector coefs(7) ;
			for(size_t j = 0 ; j < 3 ; j++)
			{
				coefficientstop[j] *= (points[i]->y-bottomY)/(topY-bottomY) ;
				coefficientsbot[j] *= (topY - points[i]->y)/(topY-bottomY) ;
				
				coefs[j] = coefficientstop[j] ;
				coefs[j+3] = coefficientsbot[j] ;
			}
			coefs[6] = -std::accumulate(&coefs[0], &coefs[6], double(0)) ;
			allids[0] = topElement->first->id*2 ;
			allids[1] = topElement->second->id*2 ;
			allids[2] = topElement->third->id*2 ;
			allids[3] = bottomElement->first->id*2 ;
			allids[4] = bottomElement->second->id*2 ;
			allids[5] = bottomElement->third->id*2 ;
			allids[6] = points[i]->id*2 ;
			a->addMultiplier(LagrangeMultiplier(allids, coefs, 0 ));
		}
		
		
	}
	
}


void applyHorizontalPlaneSections(double topX, double bottomX, Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t)
{
	std::vector<Point *> points ;
	std::vector<DelaunayTriangle *> triangles = t->getElements() ;
	
	for(size_t i  = 0 ; i < triangles.size() ; i++)
	{
		for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
		{
			if(triangles[i]->getBoundingPoint(j).x <= topX && triangles[i]->getBoundingPoint(j).x >= bottomX)
			{
				points.push_back(&triangles[i]->getBoundingPoint(j));
			}
		}
	}
	std::sort(points.begin(), points.end()) ;
	auto e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	for(size_t i = 0 ; i <  points.size() ; i++)
	{
		Point topPoint(topX, points[i]->y) ;
		Point bottomPoint(bottomX, points[i]->y) ;
		DelaunayTriangle * topElement = t->getUniqueConflictingElement(&topPoint) ;
		DelaunayTriangle * bottomElement = t->getUniqueConflictingElement(&bottomPoint) ;
		if(topElement && bottomElement)
		{
			std::vector<int> idstop ;
			std::vector<double> coefficientstop ;
			idstop.push_back(topElement->first->id*2+1);
			coefficientstop.push_back(((topElement->second->y-topElement->third->y)*(topPoint.x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topPoint.y-topElement->third->y))
			/((topElement->second->y-topElement->third->y)*(topElement->first->x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topElement->first->y-topElement->third->y)));
			
			idstop.push_back(topElement->second->id*2+1);
			coefficientstop.push_back(((topElement->third->y-topElement->first->y)*(topPoint.y-topElement->third->x)
			+(topElement->first->x-topElement->third->x)*(topPoint.y-topElement->third->y))
			/((topElement->second->y-topElement->third->y)*(topElement->first->x-topElement->third->x)
			+(topElement->third->x-topElement->second->x)*(topElement->first->y-topElement->third->y)));
			
			idstop.push_back(topElement->third->id*2+1);
			coefficientstop.push_back(1.-coefficientstop[0]-coefficientstop[1]) ;
			
			std::vector<int> idsbot ;
			std::vector<double> coefficientsbot ;
			idsbot.push_back(bottomElement->first->id*2+1);
			coefficientsbot.push_back(((bottomElement->second->y-bottomElement->third->y)*(bottomPoint.x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomPoint.y-topElement->third->y))
			/((bottomElement->second->y-bottomElement->third->y)*(bottomElement->first->x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomElement->first->y-bottomElement->third->y)));
			
			idsbot.push_back(bottomElement->second->id*2+1);
			coefficientsbot.push_back(((bottomElement->third->y-bottomElement->first->y)*(bottomPoint.x-bottomElement->third->x)
			+(bottomElement->first->x-bottomElement->third->x)*(bottomPoint.y-bottomElement->third->y))
			/((bottomElement->second->y-bottomElement->third->y)*(bottomElement->first->x-bottomElement->third->x)
			+(bottomElement->third->x-bottomElement->second->x)*(bottomElement->first->y-bottomElement->third->y)));
			
			idsbot.push_back(bottomElement->third->id*2+1);
			coefficientsbot.push_back(1.-coefficientsbot[0]-coefficientsbot[1]) ;
			
			
			std::valarray<unsigned int> allids(7) ;
			Vector coefs(7) ;
			for(size_t j = 0 ; j < 3 ; j++)
			{
				coefficientstop[j] *= (points[i]->x-bottomX)/(topX-bottomX) ;
				coefficientsbot[j] *= (topX - points[i]->x)/(topX-bottomX) ;
				
				coefs[j] = coefficientstop[j] ;
				coefs[j+3] = coefficientsbot[j] ;
			}
			coefs[6] = -std::accumulate(&coefs[0], &coefs[6], double(0)) ;
			allids[0] = topElement->first->id*2+1 ;
			allids[1] = topElement->second->id*2+1 ;
			allids[2] = topElement->third->id*2+1 ;
			allids[3] = bottomElement->first->id*2+1 ;
			allids[4] = bottomElement->second->id*2+1 ;
			allids[5] = bottomElement->third->id*2+1 ;
			allids[6] = points[i]->id*2+1 ;
			a->addMultiplier(LagrangeMultiplier(allids, coefs, 0 ));
		}
		
		
	}
	
}


void PlaneSectionsBoundaryConditions::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t)
{
	if(!isVertical)
	{
		applyHorizontalPlaneSections(uplimit, downlimit, a, t) ;
	}
	else
	{
		applyVerticalPlaneSections(uplimit, downlimit, a, t) ;
	}
}

void PlaneSectionsBoundaryConditions::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)
{
	
}


void DofDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	if ( !function )
	{
		std::vector<size_t> id_ ;
		id_.push_back( id );
		apply2DBC( surface,gp, Jinv, id_, condition, data*getScale(), a, axis ) ;
	}
	else
	{
		std::vector<Point> id_ ;

		for ( int i = 0 ; i < surface->getBoundingPoints().size() ; i++ )
		{
			if ( surface->getBoundingPoint( i ).id == id )
			{
				id_.push_back( surface->getBoundingPoint( i ) );
				apply2DBC( surface,gp,Jinv, id_, condition, dataFunction*getScale(), a, axis ) ;
				return ;
			}
		}
	}
}

void DofDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	if ( surface )
		return ;

	if ( !function )
	{
		std::vector<size_t> id_ ;
		id_.push_back( id );
		apply3DBC( volume,gp,Jinv,  id_, condition, data*getScale(),  a , axis) ;
	}
	else
	{
		std::vector<Point> id_ ;

		for ( int i = 0 ; i < volume->getBoundingPoints().size() ; i++ )
		{
			if ( volume->getBoundingPoint( i ).id == id )
			{
				id_.push_back( volume->getBoundingPoint( i ) );
				apply3DBC( volume,gp, Jinv, id_, condition, dataFunction*getScale(), a, axis ) ;
				return ;
			}
		}
	}
}

void ElementDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	if ( volume )
		return ;

	std::vector<DelaunayTriangle *> elements = t->getElements() ;

	std::set<Point *> points ;

	for ( size_t i = 0 ; i < elements.size() ; i++ )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
			continue ;

		for ( size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++ )
		{
			Point test( elements[i]->getBoundingPoint( j ) ) ;
			surface->project( &test );

			if ( dist( test, elements[i]->getBoundingPoint( j ) ) < POINT_TOLERANCE_2D )
				points.insert( &elements[i]->getBoundingPoint( j ) ) ;
		}
	}

	for ( auto i = points.begin() ; i != points.end() ; i++ )
	{
		Point local = surface->inLocalCoordinates( *( *i ) ) ;
		std::vector<Point> p ;
		p.push_back( local );

		Vector disps(0.,2) ;
		surface->getState().getField( DISPLACEMENT_FIELD, local, disps, true) ;
		a->setPointAlong( XI, disps[0], ( *i )->id ) ;
		a->setPointAlong( ETA, disps[1], ( *i )->id ) ;

	}

};

void ElementDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	if ( surface )
		return ;

	std::vector<DelaunayTetrahedron *> elements = t->getElements() ;

	std::set<Point *> points ;

	for ( size_t i = 0 ; i < elements.size() ; i++ )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
			continue ;

		for ( size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++ )
		{
			Point test( elements[i]->getBoundingPoint( j ) ) ;
			volume->project( &test );

			if ( dist( test, elements[i]->getBoundingPoint( j ) ) < POINT_TOLERANCE_3D )
				points.insert( &elements[i]->getBoundingPoint( j ) ) ;
		}
	}

	for ( auto i = points.begin() ; i != points.end() ; i++ )
	{
		Point local = volume->inLocalCoordinates( *( *i ) ) ;
		std::vector<Point> p ;
		p.push_back( local );

		Vector disps(0.,3) ;
		volume->getState().getField( DISPLACEMENT_FIELD, local, disps, true) ;
		a->setPointAlong( XI, disps[0], ( *i )->id ) ;
		a->setPointAlong( ETA, disps[1], ( *i )->id ) ;
		a->setPointAlong( ZETA, disps[2], ( *i )->id ) ;
	}
};

bool isInBoundary2D( Point test, Point min, Point max)
{
	return (test.x >= min.x && test.x <= max.x && test.y >= min.y && test.y <= max.y) ;
}

bool isInBoundary3D( Point test, Point min, Point max)
{
	return (test.x >= min.x && test.x <= max.x && test.y >= min.y && test.y <= max.y && test.z >= min.z && test.z <= max.z) ;
}


bool isOnBoundary( BoundingBoxPosition pos, Point test, Point min, Point max , double tol)
{
	switch( pos )
	{
		// 2D edges, 3D planes, 4D time planes
		case LEFT:
			return ( std::abs(test.x - min.x) < tol ) ;
		case RIGHT:
			return ( std::abs(test.x - max.x) < tol ) ;
		case BOTTOM:
			return ( std::abs(test.y - min.y) < tol ) ;
		case TOP:
			return ( std::abs(test.y - max.y) < tol ) ;
		case BACK:
			return ( std::abs(test.z - min.z) < tol ) ;
		case FRONT:
			return ( std::abs(test.z - max.z) < tol ) ;
		case BEFORE:
			return ( std::abs(test.t - min.t) < tol ) ;
		case AFTER:
			return ( std::abs(test.t - max.t) < tol ) ;
		
		// 2D corners, 3D edges, 4D planes
		case BOTTOM_LEFT:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol )) ;
		case BOTTOM_RIGHT:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol )) ;
		case TOP_LEFT:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol )) ;
		case TOP_RIGHT:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol )) ;
			
		case BACK_LEFT:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol )) ;
		case BACK_RIGHT:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol )) ;
		case FRONT_LEFT:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol )) ;
		case FRONT_RIGHT:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol )) ;

		case BOTTOM_BACK:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol )) ;
		case TOP_BACK:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol )) ;
		case FRONT_BOTTOM:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol )) ;
		case FRONT_TOP:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol )) ;

		case BOTTOM_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case BOTTOM_AFTER:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;
		case TOP_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case TOP_AFTER:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;
			
		case BACK_BEFORE:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case BACK_AFTER:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;
		case FRONT_BEFORE:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case FRONT_AFTER:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;

		case LEFT_BEFORE:
			return (isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case LEFT_AFTER:
			return (isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;
		case RIGHT_BEFORE:
			return (isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol )) ;
		case RIGHT_AFTER:
			return (isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol )) ;

		// 3D corners, 4D lines			
		case BOTTOM_LEFT_BACK:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) ) ;
		case BOTTOM_LEFT_FRONT:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) ) ;
		case BOTTOM_RIGHT_BACK:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) ) ;
		case BOTTOM_RIGHT_FRONT:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) ) ;
		case TOP_LEFT_BACK:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) ) ;
		case TOP_LEFT_FRONT:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) ) ;
		case TOP_RIGHT_BACK:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) ) ;
		case TOP_RIGHT_FRONT:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) ) ;
			
		case BOTTOM_LEFT_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case BOTTOM_RIGHT_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case TOP_LEFT_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case TOP_RIGHT_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
			
		case BACK_LEFT_BEFORE:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case BACK_RIGHT_BEFORE:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case FRONT_LEFT_BEFORE:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case FRONT_RIGHT_BEFORE:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;

		case BOTTOM_BACK_BEFORE:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case TOP_BACK_BEFORE:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case FRONT_BOTTOM_BEFORE:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case FRONT_TOP_BEFORE:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( BEFORE, test, min, max, tol) ) ;

		case BOTTOM_LEFT_AFTER:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case BOTTOM_RIGHT_AFTER:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case TOP_LEFT_AFTER:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case TOP_RIGHT_AFTER:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
			
		case BACK_LEFT_AFTER:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case BACK_RIGHT_AFTER:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case FRONT_LEFT_AFTER:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case FRONT_RIGHT_AFTER:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;

		case BOTTOM_BACK_AFTER:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case TOP_BACK_AFTER:
			return (isOnBoundary( BACK, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case FRONT_BOTTOM_AFTER:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;
		case FRONT_TOP_AFTER:
			return (isOnBoundary( FRONT, test, min, max, tol ) && isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( AFTER, test, min, max, tol) ) ;

		// 4D corners			
		case BOTTOM_LEFT_BACK_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol) ) ;
		case BOTTOM_LEFT_FRONT_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case BOTTOM_RIGHT_BACK_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case BOTTOM_RIGHT_FRONT_BEFORE:
			return (isOnBoundary( BOTTOM, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case TOP_LEFT_BACK_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case TOP_LEFT_FRONT_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( LEFT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case TOP_RIGHT_BACK_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( BACK, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
		case TOP_RIGHT_FRONT_BEFORE:
			return (isOnBoundary( TOP, test, min, max, tol ) && isOnBoundary( RIGHT, test, min, max, tol ) && isOnBoundary( FRONT, test, min, max, tol) && isOnBoundary( BEFORE, test, min, max, tol)  ) ;
			
	}
	return false ;
}




void BoundingBoxNearestNodeDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	if ( cache.empty() )
	{
		std::vector<DelaunayTriangle *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;
		
		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if ( elements[i]->getBoundingPoint( j ).x < minx )
					minx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).x > maxx )
					maxx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).y < miny )
					miny = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).y > maxy )
					maxy = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).t < mint )
					mint = elements[i]->getBoundingPoint( j ).t ;

				if ( elements[i]->getBoundingPoint( j ).t > maxt )
					maxt = elements[i]->getBoundingPoint( j ).t ;
			  
			}
		}

		double tol = std::min( maxx - minx, maxy - miny ) * .0001 ;
		
		Point pmin(minx,miny, 0., mint) ;
		Point pmax(maxx,maxy, 0., maxt) ;

		std::map<double, std::pair<Point, DelaunayTriangle*> > id  ;

		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR)
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if( isOnBoundary( pos, elements[i]->getBoundingPoint( j ), pmin, pmax, tol ) )
				{
					id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
				}
			}
		}

		std::vector<Point> target ;

		target.push_back( id.begin()->second.first ) ;
		cache2d.push_back( id.begin()->second.second ) ;
		cache.push_back( target ) ;
		GaussPointArray gp = id.begin()->second.second->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), id.begin()->second.second->getGaussPoints().gaussPoints.size() ) ;

		for ( size_t i = 0 ; i < id.begin()->second.second->getGaussPoints().gaussPoints.size() ; i++ )
		{
			id.begin()->second.second->getInverseJacobianMatrix( id.begin()->second.second->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
		}

		if ( !function )
			apply2DBC( id.begin()->second.second,gp,Jinv, target, condition, data*getScale(), a , axis ) ;
		else
			apply2DBC( id.begin()->second.second,gp,Jinv, target, condition, dataFunction*getScale(), a , axis ) ;
		
		
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			GaussPointArray gp = cache2d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache2d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}
			
			if ( !function )
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a ) ;
			else
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a ) ;
		}
	}
}

void BoundingBoxNearestNodeDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	if ( cache.empty() )
	{
		std::vector<DelaunayTetrahedron *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double minz = elements.front()->getBoundingPoint( 0 ).z ;
		double maxz = elements.front()->getBoundingPoint( 0 ).z ;
		
		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;
		
		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if ( elements[i]->getBoundingPoint( j ).x < minx )
					minx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).x > maxx )
					maxx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).y < miny )
					miny = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).y > maxy )
					maxy = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).z < minz )
					minz = elements[i]->getBoundingPoint( j ).z ;

				if ( elements[i]->getBoundingPoint( j ).z > maxz )
					maxz = elements[i]->getBoundingPoint( j ).z ;
				
				if ( elements[i]->getBoundingPoint( j ).t < mint )
					mint = elements[i]->getBoundingPoint( j ).t ;

				if ( elements[i]->getBoundingPoint( j ).t > maxt )
					maxt = elements[i]->getBoundingPoint( j ).t ;
			  
			}
		}

		double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;
		
		Point pmin( minx, miny, minz, mint) ;
		Point pmax( maxx, maxy, maxz, maxt) ;

		std::map<double, std::pair<Point, DelaunayTetrahedron*> > id  ;

		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR)
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if( isOnBoundary( pos, elements[i]->getBoundingPoint( j ), pmin, pmax, tol ) )
				{
					id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
				}
			}
		}

		std::vector<Point> target ;

		target.push_back( id.begin()->second.first ) ;
		cache3d.push_back( id.begin()->second.second ) ;
		cache.push_back( target ) ;
		GaussPointArray gp = id.begin()->second.second->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), id.begin()->second.second->getGaussPoints().gaussPoints.size() ) ;

		for ( size_t i = 0 ; i < id.begin()->second.second->getGaussPoints().gaussPoints.size() ; i++ )
		{
			id.begin()->second.second->getInverseJacobianMatrix( id.begin()->second.second->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
		}

		if ( !function )
			apply3DBC( id.begin()->second.second,gp,Jinv, target, condition, data*getScale(), a , axis) ;
		else
			apply3DBC( id.begin()->second.second,gp,Jinv, target, condition, dataFunction*getScale(), a, axis ) ;
		
		
	}
	else
	{
		for ( size_t i = 0 ; i < cache3d.size() ; ++i )
		{
			GaussPointArray gp = cache3d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache3d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}
			
			if ( !function )
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis) ;
			else
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis) ;
		}
	}
  
}

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition( LagrangeMultiplierType t, Geometry * source, double d, int a ) : BoundaryCondition( t, d, a ), domain( source ) { };

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition( LagrangeMultiplierType t, Geometry * source, const Function & d, int a ) : BoundaryCondition( t, d, a ), domain( source ) { };

void GeometryDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	std::vector<DelaunayTriangle *> elements = t->getElements() ;
	double tol = domain->getRadius() * .0001 ;


	for ( size_t i = 0 ; i < elements.size() ; ++i )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
			continue ;

		std::vector<Point> id  ;

		for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
		{
			Circle c( tol, elements[i]->getBoundingPoint( j ) ) ;

			if ( domain->intersects( &c ) || domain->in( elements[i]->getBoundingPoint( j ) ) )
			{
				id.push_back( elements[i]->getBoundingPoint( j ) ) ;
			}
		}
		GaussPointArray gp = elements[i]->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

		for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
		{
			elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}

		
		if ( !function )
			apply2DBC( elements[i],gp, Jinv, id, condition, data*getScale(), a ) ;
		else
			apply2DBC( elements[i],gp, Jinv, id, condition, dataFunction*getScale(), a ) ;
	}
}

void GeometryDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<DelaunayTetrahedron *> elements = t->getElements() ;
	double tol = domain->getRadius() * .0001 ;

	for ( size_t i = 0 ; i < elements.size() ; ++i )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
			continue ;

		std::vector<Point> id  ;

		for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
		{
			Sphere c( tol, elements[i]->getBoundingPoint( j ) ) ;

			if ( domain->intersects( &c ) || domain->in( elements[i]->getBoundingPoint( j ) ) )
			{
				id.push_back( elements[i]->getBoundingPoint( j ) ) ;
			}
		}
		
		GaussPointArray gp = elements[i]->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

		for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
		{
			elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}
		
		if ( !function )
			apply3DBC( elements[i],gp,Jinv, id, condition, data*getScale(), a ) ;
		else
			apply3DBC( elements[i],gp, Jinv, id, condition, dataFunction*getScale(), a ) ;
	}
}

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	if ( cache.empty() )
	{
		std::vector<DelaunayTriangle *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;
		
		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if ( elements[i]->getBoundingPoint( j ).x < minx )
					minx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).x > maxx )
					maxx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).y < miny )
					miny = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).y > maxy )
					maxy = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).t < mint )
					mint = elements[i]->getBoundingPoint( j ).t ;

				if ( elements[i]->getBoundingPoint( j ).t > maxt )
					maxt = elements[i]->getBoundingPoint( j ).t ;
			  
			}
		}

		double tol = std::min( maxx - minx, maxy - miny ) * .0001 ;
		
		Point pmin(minx,miny, 0., mint) ;
		Point pmax(maxx,maxy, 0., maxt) ;
		
		Point rmin(xmin, ymin) ;
		Point rmax(xmax, ymax) ;
		
		  for ( size_t i = 0 ; i < elements.size() ; ++i )
		  {
			  if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
				  continue ;

			  std::vector<Point> id  ;

			  for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			  {
				  if( isOnBoundary( pos, elements[i]->getBoundingPoint( j ), pmin, pmax, tol ) && isInBoundary2D( elements[i]->getBoundingPoint( j ), rmin, rmax ) )
				  {
					  if ( cache2d.empty() || cache2d.back() != elements[i] )
					  {
						  cache.push_back( std::vector<Point>() );
						  cache2d.push_back( elements[i] );
					  }

					  cache.back().push_back( elements[i]->getBoundingPoint( j ) ) ;
				  }
			  }


			  if ( !cache2d.empty() && cache2d.back() == elements[i] )
			  {
				  GaussPointArray gp = elements[i]->getGaussPoints() ;
				  std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

				  for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
				  {
					  elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
				  }
				  if ( !function )
					  apply2DBC( elements[i],gp,Jinv, cache.back(), condition, data*getScale(), a , axis ) ;
				  else
					  apply2DBC( elements[i],gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis) ;
			  }
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			GaussPointArray gp = cache2d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache2d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}
			
			if ( !function )
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis) ;
			else
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis) ;
		}
	}
}

void BoundingBoxDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{

	if ( cache.empty() )
	{
		std::vector<DelaunayTriangle *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;

		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
				continue ;

			if ( elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( elements[i]->getBoundingPoint( j ).x < minx )
						minx = elements[i]->getBoundingPoint( j ).x ;

					if ( elements[i]->getBoundingPoint( j ).x > maxx )
						maxx = elements[i]->getBoundingPoint( j ).x ;

					if ( elements[i]->getBoundingPoint( j ).y < miny )
						miny = elements[i]->getBoundingPoint( j ).y ;

					if ( elements[i]->getBoundingPoint( j ).y > maxy )
						maxy = elements[i]->getBoundingPoint( j ).y ;

					if ( elements[i]->getBoundingPoint( j ).t < mint )
						mint = elements[i]->getBoundingPoint( j ).t ;

					if ( elements[i]->getBoundingPoint( j ).t > maxt )
						maxt = elements[i]->getBoundingPoint( j ).t ;
				}
			}
		}
		
		Point pmin(minx, miny, 0., mint) ;
		Point pmax(maxx, maxy, 0., maxt) ;

		double tol = std::max( std::min( maxx - minx, maxy - miny ) * .001, POINT_TOLERANCE_2D ) ;

		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if( isOnBoundary( pos, elements[i]->getBoundingPoint(j), pmin, pmax, tol) )
				{
					if ( cache2d.empty() || cache2d.back() != elements[i] )
					{
						cache.push_back( std::vector<Point>() );
						cache2d.push_back( elements[i] );
					}

					cache.back().push_back( elements[i]->getBoundingPoint( j ) ) ;
				}
			}

			if ( !cache2d.empty() && cache2d.back() == elements[i] )
			{
				GaussPointArray gp = elements[i]->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
				{
					elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
				}
				
				if ( !function )
					apply2DBC( elements[i],gp,Jinv, cache.back(), condition, data*getScale(), a , axis) ;
				else
					apply2DBC( elements[i],gp,Jinv, cache.back(), condition, dataFunction*getScale(), a, axis ) ;
			}
		}
		
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			GaussPointArray gp = cache2d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache2d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}

			if ( !function )
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a, axis ) ;
			else
				apply2DBC( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a, axis ) ;
		}
	}
}

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	if ( cache.empty() )
	{
		std::vector<DelaunayTetrahedron *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double minz = elements.front()->getBoundingPoint( 0 ).z ;
		double maxz = elements.front()->getBoundingPoint( 0 ).z ;
		
		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;
		
		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if ( elements[i]->getBoundingPoint( j ).x < minx )
					minx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).x > maxx )
					maxx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).y < miny )
					miny = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).y > maxy )
					maxy = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).z < minz )
					minz = elements[i]->getBoundingPoint( j ).z ;

				if ( elements[i]->getBoundingPoint( j ).z > maxz )
					maxz = elements[i]->getBoundingPoint( j ).z ;

				if ( elements[i]->getBoundingPoint( j ).t < mint )
					mint = elements[i]->getBoundingPoint( j ).t ;

				if ( elements[i]->getBoundingPoint( j ).t > maxt )
					maxt = elements[i]->getBoundingPoint( j ).t ;
			  
			}
		}

		double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;
		
		Point pmin(minx,miny, minz, mint) ;
		Point pmax(maxx,maxy, maxz, maxt) ;
		
		Point rmin(xmin, ymin, zmin) ;
		Point rmax(xmax, ymax, zmax) ;
		
		  for ( size_t i = 0 ; i < elements.size() ; ++i )
		  {
			  if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
				  continue ;

			  std::vector<Point> id  ;

			  for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			  {
				  if( isOnBoundary( pos, elements[i]->getBoundingPoint( j ), pmin, pmax, tol ) && isInBoundary3D( elements[i]->getBoundingPoint( j ), rmin, rmax ) )
				  {
					  if ( cache3d.empty() || cache3d.back() != elements[i] )
					  {
						  cache.push_back( std::vector<Point>() );
						  cache3d.push_back( elements[i] );
					  }

					  cache.back().push_back( elements[i]->getBoundingPoint( j ) ) ;
				  }
			  }


			  if ( !cache3d.empty() && cache3d.back() == elements[i] )
			  {
				  GaussPointArray gp = elements[i]->getGaussPoints() ;
				  std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

				  for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
				  {
					  elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
				  }
				  if ( !function )
					  apply3DBC( elements[i],gp,Jinv, cache.back(), condition, data*getScale(), a , axis) ;
				  else
					  apply3DBC( elements[i],gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis ) ;
			  }
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache3d.size() ; ++i )
		{
			GaussPointArray gp = cache3d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache3d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}
			
			if ( !function )
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis) ;
			else
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis) ;
		}
	}
  
}

void BoundingBoxDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	if ( cache.empty() )
	{
		std::vector<DelaunayTetrahedron *> elements = t->getElements() ;

		if ( elements.empty() )
		{
			std::cout << "no elements in assembly" << std::endl ;
			return ;
		}

		double minx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double minz = elements.front()->getBoundingPoint( 0 ).z ;
		double maxz = elements.front()->getBoundingPoint( 0 ).z ;
		
		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;
		
		for ( size_t i = 0 ; i < elements.size() ; ++i )
		{
			if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() || elements[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				continue ;

			for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			{
				if ( elements[i]->getBoundingPoint( j ).x < minx )
					minx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).x > maxx )
					maxx = elements[i]->getBoundingPoint( j ).x ;

				if ( elements[i]->getBoundingPoint( j ).y < miny )
					miny = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).y > maxy )
					maxy = elements[i]->getBoundingPoint( j ).y ;

				if ( elements[i]->getBoundingPoint( j ).z < minz )
					minz = elements[i]->getBoundingPoint( j ).z ;

				if ( elements[i]->getBoundingPoint( j ).z > maxz )
					maxz = elements[i]->getBoundingPoint( j ).z ;

				if ( elements[i]->getBoundingPoint( j ).t < mint )
					mint = elements[i]->getBoundingPoint( j ).t ;

				if ( elements[i]->getBoundingPoint( j ).t > maxt )
					maxt = elements[i]->getBoundingPoint( j ).t ;
			  
			}
		}

		double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;
		
		Point pmin(minx,miny, minz, mint) ;
		Point pmax(maxx,maxy, maxz, maxt) ;
		
		  for ( size_t i = 0 ; i < elements.size() ; ++i )
		  {
			  if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
				  continue ;

			  std::vector<Point> id  ;

			  for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
			  {
				  if( isOnBoundary( pos, elements[i]->getBoundingPoint( j ), pmin, pmax, tol ) )
				  {
					  if ( cache3d.empty() || cache3d.back() != elements[i] )
					  {
						  cache.push_back( std::vector<Point>() );
						  cache3d.push_back( elements[i] );
					  }

					  cache.back().push_back( elements[i]->getBoundingPoint( j ) ) ;
				  }
			  }


			  if ( !cache3d.empty() && cache3d.back() == elements[i] )
			  {
				  GaussPointArray gp = elements[i]->getGaussPoints() ;
				  std::valarray<Matrix> Jinv( Matrix(), elements[i]->getGaussPoints().gaussPoints.size() ) ;

				  for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
				  {
					  elements[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
				  }
				  if ( !function )
					  apply3DBC( elements[i],gp,Jinv, cache.back(), condition, data*getScale(), a , axis) ;
				  else
					  apply3DBC( elements[i],gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis ) ;
			  }
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache3d.size() ; ++i )
		{
			GaussPointArray gp = cache3d[i]->getGaussPoints() ;
			std::valarray<Matrix> Jinv( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
			{
				cache3d[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
			}
			
			if ( !function )
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis) ;
			else
				apply3DBC( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis) ;
		}
	}  
}

BoundaryCondition::BoundaryCondition( LagrangeMultiplierType t, const double & d, int a  ) : scale( 1 ), condition( t ), data( d ), function( false ), axis(a) { } ;

BoundaryCondition::BoundaryCondition( LagrangeMultiplierType t, const Function & d, int a ) : scale( 1 ), condition( t ), dataFunction( d ), function( true ), axis(a) { } ;

void BoundaryCondition::setScale( double d )
{
	scale = d ;
}

double BoundaryCondition::getScale() const
{
	return scale ;
}

ProjectionDefinedBoundaryCondition::ProjectionDefinedBoundaryCondition( LagrangeMultiplierType t, const Point & dir, double d, int a ) : BoundaryCondition( t, d, a ), direction( dir ) { }

void ProjectionDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	std::vector<DelaunayTriangle *> tris = t->getElements() ;

	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		DelaunayTreeItem * VoidItem ;
		bool border = false ;

		for ( size_t j = 0 ; j < tris[i]->neighbour.size() ; j++ )
		{
			bool voidNeighbour = ( tris[i]->getNeighbour( j )->isTriangle
			                       && dynamic_cast<DelaunayTriangle *>( tris[i]->getNeighbour( j ) )->getBehaviour()->type == VOID_BEHAVIOUR ) ;
			border = border || tris[i]->getNeighbour( j )->isPlane
			         || voidNeighbour ;

			if ( voidNeighbour )
				VoidItem = tris[i]->getNeighbour( j ) ;

			if ( tris[i]->getNeighbour( j )->isPlane )
				VoidItem = tris[i]->getNeighbour( j ) ;
		}

		if ( tris[i]->getBehaviour()->type == VOID_BEHAVIOUR )
			border = false ;

		if ( border )
		{
			std::pair<Point *, Point*> commonSurface = tris[i]->commonEdge( VoidItem ) ;

			Segment ray(( tris[i]->getCenter() ), ( tris[i]->getCenter() ) - direction*( tris[i]->getRadius() ) ) ;
			bool isOnTheRightSide = ray.intersects( Segment( *commonSurface.first, *commonSurface.second ) ) ;

			if ( isOnTheRightSide )
			{
				std::vector<Point> id ;

				for ( size_t j = 0 ; j < tris[i]->getBoundingPoints().size() ; j++ )
				{

					Line side( tris[i]->getBoundingPoint( j ), tris[i]->getBoundingPoint( j ) - tris[i]->getBoundingPoint(( j + 1 ) % tris[i]->getBoundingPoints().size() ) ) ;

					if ( side.intersects( ray ) )
					{
						id.push_back( tris[i]->getBoundingPoint( j ) ) ;
						id.push_back( tris[i]->getBoundingPoint(( j + 1 ) % tris[i]->getBoundingPoints().size() ) ) ;
					}
				}

				GaussPointArray gp = tris[i]->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), tris[i]->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
				{
					tris[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
				}

				if ( !id.empty() )
				{
					if ( !function )
						apply2DBC( tris[i],gp,Jinv, id, condition, data*getScale(), a ) ;
					else
						apply2DBC( tris[i],gp,Jinv, id, condition, dataFunction*getScale(), a ) ;
				}
			}
		}
	}
}

void ProjectionDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<DelaunayTetrahedron *> tris = t->getElements() ;

	for ( size_t i = 0 ; i < tris.size() ; i++ )
	{
		std::vector<DelaunayDemiSpace *> space ;

		for ( size_t j = 0 ; j < tris[i]->neighbour.size() ; j++ )
		{
			if ( tris[i]->getNeighbour( j )->isSpace() )
			{
				space.push_back( static_cast<DelaunayDemiSpace *>( tris[i]->getNeighbour( j ) ) ) ;
			}
		}

		std::vector<Point> id ;

		for ( size_t s = 0 ; s < space.size() ; s++ )
		{
			Segment ray( tris[i]->getCenter(), tris[i]->getCenter() - direction*2.*tris[i]->getRadius() ) ;
			std::vector<Point *> points = space[s]->commonSurface( tris[i] ) ;
			Plane surf( *points[0], *points[1], *points[2] ) ;

			for ( size_t j = 3 ; j < tris[i]->getBoundingPoints().size() ; j++ )
			{
				if ( isCoplanar( *points[0], *points[1], *points[2], tris[i]->getBoundingPoint( j ) ) )
					points.push_back( &tris[i]->getBoundingPoint( j ) ) ;
			}

			if ( surf.intersects( ray ) && !points.empty() )
			{
				for ( size_t j = 0 ; j < points.size() ; j++ )
					id.push_back( *points[j] ) ;
			}
		}
		GaussPointArray gp = tris[i]->getGaussPoints() ;
		std::valarray<Matrix> Jinv( Matrix(), tris[i]->getGaussPoints().gaussPoints.size() ) ;

		for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
		{
			tris[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
		}
		
		if ( !function )
			apply3DBC( tris[i],gp,Jinv, id, condition, data*getScale(), a ) ;
		else
			apply3DBC( tris[i],gp,Jinv, id, condition, dataFunction*getScale(), a ) ;
	}
}


TimeContinuityBoundaryCondition::TimeContinuityBoundaryCondition() : BoundaryCondition( GENERAL, 0. ) { 
	previousDisp.resize(0) ;
} ;

void TimeContinuityBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	t->getAdditionalPoints() ;
	std::vector<DelaunayTriangle *> tri = t->getElements() ;
	std::vector<Point> id ;
	size_t timePlanes = tri[0]->timePlanes() ;

	if ( timePlanes < 2 )
		return ;

	size_t firstTimePlane = tri[0]->getBoundingPoints().size() / timePlanes ;

	size_t lastTimePlane = tri[0]->getBoundingPoints().size() * ( timePlanes - 1 ) / timePlanes ;

	size_t dof = tri[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;

	previousDisp.resize( a->getDisplacements().size()) ;
	previousDisp = a->getDisplacements() ;

	GaussPointArray gp = tri[0]->getGaussPoints() ;
	std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size() ) ;
	
	
	if(correspondance.size() == 0)
	{
		std::valarray<bool> done( false, tri.size()*tri[0]->getBoundingPoints().size() ) ;		
		for ( size_t i = 0 ; i < tri.size() ; i++ )
		{
			if(i%1000 == 0)
				std::cerr << "\r getting nodes for time continuity boundary condition 0 (triangle "  << i << "/" << tri.size() << ")" << std::flush ;
			for ( size_t k = 0 ; k < firstTimePlane ; k++ )
			{
				size_t id_ = tri[i]->getBoundingPoint(k).id ;
				if(!done[id_])
				{
					size_t corresponding = tri[i]->getBoundingPoint( lastTimePlane + k ).id ;
					correspondance.push_back(std::make_pair(id_, corresponding)) ;
					done[id_] = true ;
				}
			}
		}
		std::cerr << "... done" << std::endl ;
	}
	
	if( previousDisp.size() == 0 )
	{
		for(size_t i = 0 ; i < correspondance.size() ; i++)
		{
			if(i%1000 == 0)
				std::cerr << "\r applying time continuity boundary condition 0 (point "  << i << "/" << correspondance.size() << ")" << std::flush ;
			for(size_t n = 0 ; n < dof ; n++)
			{
				a->setPointAlongIndexedAxis( n, 0., correspondance[i].first )  ;
			}
		}
	}
	else
	{
		for(size_t i = 0 ; i < correspondance.size() ; i++)
		{
			if(i%1000 == 0)
				std::cerr << "\r applying time continuity boundary condition (point "  << i << "/" << correspondance.size() << ")" << std::flush ;
			for(size_t n = 0 ; n < dof ; n++)
			{
				a->setPointAlongIndexedAxis( n, previousDisp[ correspondance[i].second * dof + n] , correspondance[i].first )  ;
			}
		}	  
	}
	std::cerr << "... done" << std::endl ;
	return ;
/*	
	if ( previousDisp.size() == 0 )
	{
		for ( size_t i = 0 ; i < tri.size() ; i++ )
		{
			if(i%1000 == 0)
				std::cerr << "\r applying time continuity boundary condition 0 (triangle "  << i << "/" << tri.size() << ")" << std::flush ;
// 			GaussPointArray gp = tri[i]->getGaussPoints() ;
// 			std::valarray<Matrix> Jinv( Matrix(), tri[i]->getGaussPoints().gaussPoints.size() ) ;
// 			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
// 			{
// 				tri[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
// 			}
			
			for ( size_t k = 0 ; k < firstTimePlane ; k++ )
			{
				size_t id_ = tri[i]->getBoundingPoint(k).id ;
				if(!done[id_])
				{
					for(size_t n = 0 ; n < dof ; n++)
					{
						a->setPointAlongIndexedAxis( n, 0., id_ )  ;
	//					apply2DBC( tri[i], gp, Jinv, id, SET_ALONG_INDEXED_AXIS, 0., a, n) ;
					}
					done[id_] = true ;
				}
// 				id.clear() ;
// 				id.push_back( tri[i]->getBoundingPoint( k ) ) ;
// 				size_t corresponding = tri[i]->getBoundingPoint( lastTimePlane + k ).id ;
				
			}
		}
		std::cerr << " ...done" << std::endl ;
		return ;
	}

	for ( size_t i = 0 ; i < tri.size() ; i++ )
	{
		if(i%1000 == 0)
			std::cerr << "\r applying time continuity boundary condition (triangle "  << i << "/" << tri.size() << ")" << std::flush ;
		for ( size_t k = 0 ; k < firstTimePlane ; k++ )
		{
			size_t id_ = tri[i]->getBoundingPoint(k).id ;
 			size_t corresponding = tri[i]->getBoundingPoint( lastTimePlane + k ).id ;

			if(!done[id_])
			{
				for(size_t n = 0 ; n < dof ; n++)
				{
					a->setPointAlongIndexedAxis( n, previousDisp[corresponding*dof+n], id_ )  ;
//					apply2DBC( tri[i], gp, Jinv, id, SET_ALONG_INDEXED_AXIS, 0., a, n) ;
				}
				done[id_] = true ;
			}
			
// 			GaussPointArray gp = tri[i]->getGaussPoints() ;
// 			std::valarray<Matrix> Jinv( Matrix(), tri[i]->getGaussPoints().gaussPoints.size() ) ;
// 
// 			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
// 			{
// 				tri[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
// 			}
			
// 			for(size_t n = 0 ; n < dof ; n++)
// 			{
// 					a->setPointAlongIndexedAxis( n,  previousDisp[corresponding*dof+n], tri[i]->getBoundingPoint(k).id )  ;
// //				apply2DBC( tri[i], gp, Jinv, id, SET_ALONG_INDEXED_AXIS, previousDisp[corresponding*dof+n]*getScale(), a, n) ;
// 			}
/*			
			apply2DBC( tri[i],gp,Jinv, id, SET_ALONG_XI, previousDisp[corresponding*dof]*getScale(), a ) ;
			apply2DBC( tri[i],gp,Jinv, id, SET_ALONG_ETA, previousDisp[corresponding*dof+1]*getScale(), a ) ;*/
				/*				apply2DBC(tri[i], id, SET_STRESS_XI, previousStress[0], a) ;
								apply2DBC(tri[i], id, SET_STRESS_ETA, previousStress[1], a) ;
								apply2DBC(tri[i], id, SET_STRESS_XI_ETA, previousStress[2], a) ;*/
// 		}
// 	}
// 	std::cerr << " ...done" << std::endl ;
}

void TimeContinuityBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	t->getAdditionalPoints() ;
	std::vector<DelaunayTetrahedron *> tri = t->getElements() ;
	std::vector<Point> id ;
	size_t timePlanes = tri[0]->timePlanes() ;

	if ( timePlanes < 2 )
		return ;

	size_t firstTimePlane = tri[0]->getBoundingPoints().size() / timePlanes ;

	size_t lastTimePlane = tri[0]->getBoundingPoints().size() * ( timePlanes - 1 ) / timePlanes ;

	size_t dof = tri[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;

	previousDisp.resize( a->getDisplacements().size()) ;
	previousDisp = a->getDisplacements() ;
	
	GaussPointArray gp = tri[0]->getGaussPoints() ;
	std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size() ) ;
	

	if ( previousDisp.size() == 0 )
	{
		for ( size_t i = 0 ; i < tri.size() ; i++ )
		{
			
			for ( size_t k = 0 ; k < firstTimePlane ; k++ )
			{
				id.clear() ;
				id.push_back( tri[i]->getBoundingPoint( k ) ) ;
				size_t corresponding = tri[i]->getBoundingPoint( lastTimePlane + k ).id ;

// 				GaussPointArray gp = tri[i]->getGaussPoints() ;
// 				std::valarray<Matrix> Jinv( Matrix(), tri[i]->getGaussPoints().gaussPoints.size() ) ;
// 				
// 				for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
// 				{
// 					tri[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
// 				}
				
				for(size_t n = 0 ; n < dof ; n++)
				{
					apply3DBC( tri[i], gp, Jinv, id, SET_ALONG_INDEXED_AXIS, 0., a, n) ;
				}
			}
		}
		return ;
	}

	for ( size_t i = 0 ; i < tri.size() ; i++ )
	{
		for ( size_t k = 0 ; k < firstTimePlane ; k++ )
		{
			id.clear() ;
			id.push_back( tri[i]->getBoundingPoint( k ) ) ;
			size_t corresponding = tri[i]->getBoundingPoint( lastTimePlane + k ).id ;
			
// 			GaussPointArray gp = tri[i]->getGaussPoints() ;
// 			std::valarray<Matrix> Jinv( Matrix(), tri[i]->getGaussPoints().gaussPoints.size() ) ;
// 
// 			for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
// 			{
// 				tri[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
// 			}
			
			for(size_t n = 0 ; n < dof ; n++)
			{
				apply3DBC( tri[i], gp, Jinv, id, SET_ALONG_INDEXED_AXIS, previousDisp[corresponding*dof+n]*getScale(), a, n) ;
			}
/*			
			apply2DBC( tri[i],gp,Jinv, id, SET_ALONG_XI, previousDisp[corresponding*dof]*getScale(), a ) ;
			apply2DBC( tri[i],gp,Jinv, id, SET_ALONG_ETA, previousDisp[corresponding*dof+1]*getScale(), a ) ;*/
				/*				apply2DBC(tri[i], id, SET_STRESS_XI, previousStress[0], a) ;
								apply2DBC(tri[i], id, SET_STRESS_ETA, previousStress[1], a) ;
								apply2DBC(tri[i], id, SET_STRESS_XI_ETA, previousStress[2], a) ;*/
		}
	}
}


