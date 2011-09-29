// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011

#include "boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"

using namespace Mu ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double d ) : BoundaryCondition( t, d ), pos( pos ) { } ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d ) : BoundaryCondition( t, d ), pos( pos ) { } ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, double d ) : BoundaryCondition( t, d ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( zm ), zmax( zp )
{

}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, const Function & d ) : BoundaryCondition( t, d ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( zm ), zmax( zp )
{

}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, double d ) : BoundaryCondition( t, d ), pos( pos ), nearest( p ) {} ;

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, const Function & d ) : BoundaryCondition( t, d ), pos( pos ), nearest( p ) {} ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double d ): BoundaryCondition( t, d ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( 0 ), zmax( 0 )
{

}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, const Function & d ): BoundaryCondition( t, d ), pos( pos ),  xmin( xm ), xmax( xp ), ymin( ym ), ymax( yp ), zmin( 0 ), zmax( 0 )
{

}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition( ElementarySurface * surface ) : BoundaryCondition( GENERAL, 0 ), surface( surface ), volume( NULL )
{
}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition( ElementaryVolume * volume ) : BoundaryCondition( GENERAL, 0 ), surface( NULL ), volume( volume )
{
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementarySurface * surface , size_t id, double d ) : BoundaryCondition( t, d ), id( id ), surface( surface ), volume( NULL )
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementarySurface * surface , size_t id, const Function & d ) : BoundaryCondition( t, d ), id( id ), surface( surface ), volume( NULL )
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementaryVolume * volume , size_t id, double d ) : BoundaryCondition( t, d ), id( id ), surface( NULL ), volume( volume )
{

}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition( LagrangeMultiplierType t, ElementaryVolume * volume , size_t id, const Function & d ) : BoundaryCondition( t, d ), id( id ), surface( NULL ), volume( volume )
{

}

void apply2DBC( ElementarySurface *e,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a )
{
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;

	double nTimePlanes = 1 ;
	if(e->getOrder() > CONSTANT_TIME_LINEAR)
	{
		nTimePlanes = e->timePlanes() ;
	}
	
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][2][2] ;
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
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
					for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
					{
						if ( id[j] == e->getBoundingPoint( k ).id )
						{
							shapeFunctions.push_back( e->getShapeFunction( k ) ) ;
						}
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][2][2] ;
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][2][2] ;
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
				}

				return ;
			}

			default:
				break;
		}
	}
}

void apply3DBC( ElementaryVolume *e,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a )
{
// 	std::cout << "splash" << std::endl ;
// 	std::cout << (size_t)(e) << std::endl ;
// 	std::cout << (size_t)(e->getBehaviour()) << std::endl ;
// 	std::cout << (size_t)(e->getBehaviour()->type) << std::endl ;
	if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
		return ;

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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
					a->addForceOn( XI, forces[0], id[i] ) ;
					a->addForceOn( ETA, forces[1], id[i] ) ;
					a->addForceOn( ZETA, forces[2], id[i] ) ;
				}

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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
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
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( GradientDot( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) + VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					if(e->getOrder() >= CONSTANT_TIME_LINEAR)
						forces *= Jinv[0][3][3] ;
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

void apply2DBC( ElementarySurface *e,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a )
{
	std::vector<size_t> ids ;

	for ( size_t i = 0 ; i < id.size() ; i++ )
		ids.push_back( id[i].id );

	apply2DBC( e, ids, condition, data, a ) ;
}

void apply3DBC( ElementaryVolume *e,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a )
{
	std::vector<size_t> ids ;

	for ( size_t i = 0 ; i < id.size() ; i++ )
		ids.push_back( id[i].id );

	apply3DBC( e, ids, condition, data, a ) ;
}

void apply2DBC( ElementarySurface *e,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a )
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
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector imposed( 3 ) ;
					imposed[0] = vm.eval( data, id[i] ) ;
					imposed[1] = 0 ;
					imposed[2] = 0 ;
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;

				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				VirtualMachine vm ;

				Point p1( 0, 0, 0, -1 ) ;
				Point p2( 0, 0, 0, 0 ) ;
				Point p3( 0, 0, 0, 1 ) ;
				Point p4( 1, 0, 0, -1 ) ;
				Point p5( 1, 0, 0, 0 ) ;
				Point p6( 1, 0, 0, 1 ) ;
				Point p7( 0, 1, 0, -1 ) ;
				Point p8( 0, 1, 0, 0 ) ;
				Point p9( 0, 1, 0, 1 ) ;

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					imposed[0] = 0 ;
					imposed[1] = vm.eval( data, id[i] ) ;
					imposed[2] = 0 ;
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				}

				std::vector<Variable> v( 2 ) ;

				v[0] = XI ;
				v[1] = ETA ;
				if(e->getOrder() >= CONSTANT_TIME_LINEAR)
				{
					v.push_back(TIME_VARIABLE) ;
				}
				Vector imposed( 3 ) ;

				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					imposed[0] = 0 ;
					imposed[1] = 0 ;
					imposed[2] = vm.eval( data, id[i] ) ;
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i] ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
				}

				return ;
			}

			default:
				break;
		}
	}
}

void apply3DBC( ElementaryVolume *e,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a )
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

				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					imposed[0] = vm.eval( data, id[i] ) ;
					imposed[1] = 0 ;
					imposed[2] = 0 ;
					imposed[3] = 0 ;
					imposed[4] = 0 ;
					imposed[5] = 0 ;
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				imposed[1] = vm.eval( data, id[i] ) ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

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
				imposed[2] = vm.eval( data, id[i] ) ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				imposed[3] = vm.eval( data, id[i] ) ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				imposed[4] = vm.eval( data, id[i] ) ;
				imposed[5] = 0 ;
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
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
				imposed[5] = vm.eval( data, id[i] ) ;
				std::valarray<Matrix> Jinv( Matrix(), e->getGaussPoints().gaussPoints.size() ) ;

				for ( size_t i = 0 ; i < e->getGaussPoints().gaussPoints.size() ; i++ )
				{
					e->getInverseJacobianMatrix( e->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
				}

				for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
				{
					Vector forces =  VirtualMachine().ieval( Gradient( shapeFunctions[i], true ) * ( imposed ), e->getGaussPoints(), Jinv, v ) ;
					a->addForceOn( XI, forces[0], id[i].id ) ;
					a->addForceOn( ETA, forces[1], id[i].id ) ;
					a->addForceOn( ZETA, forces[2], id[i].id ) ;
				}

				return ;
			}

			default:
				break;
		}
	}
}


void DofDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
	if ( !function )
	{
		std::vector<size_t> id_ ;
		id_.push_back( id );
		apply2DBC( surface, id_, condition, data*getScale(), a ) ;
	}
	else
	{
		std::vector<Point> id_ ;

		for ( int i = 0 ; i < surface->getBoundingPoints().size() ; i++ )
		{
			if ( surface->getBoundingPoint( i ).id == id )
			{
				id_.push_back( surface->getBoundingPoint( i ) );
				apply2DBC( surface, id_, condition, dataFunction*getScale(), a ) ;
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
		apply3DBC( volume,  id_, condition, data*getScale(),  a ) ;
	}
	else
	{
		std::vector<Point> id_ ;

		for ( int i = 0 ; i < volume->getBoundingPoints().size() ; i++ )
		{
			if ( volume->getBoundingPoint( i ).id == id )
			{
				id_.push_back( volume->getBoundingPoint( i ) );
				apply3DBC( volume, id_, condition, dataFunction*getScale(), a ) ;
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

		Vector disps = surface->getState().getDisplacements( local, true ) ;
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

		Vector disps = volume->getState().getDisplacements( local, true ) ;
		a->setPointAlong( XI, disps[0], ( *i )->id ) ;
		a->setPointAlong( ETA, disps[1], ( *i )->id ) ;
		a->setPointAlong( ZETA, disps[2], ( *i )->id ) ;
	}
};

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

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

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
				}
			}
		}

		double tol = std::min( maxx - minx, maxy - miny ) * .0001 ;

		std::map<double, std::pair<Point, DelaunayTriangle*> > id  ;

		switch ( pos )
		{

			case TOP:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;

			}

			case LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case BOTTOM:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case TOP_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case TOP_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case BOTTOM_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			case BOTTOM_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
						{
							id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
						}
					}
				}

				std::vector<Point> target ;

				target.push_back( id.begin()->second.first ) ;
				cache2d.push_back( id.begin()->second.second ) ;
				cache.push_back( target ) ;

				if ( !function )
					apply2DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
				else
					apply2DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

				break ;
			}

			default:
			{
				break;
			}
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			if ( !function )
				apply2DBC( cache2d[i], cache[i], condition, data*getScale(), a ) ;
			else
				apply2DBC( cache2d[i], cache[i], condition, dataFunction*getScale(), a ) ;
		}
	}
}

void BoundingBoxNearestNodeDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<DelaunayTetrahedron *> elements = t->getElements() ;

	if ( elements.empty() )
	{
		std::cout << "no elements in assembly" << std::endl ;
		return ;
	}

	double minx = elements.front()->getBoundingPoint( 0 ).x ;

	double miny = elements.front()->getBoundingPoint( 0 ).y ;
	double minz = elements.front()->getBoundingPoint( 0 ).z ;
	double maxx = elements.front()->getBoundingPoint( 0 ).x ;
	double maxy = elements.front()->getBoundingPoint( 0 ).y ;
	double maxz = elements.front()->getBoundingPoint( 0 ).z ;

	for ( size_t i = 0 ; i < elements.size() ; ++i )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
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

			if ( elements[i]->getBoundingPoint( j ).y < minz )
				minz = elements[i]->getBoundingPoint( j ).z ;

			if ( elements[i]->getBoundingPoint( j ).y > maxz )
				maxz = elements[i]->getBoundingPoint( j ).z ;
		}
	}

	double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

	std::map<double, std::pair<Point, ElementaryVolume *> > id  ;

	switch ( pos )
	{

		case TOP:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case FRONT_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case FRONT_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;
			apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			break ;
		}

		case BACK_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BACK_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case FRONT_TOP:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;
			apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			break ;
		}

		case FRONT_BOTTOM:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_LEFT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_LEFT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_LEFT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_LEFT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_RIGHT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case TOP_RIGHT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_RIGHT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		case BOTTOM_RIGHT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id[dist( elements[i]->getBoundingPoint( j ), nearest )] = std::make_pair( elements[i]->getBoundingPoint( j ), elements[i] ) ;
					}
				}
			}

			std::vector<Point> target ;

			target.push_back( id.begin()->second.first ) ;

			if ( !function )
				apply3DBC( id.begin()->second.second, target, condition, data*getScale(), a ) ;
			else
				apply3DBC( id.begin()->second.second, target, condition, dataFunction*getScale(), a ) ;

			break ;
		}

		default:
		{
			break;
		}
	}
}

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition( LagrangeMultiplierType t, Geometry * source, double d ) : BoundaryCondition( t, d ), domain( source ) { };

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition( LagrangeMultiplierType t, Geometry * source, const Function & d ) : BoundaryCondition( t, d ), domain( source ) { };

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

		if ( !function )
			apply2DBC( elements[i], id, condition, data*getScale(), a ) ;
		else
			apply2DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
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

		if ( !function )
			apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
		else
			apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
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

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

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
				}
			}
		}

		double tol = std::min( maxx - minx, maxy - miny ) * .0001 ;

		switch ( pos )
		{

			case TOP:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;

				break ;
			}

			case TOP_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					std::vector<Point> id  ;

					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						   )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			default:
			{
				break;
			}
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			if ( !function )
				apply2DBC( cache2d[i], cache[i], condition, data*getScale(), a ) ;
			else
				apply2DBC( cache2d[i], cache[i], condition, dataFunction*getScale(), a ) ;
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

		double miny = elements.front()->getBoundingPoint( 0 ).y ;
		double maxx = elements.front()->getBoundingPoint( 0 ).x ;
		double maxy = elements.front()->getBoundingPoint( 0 ).y ;

		double mint = elements.front()->getBoundingPoint( 0 ).t ;
		double maxt = elements.front()->getBoundingPoint( 0 ).t ;

		if ( elements.front()->getOrder() >= CONSTANT_TIME_LINEAR )
		{
			for ( size_t j = 0 ; j < elements.front()->getBoundingPoints().size() ; ++j )
			{
				if ( elements.front()->getBoundingPoint( j ).t < mint )
					mint = elements.front()->getBoundingPoint( j ).t ;

				if ( elements.front()->getBoundingPoint( j ).t > maxt )
					maxt = elements.front()->getBoundingPoint( j ).t ;
			}
		}

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
				}
			}
		}

		double tol = std::max( std::min( maxx - minx, maxy - miny ) * .001, POINT_TOLERANCE_2D ) ;

		switch ( pos )
		{

			case TOP:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}
				break ;
			}

			case LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BEFORE:
			{

				if ( maxt != mint )
				{
					tol = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if ( std::abs( elements[i]->getBoundingPoint( j ).t - mint ) < tol )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}
				}

				break ;
			}

			case NOW:
			{

				if ( maxt != mint )
				{
					tol = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if ( std::abs( elements[i]->getBoundingPoint( j ).t - (mint+maxt)*0.5 ) < tol )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}
				}

				break ;
			}

			case AFTER:
			{
				if ( maxt != mint )
				{
					tol = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if ( std::abs( elements[i]->getBoundingPoint( j ).t - maxt ) < tol )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}
				}

				break ;
			}

			case TOP_BEFORE:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - mint ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

			case LEFT_BEFORE:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - mint ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

			case LEFT_NOW:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;
					
					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;
						
						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - (mint+maxt)*0.5 ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}
					
					break ;
				}
			}
			
			case LEFT_AFTER:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - maxt ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

		case TOP_NOW:
		{
			if ( maxt != mint )
			{
				double tolt = ( maxt - mint ) * 0.001 ;

				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if (( std::abs( elements[i]->getBoundingPoint( j ).t - (mint+maxt)*0.5 ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol ) )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}
		}

			case TOP_AFTER:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - maxt ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

			case BOTTOM_BEFORE:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - mint ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

		case BOTTOM_NOW:
		{
			if ( maxt != mint )
			{
				double tolt = ( maxt - mint ) * 0.001 ;

				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if (( std::abs( elements[i]->getBoundingPoint( j ).t - (mint+maxt)*0.5 ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol ) )
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
						if ( !function )
							apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}
		}

			case BOTTOM_AFTER:
			{
				if ( maxt != mint )
				{
					double tolt = ( maxt - mint ) * 0.001 ;

					for ( size_t i = 0 ; i < elements.size() ; ++i )
					{
						if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
							continue ;

						for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
						{
							if (( std::abs( elements[i]->getBoundingPoint( j ).t - maxt ) < tolt ) && ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol ) )
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
							if ( !function )
								apply2DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
							else
								apply2DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
						}
					}

					break ;
				}
			}

			default:
			{
				break;
			}
		}
	}
	else
	{
		for ( size_t i = 0 ; i < cache2d.size() ; ++i )
		{
			if ( !function )
				apply2DBC( cache2d[i], cache[i], condition, data*getScale(), a ) ;
			else
				apply2DBC( cache2d[i], cache[i], condition, dataFunction*getScale(), a ) ;
		}
	}
}

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<DelaunayTetrahedron *>  elements = t->getElements() ;
	double minx = elements.front()->getBoundingPoint( 0 ).x ;
	double miny = elements.front()->getBoundingPoint( 0 ).y ;
	double minz = elements.front()->getBoundingPoint( 0 ).z ;
	double maxx = elements.front()->getBoundingPoint( 0 ).x ;
	double maxy = elements.front()->getBoundingPoint( 0 ).y ;
	double maxz = elements.front()->getBoundingPoint( 0 ).z ;

	for ( size_t i = 0 ; i < elements.size() ; ++i )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
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

			if ( elements[i]->getBoundingPoint( j ).y < minz )
				minz = elements[i]->getBoundingPoint( j ).z ;

			if ( elements[i]->getBoundingPoint( j ).y > maxz )
				maxz = elements[i]->getBoundingPoint( j ).z ;
		}
	}

	double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

	if ( cache3d.empty() )
	{
		switch ( pos )
		{

			case TOP:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case FRONT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BACK:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case FRONT_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case FRONT_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BACK_LEFT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BACK_RIGHT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case FRONT_TOP:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case FRONT_BOTTOM:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_LEFT_FRONT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_LEFT_BACK:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_LEFT_FRONT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_LEFT_BACK:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_RIGHT_FRONT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case TOP_RIGHT_BACK:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_RIGHT_FRONT:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			case BOTTOM_RIGHT_BACK:
			{
				for ( size_t i = 0 ; i < elements.size() ; ++i )
				{
					if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
						continue ;

					std::vector<Point> id  ;

					for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
					{
						if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
						        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol
						        && elements[i]->getBoundingPoint( j ).x >= xmin
						        && elements[i]->getBoundingPoint( j ).x <= xmax
						        && elements[i]->getBoundingPoint( j ).y >= ymin
						        && elements[i]->getBoundingPoint( j ).y <= ymax
						        && elements[i]->getBoundingPoint( j ).z >= zmin
						        && elements[i]->getBoundingPoint( j ).z <= zmax
						   )
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
						if ( !function )
							apply3DBC( elements[i], cache.back(), condition, data*getScale(), a ) ;
						else
							apply3DBC( elements[i], cache.back(), condition, dataFunction*getScale(), a ) ;
					}
				}

				break ;
			}

			default:
			{
				break;
			}
		}
	}

	else
	{
		for ( size_t i = 0 ; i < cache3d.size() ; ++i )
		{
			if ( !function )
				apply3DBC( cache3d[i], cache[i], condition, data*getScale(), a ) ;
			else
				apply3DBC( cache3d[i], cache[i], condition, dataFunction*getScale(), a ) ;
		}
	}
}

void BoundingBoxDefinedBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<ElementaryVolume *> & elements = a->getElements3d() ;
	double minx = elements.front()->getBoundingPoint( 0 ).x ;
	double miny = elements.front()->getBoundingPoint( 0 ).y ;
	double minz = elements.front()->getBoundingPoint( 0 ).z ;
	double maxx = elements.front()->getBoundingPoint( 0 ).x ;
	double maxy = elements.front()->getBoundingPoint( 0 ).y ;
	double maxz = elements.front()->getBoundingPoint( 0 ).z ;

	for ( size_t i = 0 ; i < elements.size() ; ++i )
	{
		if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
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

			if ( elements[i]->getBoundingPoint( j ).y < minz )
				minz = elements[i]->getBoundingPoint( j ).z ;

			if ( elements[i]->getBoundingPoint( j ).y > maxz )
				maxz = elements[i]->getBoundingPoint( j ).z ;
		}
	}

	double tol = std::min( std::min( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

	switch ( pos )
	{

		case TOP:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case FRONT_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case FRONT_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BACK_LEFT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BACK_RIGHT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case FRONT_TOP:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case FRONT_BOTTOM:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_LEFT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_LEFT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_LEFT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_LEFT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - minx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_RIGHT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case TOP_RIGHT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - maxy ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_RIGHT_FRONT:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - minz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		case BOTTOM_RIGHT_BACK:
		{
			for ( size_t i = 0 ; i < elements.size() ; ++i )
			{
				if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
					continue ;

				std::vector<Point> id  ;

				for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
				{
					if ( std::abs( elements[i]->getBoundingPoint( j ).x - maxx ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).y - miny ) < tol
					        && std::abs( elements[i]->getBoundingPoint( j ).z - maxz ) < tol )
					{
						id.push_back( elements[i]->getBoundingPoint( j ) ) ;
					}
				}

				if ( !function )
					apply3DBC( elements[i], id, condition, data*getScale(), a ) ;
				else
					apply3DBC( elements[i], id, condition, dataFunction*getScale(), a ) ;
			}

			break ;
		}

		default:
		{
			break;
		}
	}
}

BoundaryCondition::BoundaryCondition( LagrangeMultiplierType t, const double & d ) : scale( 1 ), condition( t ), data( d ), function( false ) { } ;

BoundaryCondition::BoundaryCondition( LagrangeMultiplierType t, const Function & d ) : scale( 1 ), condition( t ), dataFunction( d ), function( true ) { } ;

void BoundaryCondition::setScale( double d )
{
	scale = d ;
}

double BoundaryCondition::getScale() const
{
	return scale ;
}

ProjectionDefinedBoundaryCondition::ProjectionDefinedBoundaryCondition( LagrangeMultiplierType t, const Point & dir, double d ) : BoundaryCondition( t, d ), direction( dir ) { }

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

				if ( !id.empty() )
				{
					if ( !function )
						apply2DBC( tris[i], id, condition, data*getScale(), a ) ;
					else
						apply2DBC( tris[i], id, condition, dataFunction*getScale(), a ) ;
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

		if ( !function )
			apply3DBC( tris[i], id, condition, data*getScale(), a ) ;
		else
			apply3DBC( tris[i], id, condition, dataFunction*getScale(), a ) ;
	}
}


TimeContinuityBoundaryCondition::TimeContinuityBoundaryCondition() : BoundaryCondition( GENERAL, 0. ) { } ;

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

	Vector previousDisp ;

	previousDisp.resize( tri[0]->getState().getDisplacements().size() ) ;

	if ( previousDisp.size() == 0 )
		return ;

	Vector previousStress ;

	previousStress.resize( tri[0]->getState().getPreviousStress( Point( 0, 0, 0, 0 ) ).size() ) ;

	if ( previousStress.size() == 0 )
		return ;


	for ( size_t i = 0 ; i < tri.size() ; i++ )
	{
		previousDisp = tri[i]->getState().getDisplacements() ;

		for ( size_t tp = 0 ; tp < 1 ; tp++ )
		{
			for ( size_t k = 0 ; k < firstTimePlane ; k++ )
			{
				id.clear() ;
				id.push_back( tri[i]->getBoundingPoint( firstTimePlane*tp + k ) ) ;
				previousStress = tri[i]->getState().getStress( id[k], false ) ;
				apply2DBC( tri[i], id, SET_ALONG_XI, previousDisp[( lastTimePlane+k )*dof+0]*getScale(), a ) ;
				apply2DBC( tri[i], id, SET_ALONG_ETA, previousDisp[( lastTimePlane+k )*dof+1]*getScale(), a ) ;
				/*				apply2DBC(tri[i], id, SET_STRESS_XI, previousStress[0], a) ;
								apply2DBC(tri[i], id, SET_STRESS_ETA, previousStress[1], a) ;
								apply2DBC(tri[i], id, SET_STRESS_XI_ETA, previousStress[2], a) ;*/
			}
		}
	}
}

void TimeContinuityBoundaryCondition::apply( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
	std::vector<DelaunayTetrahedron *> tet = t->getElements() ;
	std::vector<Point> id ;
	size_t timePlanes = tet[0]->timePlanes() ;

	if ( timePlanes < 2 )
		return ;

	size_t firstTimePlane = tet[0]->getBoundingPoints().size() / timePlanes ;

	size_t lastTimePlane = tet[0]->getBoundingPoints().size() * ( timePlanes - 1 ) / timePlanes ;

	size_t dof = tet[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;

	Vector previousDisp ;

	previousDisp.resize( tet[0]->getState().getPreviousDisplacements().size() ) ;

	if ( previousDisp.size() == 0 )
		return ;

	std::valarray<LagrangeMultiplierType> disp( SET_ALONG_XI, 3 ) ;

	disp[1] = SET_ALONG_ETA ;

	disp[2] = SET_ALONG_ZETA ;

	for ( size_t i = 0 ; i < tet.size() ; i++ )
	{
		previousDisp = tet[i]->getState().getPreviousDisplacements() ;

		for ( size_t k = 0; k < firstTimePlane ; k++ )
		{
			id.clear() ;
			id.push_back( tet[i]->getBoundingPoint( k ) ) ;

			for ( size_t j = 0 ; j < dof ; j++ )
			{
				apply3DBC( tet[i], id, disp[j], previousDisp[( lastTimePlane+k )*dof+j]*getScale(), a ) ;
			}
		}
	}
}


