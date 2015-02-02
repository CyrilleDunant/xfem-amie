// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011

#include "boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../physics/viscoelasticity.h"

using namespace Amie ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ) { } ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ) { } ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( zm ), zmax ( zp )
{
}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( zm ), zmax ( zp )
{

}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ), nearest ( p ) {} ;

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ), nearest ( p ) {} ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( 0 ), zmax ( 0 )
{
}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( 0 ), zmax ( 0 )
{

}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition ( ElementarySurface * surface ) : BoundaryCondition ( GENERAL, 0 ), surface ( surface ), volume ( nullptr )
{
}

ElementDefinedBoundaryCondition::ElementDefinedBoundaryCondition ( ElementaryVolume * volume ) : BoundaryCondition ( GENERAL, 0 ), surface ( nullptr ), volume ( volume )
{
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition ( LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & g, const std::valarray<Matrix> & J , size_t id, double d, int a ) : BoundaryCondition ( t, d, a ), id ( id ), surface ( surface ), volume ( nullptr )
{
    gp = new GaussPointArray ( g ) ;
    Jinv = new std::valarray<Matrix> ( J.size() ) ;
    std::copy ( &J[0], &J[J.size()], & ( *this->Jinv ) [0] ) ;
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition ( LagrangeMultiplierType t, ElementarySurface * surface, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a ) : BoundaryCondition ( t, d, a ), id ( id ), surface ( surface ), volume ( nullptr )
{
    this->gp = new GaussPointArray ( gp ) ;
    this->Jinv = new std::valarray<Matrix> ( Jinv.size() ) ;
    std::copy ( &Jinv[0], &Jinv[Jinv.size()], & ( *this->Jinv ) [0] ) ;
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition ( LagrangeMultiplierType t, ElementaryVolume * volume, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, double d, int a ) : BoundaryCondition ( t, d, a ), id ( id ), surface ( nullptr ), volume ( volume )
{
    this->gp = new GaussPointArray ( gp ) ;
    this->Jinv = new std::valarray<Matrix> ( Jinv.size() ) ;
    std::copy ( &Jinv[0], &Jinv[Jinv.size()], & ( *this->Jinv ) [0] ) ;
}

DofDefinedBoundaryCondition::DofDefinedBoundaryCondition ( LagrangeMultiplierType t, ElementaryVolume * volume, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv , size_t id, const Function & d, int a ) : BoundaryCondition ( t, d, a ), id ( id ), surface ( nullptr ), volume ( volume )
{
    this->gp = new GaussPointArray ( gp ) ;
    this->Jinv = new std::valarray<Matrix> ( Jinv.size() ) ;
    std::copy ( &Jinv[0], &Jinv[Jinv.size()], & ( *this->Jinv ) [0] ) ;
}

void apply2DBC ( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
    {
        return ;
    }

    double nTimePlanes = 1 ;
    if ( e->getOrder() > CONSTANT_TIME_LINEAR )
    {
        nTimePlanes = e->timePlanes() ;
    }

    VirtualMachine vm ;
    int n = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;

    for ( size_t idit = 0 ; idit < id.size() ; idit++ )
    {
        switch ( condition )
        {

        case GENERAL :
            std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
            break ;

        case FIX_ALONG_XI:
	    if(n > 2)
	    {
		for(size_t ax = 0 ; ax < n/2 ; ax++)
	            a->setPointAlongIndexedAxis ( ax*2, 0, id[idit] ) ;		
		break ;
	    }
	    else
	        a->setPointAlong ( XI, 0, id[idit] ) ;
            break ;

        case SET_ALONG_XI:
	    if(n>2)
	            a->setPointAlongIndexedAxis ( 0, data, id[idit] ) ;
	    else	
                 a->setPointAlong ( XI, data, id[idit] ) ;
            break ;

        case FIX_ALONG_ETA:
        {
	    if(n > 2)
	    {
		for(size_t ax = 0 ; ax < n/2 ; ax++)
	            a->setPointAlongIndexedAxis ( ax*2+1, 0, id[idit] ) ;	
		break ;	
	    }
	    else
	            a->setPointAlong ( ETA, 0, id[idit] ) ;
            break ;
        }

        case SET_ALONG_ETA:
	    if(n>2)
	            a->setPointAlongIndexedAxis ( 1, data, id[idit] ) ;
	    else	
	            a->setPointAlong ( ETA, data, id[idit] ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
//				std::cout << axis << "\t" << id[idit] << "\t" << data << std::endl ;
            a->setPointAlongIndexedAxis ( axis, data, id[idit] ) ;
            break ;

        case SET_FORCE_XI:
        {
            if ( !e->getBehaviour()->fractured() )
            {
                a->addForceOn ( XI, data/nTimePlanes, id[idit] ) ;
            }
            break ;
        }

        case SET_FORCE_ETA:
        {
            if ( !e->getBehaviour()->fractured() )
            {
                a->addForceOn ( ETA, data/nTimePlanes, id[idit] ) ;
            }
            break ;
        }

        case SET_VOLUMIC_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 3 ) ;
            imposed[0] = data ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;

                a->addForceOn ( XI, forces[0], id[idit] ) ;
                a->addForceOn ( ETA, forces[1], id[idit] ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model != PURE_ELASTICITY )
                {
                    for ( size_t b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 2*b, -forces[0], id[idit] ) ;
                        a->addForceOnIndexedAxis ( 2*b+1, -forces[1], id[idit] ) ;
                    }
                }

            }

            return ;
        }

        case SET_VOLUMIC_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;
            imposed[2] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;

                a->addForceOn ( XI, forces[0], id[idit] ) ;
                a->addForceOn ( ETA, forces[1], id[idit] ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model != PURE_ELASTICITY )
                {
                    for ( size_t b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 2*b, -forces[0], id[idit] ) ;
                        a->addForceOnIndexedAxis ( 2*b+1, -forces[1], id[idit] ) ;
                    }
                }

            }

            return ;
        }

        case SET_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= edge.norm() *.5 ;
                if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                    gpe.gaussPoints[i].second *= .5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }


            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 3 ) ;
            imposed[0] = data ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

                a->addForceOn ( XI, forces[0], id[j] ) ;
                a->addForceOn ( ETA, forces[1], id[j] ) ;
            }

            return ;
        }

        case SET_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= edge.norm() *.5 ;
                if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                    gpe.gaussPoints[i].second *= .5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }


            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;
            imposed[2] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

                a->addForceOn ( XI, forces[0], id[j] ) ;
                a->addForceOn ( ETA, forces[1], id[j] ) ;
            }

            return ;
        }

        case SET_NORMAL_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                                if ( dist ( first, last ) < dist ( last, &e->getBoundingPoint ( i ) ) )
                                {
                                    first = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            double enorm = edge.norm()*.5  ;
             if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                enorm *= .5 ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= enorm ;
            }
            std::valarray<Matrix> Jinve ( Jinv[0], gpe.gaussPoints.size() ) ;
            if ( e->isMoved() )
            {
                for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
                {
                    e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
                }
            }

            std::vector<Variable> v ( 2 ) ;
            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector normalv = edge.normalv ( e->getCenter() ) ;
            double nangle = atan2 ( normalv[1], normalv[0] ) ;

            Vector imposedx ( 3 ) ;
            imposedx[0] = data ;
            imposedx[1] = 0;
            imposedx[2] = 0 ;


            double c = cos ( nangle ) ;
            double s = sin ( nangle ) ;
            Matrix nrot ( 3,3 ) ;
            nrot[0][0] = c*c ;
            nrot[0][1] = s*s ;
            nrot[0][2] = -2.*s*c ;
            nrot[1][0] = s*s ;
            nrot[1][1] = c*c ;
            nrot[1][2] = 2.*s*c ;
            nrot[2][0] = s*c ;
            nrot[2][1] = -s*c ;
            nrot[2][2] = c*c-s*s ;

            Vector istr ( 0., 3 ) ;

            istr = nrot* ( imposedx ) ;
            Vector test(2) ;
            test[0] = 0 ; test[1] = 1 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[j], gpe, Jinve, v, false, normalv ) ;

                a->addForceOn ( XI,  forces[0], id[j] ) ;
                a->addForceOn ( ETA, forces[1], id[j] ) ;

            }

            return ;
        }

        case SET_TANGENT_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= edge.norm() *.5 ;
                 if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                    gpe.gaussPoints[i].second *= .5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 2 ) ;
            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Point normal = edge.normal ( e->getCenter() ) ;

// 				bool sameSide = isOnTheSameSide(e->getCenter(), edge.midPoint()+normal, edge.first(), edge.second()) ;
// 				double sign = 1. ;
// 				if(!sameSide)
// 					sign = -1 ;
            double nangle = atan2 ( normal.getY(), normal.getX() ) ;
// 				Vector normalv(2) ;
// 				normalv[0] = normal.getX()*sign ;
// 				normalv[1] = normal.getY()*sign ;

            Vector imposedx ( 3 ) ;
            imposedx[0] = 0 ;
            imposedx[1] = data;
            imposedx[2] = 0 ;


            double c = cos ( nangle ) ;
            double s = sin ( nangle ) ;
            Matrix nrot ( 3,3 ) ;
            nrot[0][0] = c*c ;
            nrot[0][1] = s*s ;
            nrot[0][2] = -2.*s*c ;
            nrot[1][0] = s*s ;
            nrot[1][1] = c*c ;
            nrot[1][2] = 2.*s*c ;
            nrot[2][0] = s*c ;
            nrot[2][1] = -s*c ;
            nrot[2][2] = c*c-s*s ;

            Vector istr ( 0., 3 ) ;

            istr = nrot* ( imposedx ) ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[j], gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

                a->addForceOn ( XI,  forces[0], id[j] ) ;
                a->addForceOn ( ETA, forces[1], id[j] ) ;
            }

            return ;
        }

        case SET_FLUX_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= edge.norm() /2 ;
                 if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                    gpe.gaussPoints[i].second *= .5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 2 ) ;
            imposed[0] = data ;
            imposed[1] = 0 ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gpe, Jinve, v ) ;
                a->addForceOn ( XI, f, id[i] ) ;
            }

            return ;
        }

        case SET_FLUX_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                gpe.gaussPoints[i].second *= edge.norm() /2 ;
                 if(e->getOrder() >= CONSTANT_TIME_LINEAR)
                    gpe.gaussPoints[i].second *= .5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 2 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gpe, Jinve, v ) ;
                a->addForceOn ( XI, f, id[i] ) ;
            }

            return ;
        }

        default:
            break;
        }
    }
}

void apply3DBC ( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<size_t> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
    {
        return ;
    }

    VirtualMachine vm ;

    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        switch ( condition )
        {

        case GENERAL :
            std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
            break ;

        case FIX_ALONG_XI:
            a->setPointAlong ( XI, 0., id[i] ) ;
            break ;

        case SET_ALONG_XI:
            a->setPointAlong ( XI, data, id[i] ) ;
            break ;

        case FIX_ALONG_ETA:
            a->setPointAlong ( ETA, 0., id[i] ) ;
            break ;

        case SET_ALONG_ETA:
            a->setPointAlong ( ETA, data, id[i] ) ;
            break ;

        case FIX_ALONG_ZETA:
            a->setPointAlong ( ZETA, 0., id[i] ) ;
            break ;

        case SET_ALONG_ZETA:
            a->setPointAlong ( ZETA, data, id[i] ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, data, id[i] ) ;
            break ;

        case SET_FORCE_XI:

            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( XI, data, id[i] ) ;

            break ;

        case SET_FORCE_ETA:
            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( ETA, data, id[i] ) ;

            break ;

        case SET_FORCE_ZETA:
            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( ZETA, data, id[i] ) ;

            break ;

        case SET_FLUX_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = data ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i] ) ;
            }

            return ;
        }

        case SET_FLUX_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;
            imposed[2] = 0 ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i] ) ;
            }

            return ;
        }

        case SET_FLUX_ZETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = 0 ;
            imposed[2] = data ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i] ) ;
            }

            return ;
        }

        case SET_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i )  ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i )  ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i )  ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 0., 6 ) ;
            imposed[0] = data ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false, std::abs ( edge.normalv() ) ) ;

                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }

            return ;
        }

        case SET_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }

                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 0.,  6 ) ;
            imposed[1] = data ;
            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false,  std::abs ( edge.normalv() ) ) ;

//					std::cout << "constant\t" << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;

                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }
//				exit(0) ;

            return ;
        }

        case SET_STRESS_ZETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 0., 6 ) ;
            imposed[2] = data ;


            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false, std::abs ( edge.normalv() ) ) ;

                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }

            return ;
        }

        case SET_NORMAL_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {
                     
                    if ( !&e->getBoundingPoint ( k ) )
                    {
                        continue ;
                    }

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( id[j] == e->getBoundingPoint ( k ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( k ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( k ) ;
                        }
                        else if ( !middle )
                        {
                            middle = &e->getBoundingPoint ( k ) ;
                        }
                        else
                        {
                            last = &e->getBoundingPoint ( k ) ;
                        }
                    }
                }
                for ( size_t k = 0 ; k < e->getEnrichmentFunctions().size() ; k++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( k ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( k ) ) ;
                    }
                }
            }
            if ( !last )
            {
                return ;
            }

            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( Matrix(3+e->getOrder() >= CONSTANT_TIME_LINEAR ,3+e->getOrder() >= CONSTANT_TIME_LINEAR),gpe.gaussPoints.size()) ;

            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector normal = edge.normalv ( e->getCenter() ) ;
            Vector imposed ( 0., 6 ) ;
            imposed[0] = data ;

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE_2D )
            {
                np0 = Point ( -normal[2], normal[0], -normal[1] ) ;
            }

            Point np1 = np^np0 ;
            np0 = np1^np ;

            double lx = normal[0] ;
            double ly = np0.getX() ;
            double lz = np1.getX() ;

            double rx = normal[1] ;
            double ry = np0.getY() ;
            double rz = np1.getY() ;

            double tx = normal[2] ;
            double ty = np0.getZ() ;
            double tz = np1.getZ() ;

            Matrix transform ( 6,6 ) ;
            transform[0][0] = lx*lx ;
            transform[0][1] = ly*ly ;
            transform[0][2] = lz*lz ;
            transform[0][3] = lx*ly ;
            transform[0][4] = lx*lz ;
            transform[0][5] = ly*lz ;
            transform[1][0] = rx*rx ;
            transform[1][1] = ry*ry ;
            transform[1][2] = rz*rz ;
            transform[1][3] = rx*ry ;
            transform[1][4] = rx*rz ;
            transform[1][5] = ry*rz ;
            transform[2][0] = tx*tx ;
            transform[2][1] = ty*ty ;
            transform[2][2] = tz*tz ;
            transform[2][3] = tx*ty ;
            transform[2][4] = tx*tz ;
            transform[2][5] = ty*tz ;
            transform[3][0] = 2.*lx*rx ;
            transform[3][1] = 2.*ly*ry ;
            transform[3][5] = 2.*lz*rz ;
            transform[3][3] = lx*ry+ly*rx ;
            transform[3][4] = lz*rx+lx*rz ;
            transform[3][5] = ly*rz+lz*ry ;
            transform[4][0] = 2.*lx*tx ;
            transform[4][1] = 2.*ly*ty ;
            transform[4][2] = 2.*lz*tz ;
            transform[4][3] = tx*ly+ly*lx ;
            transform[4][4] = tz*lx+tx*lz ;
            transform[4][5] = ty*lz+tz*ly ;
            transform[5][0] = 2.*rx*tx ;
            transform[5][1] = 2.*ry*ty ;
            transform[5][2] = 2.*rz*tz ;
            transform[5][3] = rx*ty+ry*tx ;
            transform[5][4] = rz*tx+rx*tz ;
            transform[5][5] = ry*tz+rz*ty ;

            Vector istr ( 0., 6 ) ;

            istr = transform * ( imposed ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normal ) ;
                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
//                 std::cout << forces[0] << "  " << forces[1] << "  " << forces[2] << std::endl ;
            }
            return ;
        }

        case SET_TANGENT_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t i = 0 ; i < gpe.gaussPoints.size() ; i++ )
            {
                gpe.gaussPoints[i].first = e->inLocalCoordinates ( gpe.gaussPoints[i].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[i].first, Jinve[i] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector normal = edge.normalv ( e->getCenter() ) ;
            Vector imposed ( 6 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;
            imposed[2] = 0 ;
            imposed[3] = 0 ;
            imposed[4] = 0 ;
            imposed[5] = 0 ;

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE_2D )
            {
                np0 = Point ( -normal[2], normal[0], -normal[1] ) ;
            }

            Point np1 = np^np0 ;
            np0 = np1^np ;

            double lx = normal[0] ;
            double ly = np0.getX() ;
            double lz = np1.getX() ;

            double rx = normal[1] ;
            double ry = np0.getY() ;
            double rz = np1.getY() ;

            double tx = normal[2] ;
            double ty = np0.getZ() ;
            double tz = np1.getZ() ;

            Matrix transform ( 6,6 ) ;
            transform[0][0] = lx*lx ;
            transform[0][1] = ly*ly ;
            transform[0][2] = lz*lz ;
            transform[0][3] = lx*ly ;
            transform[0][4] = lx*lz ;
            transform[0][5] = ly*lz ;
            transform[1][0] = rx*rx ;
            transform[1][1] = ry*ry ;
            transform[1][2] = rz*rz ;
            transform[1][3] = rx*ry ;
            transform[1][4] = rx*rz ;
            transform[1][5] = ry*rz ;
            transform[2][0] = tx*tx ;
            transform[2][1] = ty*ty ;
            transform[2][2] = tz*tz ;
            transform[2][3] = tx*ty ;
            transform[2][4] = tx*tz ;
            transform[2][5] = ty*tz ;
            transform[3][0] = 2.*lx*rx ;
            transform[3][1] = 2.*ly*ry ;
            transform[3][5] = 2.*lz*rz ;
            transform[3][3] = lx*ry+ly*rx ;
            transform[3][4] = lz*rx+lx*rz ;
            transform[3][5] = ly*rz+lz*ry ;
            transform[4][0] = 2.*lx*tx ;
            transform[4][1] = 2.*ly*ty ;
            transform[4][2] = 2.*lz*tz ;
            transform[4][3] = tx*ly+ly*lx ;
            transform[4][4] = tz*lx+tx*lz ;
            transform[4][5] = ty*lz+tz*ly ;
            transform[5][0] = 2.*rx*tx ;
            transform[5][1] = 2.*ry*ty ;
            transform[5][2] = 2.*rz*tz ;
            transform[5][3] = rx*ty+ry*tx ;
            transform[5][4] = rz*tx+rx*tz ;
            transform[5][5] = ry*tz+rz*ty ;

            Vector istr ( 0., 6 ) ;

            istr = transform* ( imposed ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normal ) ;

                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }

            return ;
        }

        default:
            break;
        }
    }
}

void apply2DBC ( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
    std::vector<size_t> ids ;

    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        ids.push_back ( id[i].getId() );
    }

    apply2DBC ( e,gp,Jinv, ids, condition, data, a, axis ) ;
}

void apply3DBC ( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a, int axis = 0 )
{
    std::vector<size_t> ids ;

    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        ids.push_back ( id[i].getId() );
    }

    apply3DBC ( e, gp, Jinv, ids, condition, data, a, axis ) ;
}

void apply2DBC ( ElementarySurface *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a, int axis = 0 )
{
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
    {
        return ;
    }

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
            a->setPointAlong ( XI, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_XI:
            a->setPointAlong ( XI, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ETA:
            a->setPointAlong ( ETA, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_ETA:
            a->setPointAlong ( ETA,  vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_FORCE_XI:

            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( XI,  vm.eval ( data, id[i] ), id[i].getId() ) ;

            break ;

        case SET_FORCE_ETA:
            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( ETA,  vm.eval ( data, id[i] ), id[i].getId() ) ;

            break ;

        case SET_VOLUMIC_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 0, 3, shapeFunctions[i], e, gp, Jinv, v, true ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc )
                {
                    for ( size_t b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 2*b, -forces[0], id[i].getId() ) ;
                        a->addForceOnIndexedAxis ( 2*b+1, -forces[1], id[i].getId() ) ;
                    }
                }
            }
            return ;
        }

        case SET_VOLUMIC_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 1, 3, shapeFunctions[i], e, gp, Jinv, v, true ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc )
                {
                    for ( size_t b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 2*b, -forces[0], id[i].getId() ) ;
                        a->addForceOnIndexedAxis ( 2*b+1, -forces[1], id[i].getId() ) ;
                    }
                }
            }

            return ;
        }

        case SET_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                gpe.gaussPoints[j].second *= edge.norm() *.5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }


            std::vector<Variable> v ( 2 ) ;
            Vector imposed ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 0, 3, shapeFunctions[i], e, gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
            }

            return ;
        }

        case SET_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

//				std::cout << vm.eval(data, id[i]) << std::endl ;
//				id[i].print() ;
            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                gpe.gaussPoints[j].second *= edge.norm() *.5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }



            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;

            VirtualMachine vm ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 1, 3, shapeFunctions[i], e,gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
            }

            return ;
        }

        case SET_NORMAL_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
//                 gpe.gaussPoints[j].second = edge.norm() *.5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }


            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Point normal = edge.normal() ;

            bool sameSide = isOnTheSameSide ( e->getCenter(), edge.midPoint() +normal, edge.first(), edge.second() ) ;
            double sign = 1. ;
            if ( !sameSide )
            {
                sign = -1 ;
            }
            double nangle = atan2 ( normal.getY() *sign, normal.getX() *sign ) ;
            Vector normalv ( 2 ) ;
            normalv[0] = normal.getX() *sign ;
            normalv[1] = normal.getY() *sign ;

            Vector imposedx ( 3 ) ;
            imposedx[0] = VirtualMachine().eval ( data, edge.midPoint() ) ;
            imposedx[1] = 0;
            imposedx[2] = 0 ;


            double c = cos ( nangle ) ;
            double s = sin ( nangle ) ;
            Matrix nrot ( 3,3 ) ;
            nrot[0][0] = c*c ;
            nrot[0][1] = s*s ;
            nrot[0][2] = -2.*s*c ;
            nrot[1][0] = s*s ;
            nrot[1][1] = c*c ;
            nrot[1][2] = 2.*s*c ;
            nrot[2][0] = s*c ;
            nrot[2][1] = -s*c ;
            nrot[2][2] = c*c-s*s ;

            Vector istr ( 0., 3 ) ;

            istr = nrot* ( imposedx ) ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normalv ) ;

                a->addForceOn ( XI, forces[0], id[j].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[j].getId() ) ;
            }

            return ;
        }

        case SET_TANGENT_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( i ) ) )
                                {
                                    last = &e->getBoundingPoint ( i ) ;
                                }
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *first, *last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;
            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                gpe.gaussPoints[j].second *= edge.norm() *.5 ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }


            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Point normal = edge.normal() ;

            bool sameSide = isOnTheSameSide ( e->getCenter(), edge.midPoint() +normal, edge.first(), edge.second() ) ;
            double sign = 1. ;
            if ( !sameSide )
            {
                sign = -1 ;
            }
            double nangle = atan2 ( normal.getY() *sign, normal.getX() *sign ) ;
            Vector normalv ( 2 ) ;
            normalv[0] = normal.getX() *sign ;
            normalv[1] = normal.getY() *sign ;

            Vector imposedx ( 3 ) ;
            imposedx[0] = 0;
            imposedx[1] = VirtualMachine().eval ( data, edge.midPoint() ) ;
            imposedx[2] = 0 ;


            double c = cos ( nangle ) ;
            double s = sin ( nangle ) ;
            Matrix nrot ( 3,3 ) ;
            nrot[0][0] = c*c ;
            nrot[0][1] = s*s ;
            nrot[0][2] = -2.*s*c ;
            nrot[1][0] = s*s ;
            nrot[1][1] = c*c ;
            nrot[1][2] = 2.*s*c ;
            nrot[2][0] = s*c ;
            nrot[2][1] = -s*c ;
            nrot[2][2] = c*c-s*s ;

            Vector istr ( 0., 3 ) ;

            istr = nrot* ( imposedx ) ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normalv ) ;

                a->addForceOn ( XI, forces[0], id[j].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[j].getId() ) ;
            }

            return ;
        }

        case SET_FLUX_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = vm.eval ( data, id[i] ) ; ;
            imposed[2] = 0 ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i].getId() ) ;
            }

            return ;
        }

        default:
            break;
        }
    }
}

void apply3DBC ( ElementaryVolume *e, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv,  const std::vector<Point> & id, LagrangeMultiplierType condition, const Function & data, Assembly * a, int axis = 0 )
{
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR )
    {
        return ;
    }

    VirtualMachine vm ;

    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        switch ( condition )
        {

        case GENERAL :
            std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
            break ;

        case FIX_ALONG_XI:
            a->setPointAlong ( XI, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_XI:
            a->setPointAlong ( XI, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ETA:
            a->setPointAlong ( ETA, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_ETA:
            a->setPointAlong ( ETA, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ZETA:
            a->setPointAlong ( ZETA, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_ZETA:
            a->setPointAlong ( ZETA, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_FORCE_XI:

            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( XI, vm.eval ( data, id[i] ), id[i].getId() ) ;

            break ;

        case SET_FORCE_ETA:
            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( ETA, vm.eval ( data, id[i] ), id[i].getId() ) ;

            break ;

        case SET_FORCE_ZETA:
            if ( e->getBehaviour()->fractured() )
            {
                break ;
            }

            a->setForceOn ( ZETA, vm.eval ( data, id[i] ), id[i].getId() ) ;

            break ;

        case SET_STRESS_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ) && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }

            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }
            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            Vector imposed ( 6 ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 0, 6, shapeFunctions[i], e,gpe, Jinve, v, false, edge.normalv() ) ;

//					std::cout << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;

            }

            return ;
        }

        case SET_STRESS_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ) && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }
            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 6 ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 1, 6, shapeFunctions[i], e,gpe, Jinve, v, false, edge.normalv() ) ;

//					std::cout << forces[0] << "\t" << forces[1] << "\t" << forces[2] << std::endl ;


                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;
            }
//				exit(0) ;

            return ;
        }

        case SET_STRESS_ZETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ) && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 2, 6, shapeFunctions[i], e,gpe, Jinve, v, false, edge.normalv() ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;
            }

            return ;
        }

        case SET_FLUX_XI:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = vm.eval ( data, id[i] ) ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;


            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i].getId() ) ;
            }

            return ;
        }

        case SET_TANGENT_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ) && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector normal = edge.normalv ( e->getCenter() ) ;
            Vector imposed ( 6 ) ;
            imposed[0] = 0 ;
            imposed[1] = vm.eval ( data, id[i] ) ;
            imposed[2] = 0 ;
            imposed[3] = 0 ;
            imposed[4] = 0 ;
            imposed[5] = 0 ;

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE_2D )
            {
                np0 = Point ( -normal[2], normal[0], -normal[1] ) ;
            }

            Point np1 = np^np0 ;
            np0 = np1^np ;

            double lx = normal[0] ;
            double ly = np0.getX() ;
            double lz = np1.getX() ;

            double rx = normal[1] ;
            double ry = np0.getY() ;
            double rz = np1.getY() ;

            double tx = normal[2] ;
            double ty = np0.getZ() ;
            double tz = np1.getZ() ;

            Matrix transform ( 6,6 ) ;
            transform[0][0] = lx*lx ;
            transform[0][1] = ly*ly ;
            transform[0][2] = lz*lz ;
            transform[0][3] = lx*ly ;
            transform[0][4] = lx*lz ;
            transform[0][5] = ly*lz ;
            transform[1][0] = rx*rx ;
            transform[1][1] = ry*ry ;
            transform[1][2] = rz*rz ;
            transform[1][3] = rx*ry ;
            transform[1][4] = rx*rz ;
            transform[1][5] = ry*rz ;
            transform[2][0] = tx*tx ;
            transform[2][1] = ty*ty ;
            transform[2][2] = tz*tz ;
            transform[2][3] = tx*ty ;
            transform[2][4] = tx*tz ;
            transform[2][5] = ty*tz ;
            transform[3][0] = 2.*lx*rx ;
            transform[3][1] = 2.*ly*ry ;
            transform[3][5] = 2.*lz*rz ;
            transform[3][3] = lx*ry+ly*rx ;
            transform[3][4] = lz*rx+lx*rz ;
            transform[3][5] = ly*rz+lz*ry ;
            transform[4][0] = 2.*lx*tx ;
            transform[4][1] = 2.*ly*ty ;
            transform[4][2] = 2.*lz*tz ;
            transform[4][3] = tx*ly+ly*lx ;
            transform[4][4] = tz*lx+tx*lz ;
            transform[4][5] = ty*lz+tz*ly ;
            transform[5][0] = 2.*rx*tx ;
            transform[5][1] = 2.*ry*ty ;
            transform[5][2] = 2.*rz*tz ;
            transform[5][3] = rx*ty+ry*tx ;
            transform[5][4] = rz*tx+rx*tz ;
            transform[5][5] = ry*tz+rz*ty ;

            Vector istr ( 0., 6 ) ;

            istr = transform* ( imposed ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normal ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;
            }

            return ;
        }

        case SET_NORMAL_STRESS:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j] == e->getBoundingPoint ( i ) )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j] == e->getBoundingPoint ( i ) && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE_3D  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE_3D )
                       )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( i ) ;
                        }
                        else
                        {
                            if ( !middle )
                            {
                                middle = &e->getBoundingPoint ( i ) ;
                            }
                            else
                            {
                                last = &e->getBoundingPoint ( i ) ;
                            }
                        }
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            TriPoint edge ( first, middle, last ) ;
            GaussPointArray gpe ( edge.getGaussPoints ( e->getOrder() >= CONSTANT_TIME_LINEAR ), -1 ) ;
            std::valarray<Matrix> Jinve ( gpe.gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gpe.gaussPoints.size() ; j++ )
            {
                gpe.gaussPoints[j].first = e->inLocalCoordinates ( gpe.gaussPoints[j].first ) ;
                e->getInverseJacobianMatrix ( gpe.gaussPoints[j].first, Jinve[j] ) ;
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector normal = edge.normalv ( e->getCenter() ) ;
            Vector imposed ( 6 ) ;
            imposed[0] = vm.eval ( data, id[i] ) ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;
            imposed[3] = 0 ;
            imposed[4] = 0 ;
            imposed[5] = 0 ;

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE_2D )
            {
                np0 = Point ( -normal[2], normal[0], -normal[1] ) ;
            }

            Point np1 = np^np0 ;
            np0 = np1^np ;

            double lx = normal[0] ;
            double ly = np0.getX() ;
            double lz = np1.getX() ;

            double rx = normal[1] ;
            double ry = np0.getY() ;
            double rz = np1.getY() ;

            double tx = normal[2] ;
            double ty = np0.getZ() ;
            double tz = np1.getZ() ;
            Matrix transform ( 6,6 ) ;

            transform[0][0] = lx*lx ;
            transform[0][1] = ly*ly ;
            transform[0][2] = lz*lz ;
            transform[0][3] = lx*ly ;
            transform[0][4] = lx*lz ;
            transform[0][5] = ly*lz ;
            transform[1][0] = rx*rx ;
            transform[1][1] = ry*ry ;
            transform[1][2] = rz*rz ;
            transform[1][3] = rx*ry ;
            transform[1][4] = rx*rz ;
            transform[1][5] = ry*rz ;
            transform[2][0] = tx*tx ;
            transform[2][1] = ty*ty ;
            transform[2][2] = tz*tz ;
            transform[2][3] = tx*ty ;
            transform[2][4] = tx*tz ;
            transform[2][5] = ty*tz ;
            transform[3][0] = 2.*lx*rx ;
            transform[3][1] = 2.*ly*ry ;
            transform[3][5] = 2.*lz*rz ;
            transform[3][3] = lx*ry+ly*rx ;
            transform[3][4] = lz*rx+lx*rz ;
            transform[3][5] = ly*rz+lz*ry ;
            transform[4][0] = 2.*lx*tx ;
            transform[4][1] = 2.*ly*ty ;
            transform[4][2] = 2.*lz*tz ;
            transform[4][3] = tx*ly+ly*lx ;
            transform[4][4] = tz*lx+tx*lz ;
            transform[4][5] = ty*lz+tz*ly ;
            transform[5][0] = 2.*rx*tx ;
            transform[5][1] = 2.*ry*ty ;
            transform[5][2] = 2.*rz*tz ;
            transform[5][3] = rx*ty+ry*tx ;
            transform[5][4] = rz*tx+rx*tz ;
            transform[5][5] = ry*tz+rz*ty ;

            Vector istr ( 0., 6 ) ;

            istr = transform* ( imposed ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[i], gpe, Jinve, v, false, normal ) ;

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;
            }

            return ;
        }

        case SET_FLUX_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = vm.eval ( data, id[i] ) ;
            imposed[2] = 0 ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i].getId() ) ;
            }

            return ;
        }

        case SET_FLUX_ZETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }
            }

            std::vector<Variable> v ( 3 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }
            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = 0 ;
            imposed[2] = vm.eval ( data, id[i] ) ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                double f = 0. ;
                f = vm.ieval ( VectorGradient ( shapeFunctions[i] ) * ( imposed ), gp, Jinv, v ) ;
                a->addForceOn ( XI, f, id[i].getId() ) ;
            }

            return ;
        }
        default:
            break;
        }
    }
}


void applyVerticalPlaneSections ( double topY, double bottomY, Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    std::vector<Point *> points ;

    for ( auto i  = t->begin() ; i != t->end() ; i++ )
    {
        for ( size_t j = 0 ; j <i->getBoundingPoints().size() ; j++ )
        {
            if ( i->getBoundingPoint ( j ).getY() <= topY && i->getBoundingPoint ( j ).getY() >= bottomY )
            {
                points.push_back ( &i->getBoundingPoint ( j ) );
            }
        }
    }
    std::sort ( points.begin(), points.end() ) ;
    auto e = std::unique ( points.begin(), points.end() ) ;
    points.erase ( e, points.end() ) ;

    for ( size_t i = 0 ; i <  points.size() ; i++ )
    {
        Point topPoint ( points[i]->getX(), topY ) ;
        Point bottomPoint ( points[i]->getX(), bottomY ) ;
        DelaunayTriangle * topElement = t->getUniqueConflictingElement ( &topPoint ) ;
        DelaunayTriangle * bottomElement = t->getUniqueConflictingElement ( &bottomPoint ) ;
        if ( topElement && bottomElement )
        {
            std::vector<int> idstop ;
            std::vector<double> coefficientstop ;
            idstop.push_back ( topElement->first->getId() *2 );
            coefficientstop.push_back ( ( ( topElement->second->getY()-topElement->third->getY() ) * ( topPoint.getX()-topElement->third->getX() )
                                          + ( topElement->third->getX()-topElement->second->getX() ) * ( topPoint.getY()-topElement->third->getY() ) )
                                        / ( ( topElement->second->getY()-topElement->third->getY() ) * ( topElement->first->getX()-topElement->third->getX() )
                                            + ( topElement->third->getX()-topElement->second->getX() ) * ( topElement->first->getY()-topElement->third->getY() ) ) );

            idstop.push_back ( topElement->second->getId() *2 );
            coefficientstop.push_back ( ( ( topElement->third->getY()-topElement->first->getY() ) * ( topPoint.getY()-topElement->third->getX() )
                                          + ( topElement->first->getX()-topElement->third->getX() ) * ( topPoint.getY()-topElement->third->getY() ) )
                                        / ( ( topElement->second->getY()-topElement->third->getY() ) * ( topElement->first->getX()-topElement->third->getX() )
                                            + ( topElement->third->getX()-topElement->second->getX() ) * ( topElement->first->getY()-topElement->third->getY() ) ) );

            idstop.push_back ( topElement->third->getId() *2 );
            coefficientstop.push_back ( 1.-coefficientstop[0]-coefficientstop[1] ) ;

            std::vector<int> idsbot ;
            std::vector<double> coefficientsbot ;
            idsbot.push_back ( bottomElement->first->getId() *2 );
            coefficientsbot.push_back ( ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomPoint.getX()-bottomElement->third->getX() )
                                          + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomPoint.getY()-topElement->third->getY() ) )
                                        / ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomElement->first->getX()-bottomElement->third->getX() )
                                            + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomElement->first->getY()-bottomElement->third->getY() ) ) );

            idsbot.push_back ( bottomElement->second->getId() *2 );
            coefficientsbot.push_back ( ( ( bottomElement->third->getY()-bottomElement->first->getY() ) * ( bottomPoint.getX()-bottomElement->third->getX() )
                                          + ( bottomElement->first->getX()-bottomElement->third->getX() ) * ( bottomPoint.getY()-bottomElement->third->getY() ) )
                                        / ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomElement->first->getX()-bottomElement->third->getX() )
                                            + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomElement->first->getY()-bottomElement->third->getY() ) ) );

            idsbot.push_back ( bottomElement->third->getId() *2 );
            coefficientsbot.push_back ( 1.-coefficientsbot[0]-coefficientsbot[1] ) ;

            std::valarray<unsigned int> allids ( 7 ) ;
            Vector coefs ( 7 ) ;
            for ( size_t j = 0 ; j < 3 ; j++ )
            {
                coefficientstop[j] *= ( points[i]->getY()-bottomY ) / ( topY-bottomY ) ;
                coefficientsbot[j] *= ( topY - points[i]->getY() ) / ( topY-bottomY ) ;

                coefs[j] = coefficientstop[j] ;
                coefs[j+3] = coefficientsbot[j] ;
            }
            coefs[6] = -std::accumulate ( &coefs[0], &coefs[6], double ( 0 ) ) ;
            allids[0] = topElement->first->getId() *2 ;
            allids[1] = topElement->second->getId() *2 ;
            allids[2] = topElement->third->getId() *2 ;
            allids[3] = bottomElement->first->getId() *2 ;
            allids[4] = bottomElement->second->getId() *2 ;
            allids[5] = bottomElement->third->getId() *2 ;
            allids[6] = points[i]->getId() *2 ;
            a->addMultiplier ( LagrangeMultiplier ( allids, coefs, 0 ) );
        }


    }

}


void applyHorizontalPlaneSections ( double topX, double bottomX, Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    std::vector<Point *> points ;

    for ( auto i  = t->begin() ; i != t->end() ; i++ )
    {
        for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
        {
            if ( i->getBoundingPoint ( j ).getX() <= topX && i->getBoundingPoint ( j ).getX() >= bottomX )
            {
                points.push_back ( &i->getBoundingPoint ( j ) );
            }
        }
    }
    std::sort ( points.begin(), points.end() ) ;
    auto e = std::unique ( points.begin(), points.end() ) ;
    points.erase ( e, points.end() ) ;

    for ( size_t i = 0 ; i <  points.size() ; i++ )
    {
        Point topPoint ( topX, points[i]->getY() ) ;
        Point bottomPoint ( bottomX, points[i]->getY() ) ;
        DelaunayTriangle * topElement = t->getUniqueConflictingElement ( &topPoint ) ;
        DelaunayTriangle * bottomElement = t->getUniqueConflictingElement ( &bottomPoint ) ;
        if ( topElement && bottomElement )
        {
            std::vector<int> idstop ;
            std::vector<double> coefficientstop ;
            idstop.push_back ( topElement->first->getId() *2+1 );
            coefficientstop.push_back ( ( ( topElement->second->getY()-topElement->third->getY() ) * ( topPoint.getX()-topElement->third->getX() )
                                          + ( topElement->third->getX()-topElement->second->getX() ) * ( topPoint.getY()-topElement->third->getY() ) )
                                        / ( ( topElement->second->getY()-topElement->third->getY() ) * ( topElement->first->getX()-topElement->third->getX() )
                                            + ( topElement->third->getX()-topElement->second->getX() ) * ( topElement->first->getY()-topElement->third->getY() ) ) );

            idstop.push_back ( topElement->second->getId() *2+1 );
            coefficientstop.push_back ( ( ( topElement->third->getY()-topElement->first->getY() ) * ( topPoint.getY()-topElement->third->getX() )
                                          + ( topElement->first->getX()-topElement->third->getX() ) * ( topPoint.getY()-topElement->third->getY() ) )
                                        / ( ( topElement->second->getY()-topElement->third->getY() ) * ( topElement->first->getX()-topElement->third->getX() )
                                            + ( topElement->third->getX()-topElement->second->getX() ) * ( topElement->first->getY()-topElement->third->getY() ) ) );

            idstop.push_back ( topElement->third->getId() *2+1 );
            coefficientstop.push_back ( 1.-coefficientstop[0]-coefficientstop[1] ) ;

            std::vector<int> idsbot ;
            std::vector<double> coefficientsbot ;
            idsbot.push_back ( bottomElement->first->getId() *2+1 );
            coefficientsbot.push_back ( ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomPoint.getX()-bottomElement->third->getX() )
                                          + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomPoint.getY()-topElement->third->getY() ) )
                                        / ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomElement->first->getX()-bottomElement->third->getX() )
                                            + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomElement->first->getY()-bottomElement->third->getY() ) ) );

            idsbot.push_back ( bottomElement->second->getId() *2+1 );
            coefficientsbot.push_back ( ( ( bottomElement->third->getY()-bottomElement->first->getY() ) * ( bottomPoint.getX()-bottomElement->third->getX() )
                                          + ( bottomElement->first->getX()-bottomElement->third->getX() ) * ( bottomPoint.getY()-bottomElement->third->getY() ) )
                                        / ( ( bottomElement->second->getY()-bottomElement->third->getY() ) * ( bottomElement->first->getX()-bottomElement->third->getX() )
                                            + ( bottomElement->third->getX()-bottomElement->second->getX() ) * ( bottomElement->first->getY()-bottomElement->third->getY() ) ) );

            idsbot.push_back ( bottomElement->third->getId() *2+1 );
            coefficientsbot.push_back ( 1.-coefficientsbot[0]-coefficientsbot[1] ) ;


            std::valarray<unsigned int> allids ( 7 ) ;
            Vector coefs ( 7 ) ;
            for ( size_t j = 0 ; j < 3 ; j++ )
            {
                coefficientstop[j] *= ( points[i]->getX()-bottomX ) / ( topX-bottomX ) ;
                coefficientsbot[j] *= ( topX - points[i]->getX() ) / ( topX-bottomX ) ;

                coefs[j] = coefficientstop[j] ;
                coefs[j+3] = coefficientsbot[j] ;
            }
            coefs[6] = -std::accumulate ( &coefs[0], &coefs[6], double ( 0 ) ) ;
            allids[0] = topElement->first->getId() *2+1 ;
            allids[1] = topElement->second->getId() *2+1 ;
            allids[2] = topElement->third->getId() *2+1 ;
            allids[3] = bottomElement->first->getId() *2+1 ;
            allids[4] = bottomElement->second->getId() *2+1 ;
            allids[5] = bottomElement->third->getId() *2+1 ;
            allids[6] = points[i]->getId() *2+1 ;
            a->addMultiplier ( LagrangeMultiplier ( allids, coefs, 0 ) );
        }


    }

}


void PlaneSectionsBoundaryConditions::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( !isVertical )
    {
        applyHorizontalPlaneSections ( uplimit, downlimit, a, t ) ;
    }
    else
    {
        applyVerticalPlaneSections ( uplimit, downlimit, a, t ) ;
    }
}

void PlaneSectionsBoundaryConditions::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{

}


void DofDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( !function )
    {
        std::vector<size_t> id_ ;
        id_.push_back ( id );
        apply2DBC ( surface,*gp, *Jinv, id_, condition, data*getScale(), a, axis ) ;
    }
    else
    {
        std::vector<Point> id_ ;

        for ( int i = 0 ; i < surface->getBoundingPoints().size() ; i++ )
        {
            if ( surface->getBoundingPoint ( i ).getId() == id )
            {
                id_.push_back ( surface->getBoundingPoint ( i ) );
                apply2DBC ( surface,*gp,*Jinv, id_, condition, dataFunction*getScale(), a, axis ) ;
                return ;
            }
        }
    }
}

void DofDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( surface )
    {
        return ;
    }

    if ( !function )
    {
        std::vector<size_t> id_ ;
        id_.push_back ( id );
        apply3DBC ( volume,*gp,*Jinv,  id_, condition, data*getScale(),  a , axis ) ;
    }
    else
    {
        std::vector<Point> id_ ;

        for ( int i = 0 ; i < volume->getBoundingPoints().size() ; i++ )
        {
            if ( volume->getBoundingPoint ( i ).getId() == id )
            {
                id_.push_back ( volume->getBoundingPoint ( i ) );
                apply3DBC ( volume,*gp, *Jinv, id_, condition, dataFunction*getScale(), a, axis ) ;
                return ;
            }
        }
    }
}

void ElementDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( volume )
    {
        return ;
    }

    std::set<Point *> points ;

    for ( auto i = t->begin() ; i != t->end() ; i++ )
    {
        if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
        {
            continue ;
        }

        for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
        {
            Point test ( i->getBoundingPoint ( j ) ) ;
            surface->project ( &test );

            if ( dist ( test, i->getBoundingPoint ( j ) ) < POINT_TOLERANCE_2D )
            {
                points.insert ( &i->getBoundingPoint ( j ) ) ;
            }
        }
    }

    for ( auto i = points.begin() ; i != points.end() ; i++ )
    {
        Point local = surface->inLocalCoordinates ( * ( *i ) ) ;
        std::vector<Point> p ;
        p.push_back ( local );

        Vector disps ( 0.,2 ) ;
        surface->getState().getField ( DISPLACEMENT_FIELD, local, disps, true ) ;
        a->setPointAlong ( XI, disps[0], ( *i )->getId() ) ;
        a->setPointAlong ( ETA, disps[1], ( *i )->getId() ) ;

    }

};

void ElementDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( surface )
    {
        return ;
    }

    std::set<Point *> points ;

    for ( auto i = t->begin() ; i != t->end() ; i++ )
    {
        if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
        {
            continue ;
        }

        for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
        {
            Point test ( i->getBoundingPoint ( j ) ) ;
            volume->project ( &test );

            if ( dist ( test, i->getBoundingPoint ( j ) ) < POINT_TOLERANCE_3D )
            {
                points.insert ( &i->getBoundingPoint ( j ) ) ;
            }
        }
    }

    for ( auto i = points.begin() ; i != points.end() ; i++ )
    {
        Point local = volume->inLocalCoordinates ( * ( *i ) ) ;
        std::vector<Point> p ;
        p.push_back ( local );

        Vector disps ( 0.,3 ) ;
        volume->getState().getField ( DISPLACEMENT_FIELD, local, disps, true ) ;
        a->setPointAlong ( XI, disps[0], ( *i )->getId() ) ;
        a->setPointAlong ( ETA, disps[1], ( *i )->getId() ) ;
        a->setPointAlong ( ZETA, disps[2], ( *i )->getId() ) ;
    }
};

bool isInBoundary2D ( Point test, Point min, Point max )
{
    return ( test.getX() >= min.getX() && test.getX() <= max.getX() && test.getY() >= min.getY() && test.getY() <= max.getY() ) ;
}

bool isInBoundary3D ( Point test, Point min, Point max )
{
    return ( test.getX() >= min.getX() && test.getX() <= max.getX() && test.getY() >= min.getY() && test.getY() <= max.getY() && test.getZ() >= min.getZ() && test.getZ() <= max.getZ() ) ;
}


bool isOnBoundary ( BoundingBoxPosition pos, Point & test, Point & min, Point & max , double tol )
{
    switch ( pos )
    {
        // 2D edges, 3D planes, 4D time planes
    case LEFT:
        return ( std::abs ( test.getX() - min.getX() ) < tol ) ;
    case RIGHT:
        return ( std::abs ( test.getX() - max.getX() ) < tol ) ;
    case BOTTOM:
        return ( std::abs ( test.getY() - min.getY() ) < tol ) ;
    case TOP:
        return ( std::abs ( test.getY() - max.getY() ) < tol ) ;
    case BACK:
        return ( std::abs ( test.getZ() - min.getZ() ) < tol ) ;
    case FRONT:
        return ( std::abs ( test.getZ() - max.getZ() ) < tol ) ;
    case BEFORE:
        return ( std::abs ( test.getT() - min.getT() ) < tol ) ;
    case AFTER:
        return ( std::abs ( test.getT() - max.getT() ) < tol ) ;

        // 2D corners, 3D edges, 4D planes
    case BOTTOM_LEFT:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) ) ;
    case TOP_LEFT:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) ) ;
    case TOP_RIGHT:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) ) ;

    case BACK_LEFT:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) ) ;
    case BACK_RIGHT:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) ) ;
    case FRONT_LEFT:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) ) ;
    case FRONT_RIGHT:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) ) ;

    case BOTTOM_BACK:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) ) ;
    case TOP_BACK:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) ) ;
    case FRONT_BOTTOM:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) ) ;
    case FRONT_TOP:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) ) ;

    case BOTTOM_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BOTTOM_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

    case BACK_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BACK_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

    case LEFT_BEFORE:
        return ( isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case LEFT_AFTER:
        return ( isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case RIGHT_BEFORE:
        return ( isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case RIGHT_AFTER:
        return ( isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

        // 3D corners, 4D lines
    case BOTTOM_LEFT_BACK:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) ) ;
    case BOTTOM_LEFT_FRONT:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_BACK:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_FRONT:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) ) ;
    case TOP_LEFT_BACK:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) ) ;
    case TOP_LEFT_FRONT:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) ) ;
    case TOP_RIGHT_BACK:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) ) ;
    case TOP_RIGHT_FRONT:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) ) ;

    case BOTTOM_LEFT_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_LEFT_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_RIGHT_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;

    case BACK_LEFT_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BACK_RIGHT_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_LEFT_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_RIGHT_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;

    case BOTTOM_BACK_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_BACK_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_BOTTOM_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_TOP_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;

    case BOTTOM_LEFT_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_LEFT_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_RIGHT_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

    case BACK_LEFT_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BACK_RIGHT_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_LEFT_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_RIGHT_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

    case BOTTOM_BACK_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_BACK_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_BOTTOM_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_TOP_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;

        // 4D corners
    case BOTTOM_LEFT_BACK_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BOTTOM_LEFT_FRONT_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_BACK_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_FRONT_BEFORE:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_LEFT_BACK_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_LEFT_FRONT_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_RIGHT_BACK_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_RIGHT_FRONT_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;

    }
    return false ;
}




void BoundingBoxNearestNodeDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{

    if ( cache.empty() )
    {

        if ( t->begin().size() == 0 )
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end() ; ++i )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }

            }
        }

        double tol = std::min ( maxx - minx, maxy - miny ) * .0001 ;

        Point pmin ( minx,miny, 0., mint ) ;
        Point pmax ( maxx,maxy, 0., maxt ) ;

        std::map<double, std::pair<Point, DelaunayTriangle*> > id  ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) )
                {
                    id[dist ( i->getBoundingPoint ( j ), nearest )] = std::make_pair ( i->getBoundingPoint ( j ), i ) ;
                }
            }
        }

        std::vector<Point> target ;

        target.push_back ( id.begin()->second.first ) ;
        cache2d.push_back ( id.begin()->second.second ) ;
        cache.push_back ( target ) ;
        GaussPointArray gp = id.begin()->second.second->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), id.begin()->second.second->getGaussPoints().gaussPoints.size() ) ;

        for ( size_t i = 0 ; i < id.begin()->second.second->getGaussPoints().gaussPoints.size() ; i++ )
        {
            id.begin()->second.second->getInverseJacobianMatrix ( id.begin()->second.second->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
        }

        if ( !function )
        {
            apply2DBC ( id.begin()->second.second,gp,Jinv, target, condition, data*getScale(), a , axis ) ;
        }
        else
        {
            apply2DBC ( id.begin()->second.second,gp,Jinv, target, condition, dataFunction*getScale(), a , axis ) ;
        }


    }
    else
    {
        for ( size_t i = 0 ; i < cache2d.size() ; ++i )
        {
            GaussPointArray gp = cache2d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache2d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            if ( !function )
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a, axis ) ;
            }
            else
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a, axis ) ;
            }
        }
    }
}

void BoundingBoxNearestNodeDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {

        if ( t->begin().size() == 0 )
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double minz = t->begin()->getBoundingPoint ( 0 ).getZ() ;
        double maxz = t->begin()->getBoundingPoint ( 0 ).getZ() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end() ; ++i )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() < minz )
                {
                    minz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() > maxz )
                {
                    maxz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }

            }
        }

        double tol = std::min ( std::min ( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

        Point pmin ( minx, miny, minz, mint ) ;
        Point pmax ( maxx, maxy, maxz, maxt ) ;

        std::map<double, std::pair<Point, DelaunayTetrahedron*> > id  ;

        for ( auto i = t->begin() ; i != t->end() ; ++i )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) )
                {
                    id[dist ( i->getBoundingPoint ( j ), nearest )] = std::make_pair ( i->getBoundingPoint ( j ), i ) ;
                }
            }
        }

        std::vector<Point> target ;

        target.push_back ( id.begin()->second.first ) ;
        cache3d.push_back ( id.begin()->second.second ) ;
        cache.push_back ( target ) ;
        GaussPointArray gp = id.begin()->second.second->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), id.begin()->second.second->getGaussPoints().gaussPoints.size() ) ;

        for ( size_t i = 0 ; i < id.begin()->second.second->getGaussPoints().gaussPoints.size() ; i++ )
        {
            id.begin()->second.second->getInverseJacobianMatrix ( id.begin()->second.second->getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
        }

        if ( !function )
        {
            apply3DBC ( id.begin()->second.second,gp,Jinv, target, condition, data*getScale(), a , axis ) ;
        }
        else
        {
            apply3DBC ( id.begin()->second.second,gp,Jinv, target, condition, dataFunction*getScale(), a, axis ) ;
        }


    }
    else
    {
        for ( size_t i = 0 ; i < cache3d.size() ; ++i )
        {
            GaussPointArray gp = cache3d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            if ( !function )
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }

}

GeometryDefinedSurfaceBoundaryCondition::GeometryDefinedSurfaceBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, double d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { };

GeometryDefinedSurfaceBoundaryCondition::GeometryDefinedSurfaceBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, const Function & d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { };

void GeometryDefinedSurfaceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( cache2d.empty() )
    {
        std::vector<DelaunayTriangle *> elements = t->getConflictingElements ( domain ) ;
        double tol = domain->getRadius() * .001 ;


        for (auto i = elements.begin() ; i != elements.end() ; i++ )
        {
            if ( (*i)->getBehaviour()->getDamageModel() && (*i)->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < (*i)->getBoundingPoints().size() ; ++j )
            {
                Point test = (*i)->getBoundingPoint ( j ) ;
                domain->project ( &test );

                if ( squareDist2D ( test, (*i)->getBoundingPoint ( j ) ) < tol*tol )
                {
                    id.push_back ( (*i)->getBoundingPoint ( j ) ) ;
                }
            }
            if ( !id.empty() )
            {
                cache2d.push_back ( (*i) );
                cache.push_back ( id );
            }

        }
    }

    for ( size_t i = 0 ; i < cache2d.size() ; i++ )
    {
        GaussPointArray gp = cache2d[i]->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;
        for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
        {
            cache2d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
        }

        if ( !function )
        {
            apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a ) ;
        }
        else
        {
            apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a ) ;
        }
    }
}

void GeometryDefinedSurfaceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{

    if ( cache3d.empty() )
    {
        std::vector<DelaunayTetrahedron *> elements = t->getConflictingElements ( domain ) ;
        double tol = domain->getRadius() * .001 ;
        for ( size_t i = 0 ; i < elements.size() ; ++i )
        {
            if ( elements[i]->getBehaviour()->getDamageModel() && elements[i]->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < elements[i]->getBoundingPoints().size() ; ++j )
            {
                Point test = elements[i]->getBoundingPoint ( j ) ;
                domain->project ( &test );

                if ( squareDist2D ( test, elements[i]->getBoundingPoint ( j ) ) < tol*tol )
                {
                    id.push_back ( elements[i]->getBoundingPoint ( j ) ) ;
                }
            }
            if ( !id.empty() )
            {
                cache3d.push_back ( elements[i] );
                cache.push_back ( id );
            }
        }
    }

    for ( size_t i = 0 ; i < cache3d.size() ; i++ )
    {
        GaussPointArray gp = cache3d[i]->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;
        for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
        {
            cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
        }

        if ( !function )
        {
            apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a ) ;
        }
        else
        {
            apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a ) ;
        }
    }
}


GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, double d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { };

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, const Function & d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { };

void GeometryDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( cache.empty() )
    {
        double tol = domain->getRadius() * .0001 ;


        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                Circle c ( tol, i->getBoundingPoint ( j ) ) ;

                if ( domain->intersects ( &c ) || domain->in ( i->getBoundingPoint ( j ) ) )
                {
                    id.push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }
            if ( id.empty() )
            {
                continue ;
            }
            else
            {
                cache2d.push_back ( i );
                cache.push_back ( id );
            }

            GaussPointArray gp = i->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }


            if ( !function )
            {
                apply2DBC ( i,gp, Jinv, id, condition, data*getScale(), a ) ;
            }
            else
            {
                apply2DBC ( i,gp, Jinv, id, condition, dataFunction*getScale(), a ) ;
            }
        }
    }
    else
    {
        for ( size_t i = 0 ; i < cache2d.size() ; ++i )
        {
            GaussPointArray gp = cache2d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache2d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }


            if ( !function )
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }
}

void GeometryDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {
        double tol = domain->getRadius() * .0001 ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                Sphere c ( tol, i->getBoundingPoint ( j ) ) ;

                if ( domain->intersects ( &c ) || domain->in ( i->getBoundingPoint ( j ) ) )
                {
                    id.push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }
            if ( id.empty() )
            {
                continue ;
            }
            else
            {
                cache3d.push_back ( i );
                cache.push_back ( id );
            }

            GaussPointArray gp = i->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            if ( !function )
            {
                apply3DBC ( i,gp,Jinv, id, condition, data*getScale(), a ) ;
            }
            else
            {
                apply3DBC ( i,gp, Jinv, id, condition, dataFunction*getScale(), a ) ;
            }
        }
    }
    else
    {
        for ( size_t i = 0 ; i < cache3d.size() ; ++i )
        {
            GaussPointArray gp = cache3d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }


            if ( !function )
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }
}

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( cache.empty() )
    {

        if ( t->begin().size() == 0 )
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }

            }
        }

        double tol = std::min ( maxx - minx, maxy - miny ) * .0001 ;

        Point pmin ( minx,miny, 0., mint ) ;
        Point pmax ( maxx,maxy, 0., maxt ) ;

        Point rmin ( xmin, ymin ) ;
        Point rmax ( xmax, ymax ) ;


        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {

                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) && isInBoundary2D ( i->getBoundingPoint ( j ), rmin, rmax ) )
                {
                    if ( cache2d.empty() || cache2d.back() != i )
                    {
                        cache.push_back ( std::vector<Point>() );
                        cache2d.push_back ( i );
                    }

                    cache.back().push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }

            if ( !cache2d.empty() && cache2d.back() == i )
            {
                GaussPointArray gp = i->getGaussPoints() ;
                std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

                for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
                {
                    i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
                }
                if ( !function )
                {
                    apply2DBC ( i,gp,Jinv, cache.back(), condition, data*getScale(), a , axis ) ;
                }
                else
                {
                    apply2DBC ( i,gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis ) ;
                }
            }
        }

    }
    else
    {
        for ( size_t i = 0 ; i < cache2d.size() ; ++i )
        {
            GaussPointArray gp = cache2d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache2d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }


            if ( !function )
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }
}

void BoundingBoxDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{

    if ( cache.empty() )
    {

        if ( t->begin().size() == 0 )
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }


            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }
            }
        }

        Point pmin ( minx, miny, 0., mint ) ;
        Point pmax ( maxx, maxy, 0., maxt ) ;

        double tol = std::max ( std::min ( std::min ( maxx - minx, maxy - miny ), maxt-mint ) * .001, POINT_TOLERANCE_2D ) ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) )
                {
                    if ( cache2d.empty() || cache2d.back() != (DelaunayTriangle *)i )
                    {
                        cache.push_back ( std::vector<Point>() );
                        cache2d.push_back ( i );
                    }

                    cache.back().push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }

            if ( !cache2d.empty() && cache2d.back() == (DelaunayTriangle *)i )
            {
                GaussPointArray gp = i->getGaussPoints() ;
                std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

                for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
                {
                    i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
                }

                if ( !function )
                {
                    apply2DBC ( i,gp,Jinv, cache.back(), condition, data*getScale(), a , axis ) ;
                }
                else
                {
                    apply2DBC ( i,gp,Jinv, cache.back(), condition, dataFunction*getScale(), a, axis ) ;
                }
            }
        }

    }
    else
    {
        for ( size_t i = 0 ; i < cache2d.size() ; ++i )
        {
            GaussPointArray gp = cache2d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache2d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache2d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            for ( size_t j = 0 ; j < cache[i].size() ; j++ )
            {
                cache[i][j].setT ( cache[i][j].getT() + cache2d[i]->getState().getDeltaTime() ) ;
            }

            if ( !function )
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, data*getScale(), a, axis ) ;
            }
            else
            {
                apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a, axis ) ;
            }
        }
    }
}

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {

        if ( t->begin().size() == 0)
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double minz = t->begin()->getBoundingPoint ( 0 ).getZ() ;
        double maxz = t->begin()->getBoundingPoint ( 0 ).getZ() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end()  ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() < minz )
                {
                    minz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() > maxz )
                {
                    maxz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }

            }
        }

        double tol = std::min ( std::min ( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

        Point pmin ( minx,miny, minz, mint ) ;
        Point pmax ( maxx,maxy, maxz, maxt ) ;

        Point rmin ( xmin, ymin, zmin ) ;
        Point rmax ( xmax, ymax, zmax ) ;

        for ( auto i = t->begin() ; i != t->end()  ; i++  )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) && isInBoundary3D ( i->getBoundingPoint ( j ), rmin, rmax ) )
                {
                    if ( cache3d.empty() || cache3d.back() != i )
                    {
                        cache.push_back ( std::vector<Point>() );
                        cache3d.push_back ( i );
                    }

                    cache.back().push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }


            if ( !cache3d.empty() && cache3d.back() == i )
            {
                GaussPointArray gp = i->getGaussPoints() ;
                std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

                for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
                {
                    i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
                }
                if ( !function )
                {
                    apply3DBC ( i,gp,Jinv, cache.back(), condition, data*getScale(), a , axis ) ;
                }
                else
                {
                    apply3DBC ( i,gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis ) ;
                }
            }
        }
    }
    else
    {
        for ( size_t i = 0 ; i < cache3d.size() ; ++i )
        {
            GaussPointArray gp = cache3d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            if ( !function )
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }

}

void BoundingBoxDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {

        if ( t->begin().size() == 0 )
        {
            std::cout << "no elements in assembly" << std::endl ;
            return ;
        }

        double minx = t->begin()->getBoundingPoint ( 0 ).getX() ;
        double maxx = t->begin()->getBoundingPoint ( 0 ).getX() ;

        double miny = t->begin()->getBoundingPoint ( 0 ).getY() ;
        double maxy = t->begin()->getBoundingPoint ( 0 ).getY() ;

        double minz = t->begin()->getBoundingPoint ( 0 ).getZ() ;
        double maxz = t->begin()->getBoundingPoint ( 0 ).getZ() ;

        double mint = t->begin()->getBoundingPoint ( 0 ).getT() ;
        double maxt = t->begin()->getBoundingPoint ( 0 ).getT() ;

        for ( auto i = t->begin() ; i != t->end()  ; i++  )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                continue ;
            }

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( i->getBoundingPoint ( j ).getX() < minx )
                {
                    minx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getX() > maxx )
                {
                    maxx = i->getBoundingPoint ( j ).getX() ;
                }

                if ( i->getBoundingPoint ( j ).getY() < miny )
                {
                    miny = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getY() > maxy )
                {
                    maxy = i->getBoundingPoint ( j ).getY() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() < minz )
                {
                    minz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getZ() > maxz )
                {
                    maxz = i->getBoundingPoint ( j ).getZ() ;
                }

                if ( i->getBoundingPoint ( j ).getT() < mint )
                {
                    mint = i->getBoundingPoint ( j ).getT() ;
                }

                if ( i->getBoundingPoint ( j ).getT() > maxt )
                {
                    maxt = i->getBoundingPoint ( j ).getT() ;
                }

            }
        }

        double tol = std::min ( std::min ( maxx - minx, maxy - miny ), maxz - minz ) * .0001 ;

        Point pmin ( minx,miny, minz, mint ) ;
        Point pmax ( maxx,maxy, maxz, maxt ) ;

        for ( auto i = t->begin() ; i != t->end()  ; i++  )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                if ( isOnBoundary ( pos, i->getBoundingPoint ( j ), pmin, pmax, tol ) )
                {
                    if ( cache3d.empty() || cache3d.back() != i )
                    {
                        cache.push_back ( std::vector<Point>() );
                        cache3d.push_back ( i );
                    }

                    cache.back().push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }


            if ( !cache3d.empty() && cache3d.back() == i )
            {
                GaussPointArray gp = i->getGaussPoints() ;
                std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

                for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
                {
                    i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
                }
                if ( !function )
                {
                    apply3DBC ( i,gp,Jinv, cache.back(), condition, data*getScale(), a , axis ) ;
                }
                else
                {
                    apply3DBC ( i,gp,Jinv, cache.back(), condition, dataFunction*getScale(), a , axis ) ;
                }
            }
        }
    }
    else
    {
        for ( size_t i = 0 ; i < cache3d.size() ; ++i )
        {
            GaussPointArray gp = cache3d[i]->getGaussPoints() ;
            std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

            for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
            {
                cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
            }

            if ( !function )
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a , axis ) ;
            }
            else
            {
                apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, dataFunction*getScale(), a , axis ) ;
            }
        }
    }
}

BoundaryCondition::BoundaryCondition ( LagrangeMultiplierType t, const double & d, int a ) : scale ( 1 ), condition ( t ), data ( d ), function ( false ), axis ( a ) { } ;

BoundaryCondition::BoundaryCondition ( LagrangeMultiplierType t, const Function & d, int a ) : scale ( 1 ), condition ( t ), dataFunction ( d ), function ( true ), axis ( a ) { } ;

void BoundaryCondition::setScale ( double d )
{
    scale = d ;
}

double BoundaryCondition::getScale() const
{
    return scale ;
}

ProjectionDefinedBoundaryCondition::ProjectionDefinedBoundaryCondition ( LagrangeMultiplierType t, const Point & dir, double d, int a ) : BoundaryCondition ( t, d, a ), direction ( dir ) { }

void ProjectionDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{

    for ( auto i = t->begin() ; i != t->end()  ; i++  )
    {
        DelaunayTreeItem * VoidItem ;
        bool border = false ;

        for ( size_t j = 0 ; j < i->neighbour.size() ; j++ )
        {
            bool voidNeighbour = ( i->getNeighbour ( j )->isTriangle
                                   && dynamic_cast<DelaunayTriangle *> ( i->getNeighbour ( j ) )->getBehaviour()->type == VOID_BEHAVIOUR ) ;
            border = border || i->getNeighbour ( j )->isPlane
                     || voidNeighbour ;

            if ( voidNeighbour )
            {
                VoidItem = i->getNeighbour ( j ) ;
            }

            if ( i->getNeighbour ( j )->isPlane )
            {
                VoidItem = i->getNeighbour ( j ) ;
            }
        }

        if ( i->getBehaviour()->type == VOID_BEHAVIOUR )
        {
            border = false ;
        }

        if ( border )
        {
            std::pair<Point *, Point*> commonSurface = i->commonEdge ( VoidItem ) ;

            Segment ray ( ( i->getCenter() ), ( i->getCenter() ) - direction* ( i->getRadius() ) ) ;
            bool isOnTheRightSide = ray.intersects ( Segment ( *commonSurface.first, *commonSurface.second ) ) ;

            if ( isOnTheRightSide )
            {
                std::vector<Point> id ;

                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {

                    Line side ( i->getBoundingPoint ( j ), i->getBoundingPoint ( j ) - i->getBoundingPoint ( ( j + 1 ) % i->getBoundingPoints().size() ) ) ;

                    if ( side.intersects ( ray ) )
                    {
                        id.push_back ( i->getBoundingPoint ( j ) ) ;
                        id.push_back ( i->getBoundingPoint ( ( j + 1 ) % i->getBoundingPoints().size() ) ) ;
                    }
                }

                GaussPointArray gp = i->getGaussPoints() ;
                std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

                for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
                {
                    i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
                }

                if ( !id.empty() )
                {
                    if ( !function )
                    {
                        apply2DBC ( i,gp,Jinv, id, condition, data*getScale(), a ) ;
                    }
                    else
                    {
                        apply2DBC ( i,gp,Jinv, id, condition, dataFunction*getScale(), a ) ;
                    }
                }
            }
        }
    }
}

void ProjectionDefinedBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{

    for ( auto i = t->begin() ; i != t->end()  ; i++  )
    {
        std::vector<DelaunayDemiSpace *> space ;

        for ( size_t j = 0 ; j < i->neighbour.size() ; j++ )
        {
            if (i->getNeighbour ( j )->isSpace() )
            {
                space.push_back ( static_cast<DelaunayDemiSpace *> ( i->getNeighbour ( j ) ) ) ;
            }
        }

        std::vector<Point> id ;

        for ( size_t s = 0 ; s < space.size() ; s++ )
        {
            Segment ray ( i->getCenter(), i->getCenter() - direction*2.*i->getRadius() ) ;
            std::vector<Point *> points = space[s]->commonSurface ( i ) ;
            Plane surf ( *points[0], *points[1], *points[2] ) ;

            for ( size_t j = 3 ; j < i->getBoundingPoints().size() ; j++ )
            {
                if ( isCoplanar ( *points[0], *points[1], *points[2], i->getBoundingPoint ( j ) ) )
                {
                    points.push_back ( &i->getBoundingPoint ( j ) ) ;
                }
            }

            if ( surf.intersects ( ray ) && !points.empty() )
            {
                for ( size_t j = 0 ; j < points.size() ; j++ )
                {
                    id.push_back ( *points[j] ) ;
                }
            }
        }
        GaussPointArray gp = i->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), i->getGaussPoints().gaussPoints.size() ) ;

        for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
        {
            i->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
        }

        if ( !function )
        {
            apply3DBC ( i,gp,Jinv, id, condition, data*getScale(), a ) ;
        }
        else
        {
            apply3DBC ( i,gp,Jinv, id, condition, dataFunction*getScale(), a ) ;
        }
    }
}


TimeContinuityBoundaryCondition::TimeContinuityBoundaryCondition ( double i ) : BoundaryCondition ( GENERAL, 0. ), initialValue ( i )
{
    goToNext = true ;
    previousDisp.resize ( 0 ) ;
} ;

void TimeContinuityBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
std::cout << "boum" << std::endl ;
//     t->getAdditionalPoints() ;
    auto j = t->begin() ;
    size_t dof = 0 ;
    while ( dof == 0 )
    {
        
        if ( j == t->end() )
        {
            return ;
        }
        dof = j->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        j++ ;
    }

    size_t timePlanes = j->timePlanes() ;

    if ( timePlanes < 2 )
    {
        return ;
    }

    size_t ndofmax = a->getMaxDofID() ;
    size_t dofPerPlane = ndofmax / timePlanes ;

    previousDisp.resize ( a->getDisplacements().size() ) ;
    previousDisp = a->getDisplacements() ;

    if ( previousDisp.size() == 0 )
    {
        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
//			#pragma omp for
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, initialValue, dofPerPlane*i + j, true )  ;
                }
            }
        }
    }
    else
    {
        size_t extradof = previousDisp.size() - ndofmax*dof ;
        size_t extradofPerPlane = extradof / timePlanes ;

        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, previousDisp[ dofPerPlane* ( i+ ( int ) goToNext ) *dof + j*dof + n], dofPerPlane*i + j, true )  ;
                }
            }

        }
    }
}

void TimeContinuityBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
//     t->getAdditionalPoints() ;
    auto j = t->begin() ;
    size_t dof = 0 ;
    while ( dof == 0 )
    {
        
        if ( j == t->end() )
        {
            return ;
        }
        dof = j->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        j++ ;
    }

    size_t timePlanes = j->timePlanes() ;

    if ( timePlanes < 2 )
    {
        return ;
    }

    size_t ndofmax = a->getMaxDofID() ;
    size_t dofPerPlane = ndofmax / timePlanes ;

    previousDisp.resize ( a->getDisplacements().size() ) ;
    previousDisp = a->getDisplacements() ;

    if ( previousDisp.size() == 0 )
    {
        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
//          #pragma omp for
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, initialValue, dofPerPlane*i + j, true )  ;
                }
            }
        }
    }
    else
    {
        size_t extradof = previousDisp.size() - ndofmax*dof ;
        size_t extradofPerPlane = extradof / timePlanes ;

        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, previousDisp[ dofPerPlane* ( i+ ( int ) goToNext ) *dof + j*dof + n], dofPerPlane*i + j, true )  ;
                }
            }

        }
    }
}


GlobalForceBoundaryCondition::GlobalForceBoundaryCondition ( Vector & data ) : BoundaryCondition ( SET_FORCE_XI, 0., 0 )
{
    dataVector = data ;
}

void GlobalForceBoundaryCondition::setDataVector ( Vector & data )
{
    dataVector = data ;
}

void GlobalForceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    a->addForceVector ( dataVector ) ;
}

void GlobalForceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    a->addForceVector ( dataVector ) ;
}

