// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011

#include "boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../physics/viscoelasticity.h"
#include "../features/features.h"
#include "../physics/material_laws/material_laws.h"

namespace Amie {

void BoundaryCondition::setInterpolation( LinearInterpolatedMaterialLaw * inter ) 
{
    dataInterpolation = inter ;
}

void BoundaryCondition::setInterpolation( std::string f ) 
{
    dataInterpolation = new LinearInterpolatedMaterialLaw( std::make_pair("data","t" ), f) ;
}


void BoundaryCondition::step(double realtime, double dt) 
{
/*    if(dataInterpolation == nullptr)
        std::cout << "dafuq" << std::endl ;*/
    if(dataInterpolation != nullptr)
        data = dataInterpolation->get( realtime ) ;
}

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ) { }

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ) { }

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( zm ), zmax ( zp )
{
}

BoundingBoxCycleDefinedBoundaryCondition::BoundingBoxCycleDefinedBoundaryCondition(std::vector<LoadingCycle> cycles, const std::vector<LagrangeMultiplierType> t, const std::vector<BoundingBoxPosition> & pos) : BoundaryCondition ( t.front(), 0 ), positions(pos), types(t), cycles(cycles), currentCycle(-1)
{ 
    currentBC = new BoundingBoxDefinedBoundaryCondition(t.front(), pos.front(), 0.) ;
}

void BoundingBoxCycleDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t)
{
    if(currentCycle < 0)
    {
        currentBC->apply(a, t);
        currentCycle++ ;
        return ;
    }
    currentBC->setData(cycles[currentCycle].getValue());
    currentBC->apply(a, t);
    if(cycles[currentCycle].isAtEnd())
    {
        delete currentBC ;
        currentCycle++ ;
        currentBC = new BoundingBoxDefinedBoundaryCondition(types[currentCycle], positions[currentCycle], 0.) ;
    }
}

void BoundingBoxCycleDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)
{
    if(currentCycle < 0)
    {
        currentBC->apply(a, t);
        currentCycle++ ;
        return ;
    }
    
    currentBC->setData(cycles[currentCycle].getValue());
    currentBC->apply(a, t);
    if(cycles[currentCycle].isAtEnd())
    {
        delete currentBC ;
        currentCycle++ ;
        currentBC = new BoundingBoxDefinedBoundaryCondition(types[currentCycle], positions[currentCycle], 0.) ;
    }
}

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ),  xmin ( xm ), xmax ( xp ), ymin ( ym ), ymax ( yp ), zmin ( zm ), zmax ( zp )
{

}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, double d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ), nearest ( p ) {}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition ( LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, const Function & d, int a ) : BoundaryCondition ( t, d, a ), pos ( pos ), nearest ( p ) {}

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
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR || e->getBehaviour()->fractured())
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
                for(int ax = 0 ; ax < n/2 ; ax++)
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

        case INCREMENT_ALONG_XI:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n>2)
                a->setPointAlongIndexedAxis ( 0, a->getPreviousDisplacements()[ id[idit]*n ] + data, id[idit] ) ;
            else
                a->setPointAlong ( XI, a->getPreviousDisplacements()[ id[idit]*n ] + data, id[idit] ) ;
            break ;

        case FIX_ALONG_ETA:
        {
            if(n > 2)
            {
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*2+1, 0, id[idit] ) ;
                break ;
            }
            else
                a->setPointAlong ( ETA, 0, id[idit] ) ;
            break ;
        }
        case FIX_ALONG_ALL:
        {
            if(n > 2)
            {
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*2, 0, id[idit] ) ;
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*2+1, 0, id[idit] ) ;
            }
            else
            {
                a->setPointAlong ( XI, 0, id[idit] ) ;
                a->setPointAlong ( ETA, 0, id[idit] ) ;
            }
            break ;
        }

        case SET_ALONG_ETA:
            if(n>2)
                a->setPointAlongIndexedAxis ( 1, data, id[idit] ) ;
            else
                a->setPointAlong ( ETA, data, id[idit] ) ;
            break ;

        case INCREMENT_ALONG_ETA:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n>2)
                a->setPointAlongIndexedAxis ( 1, a->getPreviousDisplacements()[ id[idit]*n + 1 ] + data, id[idit] ) ;
            else
                a->setPointAlong ( ETA, a->getPreviousDisplacements()[ id[idit]*n + 1 ] + data, id[idit] ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, data, id[idit] ) ;
            break ;

        case INCREMENT_ALONG_INDEXED_AXIS:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            a->setPointAlongIndexedAxis ( axis, a->getPreviousDisplacements()[ id[idit]*n + axis ] + data, id[idit] ) ;
            break ;

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ETA:
        {
            a->setPointProportional( XI, ETA, data, 0., id[idit] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_XI:
        {
            a->setPointProportional( ETA, XI, data, 0., id[idit] ) ;
            break ;
        }

        case FIX_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
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
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            if(vec.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( XI, 0., id[idit] ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, 0., id[idit] ) ;
            else
               a->setPointProportional( ETA, XI, vec.getY()/vec.getX(), 0., id[idit] ) ;
            break ;
        }

        case SET_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
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
            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( XI, data, id[idit] ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, data, id[idit] ) ;
            else
            {
               double slope = vec.getY()/vec.getX() ;
               double normx = data/std::sqrt( 1.+(1./slope)*(1./slope)) ;
//               double normalSlope = norm.getY()/norm.getX() ;
               a->setPointProportional( ETA, XI, slope, normx*(-slope-1./slope), id[idit] ) ;
            }
            break ;
        }

        case FIX_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
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
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, 0., id[idit] ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( XI, 0., id[idit] ) ;
            else
               a->setPointProportional( ETA, XI, -vec.getX()/vec.getY(), 0., id[idit] ) ;
            break ;
        }

        case SET_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
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
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, data, id[idit] ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( XI, data, id[idit] ) ;
            else
            {
               double slope = norm.getY()/norm.getX() ;
               double tanx = data/std::sqrt( 1.+(1./slope)*(1./slope)) ;
               a->setPointProportional( ETA, XI, slope, tanx*(-slope-1./slope), id[idit] ) ;
            }
            break ;
        }

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

//             for ( size_t j = 0 ; j < id.size() ; j++ )
//             {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[idit] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)id[idit] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

//             }

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

                if(v.size() == 3)
                {
                    Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                    if ( visc && visc->model > KELVIN_VOIGT )
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                        for ( int b = 1 ; b < visc->blocks ;  b++ )
                        {
                            a->addForceToExternalForces ( 2*b, -forces[0], id[idit] ) ;
                            a->addForceToExternalForces ( 2*b+1, -forces[1], id[idit] ) ;
                        }
                    } 
                    else
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[idit] ) ;
                    a->addForceOn ( ETA, forces[1], id[idit] ) ;
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
                    if ( (int)id[idit] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)id[idit] == e->getEnrichmentFunction ( i ).getDofID() )
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

                if(v.size() == 3)
                {
                    Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                    if ( visc && visc->model > KELVIN_VOIGT )
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                        for ( int b = 1 ; b < visc->blocks ;  b++ )
                        {
                            a->addForceToExternalForces ( 2*b, -forces[0], id[idit] ) ;
                            a->addForceToExternalForces ( 2*b+1, -forces[1], id[idit] ) ;
                        }
                    } 
                    else
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[idit] ) ;
                    a->addForceOn ( ETA, forces[1], id[idit] ) ;
                }

            }

            return ;
        }

        case SET_VOLUMIC_STRESS_XI_ETA:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            std::vector<Function> shapeFunctions ;

//             for ( size_t j = 0 ; j < id.size() ; j++ )
//             {
                for ( size_t i = 0 ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[idit] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)id[idit] == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

//             }

            std::vector<Variable> v ( 2 ) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }

            Vector imposed ( 3 ) ;
            imposed[0] = 0 ;
            imposed[1] = 0 ;
            imposed[2] = data ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;

                if(v.size() == 3)
                {
                    Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                    if ( visc && visc->model > KELVIN_VOIGT )
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                        for ( int b = 1 ; b < visc->blocks ;  b++ )
                        {
                            a->addForceToExternalForces ( 2*b, -forces[0], id[idit] ) ;
                            a->addForceToExternalForces ( 2*b+1, -forces[1], id[idit] ) ;
                        }
                    } 
                    else
                    {
                        a->addForceToExternalForces ( 0, forces[0], id[idit] ) ;
                        a->addForceToExternalForces ( 1, forces[1], id[idit] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[idit] ) ;
                    a->addForceOn ( ETA, forces[1], id[idit] ) ;
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                for ( size_t i = (nTimePlanes-1)*e->getBoundingPoints().size()/nTimePlanes ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

            Vector imposed (0., 3 ) ;
            imposed[1] = data ;

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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
            test[0] = 0 ;
            test[1] = 1 ;

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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR || e->getBehaviour()->fractured())
    {
        return ;
    }

    VirtualMachine vm ;
    int n = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
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

        case INCREMENT_ALONG_XI:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n>3)
                a->setPointAlongIndexedAxis ( 0, a->getPreviousDisplacements()[ id[i]*n ] + data, id[i] ) ;
            else
                a->setPointAlong ( XI, a->getPreviousDisplacements()[ id[i]*n ] + data, id[i] ) ;
            break ;

        case FIX_ALONG_ETA:
            a->setPointAlong ( ETA, 0., id[i] ) ;
            break ;

        case SET_ALONG_ETA:
            a->setPointAlong ( ETA, data, id[i] ) ;
            break ;

        case INCREMENT_ALONG_ETA:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n>3)
                a->setPointAlongIndexedAxis ( 1, a->getPreviousDisplacements()[ id[i]*n + 1 ] + data, id[i] ) ;
            else
                a->setPointAlong ( ETA, a->getPreviousDisplacements()[ id[i]*n + 1 ] + data, id[i] ) ;
            break ;


        case FIX_ALONG_ZETA:
            a->setPointAlong ( ZETA, 0., id[i] ) ;
            break ;

        case FIX_ALONG_ALL:
        {
            if(n > 3)
            {
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3, 0, id[i] ) ;
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3+1, 0, id[i] ) ;
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3+2, 0, id[i] ) ;
            }
            else
            {
                a->setPointAlong ( XI, 0, id[i] ) ;
                a->setPointAlong ( ETA, 0, id[i] ) ;
                a->setPointAlong ( ZETA, 0, id[i] ) ;
            }
            break ;
        }

        case SET_ALONG_ZETA:
            a->setPointAlong ( ZETA, data, id[i] ) ;
            break ;

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ETA:
        {
            a->setPointProportional( XI, ETA, data, 0., id[i] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_XI:
        {
            a->setPointProportional( ETA, XI, data, 0., id[i] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA:
        {
            a->setPointProportional( XI, ZETA, data, 0., id[i] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI:
        {
            a->setPointProportional( ZETA, XI, data, 0., id[i] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA:
        {
            a->setPointProportional( ETA, ZETA, data, 0., id[i] ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA:
        {
            a->setPointProportional( ZETA, XI, data, 0., id[i] ) ;
            break ;
        }

        case INCREMENT_ALONG_ZETA:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n>3)
                a->setPointAlongIndexedAxis ( 2, a->getPreviousDisplacements()[ id[i]*n + 2 ] + data, id[i] ) ;
            else
                a->setPointAlong ( ZETA, a->getPreviousDisplacements()[ id[i]*n + 2] + data, id[i] ) ;
            break ;


        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, data, id[i] ) ;
            break ;

        case INCREMENT_ALONG_INDEXED_AXIS:
            if( std::abs( data ) < POINT_TOLERANCE || a->getPreviousDisplacements().size() == 0 )
                break ;

            a->setPointAlongIndexedAxis ( axis, a->getPreviousDisplacements()[ id[i]*n + axis ] + data, id[i] ) ;
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }

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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false, std::abs (  edge.normalv ( e->getCenter() ) ) ) ;
                a->addForceOn ( XI, forces[0], id[i] ) ;
                a->addForceOn ( ETA, forces[1], id[i] ) ;
                a->addForceOn ( ZETA, forces[2], id[i] ) ;
            }

            return ;
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

            Vector imposed ( 6 ) ;
            imposed[0] = data ;
            imposed[1] = 0 ;
            imposed[2] = 0 ;
            imposed[3] = 0 ;
            imposed[4] = 0 ;
            imposed[5] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;


                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT )
                {
                    a->addForceToExternalForces ( 0, forces[0], id[j] ) ;
                    a->addForceToExternalForces ( 1, forces[1], id[j] ) ;
                    a->addForceToExternalForces ( 2, forces[2], id[j] ) ;
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceToExternalForces ( 3*b, -forces[0], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+1, -forces[1], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+2, -forces[2], id[j] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[j] ) ;
                    a->addForceOn ( ETA, forces[1], id[j] ) ;
                    a->addForceOn ( ZETA, forces[1], id[j] ) ;
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

            Vector imposed ( 6 ) ;
            imposed[0] = 0 ;
            imposed[1] = data ;
            imposed[2] = 0 ;
            imposed[3] = 0 ;
            imposed[4] = 6 ;
            imposed[5] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT )
                {
                    a->addForceToExternalForces ( 0, forces[0], id[j] ) ;
                    a->addForceToExternalForces ( 1, forces[1], id[j] ) ;
                    a->addForceToExternalForces ( 2, forces[2], id[j] ) ;
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceToExternalForces ( 3*b, -forces[0], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+1, -forces[1], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+2, -forces[2], id[j] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[j] ) ;
                    a->addForceOn ( ETA, forces[1], id[j] ) ;
                    a->addForceOn ( ZETA, forces[2], id[j] ) ;
                }

            }

            return ;
        }
        case SET_VOLUMIC_STRESS_ZETA:
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
                    }
                }
                for ( size_t i = 0 ; i < e->getEnrichmentFunctions().size() ; i++ )
                {
                    if ( (int)(int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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

            Vector imposed ( 6 ) ;
            imposed[0] = 0 ;
            imposed[1] = 0 ;
            imposed[2] = data ;
            imposed[3] = 0 ;
            imposed[4] = 6 ;
            imposed[5] = 0 ;

            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[j], gp, Jinv, v, true, Vector() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT )
                {
                    a->addForceToExternalForces ( 0, forces[0], id[j] ) ;
                    a->addForceToExternalForces ( 1, forces[1], id[j] ) ;
                    a->addForceToExternalForces ( 2, forces[2], id[j] ) ;
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceToExternalForces ( 3*b, -forces[0], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+1, -forces[1], id[j] ) ;
                        a->addForceToExternalForces ( 3*b+2, -forces[2], id[j] ) ;
                    }
                } else {

                    a->addForceOn ( XI, forces[0], id[j] ) ;
                    a->addForceOn ( ETA, forces[1], id[j] ) ;
                    a->addForceOn ( ZETA, forces[2], id[j] ) ;
                }

            }

            return ;
        }

        case FIX_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][1]*invTransform[1][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][1]*invTransform[1][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][1]*invTransform[1][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, 0., id[i] ) ;

            return ;
        }

        case SET_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][1]*invTransform[1][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][1]*invTransform[1][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][1]*invTransform[1][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, transform[0][0]*data/coef, id[i] ) ;

            return ;
        }

        case FIX_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][0]*invTransform[0][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][0]*invTransform[0][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][0]*invTransform[0][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, 0., id[i] ) ;

            return ;
        }

        case SET_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][0]*invTransform[0][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][0]*invTransform[0][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][0]*invTransform[0][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, transform[0][1]*data/coef, id[i] ) ;

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

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( k ) ) ;
                    }
                    if ( (int)id[j] == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( k ).getDofID() )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( (int)id[j] == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                    if ( (int)id[j] == e->getEnrichmentFunction ( i ).getDofID() )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

    if ( e->getBehaviour()->type == VOID_BEHAVIOUR || e->getBehaviour()->fractured())
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
    

//	std::cerr << id.size() << std::endl ;

    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        switch ( condition )
        {

        case GENERAL :
            std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
            break ;

        case FIX_ALONG_XI:
            if(n > 2)
            {
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis( ax*2, 0, id[i].getId() ) ;
            }
            else
                a->setPointAlong ( XI, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_XI:
            if(n > 2)
                a->setPointAlongIndexedAxis( 0, vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( XI, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_XI:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n > 2)
                a->setPointAlongIndexedAxis( 0, a->getPreviousDisplacements()[ id[i].getId()*n ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( XI, a->getPreviousDisplacements()[ id[i].getId()*n ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ETA:
            if(n > 2)
            {
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis( ax*2+1, 0, id[i].getId() ) ;
            }
            else
                a->setPointAlong ( ETA, 0, id[i].getId() ) ;
            break ;

        case FIX_ALONG_ALL:
        {
            if(n > 2)
            {
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*2, 0, id[i].getId() ) ;
                for(int ax = 0 ; ax < n/2 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*2+1, 0, id[i].getId() ) ;
            }
            else
            {
                a->setPointAlong ( XI, 0, id[i].getId() ) ;
                a->setPointAlong ( ETA, 0, id[i].getId() ) ;
            }
            break ;
        }

        case SET_ALONG_ETA:
            if(n > 2)
                a->setPointAlongIndexedAxis( 1, vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ETA, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ETA:
        {
            a->setPointProportional( XI, ETA, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_XI:
        {
            a->setPointProportional( ETA, XI, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case FIX_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( k ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( k ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( k ) ) )
                                {
                                    last = &e->getBoundingPoint ( k ) ;
                                }
                                if ( dist ( first, last ) < dist ( last, &e->getBoundingPoint ( k ) ) )
                                {
                                    first = &e->getBoundingPoint ( k ) ;
                                }
                            }
                        }
                    }
                }
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            if(vec.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( XI, 0., id[i].getId() ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, 0., id[i].getId() ) ;
            else
               a->setPointProportional( ETA, XI, vec.getY()/vec.getX(), 0., id[i].getId() ) ;
            break ;
        }

        case SET_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( k ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( k ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( k ) ) )
                                {
                                    last = &e->getBoundingPoint ( k ) ;
                                }
                                if ( dist ( first, last ) < dist ( last, &e->getBoundingPoint ( k ) ) )
                                {
                                    first = &e->getBoundingPoint ( k ) ;
                                }
                            }
                        }
                    }
                }
            }

            if ( !last )
            {
                return ;
            }


            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( XI, vm.eval(data, id[i]), id[i].getId() ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, vm.eval(data, id[i]), id[i].getId() ) ;
            else
            {
               double slope = vec.getY()/vec.getX() ;
               double normx = vm.eval(data, id[i])/std::sqrt( 1.+(1./slope)*(1./slope)) ;
//               double normalSlope = norm.getY()/norm.getX() ;
               a->setPointProportional( ETA, XI, slope, normx*(-slope-1./slope), id[i].getId() ) ;
            }
            break ;
        }

        case FIX_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( k ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( k ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( k ) ) )
                                {
                                    last = &e->getBoundingPoint ( k ) ;
                                }
                                if ( dist ( first, last ) < dist ( last, &e->getBoundingPoint ( k ) ) )
                                {
                                    first = &e->getBoundingPoint ( k ) ;
                                }
                            }
                        }
                    }
                }
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, 0., id[i].getId() ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( XI, 0., id[i].getId() ) ;
            else
               a->setPointProportional( ETA, XI, -vec.getX()/vec.getY(), 0., id[i].getId() ) ;
            break ;
        }

        case SET_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() )
                    {
                        if ( !first )
                        {
                            first = &e->getBoundingPoint ( k ) ;
                        }
                        else
                        {
                            if ( !last )
                            {
                                last = &e->getBoundingPoint ( k ) ;
                            }
                            else
                            {
                                if ( dist ( first, last ) < dist ( first, &e->getBoundingPoint ( k ) ) )
                                {
                                    last = &e->getBoundingPoint ( k ) ;
                                }
                                if ( dist ( first, last ) < dist ( last, &e->getBoundingPoint ( k ) ) )
                                {
                                    first = &e->getBoundingPoint ( k ) ;
                                }
                            }
                        }
                    }
                }
            }

            if ( !last )
            {
                return ;
            }

            Segment edge ( *last, *first ) ;
            Point vec = edge.vector() ;
            Point norm = edge.normal() ;
            if(vec.norm() < POINT_TOLERANCE || norm.norm() < POINT_TOLERANCE)
               return ;

            if(std::abs(vec.getX()) < POINT_TOLERANCE)
               a->setPointAlong( ETA, vm.eval(data, id[i]), id[i].getId() ) ;
            else if(std::abs(vec.getY()) < POINT_TOLERANCE)
               a->setPointAlong( XI, vm.eval(data, id[i]), id[i].getId() ) ;
            else
            {
               double slope = norm.getY()/norm.getX() ;
               double tanx = vm.eval(data, id[i])/std::sqrt( 1.+(1./slope)*(1./slope)) ;
               a->setPointProportional( ETA, XI, slope, tanx*(-slope-1./slope), id[i].getId() ) ;
            }
            break ;
        }

        case INCREMENT_ALONG_ETA:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n > 2)
                a->setPointAlongIndexedAxis( 1, a->getPreviousDisplacements()[ id[i].getId()*n+1 ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ETA, a->getPreviousDisplacements()[ id[i].getId()*n+1 ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_INDEXED_AXIS:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            a->setPointAlongIndexedAxis( axis, a->getPreviousDisplacements()[ id[i].getId()*n + axis ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
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
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
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
                if ( visc && visc->model > KELVIN_VOIGT)
                {
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
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
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
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
                if ( visc && visc->model > KELVIN_VOIGT)
                {
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
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
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = (nTimePlanes-1)*e->getBoundingPoints().size()/nTimePlanes ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
                    {
                        shapeFunctions.push_back ( e->getEnrichmentFunction ( i ) ) ;
                    }
                }

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
            Vector imposed ( 3 ) ;
            imposed[0] = VirtualMachine().eval(data, edge.midPoint()) ;

            v[0] = XI ;
            v[1] = ETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

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

            std::vector<Function> shapeFunctions ;
            Point* first = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t i = (nTimePlanes-1)*e->getBoundingPoints().size()/nTimePlanes ; i < e->getBoundingPoints().size() ; i++ )
                {
                    if ( (int)id[j].getId() == e->getBoundingPoint ( i ).getId() )
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
                    if ( (int)id[j].getId() == e->getEnrichmentFunction ( i ).getDofID() )
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
            Vector imposed (0.,  3 ) ;
            imposed[1] = VirtualMachine().eval(data, 0,0,0,  e->getBoundingPoint(0).getT()) ;
            
            VirtualMachine vm ;

            for ( size_t i = 0 ; i < shapeFunctions.size() ; ++i )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( imposed, shapeFunctions[i], gpe, Jinve, v, false, edge.normalv ( e->getCenter() ) ) ;

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
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
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
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
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
    if ( e->getBehaviour()->type == VOID_BEHAVIOUR || e->getBehaviour()->fractured())
    {
        return ;
    }

    VirtualMachine vm ;
    int n = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    for ( size_t i = 0 ; i < id.size() ; i++ )
    {
        switch ( condition )
        {

        case GENERAL :
            std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
            break ;

        case FIX_ALONG_XI:
            if(n > 3)
                a->setPointAlongIndexedAxis( 0, 0, id[i].getId() ) ;
            else
                a->setPointAlong ( XI, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_XI:
            if(n > 3)
                a->setPointAlongIndexedAxis( 0, vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( XI, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_XI:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n > 3)
                a->setPointAlongIndexedAxis( 0, a->getPreviousDisplacements()[ id[i].getId()*n ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( XI, a->getPreviousDisplacements()[ id[i].getId()*n ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ETA:
            if(n > 3)
                a->setPointAlongIndexedAxis( 1, 0, id[i].getId() ) ;
            else
                a->setPointAlong ( XI, 0, id[i].getId() ) ;
            break ;

        case FIX_ALONG_ALL:
        {
            if(n > 3)
            {
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3, 0, id[i].getId() ) ;
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3+1, 0, id[i].getId() ) ;
                for(int ax = 0 ; ax < n/3 ; ax++)
                    a->setPointAlongIndexedAxis ( ax*3+2, 0, id[i].getId() ) ;
            }
            else
            {
                a->setPointAlong ( XI, 0, id[i].getId() ) ;
                a->setPointAlong ( ETA, 0, id[i].getId() ) ;
                a->setPointAlong ( ZETA, 0, id[i].getId() ) ;
            }
            break ;
        }

        case SET_ALONG_ETA:
            if(n > 3)
                a->setPointAlongIndexedAxis( 1, vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ETA, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_ETA:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n > 3)
                a->setPointAlongIndexedAxis( 1, a->getPreviousDisplacements()[ id[i].getId()*n + 1 ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ETA, a->getPreviousDisplacements()[ id[i].getId()*n + 1] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case FIX_ALONG_ZETA:
            if(n > 3)
                a->setPointAlongIndexedAxis( 2, 0, id[i].getId() ) ;
            else
                a->setPointAlong ( ZETA, 0, id[i].getId() ) ;
            break ;

        case SET_ALONG_ZETA:
            if(n > 3)
                a->setPointAlongIndexedAxis( 2, vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ZETA, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_ZETA:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            if(n > 3)
                a->setPointAlongIndexedAxis( 2, a->getPreviousDisplacements()[ id[i].getId()*n + 2 ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            else
                a->setPointAlong ( ZETA, a->getPreviousDisplacements()[ id[i].getId()*n + 2] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_ALONG_INDEXED_AXIS:
            a->setPointAlongIndexedAxis ( axis, vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case INCREMENT_ALONG_INDEXED_AXIS:
            if( a->getPreviousDisplacements().size() == 0 )
                break ;

            a->setPointAlongIndexedAxis( axis , a->getPreviousDisplacements()[ id[i].getId()*n + axis ] + vm.eval ( data, id[i] ), id[i].getId() ) ;
            break ;

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ETA:
        {
            a->setPointProportional( XI, ETA, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_XI:
        {
            a->setPointProportional( ETA, XI, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA:
        {
            a->setPointProportional( XI, ZETA, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI:
        {
            a->setPointProportional( ZETA, XI, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA:
        {
            a->setPointProportional( ZETA, ETA, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

        case SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA:
        {
            a->setPointProportional( ETA, ZETA, vm.eval(data, id[i]), 0., id[i].getId() ) ;
            break ;
        }

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
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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

                a->addForceOn ( XI, forces[0], id[i].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[i].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[i].getId() ) ;
            }

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
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
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


            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 0, 3, shapeFunctions[j], e, gp, Jinv, v, true ) ;

                a->addForceOn ( XI, forces[0], id[j].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[j].getId() ) ;
                a->addForceOn ( ZETA, forces[1], id[j].getId() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT)
                {
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 3*b, -forces[0], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+1, -forces[1], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+2, -forces[2], id[j].getId() ) ;
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
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
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
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 1, 3, shapeFunctions[j], e, gp, Jinv, v, true ) ;

                a->addForceOn ( XI, forces[0], id[j].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[j].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[j].getId() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT)
                {
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 3*b, -forces[0], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+1, -forces[1], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+2, -forces[2], id[j].getId() ) ;
                    }
                }
            }

            return ;
        }
        
        case SET_VOLUMIC_STRESS_ZETA:
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
                        shapeFunctions.push_back ( Function(e->getShapeFunction ( i ), 2) ) ;
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
            v[2] = ZETA ;
            if ( e->getOrder() >= CONSTANT_TIME_LINEAR )
            {
                v.push_back ( TIME_VARIABLE ) ;
            }


            for ( size_t j = 0 ; j < shapeFunctions.size() ; ++j )
            {
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( data, 2, 3, shapeFunctions[j], e, gp, Jinv, v, true ) ;

                a->addForceOn ( XI, forces[0], id[j].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[j].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[j].getId() ) ;

                Viscoelasticity * visc = dynamic_cast<Viscoelasticity *> ( e->getBehaviour() ) ;
                if ( visc && visc->model > KELVIN_VOIGT)
                {
                    for ( int b = 1 ; b < visc->blocks ;  b++ )
                    {
                        a->addForceOnIndexedAxis ( 3*b, -forces[0], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+1, -forces[1], id[j].getId() ) ;
                        a->addForceOnIndexedAxis ( 3*b+2, -forces[2], id[j].getId() ) ;
                    }
                }
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

        case FIX_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][1]*invTransform[1][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][1]*invTransform[1][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][1]*invTransform[1][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, 0., id[i].getId() ) ;

            return ;
        }

        case SET_NORMAL_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][1]*invTransform[1][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][1]*invTransform[1][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][1]*invTransform[1][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, transform[0][0]*vm.eval(data, id[i])/coef, id[i].getId() ) ;

            return ;
        }

        case FIX_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][0]*invTransform[0][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][0]*invTransform[0][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][0]*invTransform[0][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, 0., id[i].getId() ) ;

            return ;
        }

        case SET_TANGENT_DISPLACEMENT:
        {
            if ( e->getBehaviour()->fractured() )
            {
                return ;
            }

            Point* first = nullptr ;
            Point* middle = nullptr ;
            Point* last = nullptr ;

            for ( size_t j = 0 ; j < id.size() ; j++ )
            {
                for ( size_t k = 0 ; k < e->getBoundingPoints().size() ; k++ )
                {

                    DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *> ( e ) ;

                    if ( (int)id[j].getId() == e->getBoundingPoint ( k ).getId() && (
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( e->getBoundingPoint ( k ), *tet->fourth ) < POINT_TOLERANCE )
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

            Point np ( normal[0], normal[1], normal[2] ) ;
            Point np0 ( -normal[1], normal[0], -normal[2] ) ;
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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

            Matrix transform(3,3) ;
            transform[0][0] = lx ;
            transform[0][1] = ly ;
            transform[0][2] = lz ;
            transform[1][0] = rx ;
            transform[1][1] = ry ;
            transform[1][2] = rz ;
            transform[2][0] = tx ;
            transform[2][1] = ty ;
            transform[2][2] = tz ;

            Matrix invTransform = inverse3x3Matrix( transform ) ;

            std::vector< std::pair< Variable, double > > coefs ;
            double coef = 1.-transform[0][0]*invTransform[0][0]+transform[0][2]*invTransform[2][0] ;
            coefs.push_back( std::make_pair( ETA, (transform[0][0]*invTransform[0][1]+transform[0][2]*invTransform[2][1])/coef ) ) ;
            coefs.push_back( std::make_pair( ZETA, (transform[0][0]*invTransform[0][2]+transform[0][2]*invTransform[2][2])/coef ) ) ;

            a->setPointProportional( XI, coefs, transform[0][1]*vm.eval(data, id[i])/coef, id[i].getId() ) ;

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
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() )
                    {
                        shapeFunctions.push_back ( e->getShapeFunction ( i ) ) ;
                    }
                    if ( id[j].getId() == e->getBoundingPoint ( i ).getId() && (
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->first ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->second ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->third ) < POINT_TOLERANCE  ||
                                squareDist3D ( &e->getBoundingPoint ( i ), dynamic_cast<DelaunayTetrahedron *> ( e )->fourth ) < POINT_TOLERANCE )
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
            if ( std::abs ( normal[1] - normal[0] ) < POINT_TOLERANCE )
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


            for ( size_t k = 0 ; k < shapeFunctions.size() ; ++k )
            {
                imposed[0] = vm.eval ( data, id[k] ) ;
                istr = transform* ( imposed ) ;
                Vector forces = e->getBehaviour()->getForcesFromAppliedStress ( istr, shapeFunctions[k], gpe, Jinve, v, false, normal ) ;

                a->addForceOn ( XI, forces[0], id[k].getId() ) ;
                a->addForceOn ( ETA, forces[1], id[k].getId() ) ;
                a->addForceOn ( ZETA, forces[2], id[k].getId() ) ;
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

        for ( size_t i = 0 ; i < surface->getBoundingPoints().size() ; i++ )
        {
            if ( surface->getBoundingPoint ( i ).getId() == (int)id )
            {
                id_.push_back ( surface->getBoundingPoint ( i ) );
                Function effective = dataFunction*getScale() ;
                apply2DBC ( surface,*gp,*Jinv, id_, condition, effective, a, axis ) ;
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

        for ( size_t i = 0 ; i < volume->getBoundingPoints().size() ; i++ )
        {
            if ( volume->getBoundingPoint ( i ).getId() == (int)id )
            {
                id_.push_back ( volume->getBoundingPoint ( i ) );
                Function effective = dataFunction*getScale() ;
                apply3DBC ( volume,*gp, *Jinv, id_, condition, effective, a, axis ) ;
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

            if ( dist ( test, i->getBoundingPoint ( j ) ) < POINT_TOLERANCE )
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

}

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

            if ( dist ( test, i->getBoundingPoint ( j ) ) < POINT_TOLERANCE )
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
}

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
    case NOW:
        return ( std::abs ( test.getT() - (min.getT()+max.getT())*.5 ) < tol ) ;

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
    case BOTTOM_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_BEFORE:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case TOP_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

    case BACK_BEFORE:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case BACK_AFTER:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BACK_NOW:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case FRONT_BEFORE:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case FRONT_AFTER:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case FRONT_NOW:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

    case LEFT_BEFORE:
        return ( isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case LEFT_AFTER:
        return ( isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case LEFT_NOW:
        return ( isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case RIGHT_BEFORE:
        return ( isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BEFORE, test, min, max, tol ) ) ;
    case RIGHT_AFTER:
        return ( isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case RIGHT_NOW:
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

    case BOTTOM_LEFT_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_LEFT_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_RIGHT_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

    case BACK_LEFT_NOW:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case BACK_RIGHT_NOW:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case FRONT_LEFT_NOW:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case FRONT_RIGHT_NOW:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

    case BOTTOM_BACK_NOW:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_BACK_NOW:
        return ( isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case FRONT_BOTTOM_NOW:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case FRONT_TOP_NOW:
        return ( isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

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
    case BOTTOM_LEFT_BACK_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BOTTOM_LEFT_FRONT_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_BACK_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_FRONT_AFTER:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_LEFT_BACK_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_LEFT_FRONT_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_RIGHT_BACK_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case TOP_RIGHT_FRONT_AFTER:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( AFTER, test, min, max, tol ) ) ;
    case BOTTOM_LEFT_BACK_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case BOTTOM_LEFT_FRONT_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_BACK_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case BOTTOM_RIGHT_FRONT_NOW:
        return ( isOnBoundary ( BOTTOM, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_LEFT_BACK_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_LEFT_FRONT_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( LEFT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_RIGHT_BACK_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( BACK, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;
    case TOP_RIGHT_FRONT_NOW:
        return ( isOnBoundary ( TOP, test, min, max, tol ) && isOnBoundary ( RIGHT, test, min, max, tol ) && isOnBoundary ( FRONT, test, min, max, tol ) && isOnBoundary ( NOW, test, min, max, tol ) ) ;

    }
    return false ;
}


void BoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t)
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
            Function effective = dataFunction*getScale() ;
            apply2DBC ( cache2d[i],gp,Jinv, cache[i], condition, effective, a, axis ) ;
        }
    }
}

void BoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)
{
    for ( size_t i = 0 ; i < cache3d.size() ; ++i )
    {
        GaussPointArray gp = cache3d[i]->getGaussPoints() ;
        std::valarray<Matrix> Jinv ( Matrix(), cache3d[i]->getGaussPoints().gaussPoints.size() ) ;

        for ( size_t j = 0 ; j < gp.gaussPoints.size() ; j++ )
        {
            cache3d[i]->getInverseJacobianMatrix ( gp.gaussPoints[j].first, Jinv[j] ) ;
        }

        for ( size_t j = 0 ; j < cache[i].size() ; j++ )
        {
            cache[i][j].setT ( cache[i][j].getT() + cache3d[i]->getState().getDeltaTime() ) ;
        }

        if ( !function )
        {
            apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, data*getScale(), a, axis ) ;
        }
        else
        {
            Function effective = dataFunction*getScale() ;
            apply3DBC ( cache3d[i],gp,Jinv, cache[i], condition, effective, a, axis ) ;
        }
    }
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
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
            Function effective = dataFunction*getScale() ;
            apply2DBC ( id.begin()->second.second,gp,Jinv, target, condition, effective, a , axis ) ;
        }


    }

    BoundaryCondition::apply(a,t);
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
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
            Function effective = dataFunction*getScale() ;
            apply3DBC ( id.begin()->second.second,gp,Jinv, target, condition, effective, a, axis ) ;
        }


    }
    
    BoundaryCondition::apply(a,t);

}


GeometryAndFaceDefinedSurfaceBoundaryCondition::GeometryAndFaceDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & normal, double d , int a ) : BoundaryCondition ( t, d, a ), domain ( source ),faceNormal(normal/normal.norm()) { }
GeometryAndFaceDefinedSurfaceBoundaryCondition::GeometryAndFaceDefinedSurfaceBoundaryCondition(LagrangeMultiplierType t, Geometry * source, const Point & normal, const Function & d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ),faceNormal(normal/normal.norm()) { }

void GeometryAndFaceDefinedSurfaceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( cache2d.empty() )
    {
        std::vector<DelaunayTriangle *> elements = t->getConflictingElements ( domain ) ;
        double tol = domain->getRadius() * .001 ;


        for (const auto & i : elements )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }

            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;

            if(domain->in(i->getCircumCenter()+faceNormal*i->getRadius()*1.2))
                continue ;

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                Point test = i->getBoundingPoint ( j ) +faceNormal*i->getRadius() ;
                if(domain->in(test))
                    continue ;
                domain->project ( &test );

                if ( squareDist2D ( test, i->getBoundingPoint ( j ) ) < tol*tol )
                {
                    id.push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }
            if ( !id.empty() )
            {
                cache2d.push_back ( i );
                cache.push_back ( id );
            }

        }
    }

    BoundaryCondition::apply(a,t);
}

void GeometryAndFaceDefinedSurfaceBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {
//         std::vector<DelaunayTetrahedron *> elements = t->getConflictingElements ( domain ) ;
        unsigned int cacheid = t->generateCache(domain, nullptr) ;
        double tol = domain->getRadius() * .001 ;
        for ( auto i = t->begin(cacheid) ; i != t->end(cacheid) ; i++ )
        {
//             std::cerr << i.getPosition() << std::endl ;
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
                continue ;

            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;

            if(domain->in(i->getCircumCenter()+faceNormal*i->getRadius()*1.2))
                continue ;

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {
                Point test = i->getBoundingPoint ( j ) ;

                domain->project ( &test );

                if ( squareDist3D ( test, i->getBoundingPoint ( j ) ) < tol*tol )
                {
                    id.push_back ( i->getBoundingPoint ( j ) ) ;
                }
            }
            if ( !id.empty() )
            {
                cache3d.push_back ( i );
                cache.push_back ( id );
            }
        }
    }
    
    BoundaryCondition::apply(a,t);
}

GeometryDefinedSurfaceBoundaryCondition::GeometryDefinedSurfaceBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, double d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { }

GeometryDefinedSurfaceBoundaryCondition::GeometryDefinedSurfaceBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, const Function & d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { }

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

    BoundaryCondition::apply(a,t);
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

    BoundaryCondition::apply(a,t);
}



GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, double d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { }

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition ( LagrangeMultiplierType t, Geometry * source, const Function & d, int a ) : BoundaryCondition ( t, d, a ), domain ( source ) { }

GlobalBoundaryCondition::GlobalBoundaryCondition ( LagrangeMultiplierType t, double d, int a ) : BoundaryCondition ( t, d, a ) { }

GlobalBoundaryCondition::GlobalBoundaryCondition ( LagrangeMultiplierType t, const Function & d, int a ) : BoundaryCondition ( t, d, a ) { }

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
            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;

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
                Function effective = dataFunction*getScale() ;
                apply2DBC ( i,gp, Jinv, id, condition, effective, a ) ;
            }
        }
    }

    BoundaryCondition::apply(a,t);
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
            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;


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
                Function effective = dataFunction*getScale() ;
                apply3DBC ( i,gp, Jinv, id, condition, effective, a ) ;
            }
        }
    }

    BoundaryCondition::apply(a,t);
}

void GlobalBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
    if ( cache.empty() )
    {

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }
            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;

            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {

                id.push_back ( i->getBoundingPoint ( j ) ) ;
            }

            cache2d.push_back ( i );
            cache.push_back ( id );

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
                Function effective = dataFunction*getScale() ;
                apply2DBC ( i,gp, Jinv, id, condition, effective, a ) ;
            }
        }
    }

    BoundaryCondition::apply(a,t);
}

void GlobalBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
    if ( cache.empty() )
    {

        
        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if(i.getPosition()%1000 == 0)
                std::cerr << "\r dof " <<i.getPosition() << "/" << i.size() << std::flush ;
            if ( i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured() )
            {
                continue ;
            }
            if(i->getBehaviour()->type == VOID_BEHAVIOUR)
                continue ;


            std::vector<Point> id  ;

            for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
            {

                id.push_back ( i->getBoundingPoint ( j ) ) ;
            }


            cache3d.push_back ( i );
            cache.push_back ( id );
        }
        std::cerr << "\r dof " <<t->begin().size() << "/" << t->begin().size() << std::endl ;
    }

    BoundaryCondition::apply(a,t);
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
        auto i = t->begin() ;
        while((i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR)
        {
            i++ ;
            if(i == t->end())
                return ;
        }

        double minx = i->getBoundingPoint ( 0 ).getX() ;
        double maxx = i->getBoundingPoint ( 0 ).getX() ;

        double miny = i->getBoundingPoint ( 0 ).getY() ;
        double maxy = i->getBoundingPoint ( 0 ).getY() ;

        double mint = i->getBoundingPoint ( 0 ).getT() ;
        double maxt = i->getBoundingPoint ( 0 ).getT() ;

        for (  ; i != t->end() ; i++ )
        {
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
            
        }

    }

    BoundaryCondition::apply(a,t);
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
        auto i = t->begin() ;
        while((i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR)
        {
            i++ ;
            if(i == t->end())
                return ;
        }

        double minx = i->getBoundingPoint ( 0 ).getX() ;
        double maxx = i->getBoundingPoint ( 0 ).getX() ;

        double miny = i->getBoundingPoint ( 0 ).getY() ;
        double maxy = i->getBoundingPoint ( 0 ).getY() ;

        double mint = i->getBoundingPoint ( 0 ).getT() ;
        double maxt = i->getBoundingPoint ( 0 ).getT() ;

        for (  ; i != t->end() ; i++ )
        {
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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

        double tol = std::max ( std::min ( std::min ( maxx - minx, maxy - miny ), maxt-mint ) * .001, POINT_TOLERANCE ) ;

        for ( auto i = t->begin() ; i != t->end() ; i++ )
        {
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured())  || i->getBehaviour()->type == VOID_BEHAVIOUR)
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
        }

    }

    BoundaryCondition::apply(a,t);

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

        auto i = t->begin() ;
        while((i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR)
        {
            i++ ;
            if(i == t->end())
                return ;
        }

        double minx = i->getBoundingPoint ( 0 ).getX() ;
        double maxx = i->getBoundingPoint ( 0 ).getX() ;

        double miny = i->getBoundingPoint ( 0 ).getY() ;
        double maxy = i->getBoundingPoint ( 0 ).getY() ;
        
        double minz = i->getBoundingPoint ( 0 ).getZ() ;
        double maxz = i->getBoundingPoint ( 0 ).getZ() ; 

        double mint = i->getBoundingPoint ( 0 ).getT() ;
        double maxt = i->getBoundingPoint ( 0 ).getT() ;

        for (  ; i != t->end() ; i++ )
        {
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR )
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
        }
    }
    
    BoundaryCondition::apply(a,t);

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

        auto i = t->begin() ;
        while((i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR)
        {
            i++ ;
            if(i == t->end())
                return ;
        }

        double minx = i->getBoundingPoint ( 0 ).getX() ;
        double maxx = i->getBoundingPoint ( 0 ).getX() ;

        double miny = i->getBoundingPoint ( 0 ).getY() ;
        double maxy = i->getBoundingPoint ( 0 ).getY() ;
        
        double minz = i->getBoundingPoint ( 0 ).getZ() ;
        double maxz = i->getBoundingPoint ( 0 ).getZ() ; 

        double mint = i->getBoundingPoint ( 0 ).getT() ;
        double maxt = i->getBoundingPoint ( 0 ).getT() ;

        for (  ; i != t->end() ; i++ )
        {
            if ( (i->getBehaviour()->getDamageModel() && i->getBehaviour()->getDamageModel()->fractured()) || i->getBehaviour()->type == VOID_BEHAVIOUR  )
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
                    Function effective = dataFunction*getScale() ;
                    apply3DBC ( i,gp,Jinv, cache.back(), condition, effective, a , axis ) ;
                }
            }
        }
    }

    BoundaryCondition::apply(a,t);
}

BoundaryCondition::BoundaryCondition ( LagrangeMultiplierType t, const double & d, int a ) :  condition ( t ),data ( d ), scale ( 1 ), active(false), dataInterpolation(nullptr), axis ( a ), function ( false ) { }

BoundaryCondition::BoundaryCondition ( LagrangeMultiplierType t, const Function & d, int a ) :  condition ( t ), scale ( 1 ), active(false), dataFunction ( d ), dataInterpolation(nullptr), axis ( a ), function ( true ) { }

void BoundaryCondition::setScale ( double d )
{
    if(active)
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
        DelaunayTreeItem * VoidItem = nullptr;
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
                        Function effective = dataFunction*getScale() ;
                        apply2DBC ( i,gp,Jinv, id, condition, effective, a ) ;
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
            if (i->getNeighbour ( j )->isSpace )
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
            Function effective = dataFunction*getScale() ;
            apply3DBC ( i,gp,Jinv, id, condition, effective, a ) ;
        }
    }
}


TimeContinuityBoundaryCondition::TimeContinuityBoundaryCondition ( double i ) : BoundaryCondition ( GENERAL, 0. ), initialValue ( i ), instant(1.), minDeltaTime(POINT_TOLERANCE)
{
    goToNext = true ;
    previousDisp.resize ( 0 ) ;
}

void TimeContinuityBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t )
{
//     t->getAdditionalPoints() ;
    auto j = t->begin() ;
    size_t dof = 0 ;
    double dt = 0 ;
    size_t timePlanes = 0;
    while ( dof == 0 )
    {

        if ( j == t->end() )
        {
            return ;
        }
        dof = j->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        dt = j->getState().getNodalDeltaTime() ;
        timePlanes = j->timePlanes() ;
        j++ ;
    }

//    std::cout << "dt=" << dt << std::endl ;


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
/*        size_t extraDofPerPlane = 0 ;
        if( ndofmax <  previousDisp.size()/dof )
            extraDofPerPlane = (previousDisp.size()/dof-ndofmax)/timePlanes ;*/

        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, previousDisp[ dofPerPlane* ( i ) *dof + j*dof + n]*(1.-instant) + previousDisp[ dofPerPlane* ( i + 1 ) *dof + j*dof + n]*instant, dofPerPlane*i + j, true )  ;
                }
            }

/*            for(size_t j = 0 ; j < extraDofPerPlane ; j++)
            {
               for ( size_t n = 0 ; n < dof ; n++ )
               {
                    a->setPointAlongIndexedAxis ( n, previousDisp[ ndofmax + extraDofPerPlane* ( i ) *dof + j*dof + n]*(1.-instant) + previousDisp[ ndofmax + extraDofPerPlane* ( i + 1 ) *dof + j*dof + n]*instant, ndofmax + extraDofPerPlane*i + j, true )  ;
               }
            }*/
        }

        if(dt < minDeltaTime)
        {
            size_t i = timePlanes-1 ;
            for(size_t j = 0 ; j < dofPerPlane ; j++)
            {
                for(size_t n = 2 ; n < dof ; n++)
                    a->setPointAlongIndexedAxis( n, previousDisp[ dofPerPlane*i*dof+j*dof+n ], dofPerPlane*i+j, true ) ;
            }
        }
    }



}

void TimeContinuityBoundaryCondition::apply ( Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t )
{
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

        for ( size_t i = 0 ; i < timePlanes-1 ; i++ )
        {
            for ( size_t j = 0 ; j < dofPerPlane ; j++ )
            {
                for ( size_t n = 0 ; n < dof ; n++ )
                {
                    a->setPointAlongIndexedAxis ( n, previousDisp[ dofPerPlane* ( i ) *dof + j*dof + n]*(1.-instant) + previousDisp[ dofPerPlane* ( i + 1 ) *dof + j*dof + n]*instant, dofPerPlane*i + j, true )  ;
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


double LoadingCycle::getValue() 
{
    double teststate = (condition == ULTIMATE_STRESS) ? (ft->getAverageField(REAL_STRESS_FIELD, 1., -1))[axisIndex] : (ft->getAverageField(STRAIN_FIELD, 1., -1))[axisIndex];
    if(!cycleStarted)
    {
       if(teststate > ultimate) 
           type = UNLOADING ;
       else
         type = LOADING ;
       
       if(chainedCycle)
           currentState = chainedCycle->getValue() ;
       cycleStarted = true ;
    }
    
    double dt = ft->getDeltaTime() ;
    
    if(type == LOADING)
    {
        if(teststate < ultimate)
        {
            currentState = std::min(currentState+rate*dt, ultimate) ;
            return currentState;
        }
        cycleAtEnd = true ;
        return currentState ;
    }
    
    if(teststate > ultimate)
    {
        currentState = std::min(currentState+rate*dt, ultimate) ;
        return currentState ;
    }
    cycleAtEnd = true ;
    return currentState ;
}
double LoadingCycle::isAtEnd() const 
{
    return cycleAtEnd ;
    
}
}
