// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include <typeinfo>
#include "expansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/damagemodels/fractiondamage.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/materials/aggregate_behaviour.h"

using namespace Amie ;

ExpansiveZone::ExpansiveZone( Feature *father, double radius, double x, double y, const Matrix &tensor, Vector def ) : EnrichmentInclusion( father, radius, x, y ),  imposedDef( def ), cgTensor( tensor )
{
    setBehaviour( new StiffnessWithImposedDeformation( cgTensor, imposedDef ) ) ;
    homogeneized = false ;
}

ExpansiveZone::ExpansiveZone( Feature *father, double radius, double x, double y, StiffnessWithImposedDeformation * gel ) : EnrichmentInclusion( father, radius, x, y ),  imposedDef( gel->imposed ), cgTensor( gel->param )
{
    setBehaviour( gel->getCopy() ) ;
    homogeneized = false ;
}

ExpansiveZone::~ExpansiveZone() {}


void ExpansiveZone::reset()
{
    cache.clear() ;
    updated = true ;
}

void ExpansiveZone::addMeshPointsInFather()
{
    /*	PointArray back = getFather()->getInPoints() ;
    	if(back.size() < 2)
    		return ;
    	PointArray next = back ;
    	back.resize(back.size()+1) ;
    	for(size_t i = 0 ; i < next.size() ; i++)
    		back[i]=next[i] ;
    	back[next.size()] = new Point(this->getCenter()) ;
    	getFather()->getInPoints() = back ;*/
}


void ExpansiveZone::enrich( size_t &lastId , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
// 	this->setBehaviour( new StiffnessWithImposedDeformation( cgTensor, imposedDef ) ) ;
    EnrichmentInclusion::enrich( lastId, dtree) ;
    //first we get All the triangles affected
    std::vector<DelaunayTriangle *> & disc = EnrichmentInclusion::cache ;//dtree->getConflictingElements(getPrimitive()) ;

    if( disc.size() == 1 )
    {
        homogeneized = true ;
        return ;
    }
    homogeneized = false ;
    //then we select those that are cut by the circle
    std::vector<DelaunayTriangle *> ring ;
    std::vector<DelaunayTriangle *> inDisc ;

    for( size_t i = 0 ; i < disc.size() ; i++ )
    {
        Segment s0( *disc[i]->first, *disc[i]->second ) ;
        Segment s1( *disc[i]->second, *disc[i]->third ) ;
        Segment s2( *disc[i]->third, *disc[i]->first ) ;

        if( !( s0.intersection( getPrimitive() ).empty() && s1.intersection( getPrimitive() ).empty() && s2.intersection( getPrimitive() ).empty() ) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR )
            ring.push_back( disc[i] ) ;
        else if( in( disc[i]->getCenter() ) )
            inDisc.push_back( disc[i] ) ;
    }

    std::set<DelaunayTriangle *> newInterface ;

    for( size_t i = 0 ; i < ring.size() ; i++ )
    {
        if( bimateralInterfaced.find( ring[i] ) == bimateralInterfaced.end() )
        {
            BimaterialInterface *bi = nullptr ;

            if( dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() ) )
            {
                Matrix p = dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() )->getOriginalBehaviour()->getTensor(Point(1./3,1./3)) ;
                bi = new BimaterialInterface( getPrimitive(),
                                              new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
                                              dynamic_cast<HomogeneisedBehaviour *>( ring[i]->getBehaviour() )->getOriginalBehaviour()->getCopy() ) ;
            }
            else
            {
                bi = new BimaterialInterface( getPrimitive(),
                                              new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
                                              ring[i]->getBehaviour()->getCopy() ) ;
            }

            const Geometry * src =  ring[i]->getBehaviour()->getSource() ;
//			delete ring[i]->getBehaviour() ;
            ring[i]->setBehaviour( dtree, bi ) ;
            bi->transform( ring[i]) ;
            bi->setSource( src );
        }

        newInterface.insert( ring[i] ) ;
    }

    std::set<DelaunayTriangle *> newExpansive ;

    for( size_t i = 0 ; i < inDisc.size() ; i++ )
    {
        if( expansive.find( inDisc[i] ) == expansive.end() )
        {
            StiffnessWithImposedDeformation * bi = new StiffnessWithImposedDeformation( cgTensor, imposedDef ) ;
//			delete inDisc[i]->getBehaviour() ;
            inDisc[i]->setBehaviour(dtree, bi) ;
            inDisc[i]->getBehaviour()->setSource( getPrimitive() );
        }

        newExpansive.insert( inDisc[i] ) ;
    }

    expansive = newExpansive ;

    if( disc.size() == 1 )
    {
        std::cout << "SHOULD NEVER HAPPEN" << std::endl ;
        if( bimateralInterfaced.find( disc[0] ) == bimateralInterfaced.end() )
        {
            if( dynamic_cast<HomogeneisedBehaviour *>( disc[0]->getBehaviour() ) )
            {
                std::cout << "get original" << std::endl ;
                BimaterialInterface * bi = new BimaterialInterface( getPrimitive(),
                        new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
                        dynamic_cast<HomogeneisedBehaviour *>( disc[0]->getBehaviour() )->original->getCopy()
                                                                  ) ;
                delete disc[0]->getBehaviour() ;
                disc[0]->setBehaviour(dtree, bi) ;
            }
            else
            {
                BimaterialInterface * bi = new BimaterialInterface( getPrimitive(),
                        new StiffnessWithImposedDeformation( cgTensor, imposedDef ),
                        disc[0]->getBehaviour()->getCopy()
                                                                  ) ;
                delete disc[0]->getBehaviour() ;
                disc[0]->setBehaviour(dtree, bi) ;
            }

            disc[0]->getBehaviour()->transform( disc[0]) ;
            disc[0]->getBehaviour()->setSource( getPrimitive());
        }

        newInterface.insert( disc[0] ) ;
    }

    bimateralInterfaced = newInterface ;
}


void ExpansiveZone::setExpansion( Vector a )
{
    imposedDef = a ;
    updated = true ;

    for( auto i = bimateralInterfaced.begin() ; i != bimateralInterfaced.end() ; i++ )
    {
        BimaterialInterface *interface = static_cast<BimaterialInterface *>( ( *i )->getBehaviour() ) ;
        static_cast<StiffnessWithImposedDeformation *>( interface->inBehaviour )->imposed = a ;
    }

    for( auto i = expansive.begin() ; i != expansive.end() ; i++ )
    {
        static_cast<StiffnessWithImposedDeformation *>( ( *i )->getBehaviour() )->imposed = a ;
    }
}


MaterialInclusion::MaterialInclusion( Feature *father, double radius, double x, double y, LinearForm *inclusionBehaviour ) : EnrichmentInclusion( father, radius, x, y ),  inclusionBehaviour( inclusionBehaviour )
{

}

MaterialInclusion::~MaterialInclusion() {}

void MaterialInclusion::reset()
{
    cache.clear() ;
    updated = true ;
}

void MaterialInclusion::enrich( size_t &counter , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{
    EnrichmentInclusion::enrich( counter, dtree) ;
    //first we get All the triangles affected
    std::vector<DelaunayTriangle *> disc = cache ;

    //then we select those that are cut by the circle
    std::vector<DelaunayTriangle *> ring ;
    std::vector<DelaunayTriangle *> inDisc ;

    if( disc.size() < 2 )
        return ;


    for( size_t i = 0 ; i < disc.size() ; i++ )
    {
        if( this->intersects( static_cast<Triangle *>( disc[i] ) ) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR )
            ring.push_back( disc[i] ) ;
        else if( this->in( disc[i]->getCenter() ) && disc[i]->getBehaviour()->type != VOID_BEHAVIOUR )
            inDisc.push_back( disc[i] ) ;
    }

    std::set<DelaunayTriangle *> newInterface ;

    for( size_t i = 0 ; i < ring.size() ; i++ )
    {
        if( bimateralInterfaced.find( ring[i] ) == bimateralInterfaced.end() )
        {
            BimaterialInterface * bi = new BimaterialInterface( getPrimitive(),
                    inclusionBehaviour->getCopy(),
                    ring[i]->getBehaviour()->getCopy()) ;
            delete ring[i]->getBehaviour() ;
            ring[i]->setBehaviour(dtree, bi) ;
            bi->transform( ring[i]) ;
            bi->setSource( getPrimitive() );
        }

        newInterface.insert( ring[i] ) ;
    }

    std::set<DelaunayTriangle *> newExpansive ;

    for( size_t i = 0 ; i < inDisc.size() ; i++ )
    {
        if( internal.find( inDisc[i] ) == internal.end() )
        {
            delete inDisc[i]->getBehaviour() ;
            inDisc[i]->setBehaviour(dtree, inclusionBehaviour->getCopy()) ;
        }

        newExpansive.insert( inDisc[i] ) ;
    }

    internal = newExpansive ;

    if( disc.size() == 1 )
    {
        if( bimateralInterfaced.find( disc[0] ) == bimateralInterfaced.end() )
        {
            BimaterialInterface * bi = new BimaterialInterface( getPrimitive(),
                    inclusionBehaviour->getCopy(),
                    disc[0]->getBehaviour()->getCopy()
                                                              ) ;
            delete disc[0]->getBehaviour() ;
            disc[0]->setBehaviour(dtree, bi) ;
            disc[0]->getBehaviour()->transform( disc[0]) ;
        }

        newInterface.insert( disc[0] ) ;
    }

    bimateralInterfaced = newInterface ;

}



