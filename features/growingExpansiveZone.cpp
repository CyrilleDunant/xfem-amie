// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "growingExpansiveZone.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/dual_behaviour.h"
#include "../physics/homogeneised_behaviour.h"
#include "../physics/damagemodels/fractiondamage.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/materials/aggregate_behaviour.h"

using namespace Amie ;

GrowingExpansiveZone::GrowingExpansiveZone(Feature* father, const Function & g, double x, double y, Form* i, double start) : TimeDependentEnrichmentInclusion(father,g,x,y), imp(i)
{
    Feature::setBehaviour(i) ;
    changed = true ;
    this->TimeDependentCircle::center.getT() = start ;
}

GrowingExpansiveZone::GrowingExpansiveZone(Feature* father, const Function & g, double x, double y) : TimeDependentEnrichmentInclusion(father,g,x,y)
{
    changed = true ;
}

GrowingExpansiveZone::~GrowingExpansiveZone() {
//  delete imp ;
}

void GrowingExpansiveZone::enrich(size_t & counter, Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree)
{


    TimeDependentEnrichmentInclusion::enrich(counter,dtree) ;

    std::vector<DelaunayTriangle *> & disc = TimeDependentEnrichmentInclusion::cache ;

    pointsAndValues.clear() ;
    for(size_t i = 0 ; i < TimeDependentEnrichmentInclusion::cache.size() ; i++)
    {
        for(  size_t j = 0 ; j < TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunctions().size() ; j++)
        {
            if(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentSource(j) == getPrimitive())
            {
                std::vector<double> disp ;//(imp->getNumberOfDegreesOfFreedom()) ;
                if(TimeDependentEnrichmentInclusion::cache[i]->getState().getEnrichedDisplacements().size() && TimeDependentEnrichmentInclusion::cache[i]->getState().getEnrichedDisplacements().size() >= j*disp.size())
                {
                    for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
                    {
                        disp.push_back(TimeDependentEnrichmentInclusion::cache[i]->getState().getEnrichedDisplacements()[j*imp->getNumberOfDegreesOfFreedom()+n]);
                    }
                }
                else
                {
                    for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
                    {
                        disp.push_back(0.) ;
                    }
                }

                for(size_t k = 0 ; k < TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoints().size() ; k++)
                {

                    if(squareDist2D(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint(), &TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoint(k)) < POINT_TOLERANCE
                            && std::abs(TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint()->getT() - TimeDependentEnrichmentInclusion::cache[i]->getBoundingPoint(k).getT()) < POINT_TOLERANCE )
                    {
                        pointsAndValues[TimeDependentEnrichmentInclusion::cache[i]->getEnrichmentFunction(j).getPoint()] = disp ;
                        break ;
                    }
                }
            }
        }
    }

    std::set<Point *> donePoints ;
    for(size_t i = 0 ; i < disc.size() ; i++)
    {
        disc[i]->clearBoundaryConditions() ;
        if( disc[i]->getEnrichmentFunctions().size() == 0)
        {
            continue ;
        }


        GaussPointArray gp = disc[i]->getGaussPoints() ;
        std::valarray< Matrix > Jinv(Matrix(imp->getNumberOfDegreesOfFreedom(), imp->getNumberOfDegreesOfFreedom()), gp.gaussPoints.size()) ;

        for( size_t j = 0 ; j < gp.gaussPoints.size() ;  j++ )
        {
            disc[i]->getInverseJacobianMatrix( gp.gaussPoints[j].first, Jinv[j] ) ;
        }

        size_t enrichedDofPerTimePlanes = disc[i]->getEnrichmentFunctions().size()/disc[i]->timePlanes() ;

        for(size_t j = 0 ; j < disc[i]->getEnrichmentFunctions().size() - enrichedDofPerTimePlanes ; j++)
        {
            if(donePoints.find(disc[i]->getEnrichmentFunction(j).getPoint()) == donePoints.end())
            {
                donePoints.insert(disc[i]->getEnrichmentFunction(j).getPoint() ) ;
                if( true ) //pointsAndValues.find( disc[i]->getEnrichmentFunction(j+enrichedDofPerTimePlanes).getPoint() ) == pointsAndValues.end() )
                {
                    for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
                    {
                        disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, disc[i]->getEnrichmentFunction(j).getDofID(), 0., n) ) ;
                    }
                }
                else
                {
                    for(size_t n = 0 ; n < imp->getNumberOfDegreesOfFreedom() ; n++)
                    {
                        disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, disc[i]->getEnrichmentFunction(j).getDofID(), pointsAndValues[disc[i]->getEnrichmentFunction(j+enrichedDofPerTimePlanes).getPoint()][n], n) ) ;
                    }
                }
            }

//disc[i]->addBoundaryCondition( new DofDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, disc[i], gp, Jinv, dofIdCurrent[disc[i]->getEnrichmentFunction(j).getPoint()], 0., n) ) ;
        }
    }



// 	if(disc.size() < 6)
// 		return ;

    std::vector<DelaunayTriangle *> ring ;
    std::vector<DelaunayTriangle *> inDisc ;


    std::vector<std::vector<size_t > > nodesIterator ;
    size_t nodesPerPlane = disc[0]->getBoundingPoints().size() / disc[0]->timePlanes() ;
    size_t nodesPerEdge = nodesPerPlane/3 ;
    for(size_t i = 0 ; i < disc[0]->timePlanes() ; i++)
    {
        std::vector<size_t> tmp ;
        tmp.push_back(nodesPerPlane*i + nodesPerEdge*0);
        tmp.push_back(nodesPerPlane*i + nodesPerEdge*1);
        tmp.push_back(nodesPerPlane*i + nodesPerEdge*2);
        nodesIterator.push_back(tmp) ;
    }

    for(size_t i = 0 ; i < disc.size() ; i++)
    {
        bool added = false ;

        bool bin = false ;
        bool bout = false ;
        for(size_t j = 0 ; j < nodesIterator.size() && !added ; j++)
        {
            Point A = disc[i]->getBoundingPoint( nodesIterator[j][0] ) ;
            if(getPrimitive()->in(A))
                bin= true ;
            else
                bout = true ;
            Point B = disc[i]->getBoundingPoint( nodesIterator[j][1] ) ;
            if(getPrimitive()->in(B))
                bin= true ;
            else
                bout = true ;
            Point C = disc[i]->getBoundingPoint( nodesIterator[j][2] ) ;
            if(getPrimitive()->in(C))
                bin= true ;
            else
                bout = true ;

            if(bin && bout)
            {
                ring.push_back(disc[i]) ;
                added = true ;
                break ;
            }

        }
        if(!added)
        {
            Point p = disc[i]->getCenter() ;
            p.getT() = VirtualMachine().eval( disc[i]->getTTransform(), 0,0,0,0) ;
            if(getPrimitive()->in(p))
                inDisc.push_back( disc[i] ) ;
        }
    }

    std::set<DelaunayTriangle *> newInterface ;

    for( size_t i = 0 ; i < ring.size() ; i++ )
    {
        if( bimateralInterfaced.find( ring[i] ) == bimateralInterfaced.end() )
        {
            BimaterialInterface *bi = new BimaterialInterface( getPrimitive(),
                    imp->getCopy(),
                    ring[i]->getBehaviour()->getCopy() );

            const Geometry * src =  ring[i]->getBehaviour()->getSource() ;
//			delete ring[i]->getBehaviour() ;
            ring[i]->setBehaviour(dtree, bi ) ;
            bi->transform( ring[i]) ;
            bi->setSource( src );
        }
        else
            dynamic_cast<BimaterialInterface *>(ring[i]->getBehaviour())->transform( ring[i] ) ;


        newInterface.insert( ring[i] ) ;
    }


    std::set<DelaunayTriangle *> newExpansive ;

    for( size_t i = 0 ; i < inDisc.size() ; i++ )
    {
        if( expansive.find( inDisc[i] ) == expansive.end() )
        {
            inDisc[i]->setBehaviour(dtree, imp->getCopy() ) ;
            inDisc[i]->getBehaviour()->setSource( getPrimitive() );
        }

        newExpansive.insert( inDisc[i] ) ;
    }
    expansive = newExpansive ;
    bimateralInterfaced = newInterface ;

    std::vector<DelaunayTriangle *> enriched(enrichedElem.begin(), enrichedElem.end());

//    std::cout << enriched.size() << " " << expansive.size() << " " << bimateralInterfaced.size() << std::endl ;


    dofIdPrev = dofIdCurrent ;

}

void GrowingExpansiveZone::step(double dt, std::valarray< double >* results, Amie::Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree)
{
    TimeDependentEnrichmentInclusion::step(dt, results, dtree) ;

    if(dt > 0 )
    {

    }
}



