

// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "features.h"
#include "crackinitiation.h"
#include "layeredinclusion.h"
#include "sample.h"
#include "sample3d.h"
#include "../utilities/configuration.h"
#include "../utilities/enumeration_translator.h"
#include "../physics/void_form.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "../physics/homogeneised_behaviour.h"
#include "../physics/diffusion.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../solvers/multigrid.h"
#include "../solvers/multigridstep.h"
#include "../mesher/parallel_delaunay.h"
#include "../utilities/behaviour_converter.h"
#include <time.h>
#include <sys/time.h>



namespace Amie {



Mesh<DelaunayTriangle, DelaunayTreeItem> * FeatureTree::get2DMesh ( int g )
{
//     state.setStateTo ( MESHED, false ) ;

    if ( g == -1 )
    {
        return dtree ;
    }
    else
    {
        return layer2d[g] ;
    }
}

Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * FeatureTree::get3DMesh ()
{
//     state.setStateTo ( MESHED, false ) ;

    return dtree3D ;

}

std::vector<DelaunayTriangle *> FeatureTree::getBoundingTriangles ( const Feature *f )
{
    state.setStateTo ( MESHED, false ) ;

    if ( f == nullptr )
    {
        return tree[0]->getBoundingElements2D ( this ) ;
    }
    else
    {
        return f->getBoundingElements2D ( this ) ;
    }
}

bool FeatureTree::getStateConvergence() const {
    return stateConverged ;
}

FeatureTree::FeatureTree ( Feature *first, int layer, double fraction, size_t gridsize ) :  state ( this ), nodes ( 0 ), grid ( nullptr ),grid3d ( nullptr )
{
    initialValue = 0 ;
    previousDeltaTime = 0 ;
    minDeltaTime = 0.001 ;
    reuseDisplacements = false ;
    foundCheckPoint = true ;
    averageDamage = 0 ;
    behaviourSet =false ;
    damageConverged = false ;
    stateConverged = false ;
    dtree = nullptr ;
    dtree3D = nullptr ;
    structuredMesh = false ;
    alternating = false ;
    maxScore = -1 ;

    std::vector<Point> bbox = first->getBoundingBox() ;
    double min_x = 0, min_y = 0, max_x = 0, max_y = 0, max_z = 0, min_z = 0;

    for ( size_t j  =  0 ; j <  bbox.size() ; j++ )
    {
        if ( bbox[j].getY() < min_y )
        {
            min_y = bbox[j].getY() ;
        }

        if ( bbox[j].getY() > max_y )
        {
            max_y = bbox[j].getY() ;
        }

        if ( bbox[j].getX() < min_x )
        {
            min_x = bbox[j].getX() ;
        }

        if ( bbox[j].getX() > max_x )
        {
            max_x = bbox[j].getX() ;
        }

        if ( bbox[j].getZ() < min_z )
        {
            min_z = bbox[j].getZ() ;
        }

        if ( bbox[j].getZ() > max_z )
        {
            max_z = bbox[j].getZ() ;
        }
    }

    samplingRestriction = 0 ;

    if ( first )
    {
        addFeature ( nullptr, first, layer, fraction ) ;
    }

    if ( is2D() )
    {
        grid = new Grid ( ( bbox[1].getX() - bbox[0].getX() ) * 1.1,
                          ( bbox[1].getY() - bbox[2].getY() ) * 1.1, gridsize,
                          Point ( ( bbox[1].getX() + bbox[0].getX() ) *.5,
                                  ( bbox[1].getY() + bbox[2].getY() ) *.5
                                ) ) ;
        domains.push_back ( new Rectangle ( max_x-min_x, max_y-min_y, ( max_x+min_x ) *.5, ( max_y+min_y ) *.5 ) );
    }

    if ( is3D() )
    {
        grid3d = new Grid3D ( ( bbox[7].getX() - bbox[0].getX() ) * 1.1,
                              ( bbox[7].getY() - bbox[0].getY() ) * 1.1,
                              ( bbox[7].getZ() - bbox[0].getZ() ) * 1.1, gridsize / 5, ( bbox[7] + bbox[0] ) *.5 );
        domains.push_back ( new Hexahedron ( max_x-min_x, max_y-min_y, max_z-min_y , ( max_x+min_x ) *.5, ( max_y+min_y ) *.5, ( max_z+min_z ) *.5 ) );
    }

    father3D = nullptr;
    father2D = nullptr ;
    elemOrder = LINEAR ;
    renumbered = false ;
    setBehaviours = false ;
    behaviourChange = true ;
    solverConvergence = false ;
    enrichmentChange = true ;
    needMeshing = true ;

    elastic = false ;
    projectOnBoundaries = true ;

    K = new Assembly() ;

    if ( is2D() )
    {
        K->setSpaceDimension ( SPACE_TWO_DIMENSIONAL ) ;
    }
    else
    {
        K->setSpaceDimension ( SPACE_THREE_DIMENSIONAL ) ;
    }

    crackedVolume = 0 ;
    damagedVolume = 0 ;
    residualError = 1e9 ;
    samplingNumber = 0 ;
    previousSamplingNumber = 0 ;


    lastNodeId = 0;
    lastEnrichmentId = 0;
    maxitPerStep = 200 ;
    deltaTime = .1 ;
    realDeltaTime = deltaTime ;
    now = 0 ;

    setElementGenerationMethod() ;

}

FeatureTree::FeatureTree ( const char * voxelSource, std::map<unsigned char,Form *> behaviourMap , const std::vector<double> & times ) :  state ( this ), nodes ( 0 ), grid ( nullptr ),grid3d ( nullptr )
{
    initialValue = 0 ;
    previousDeltaTime = 0 ;
    minDeltaTime = 0.001 ;
    reuseDisplacements = false ;
    foundCheckPoint = true ;
    averageDamage = 0 ;
    maxScore = -1 ;
    behaviourSet =true ;
    damageConverged = true ;
    stateConverged = false ;
    dtree = nullptr ;
    alternating = false ;
    if(times.empty())
        dtree3D = new MicDerivedMesh(voxelSource, behaviourMap) ;
    else
        dtree3D = new MicDerivedMesh(voxelSource, behaviourMap,times) ;
    structuredMesh = true ;

    samplingRestriction = 0 ;

    father3D = nullptr;
    father2D = nullptr ;
    elemOrder = LINEAR ;
    renumbered = false ;
    setBehaviours = true ;
    behaviourChange = false ;
    solverConvergence = false ;
    enrichmentChange = false ;
    needMeshing = true ;

    elastic = false ;
    projectOnBoundaries = false ;

    K = new Assembly() ;
    K->setSpaceDimension ( SPACE_THREE_DIMENSIONAL ) ;


    crackedVolume = 0 ;
    damagedVolume = 0 ;
    residualError = 1e9 ;
    samplingNumber = 0 ;
    previousSamplingNumber = 0 ;


    lastNodeId = 0;
    lastEnrichmentId = 0;
    maxitPerStep = 200 ;
    deltaTime = .1 ;
    realDeltaTime = deltaTime ;
    now = 0 ;

    setElementGenerationMethod() ;

}

void FeatureTree::setPartition ( size_t partitionNumber )
{
    if(structuredMesh)
        return ;

    for ( size_t i = 0 ; i < domains.size() ; i++ )
    {
        delete domains[i] ;
    }

    domains.clear() ;
    std::vector<Point> bbox = tree[0]->getBoundingBox() ;
    double min_x = bbox[0].getX(), min_y = bbox[0].getX(), max_x = bbox[0].getY(), max_y = bbox[0].getY(), max_z = bbox[0].getZ(), min_z = bbox[0].getZ();

    for ( size_t j  =  0 ; j <  bbox.size() ; j++ )
    {
        if ( bbox[j].getY() < min_y )
        {
            min_y = bbox[j].getY() ;
        }

        if ( bbox[j].getY() > max_y )
        {
            max_y = bbox[j].getY() ;
        }

        if ( bbox[j].getX() < min_x )
        {
            min_x = bbox[j].getX() ;
        }

        if ( bbox[j].getX() > max_x )
        {
            max_x = bbox[j].getX() ;
        }

        if ( bbox[j].getZ() < min_z )
        {
            min_z = bbox[j].getZ() ;
        }

        if ( bbox[j].getZ() > max_z )
        {
            max_z = bbox[j].getZ() ;
        }
    }

    if ( is3D() )
    {
        std::map<double, std::vector<int>> triplet ;
        for ( int i = 0 ; i < (int)partitionNumber ; i++ )
        {
            for ( int j = 0 ; j < (int)partitionNumber ; j++ )
            {
                for ( int k = 0 ; k < (int)partitionNumber ; k++ )
                {
                    double avg = ( i+1+j+1+k+1 ) /3. ;
                    double sigma = ( i+1-avg ) * ( i+1-avg ) + ( j+1-avg ) * ( j+1-avg ) + ( k+1-avg ) * ( k+1-avg ) ;
                    if (  ( i+1 ) * ( j+1 ) * ( k+1 ) == (int)partitionNumber )
                        triplet[sigma] = {i+1,j+1,k+1} ;
                }
            }
        }


        double lx = max_x-min_x ;
        double ly = max_y-min_y ;
        double lz = max_z-min_z ;
        double di = triplet.begin()->second[0] ;
        double dj = triplet.begin()->second[1] ;
        double dk = triplet.begin()->second[2] ;
        std::cerr << "New partition : " << di << "x" << dj << "x" << dk << ", " << di*dj*dk << " domains. ...done." << std::endl ;
        for ( size_t i = 0 ; i < di ; i++ )
        {
            for ( size_t j = 0 ; j < dj ; j++ )
            {
                for ( size_t k = 0 ; k < dk ; k++ )
                {
                    domains.push_back ( new Hexahedron ( lx/di, ly/dj, lz/dk, lx/ ( 2.*di ) +lx/di*i+min_x, ly/ ( 2.*dj ) +ly/dj*j+min_y, lz/ ( 2.*dk ) +lz/dk*k+min_z ) );
                }
            }
        }
    }
    else
    {
        std::map<double, std::vector<int>> triplet ;
        for ( int i = 0 ; i < (int)partitionNumber ; i++ )
        {
            for ( int j = 0 ; j < (int)partitionNumber ; j++ )
            {
                double avg = ( i+1+j+1 ) /2. ;
                double sigma = ( i+1-avg ) * ( i+1-avg ) + ( j+1-avg ) * ( j+1-avg ) ;
                if ( ( i+1 ) * ( j+1 ) == (int)partitionNumber )
                    triplet[sigma] = {i+1,j+1} ;
            }
        }


        double lx = max_x-min_x ;
        double ly = max_y-min_y ;
        double di = triplet.begin()->second[0] ;
        double dj = triplet.begin()->second[1] ;
        std::cerr << "New partition : " << di << "x" << dj << ", " << di*dj << " domains. ...done." << std::endl ;
        for ( int i = 0 ; i < di ; i++ )
        {
            for ( int j = 0 ; j < dj ; j++ )
            {
                domains.push_back ( new Rectangle ( lx/di, ly/dj, lx/ ( 2.*di ) +lx/di*i+min_x, ly/ ( 2.*dj ) +ly/dj*j+min_y ) );
            }
        }
    }
}

void FeatureTree::twineFeature ( CompositeFeature *father, CompositeFeature *f )
{
    std::vector<Feature *> chain = father->getDescendants() ;
    std::vector<VirtualFeature *> fatherComponents = father->getComponents() ;
    std::vector<VirtualFeature *> childComponents = f->getComponents() ;

    //we look for the father of the second descendant of the "father" feature
    Feature *headerEnd = ( *std::find ( chain.begin(), chain.end(), fatherComponents[0] ) ) ;
    Feature *headerEndChild = ( *headerEnd->getChildren().rbegin() ) ;

    //we attach the first component of the twinee to the twinee
    addFeature ( f, childComponents[0] ) ;

    //the header end loses its child and gets a new one
    headerEnd->removeChild ( *headerEnd->getChildren().rbegin() );
    addFeature ( headerEnd, f ) ;

    //we re-attach the new header to the rest
    childComponents[0]->addChild ( headerEndChild ) ;
    headerEndChild->setFather ( childComponents[0] ) ;

    //now that the header is complete, we find each father component in the chain, and insert at that
    //point a new link: the corresponding child component

    for ( size_t i = 1 ; i < fatherComponents.size() ; i++ )
    {
        Feature *attachPoint = ( *std::find ( chain.begin(), chain.end(), fatherComponents[i] ) ) ;
        Feature *attachPointChild = nullptr;

        if ( !attachPoint->getChildren().empty() )
        {
            attachPointChild = ( *attachPoint->getChildren().rbegin() ) ;
        }

        //we detach the attach point from its father
        attachPoint->removeChild ( attachPointChild ) ;

        //we attach the corresponding child component to the loose end
        addFeature ( attachPoint, childComponents[i] ) ;

        //we re-attach the attach point to the new feature
        if ( attachPointChild )
        {
            childComponents[i]->addChild ( attachPointChild ) ;
            attachPointChild->setFather ( childComponents[i] ) ;

        }

    }
}

void FeatureTree::addPoint ( Point *p )
{
    extraPoints.push_back ( p );
}

void FeatureTree::addFeature ( Feature *father, Feature *f, int layer, double fraction )
{
    f->setLayer ( layer ) ;
    f->setFraction ( fraction ) ;
    scalingFactors[layer] = fraction ;
    if ( !f->isEnrichmentFeature )
    {
        needMeshing = true ;
    }

    if ( !tree.empty() && f->spaceDimensions() == SPACE_TWO_DIMENSIONAL && !f->isEnrichmentFeature )
    {
        grid->forceAdd ( f ) ;
    }
    else if ( !tree.empty() && !f->isEnrichmentFeature )
    {
        grid3d->forceAdd ( f ) ;
    }

    if ( f->isCompositeFeature && father && !father->isCompositeFeature )
    {
        std::vector<VirtualFeature *> pile = dynamic_cast<CompositeFeature *> ( f )->getComponents();
        f->setFather ( father ) ;

        if ( father != nullptr )
        {
            father->addChild ( f ) ;
        }

        this->tree.push_back ( f ) ;
        addFeature ( f, pile[0] ) ;

        for ( size_t i = 0 ; i < pile.size() - 1 ; i++ )
        {
            addFeature ( pile[i], pile[i + 1] ) ;
        }

        return ;

    }

    f->setFather ( father ) ;


    if ( father != nullptr )
    {
        father->addChild ( f ) ;
    }

    this->tree.push_back ( f ) ;


}

void FeatureTree::addManager ( EnrichmentManager * fm)
{
    enrichmentManagers.push_back(fm);
}

FeatureTree::~FeatureTree()
{
    delete father3D ;
    delete father2D ;
    delete grid ;
    delete grid3d ;
    for ( auto j = layer2d.begin() ; j!=layer2d.end() ; ++j )
    {
        delete j->second ;
    }

    delete this->dtree3D ;
    delete this->K ;

    for ( size_t i = 0 ; i < additionalPoints.size() ; i++ )
    {
        delete additionalPoints[i] ;
    }

    for ( size_t i = 0 ; i < boundaryCondition.size() ; ++i )
    {
        delete boundaryCondition[i] ;
    }
    for ( size_t i = 0 ; i < extraPoints.size() ; ++i )
    {
        delete extraPoints[i] ;
    }

}

void FeatureTree::scaleBoundaryConditions ( double scale )
{
    for ( size_t i = 0 ; i < boundaryCondition.size() ; i++ )
    {
        if(dynamic_cast<TimeContinuityBoundaryCondition *>( boundaryCondition[i]) == nullptr)
            boundaryCondition[i]->setScale ( scale );
    }
}

void FeatureTree::addBoundaryCondition ( BoundaryCondition *bc )
{
    boundaryCondition.push_back ( bc ) ;
}

void FeatureTree::removeBoundaryCondition ( BoundaryCondition *bc )
{
    std::vector<BoundaryCondition *>::iterator toDelete = std::find ( boundaryCondition.begin(), boundaryCondition.end(), bc ) ;
    if(toDelete != boundaryCondition.end())
        boundaryCondition.erase ( toDelete ) ;
}

void FeatureTree::setOrder ( Order ord )
{

    if(ord == elemOrder)
        return ;



    state.stitched = false ;
    state.renumbered = false ;
    state.initialised = false ;

    elemOrder = ord ;

    if ( father3D )
    {
        delete father3D ;
    }

    father3D = new TetrahedralElement ( elemOrder ) ;
//     father3D->compileAndPrecalculate() ;


    if ( father2D )
    {
        delete father2D ;
    }

    father2D = new TriElement ( elemOrder ) ;
//     father2D->compileAndPrecalculate() ;

    if ( ord >= CONSTANT_TIME_LINEAR )
    {
        addBoundaryCondition ( new TimeContinuityBoundaryCondition ( initialValue ) ) ;
    }

}

void FeatureTree::renumber()
{
    for(auto & bc : boundaryCondition)
        bc->clearCache() ;
    if ( is2D() )
    {
        size_t count = 0 ;
        std::cerr << " renumbering... " << std::flush ;

        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                i->getBoundingPoint ( j ).setId ( -1 ) ;
            }
        }

        for (auto i = dtree->begin() ; i != dtree->end() ; i++  )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() / i->timePlanes() ; j++ )
            {
                if ( i->getBoundingPoint ( j ).getId() == -1 )
                {
                    i->getBoundingPoint ( j ).setId ( count++ ) ;
                }
            }
        }


        lastNodeId = count ;

        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            for ( size_t k = 1 ; k < i->timePlanes() ; k++ )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() / i->timePlanes() ; j++ )
                {
                    if ( i->getBoundingPoint ( j + k* i->getBoundingPoints().size() / i->timePlanes() ).getId() == -1 )
                    {
                        i->getBoundingPoint ( j + k* i->getBoundingPoints().size() / i->timePlanes() ).setId ( i->getBoundingPoint ( j ).getId() + lastNodeId*k ) ;
                        count++ ;
                    }
                }
            }
        }

        lastNodeId = count ;

        nodes.resize ( count ) ;
        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {

            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                nodes[ i->getBoundingPoint ( j ).getId() ] = & i->getBoundingPoint ( j ) ;
            }


        }

        std::cerr << count * 2 << " ...done " << std::endl ;

    }
    else if ( is3D() )
    {
        size_t count = 0 ;
        std::cerr << " renumbering... " << std::flush ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                i->getBoundingPoint ( j ).setId ( -1 ) ;
            }
        }

        for (auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() / i->timePlanes() ; j++ )
            {
                if ( i->getBoundingPoint ( j ).getId() == -1 )
                {
                    i->getBoundingPoint ( j ).setId ( count++ ) ;
                }
            }
        }


        lastNodeId = count ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            for ( size_t k = 1 ; k < i->timePlanes() ; k++ )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() / i->timePlanes() ; j++ )
                {
                    if ( i->getBoundingPoint ( j + k* i->getBoundingPoints().size() / i->timePlanes() ).getId() == -1 )
                    {
                        i->getBoundingPoint ( j + k* i->getBoundingPoints().size() / i->timePlanes() ).setId ( i->getBoundingPoint ( j ).getId() + lastNodeId*k ) ;
                        count++ ;
                    }
                }
            }
        }

        lastNodeId = count ;

        nodes.resize ( count ) ;
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {

            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                nodes[ i->getBoundingPoint ( j ).getId() ] = & i->getBoundingPoint ( j ) ;
            }


        }

        std::cerr << count * 3 << " ...done " << std::endl ;


    }

    renumbered = true ;

}

bool FeatureTree::inRoot ( const Point &p ) const
{
    if(structuredMesh)
        return true ;

    if ( is2D() )
    {
        Point p0 ( p.getX(), p.getY() + POINT_TOLERANCE ) ;
        Point p1 ( p.getX(), p.getY() - POINT_TOLERANCE ) ;
        Point p2 ( p.getX() + POINT_TOLERANCE, p.getY() ) ;
        Point p3 ( p.getX() - POINT_TOLERANCE, p.getY() ) ;
        return ( tree[0]->in ( p ) || tree[0]->in ( p0 ) || tree[0]->in ( p1 ) || tree[0]->in ( p2 ) || tree[0]->in ( p3 ) ) ;
    }
    else
    {
        Point p0 ( p.getX(), p.getY() + POINT_TOLERANCE, p.getZ() ) ;
        Point p1 ( p.getX(), p.getY() - POINT_TOLERANCE, p.getZ() ) ;
        Point p2 ( p.getX() + POINT_TOLERANCE, p.getY(), p.getZ() ) ;
        Point p3 ( p.getX() - POINT_TOLERANCE, p.getY(), p.getZ() ) ;
        Point p4 ( p.getX(), p.getY(), p.getZ() + POINT_TOLERANCE ) ;
        Point p5 ( p.getX(), p.getY(), p.getZ() - POINT_TOLERANCE ) ;
        return ( tree[0]->in ( p ) || tree[0]->in ( p0 ) || tree[0]->in ( p1 ) || tree[0]->in ( p2 ) || tree[0]->in ( p3 ) || tree[0]->in ( p4 ) || tree[0]->in ( p5 ) ) ;
    }
}

void FeatureTree::projectTetrahedronsOnBoundaries ( size_t edge, size_t time )
{
    if(structuredMesh)
        return ;

    if ( edge + time == 0 )
    {
        return ;
    }

    std::vector<size_t> indexes ( edge * ( time + 1 ) * 6 ) ;



    size_t count = 0 ;


    Point a ( 0.25, 0.25, 0.25 ) ;
    Point b ( 0.166666666666667, 0.166666666666667, 0.166666666666667 ) ;
    Point c ( 0.5, 0.166666666666667, 0.166666666666667 ) ;
    Point d ( 0.166666666666667, 0.5, 0.166666666666667 ) ;
    Point e ( 0.166666666666667, 0.166666666666667, 0.5 ) ;

    for ( size_t j = 1 ; j < tree.size() ; j++ )
    {
        if ( !tree[j]->isEnrichmentFeature )
        {
            //In two pass
            std::vector<DelaunayTetrahedron *> tets = tree[j]->getElements3D ( this ) ;
            std::valarray<Point> originalPoints ( tets[0]->getBoundingPoints().size() ) ;

            for ( size_t i = 0 ; i < tets.size() ; i++ )
            {

                tets[i]->refresh ( father3D ) ;
                Point proj_0 ( *tets[i]->first ) ;
                tree[j]->project ( &proj_0 ) ;
                Point proj_1 ( *tets[i]->second ) ;
                tree[j]->project ( &proj_1 ) ;
                Point proj_2 ( *tets[i]->third ) ;
                tree[j]->project ( &proj_2 ) ;
                Point proj_3 ( *tets[i]->fourth ) ;
                tree[j]->project ( &proj_3 ) ;

                if (
                    squareDist3D ( proj_0 , *tets[i]->first ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_1 , *tets[i]->second ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;
                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->first->getX() +tets[i]->second->getX() ),0.5* ( tets[i]->first->getY() +tets[i]->second->getY() ),0.5* ( tets[i]->first->getZ() +tets[i]->second->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if (
                    squareDist3D ( proj_1 , *tets[i]->second ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_2 , *tets[i]->third ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;

                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->third->getX() +tets[i]->second->getX() ),0.5* ( tets[i]->third->getY() +tets[i]->second->getY() ),0.5* ( tets[i]->third->getZ() +tets[i]->second->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if (
                    squareDist3D ( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_2 , *tets[i]->third ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;
                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->third->getX() +tets[i]->fourth->getX() ),0.5* ( tets[i]->third->getY() +tets[i]->fourth->getY() ),0.5* ( tets[i]->third->getZ() +tets[i]->fourth->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if (
                    squareDist3D ( proj_0 , *tets[i]->first ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;
                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->first->getX() +tets[i]->fourth->getX() ),0.5* ( tets[i]->first->getY() +tets[i]->fourth->getY() ),0.5* ( tets[i]->first->getZ() +tets[i]->fourth->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if (
                    squareDist3D ( proj_1 , *tets[i]->second ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;
                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->second->getX() +tets[i]->fourth->getX() ),0.5* ( tets[i]->second->getY() +tets[i]->fourth->getY() ),0.5* ( tets[i]->second->getZ() +tets[i]->fourth->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if (
                    squareDist3D ( proj_0 , *tets[i]->first ) < POINT_TOLERANCE * POINT_TOLERANCE &&
                    squareDist3D ( proj_2 , *tets[i]->third ) < POINT_TOLERANCE * POINT_TOLERANCE
                )
                {
                    count++;
                    indexes.clear() ;
                    for ( size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++ )
                    {
                        if ( tets[i]->getBoundingPoint ( k ) == Point ( 0.5* ( tets[i]->first->getX() +tets[i]->third->getX() ),0.5* ( tets[i]->first->getY() +tets[i]->third->getY() ),0.5* ( tets[i]->first->getZ() +tets[i]->third->getZ() ),tets[i]->getBoundingPoint ( k ).getT() ) )
                        {
                            indexes.push_back ( k ) ;
                        }
                    }
                    Point test = tets[i]->getBoundingPoint ( indexes[0] ) ;
                    tree[j]->project ( &test ) ;

                    if ( inRoot ( test ) )
                    {
                        double dtest = 0 ;
                        for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                        {
                            size_t k = indexes[ni] ;
                            originalPoints[ni] = tets[i]->getBoundingPoint ( k ) ;
                            tree[j]->project ( &tets[i]->getBoundingPoint ( k ) ) ;
                            dtest = std::max(dtest, dist(originalPoints[ni],tets[i]->getBoundingPoint ( k ))) ;
                        }

                        if ( tets[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE &&
                                tets[i]->jacobianAtPoint ( e ) > POINT_TOLERANCE &&
                                dtest < 0.1*tets[i]->getRadius()
                           )
                        {
                            tets[i]->moved = true ;

                            for ( size_t j = 0 ; j < 4 ; j++ )
                            {
                                if ( tets[i]->getNeighbour ( j )->isTetrahedron )
                                {
                                    dynamic_cast<DelaunayTetrahedron *> ( tets[i]->getNeighbour ( j ) )->moved = true ;
                                }
                            }
                        }
                        else
                        {
                            for ( size_t ni = 0 ; ni < indexes.size() ; ni++ )
                            {
                                size_t k = indexes[ni] ;
                                tets[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                            }
                        }
                    }
                }

                if ( count % 1000 == 0 )
                {
                    std::cerr << "\r projecting points on boundaries... point " << count << "/xx" << " feature " << j << std::flush ;
                }
            }


        }
    }

    std::cerr << "\r projecting points on boundaries... point " << count << "/xx" << " ...done." << std::endl ;

}

void FeatureTree::projectTrianglesOnBoundaries ( size_t edge, size_t time )
{
    if(structuredMesh)
        return ;

    if ( edge + time == 0 )
    {
        return ;
    }

    std::valarray<size_t> indexes ( edge * ( time + 1 ) * 3 ) ;

    for ( size_t j = 0 ; j < 3 ; j++ )
    {
        for ( size_t e = 0 ; e < edge ; e++ )
        {
            indexes[3 * e + j] = j * ( edge + 1 ) + e + 1 ;

            for ( size_t t = 0 ; t < time ; t++ )
            {
                indexes[3 * ( edge * ( t + 1 ) + e ) + j] = j * ( edge + 1 ) + e + 1 + 3 * ( edge + 1 ) * ( t + 1 ) ;
            }
        }
    }

    size_t count = 0 ;
    size_t pd = 0 ;
    size_t k = 0 ;
    size_t n = indexes.size() / 3 ;

    std::valarray<Point> originalPoints ( n ) ;

    Point a ( 0.2, 0.2 ) ;
    Point b ( 0.6, 0.2 ) ;
    Point c ( 0.2, 0.6 ) ;
    Point d ( 1. / 3., 1. / 3. ) ;

    for ( size_t j = 1 ; j < this->tree.size() ; j++ )
    {
        if ( !tree[j]->isEnrichmentFeature && tree[j]->getGeometryType() != TRIANGLE && tree[j]->getGeometryType() != RECTANGLE )
        {

            std::vector<DelaunayTriangle *> triangles = this->tree[j]->getElements2D ( this ) ;

            for ( size_t i = 0 ; i < triangles.size() ; i++ )
            {
                triangles[i]->refresh ( father2D ) ;

                if ( triangles[i]->getPrimitive()->intersects ( tree[j] ) )
                {

                    Point proj_0 ( *triangles[i]->first ) ;
                    tree[j]->project ( &proj_0 ) ;
                    Point proj_1 ( *triangles[i]->second ) ;
                    tree[j]->project ( &proj_1 ) ;
                    Point proj_2 ( *triangles[i]->third ) ;
                    tree[j]->project ( &proj_2 ) ;
                    bool changed  = true;

                    if ( squareDist2D ( &proj_0 , triangles[i]->first ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_1 , triangles[i]->second ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_2 , triangles[i]->third ) > 10.*POINT_TOLERANCE )
                    {
                        count += changed ;
                        changed = false ;
                        Point test = triangles[i]->getBoundingPoint ( indexes[0] ) ;
                        tree[j]->project ( &test ) ;

                        if ( inRoot ( test ) )
                        {
                            for ( size_t ni = 0 ; ni < n ; ni++ )
                            {
                                k = indexes[3 * ni + 0] ;
                                originalPoints[ni] = triangles[i]->getBoundingPoint ( k ) ;
                                tree[j]->project ( &triangles[i]->getBoundingPoint ( k ) ) ;
                            }

                            if ( triangles[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE
                               )
                            {
                                triangles[i]->moved = true ;

                                for ( size_t j = 0 ; j < 3 ; j++ )
                                {
                                    if ( triangles[i]->getNeighbour ( j )->isTriangle )
                                    {
                                        dynamic_cast<DelaunayTriangle *> ( triangles[i]->getNeighbour ( j ) )->moved = true ;
                                    }
                                }
                            }
                            else
                            {
                                for ( size_t ni = 0 ; ni < n ; ni++ )
                                {
                                    k = indexes[3 * ni + 0] ;
                                    triangles[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                                }
                            }
                        }

// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(1)->getX() << ", " << (*triangles)[i]->getBoundingPoint(1)->getY() << std::endl ;
                    }

                    if ( squareDist2D ( &proj_1 , triangles[i]->second ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_2 , triangles[i]->third ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_0 , triangles[i]->first ) > 10.*POINT_TOLERANCE
                       )
                    {
                        count += changed ;
                        changed = false ;
                        Point test = triangles[i]->getBoundingPoint ( indexes[1] ) ;
                        tree[j]->project ( &test ) ;

                        if ( inRoot ( test ) )
                        {
                            for ( size_t ni = 0 ; ni < n ; ni++ )
                            {
                                k = indexes[3 * ni + 1] ;
                                originalPoints[ni] = triangles[i]->getBoundingPoint ( k ) ;
                                tree[j]->project ( &triangles[i]->getBoundingPoint ( k ) ) ;
                            }

                            if ( triangles[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE
                               )
                            {
                                triangles[i]->moved = true ;

                                for ( size_t j = 0 ; j < 3 ; j++ )
                                {
                                    if ( triangles[i]->getNeighbour ( j )->isTriangle )
                                    {
                                        dynamic_cast<DelaunayTriangle *> ( triangles[i]->getNeighbour ( j ) )->moved = true ;
                                    }
                                }
                            }
                            else
                            {
                                for ( size_t ni = 0 ; ni < n ; ni++ )
                                {
                                    k = indexes[3 * ni + 1] ;
                                    triangles[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                                }
                            }
                        }

// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(3)->getX() << ", " << (*triangles)[i]->getBoundingPoint(3)->getY() << std::endl ;
                    }

                    if ( squareDist2D ( &proj_2 , triangles[i]->third ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_0, triangles[i]->first ) < POINT_TOLERANCE &&
                            squareDist2D ( &proj_1, triangles[i]->second ) > 10.*POINT_TOLERANCE
                       )
                    {
                        count += changed ;
                        changed = false ;
                        Point test = triangles[i]->getBoundingPoint ( indexes[2] ) ;
                        tree[j]->project ( &test ) ;

                        if ( inRoot ( test ) )
                        {
                            for ( size_t ni = 0 ; ni < n ; ni++ )
                            {
                                k = indexes[3 * ni + 2] ;
                                originalPoints[ni] = triangles[i]->getBoundingPoint ( k ) ;
                                tree[j]->project ( &triangles[i]->getBoundingPoint ( k ) ) ;
                            }

                            if ( triangles[i]->jacobianAtPoint ( a ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( b ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( c ) > POINT_TOLERANCE &&
                                    triangles[i]->jacobianAtPoint ( d ) > POINT_TOLERANCE
                               )
                            {
                                triangles[i]->moved = true ;

                                for ( size_t j = 0 ; j < 3 ; j++ )
                                {
                                    if ( triangles[i]->getNeighbour ( j )->isTriangle )
                                    {
                                        dynamic_cast<DelaunayTriangle *> ( triangles[i]->getNeighbour ( j ) )->moved = true ;
                                    }
                                }
                            }
                            else
                            {
                                for ( size_t ni = 0 ; ni < n ; ni++ )
                                {
                                    k = indexes[3 * ni + 2] ;
                                    triangles[i]->getBoundingPoint ( k ) = originalPoints[ni] ;
                                }
                            }
                        }
                    }

                }

                if ( count % 1000 == 0 )
                {
                    std::cerr << "\r projecting points on boundaries... triangle " << count << "/" << triangles.size() << " feature " << i << std::flush ;
                }

            }
        }
    }

    std::cerr << "\r projecting points on boundaries... point " << count << "/" << pd << " ...done." << std::endl ;
}

void FeatureTree::stitch()
{

    for(auto & bc : boundaryCondition)
        bc->clearCache() ;

    if ( is2D() )
    {
        if ( elemOrder >= QUADRATIC )
        {

            layer2d.begin()->second->setElementOrder ( elemOrder, realDeltaTime ) ;
            for ( auto i = ++layer2d.begin() ; i != layer2d.end() ; i++ )
            {
                i->second->clearCaches();
                dynamic_cast<DelaunayTree *> ( i->second )->addSharedNodes ( dynamic_cast<DelaunayTree *> ( layer2d.begin()->second ) ) ;
            }

            if ( projectOnBoundaries && !structuredMesh)
            {
                switch ( elemOrder )
                {
                case QUADRATIC:
                    projectTrianglesOnBoundaries ( 1, 0 ) ;
                    break ;
                case CUBIC:
                    projectTrianglesOnBoundaries ( 2, 0 ) ;
                    break ;
                case QUADRIC:
                case QUINTIC:
                    projectTrianglesOnBoundaries ( 3, 0 ) ;
                    break ;
                case CONSTANT_TIME_LINEAR:
                case LINEAR_TIME_LINEAR:
                    break ;
                case QUADRATIC_TIME_LINEAR:
                    projectTrianglesOnBoundaries ( 1, 1 ) ;
                    break ;
                case QUADRATIC_TIME_QUADRATIC:
                    projectTrianglesOnBoundaries ( 1, 2 ) ;
                    break ;
                case CUBIC_TIME_LINEAR:
                    projectTrianglesOnBoundaries ( 2, 1 ) ;
                    break ;
                case CUBIC_TIME_QUADRATIC:
                    projectTrianglesOnBoundaries ( 2, 2 ) ;
                    break ;
                case QUADRIC_TIME_LINEAR:
                case QUINTIC_TIME_LINEAR:
                    projectTrianglesOnBoundaries ( 3, 1 ) ;
                    break ;
                case QUADRIC_TIME_QUADRATIC:
                case QUINTIC_TIME_QUADRATIC:
                    projectTrianglesOnBoundaries ( 3, 2 ) ;
                    break ;
                default:
                    std::cerr << "unsupported order" << std::endl ;
                    exit(0) ;
                    break ;
                }
            }
        }
    }
    else if ( is3D() )
    {
        if ( elemOrder >= QUADRATIC )
        {

            dtree3D->setElementOrder ( elemOrder, realDeltaTime ) ;
            dtree3D->clearCaches();

            if ( projectOnBoundaries )
            {

                switch ( elemOrder )
                {
                case QUADRATIC:
                    projectTetrahedronsOnBoundaries ( 1, 0 ) ;
                    break ;
                case CUBIC:
                    projectTetrahedronsOnBoundaries ( 2, 0 ) ;
                    break ;
                case QUADRIC:
                case QUINTIC:
                    projectTetrahedronsOnBoundaries ( 3, 0 ) ;
                    break ;
                case CONSTANT_TIME_LINEAR:
                case LINEAR_TIME_LINEAR:
                    break ;
                case QUADRATIC_TIME_LINEAR:
                    projectTetrahedronsOnBoundaries ( 1, 1 ) ;
                    break ;
                case QUADRATIC_TIME_QUADRATIC:
                    projectTetrahedronsOnBoundaries ( 1, 2 ) ;
                    break ;
                case CUBIC_TIME_LINEAR:
                    projectTetrahedronsOnBoundaries ( 2, 1 ) ;
                    break ;
                case CUBIC_TIME_QUADRATIC:
                    projectTetrahedronsOnBoundaries ( 2, 2 ) ;
                    break ;
                case QUADRIC_TIME_LINEAR:
                case QUINTIC_TIME_LINEAR:
                    projectTetrahedronsOnBoundaries ( 3, 1 ) ;
                    break ;
                case QUADRIC_TIME_QUADRATIC:
                case QUINTIC_TIME_QUADRATIC:
                    projectTetrahedronsOnBoundaries ( 3, 2 ) ;
                    break ;
                default:
                    std::cerr << "unsupported order" << std::endl ;
                    exit(0) ;
                    break ;
                }
            }

        }
    }

    if ( instants.size() > 2 && elemOrder >= CONSTANT_TIME_LINEAR && is2D() )
        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            i->second->extrude ( instants ) ;
        }
    else if ( instants.size() > 2 && elemOrder >= CONSTANT_TIME_LINEAR && is3D() )
        dtree3D->extrude(instants);

}

void FeatureTree::setSamplingNumber ( size_t news )
{
    samplingNumber = news ;
    needMeshing = true ;
    state.enriched = false ;
}

void FeatureTree::quadTreeRefine ( const Geometry * location )
{
    if(structuredMesh)
        return ;

    if ( location->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
    {
        std::cerr << "quadtree refine... " << std::flush ;
        int cacheID = dtree->begin().getId() ;
        if(location)
        {
            cacheID = dtree->generateCache(location) ;
        }

        std::cerr << dtree->begin(cacheID).size() << " elements... " << std::flush ;
        std::vector<Point> pointsToAdd ;
        std::vector<Point> illegalPoints ;
        for ( auto i = dtree->begin(cacheID) ; i != dtree->end(cacheID) ; i++ )
        {
            if ( !MinimumAngle ( M_PI/7. ).meetsCriterion ( i ) )
            {
                bool inrefinedFeature= false ;
                for ( size_t j = 0 ; j < refinedFeatures.size() ; j++ )
                {
                    if ( !refinedFeatures[j]->isVirtualFeature && refinedFeatures[j]->in ( i->getCenter() ) )
                    {
                        inrefinedFeature = true ;
                        break ;
                    }
                }
                if ( inrefinedFeature )
                {
                    continue ;
                }

                // 			i->print() ;
                Point a = *i->first*.5+*i->second*.5 ;
                Point b = *i->first*.5+*i->third*.5 ;
                Point c = ( a+b+*i->second+*i->third ) *.25 ;
                double d0 = dist ( i->first, i->second ) ;
                double d1 = dist ( i->first, i->third ) ;
                double d2 = dist ( i->third, i->second ) ;
                if ( d0 < d1 && d0 < d2 )
                {
                    illegalPoints.push_back ( a );
                    a = *i->second*.5 + *i->third*.5 ;
                    c = ( a+b+*i->second+*i->first ) *.25 ;
                }
                if ( d1 < d0 && d1 < d2 )
                {
                    illegalPoints.push_back ( b );
                    b = *i->third*.5 + *i->second*.5 ;
                    c = ( a+b+*i->third+*i->first ) *.25 ;
                }

                bool uniquea = true ;
                bool uniqueb = true ;
                bool uniquec = true ;

                for ( size_t j = 0 ; j < pointsToAdd.size() ; j++ )
                {
                    if ( uniquea && dist ( a, pointsToAdd[j] ) < POINT_TOLERANCE )
                    {
                        uniquea = false ;
                    }
                    if ( uniqueb && dist ( b, pointsToAdd[j] ) < POINT_TOLERANCE )
                    {
                        uniqueb = false ;
                    }
                    if ( uniquec && dist ( c, pointsToAdd[j] ) < POINT_TOLERANCE )
                    {
                        uniquec = false ;
                    }
                    if ( !uniquea && !uniqueb && !uniquec )
                    {
                        break ;
                    }
                }

                for ( size_t j = 0 ; j < illegalPoints.size() ; j++ )
                {
                    if ( uniquea && dist ( a, illegalPoints[j] ) < POINT_TOLERANCE )
                    {
                        uniquea = false ;
                    }
                    if ( uniqueb && dist ( b, illegalPoints[j] ) < POINT_TOLERANCE )
                    {
                        uniqueb = false ;
                    }
                    if ( uniquec && dist ( c, illegalPoints[j] ) < POINT_TOLERANCE )
                    {
                        uniquec = false ;
                    }
                    if ( !uniquea && !uniqueb && !uniquec )
                    {
                        break ;
                    }
                }

                if ( uniquea )
                {
                    pointsToAdd.push_back ( a );
                }
                if ( uniqueb )
                {
                    pointsToAdd.push_back ( b );
                }
            }

            bool inrefinedFeature= false ;
            for ( size_t j = 0 ; j < refinedFeatures.size() ; j++ )
            {
                if ( refinedFeatures[j]->in ( i->getCenter() ) )
                {
                    Point proj = i->getCenter() ;
                    tree[0]->project ( &proj ) ;
                    if ( dist ( proj,i->getCenter() ) > 2.*i->getRadius() )
                    {
                        inrefinedFeature = true ;
                        break ;
                    }
                }
            }
            if ( inrefinedFeature )
            {
                continue ;
            }

// 			i->print() ;
            Point a = *i->first*.5+*i->second*.5 ;
            Point b = *i->first*.5+*i->third*.5 ;
            Point c = *i->third*.5+*i->second*.5 ;

            bool uniquea = true ;
            bool uniqueb = true ;
            bool uniquec = true ;

            for ( size_t j = 0 ; j < pointsToAdd.size() ; j++ )
            {
                if ( uniquea && dist ( a, pointsToAdd[j] ) < POINT_TOLERANCE )
                {
                    uniquea = false ;
                }
                if ( uniqueb && dist ( b, pointsToAdd[j] ) < POINT_TOLERANCE )
                {
                    uniqueb = false ;
                }
                if ( uniquec && dist ( c, pointsToAdd[j] ) < POINT_TOLERANCE )
                {
                    uniquec = false ;
                }
                if ( !uniquea && !uniqueb && !uniquec )
                {
                    break ;
                }
            }

            for ( size_t j = 0 ; j < illegalPoints.size() ; j++ )
            {
                if ( uniquea && dist ( a, illegalPoints[j] ) < POINT_TOLERANCE )
                {
                    uniquea = false ;
                }
                if ( uniqueb && dist ( b, illegalPoints[j] ) < POINT_TOLERANCE )
                {
                    uniqueb = false ;
                }
                if ( uniquec && dist ( c, illegalPoints[j] ) < POINT_TOLERANCE )
                {
                    uniquec = false ;
                }
                if ( !uniquea && !uniqueb && !uniquec )
                {
                    break ;
                }
            }

            if ( uniquea )
            {
                pointsToAdd.push_back ( a );
            }
            if ( uniqueb )
            {
                pointsToAdd.push_back ( b );
            }
            if ( uniquec )
            {
                pointsToAdd.push_back ( c );
            }
        }

        for ( size_t i = 0 ; i < pointsToAdd.size() ; i++ )
        {
            additionalPoints.push_back ( new Point ( pointsToAdd[i] ) );
            for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            {
                j->second->insert ( additionalPoints.back() ) ;
            }
        }
        std::cerr <<  " ...done. " << std::endl ;

        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            j->second->clearCaches();
    }
}

void FeatureTree::sample()
{
    if(structuredMesh)
        return ;

    if ( samplingNumber != previousSamplingNumber )
    {
        meshPoints.clear();
        previousSamplingNumber = samplingNumber ;

        if ( is2D() )
        {
            std::cerr << "2D features " << tree.size() << std::endl ;
            double total_area = tree[0]->area() ;

            double correctionfactor = 1. ;
            if ( samplingFactors.find ( tree[0] ) != samplingFactors.end() )
            {
                correctionfactor =  samplingFactors[tree[0]] ;
            }
            tree[0]->sample ( correctionfactor * samplingNumber * 4 ) ;
            int count = 0 ;

//             #pragma omp parallel for reduction(+:count) schedule(auto)
            for ( size_t i  = 1 ; i < tree.size() ; i++ )
            {
                double shape_factor = ( sqrt ( tree[0]->area() ) / ( 2.*M_PI * tree[0]->getRadius() ) ) / ( sqrt ( tree[i]->area() ) / ( 2.*M_PI * tree[i]->getRadius() ) );
                if(tree[i]->getGeometryType() == POLYGON)
                    shape_factor = 1. ;

                if ( !tree[i]->isVirtualFeature )
                {
                    tree[i]->isUpdated = false ;
                }


                if ( shape_factor < POINT_TOLERANCE )
                {
                    continue ;
                }

                size_t npoints = ( size_t ) round ( sqrt ( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ) ;
                if ( npoints < 8 && npoints >= 5 )
                {
                    npoints = 8 ;
                }
// 				size_t npoints = std::max( ( size_t )round( sqrt( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ), (size_t) 8 ) ;
                correctionfactor = 1. ;
                if ( samplingFactors.find ( tree[i] ) != samplingFactors.end() )
                {
                    correctionfactor = samplingFactors[tree[i]] ;
                    npoints = ( size_t ) round ( correctionfactor*npoints ) ;
                }

                if ( !tree[i]->isVirtualFeature )
                {
                    for ( size_t j = 0 ;  j < refinementZones.size() ; j++ )
                    {
                        if ( tree[i]->intersects ( refinementZones[j] ) || refinementZones[j]->in ( tree[i]->getCenter() ) )
                        {
                            npoints *= 2 ;

                            refinedFeatures.push_back ( tree[i] );

                        }
                    }
                }

                if( npoints >= samplingRestriction && !tree[i]->isVirtualFeature)
                {
                    count++ ;
                    tree[i]->sample ( npoints ) ;
                }

                if ( !tree[i]->isVirtualFeature )
                {
                    tree[i]->isUpdated = false ;
                }
            }

            std::cerr << count << " particles meshed" << std::endl ;
        }
        else if ( is3D() )
        {
            std::cerr << "\r 3D features... sampling feature 0/" << this->tree.size() << "          " << std::flush ;
            if ( samplingFactors.find ( tree[0] ) != samplingFactors.end() )
                tree[0]->sample ( samplingFactors[tree[0]]*samplingNumber ) ;
            else
                tree[0]->sample ( samplingNumber ) ;
            double total_area = tree[0]->area()  ;
            int count = 0 ;
            double base_shape_factor = tree[0]->area()/pow(tree[0]->volume(), .666666666) ;

            #pragma omp parallel for reduction(+:count) schedule(runtime)
            for ( int i  = 1 ; i < ( int ) tree.size() ; i++ )
            {
                std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << "          " << std::flush ;

                double shape_factor = tree[i]->area() /pow(tree[i]->volume() , .666666666);
                size_t npoints = ( size_t ) round ( ( 1.5 * samplingNumber * tree[i]->area() * shape_factor/base_shape_factor ) / ( total_area ) ) ;
                if ( samplingFactors.find ( tree[i] ) != samplingFactors.end() )
                {
                    npoints = ( size_t ) round ( samplingFactors[tree[i]]*npoints ) ;
                }

                if( npoints >= samplingRestriction && !tree[i]->isVirtualFeature)
                {
                    count++ ;
                    tree[i]->sample ( npoints ) ;
                }


                if ( !tree[i]->isVirtualFeature )
                {
                    tree[i]->isUpdated = false ;
                }

            }

            std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << " ...done" << std::endl ;
        }
    }
    else
    {

        for ( size_t i = 0 ; i < additionalPoints.size() ; i++ )
        {
            delete additionalPoints[i] ;
        }

        additionalPoints.clear();

        meshPoints.clear();
        previousSamplingNumber = samplingNumber ;

        if ( is2D() )
        {
            std::cerr << "2D features (updating sampling)" << std::endl ;
            double total_area = tree[0]->area() ;

            if ( tree[0]->isUpdated )
            {
                tree[0]->sample ( samplingNumber ) ;
            }

            #pragma omp parallel for schedule(runtime)
            for ( size_t i  = 1 ; i < tree.size() ; i++ )
            {
                if ( tree[i]->isUpdated )
                {
                    #pragma omp critical
                    {
                        grid->remove ( tree[i] ) ;
                        grid->forceAdd ( tree[i] );
                    }
                    double shape_factor = ( sqrt ( tree[0]->area() ) / ( 2.*M_PI * tree[0]->getRadius() ) ) / ( sqrt ( tree[i]->area() ) / ( 2.*M_PI * tree[i]->getRadius() ) );

                    if ( shape_factor < POINT_TOLERANCE )
                    {
                        continue ;
                    }

                    size_t npoints = std::max ( ( size_t ) round ( sqrt ( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ), ( size_t ) 8 ) ;
                    if ( samplingFactors.find ( tree[i] ) != samplingFactors.end() )
                    {
                        npoints = ( size_t ) round ( samplingFactors[tree[i]]*npoints ) ;
                    }
                    if ( npoints >= samplingRestriction && !tree[i]->isVirtualFeature && npoints < samplingNumber )
                    {
                        tree[i]->sample ( npoints ) ;
                    }
                }
            }
        }
        else if ( is3D() )
        {
            std::cerr << "\r 3D features... sampling feature 0/" << this->tree.size() << "          " << std::flush ;

            if ( tree[0]->isUpdated )
            {
                tree[0]->sample ( 2.5 * samplingNumber ) ;
            }

            double total_area = tree[0]->area()  ;
            int count = 0 ;
            double base_shape_factor = tree[0]->area()/pow(tree[0]->volume() , .666666666);

            #pragma omp parallel for schedule(runtime)
            for ( int i  = 1 ; i < ( int ) tree.size() ; i++ )
            {
                if ( tree[i]->isUpdated )
                {
                    #pragma omp critical
                    {
                        grid3d->remove ( tree[i] ) ;
                        grid3d->forceAdd ( tree[i] );
                    }
                    std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << "          " << std::flush ;

                    double shape_factor = tree[i]->area() / pow(tree[i]->volume(), .666666666) ;
                    size_t npoints = ( size_t ) round ( ( 1.5 * samplingNumber * tree[i]->area() * shape_factor/base_shape_factor ) / ( total_area ) ) ;
                    if ( samplingFactors.find ( tree[i] ) != samplingFactors.end() )
                    {
                        npoints = ( size_t ) round ( samplingFactors[tree[i]]*npoints ) ;
                    }
                    if ( npoints >= samplingRestriction && !tree[i]->isVirtualFeature )
                    {
                        tree[i]->sample ( npoints ) ;
                    }

                    count++ ;
                }
            }

            std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << " ...done" << std::endl ;
        }
    }

    for(size_t i = 0 ; i < tree.size() ; i++)
        tree[i]->makeMaskBoundary() ;


}


void FeatureTree::refine ( size_t nit, SamplingCriterion *cri )
{
    if(structuredMesh)
        return ;
    for ( size_t t = 0 ; t < nit ; t++ )
    {
        bool corrected = false ;

        int count = 0 ;

        for (auto j = dtree->begin() ; j != dtree->end()  ; j++ )
        {
            if ( !cri->meetsCriterion ( (DelaunayTriangle *)j ) )
            {
                count++ ;
            }
        }

        std::cerr << count << " non-conformant triangles " << std::endl ;

        for ( auto j = dtree->begin() ; j != dtree->end()  ; j++ )
        {
            if ( !cri->meetsCriterion (j ) )
            {
                std::vector<Point> temp = cri->suggest ( j ) ;

                if ( !temp.empty() )
                {
                    std::random_shuffle ( temp.begin(), temp.end() ) ;
                    std::cerr << "inserting " << temp.size() << " points" << std::endl ;

                    for ( size_t i = 0 ; i < temp.size() ; i++ )
                    {
                        dtree->insert ( new Point ( temp[i] ) ) ;
                        corrected = true ;
                    }

                    break ;
                }
            }
        }

        if ( !corrected )
        {
            break ;
        }
    }
}

void FeatureTree::refine ( size_t level )
{
    if(structuredMesh)
        return ;
    state.setStateTo ( RENUMBERED, false ) ;
    double pointDensity = 0 ;

    if ( is2D() )
    {
        pointDensity = .5 * sqrt ( tree[0]->area() ) / meshPoints.size() ;
    }
    else
    {
        pointDensity = .5 * pow ( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;
    }

    if ( level < 1 )
    {
        return ;
    }

    if ( this->dtree == nullptr && this->dtree3D != nullptr )
    {
        std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;

        for ( size_t j = 1;  j < this->tree.size() ; j++ )
        {
            zonesVec.push_back ( std::pair<std::vector<Geometry *>, Feature *> ( this->tree[j]->getRefinementZones ( level ),  this->tree[j] ) ) ;
        }

        std::vector<Feature *> enrichmentFeature ;

        for ( size_t i  = 0 ; i < this->tree.size() ; i++ )
        {
            if ( tree[i]->isEnrichmentFeature )
            {
                enrichmentFeature.push_back ( tree[i] ) ;
            }
        }

        size_t points_added = 0 ;

        for ( size_t i = 0 ; i < zonesVec.size() ; i++ )
        {
            for ( size_t j = 0 ; j < zonesVec[i].first.size() ; j++ )
            {
                if ( !zonesVec[i].second->isEnrichmentFeature )
                {


                    std::vector<Point> toAdd ;
                    std::vector<DelaunayTetrahedron *>  tet = this->dtree3D->getConflictingElements ( zonesVec[i].first[j] ) ;
                    // 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;

                    for ( size_t k = 0 ; k < tet.size() ; k++ )
                    {
                        size_t count_0 = 0 ;
                        size_t count_1 = 0 ;
                        size_t count_2 = 0 ;
                        size_t count_3 = 0 ;
                        size_t count_4 = 0 ;
                        size_t count_5 = 0 ;

                        Point p0  = *tet[k]->first * ( 0.5 ) + *tet[k]->second * ( 0.5 ) ;
                        p0.setId ( -1 ) ;
                        Point p1  = *tet[k]->first * ( 0.5 ) + *tet[k]->third * ( 0.5 ) ;
                        p1.setId ( -1 ) ;
                        Point p2  = *tet[k]->first * ( 0.5 ) + *tet[k]->fourth * ( 0.5 ) ;
                        p2.setId ( -1 ) ;
                        Point p3  = *tet[k]->third * ( 0.5 ) + *tet[k]->second * ( 0.5 ) ;
                        p3.setId ( -1 ) ;
                        Point p4  = *tet[k]->fourth * ( 0.5 ) + *tet[k]->third * ( 0.5 ) ;
                        p4.setId ( -1 ) ;
                        Point p5  = *tet[k]->second * ( 0.5 ) + *tet[k]->fourth * ( 0.5 ) ;
                        p5.setId ( -1 ) ;


                        for ( size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++ )
                        {
                            if ( !zonesVec[i].second->getFather()->getChild ( l )->isEnrichmentFeature )
                            {
                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p0, pointDensity ) )
                                {
                                    count_0++ ;
                                }

                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p1, pointDensity ) )
                                {
                                    count_1++ ;
                                }

                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p1, pointDensity ) )
                                {
                                    count_2++ ;
                                }

                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p2, pointDensity ) )
                                {
                                    count_2++ ;
                                }

                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p3, pointDensity ) )
                                {
                                    count_3++ ;
                                }

                                if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p4, pointDensity ) )
                                {
                                    count_4++ ;
                                }
                            }
                        }


                        // 					for(size_t m = 0 ; m < enrichmentFeature.size() ; m++)
                        // 					{
                        // 						if(enrichmentFeature[m]->inBoundary(p0))
                        // 							count_0++ ;
                        // 						if(enrichmentFeature[m]->inBoundary(p1))
                        // 							count_1++ ;
                        // 						if(enrichmentFeature[m]->inBoundary(p2))
                        // 							count_2++ ;
                        // 						if(enrichmentFeature[m]->inBoundary(p3))
                        // 							count_3++ ;
                        // 						if(enrichmentFeature[m]->inBoundary(p4))
                        // 							count_4++ ;
                        // 						if(enrichmentFeature[m]->inBoundary(p5))
                        // 							count_5++ ;
                        // 					}

                        if ( count_0 == 0 )
                        {
                            toAdd.push_back ( p0 ) ;
                        }

                        if ( count_1 == 0 )
                        {
                            toAdd.push_back ( p1 ) ;
                        }

                        if ( count_2 == 0 )
                        {
                            toAdd.push_back ( p2 ) ;
                        }

                        if ( count_3 == 0 )
                        {
                            toAdd.push_back ( p3 ) ;
                        }

                        if ( count_4 == 0 )
                        {
                            toAdd.push_back ( p4 ) ;
                        }

                        if ( count_5 == 0 )
                        {
                            toAdd.push_back ( p5 ) ;
                        }

                    }


                    std::sort ( toAdd.begin(), toAdd.end() ) ;
                    auto e = std::unique ( toAdd.begin(), toAdd.end() ) ;
                    toAdd.erase ( e, toAdd.end() ) ;

                    // 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;

                    std::random_shuffle ( toAdd.begin(), toAdd.end() ) ;

                    for ( size_t k = 0 ; k < toAdd.size() ; k++ )
                    {
                        std::cerr << "\r refining feature " << i + 1 << "/" << zonesVec.size() << "...added " << points_added << " points" << std::flush ;
                        Point *p = new Point ( toAdd[k] ) ;
                        points_added++ ;
                        meshPoints.push_back ( std::pair<Point *, const Feature *> ( p, zonesVec[i].second->getFather() ) ) ;
                        dtree3D->insert ( p ) ;
                    }
                }
            }
        }

        std::cerr << "...done" << std::endl ;
    }
    else
    {
        std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;

        for ( size_t j = 1;  j < this->tree.size() ; j++ )
        {
            zonesVec.push_back ( std::pair<std::vector<Geometry *>, Feature *> ( this->tree[j]->getRefinementZones ( level ),  this->tree[j] ) ) ;
        }

        std::vector<Feature *> enrichmentFeature ;

        for ( size_t i  = 0 ; i < this->tree.size() ; i++ )
        {
            if ( tree[i]->isEnrichmentFeature )
            {
                enrichmentFeature.push_back ( tree[i] ) ;
            }
        }

        size_t points_added = 0 ;

        for ( size_t i = 0 ; i < zonesVec.size() ; i++ )
        {

            for ( size_t j = 0 ; j < zonesVec[i].first.size() ; j++ )
            {

                if ( !zonesVec[i].second->isEnrichmentFeature )
                {
                    std::vector<Point *> sample = zonesVec[i].second->doubleSurfaceSampling() ;

                    for ( size_t k = 0 ; k < sample.size() ; k++ )
                    {
                        if ( tree[0]->in ( *sample[k] ) )
                        {
                            bool yes = true ;

                            for ( size_t l = 0 ; l < enrichmentFeature.size() ; l++ )
                            {
                                if ( enrichmentFeature[l]->inBoundary ( *sample[k], pointDensity ) )
                                {
                                    yes = false ;
                                    break ;
                                }
                            }

                            for ( size_t l = 0 ; l < zonesVec[i].second->getChildren().size() ; l++ )
                            {
                                if ( zonesVec[i].second->getChild ( l )->inBoundary ( *sample[k], pointDensity ) )
                                {
                                    yes = false ;
                                    break ;
                                }
                            }

                            if ( yes )
                            {
                                dtree->insert ( sample[k] ) ;
                                this->meshPoints.push_back ( std::make_pair ( sample[k], zonesVec[i].second ) ) ;
                            }
                        }
                    }
                }

                std::vector<Point> toAdd ;

                std::vector<DelaunayTriangle *>  tri = this->dtree->getConflictingElements ( zonesVec[i].first[j] ) ;
// 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;

                for ( size_t k = 0 ; k < tri.size() ; k++ )
                {
// 				if((*tri)[k]->area() > 4e-4)
// 				{
                    size_t count_0 = 0 ;
                    size_t count_1 = 0 ;
                    size_t count_2 = 0 ;
                    double rand0 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
                    double rand1 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
                    double rand2 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/

                    Point p0  = *tri[k]->first * ( 0.5 + rand0 ) + *tri[k]->second * ( 0.5 - rand0 ) ;
                    p0.setId ( -1 ) ;
                    Point p1  = *tri[k]->first * ( 0.5 + rand1 ) + *tri[k]->third * ( 0.5 - rand1 ) ;
                    p1.setId ( -1 ) ;
                    Point p2  = *tri[k]->second * ( 0.5 + rand2 ) + *tri[k]->third * ( 0.5 - rand2 ) ;
                    p2.setId ( -1 ) ;
//

                    for ( size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++ )
                    {
                        if ( !zonesVec[i].second->getFather()->getChild ( l )->isEnrichmentFeature )
                        {
                            if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p0, pointDensity ) )
                            {
                                count_0++ ;
                            }

                            if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p1, pointDensity ) )
                            {
                                count_1++ ;
                            }

                            if ( zonesVec[i].second->getFather()->getChild ( l )->inBoundary ( p1, pointDensity ) )
                            {
                                count_2++ ;
                            }
                        }
                    }



                    for ( size_t m = 0 ; m < enrichmentFeature.size() ; m++ )
                    {
                        if ( enrichmentFeature[m]->inBoundary ( p0, pointDensity ) )
                        {
                            count_0++ ;
                        }

                        if ( enrichmentFeature[m]->inBoundary ( p1, pointDensity ) )
                        {
                            count_1++ ;
                        }

                        if ( enrichmentFeature[m]->inBoundary ( p2, pointDensity ) )
                        {
                            count_2++ ;
                        }
                    }

                    if ( count_0 == 0 && zonesVec[i].first[j]->in ( *tri[k]->first ) && zonesVec[i].first[j]->in ( *tri[k]->second ) )
                    {
                        toAdd.push_back ( p0 ) ;
                    }

                    if ( count_1 == 0 && zonesVec[i].first[j]->in ( *tri[k]->first ) && zonesVec[i].first[j]->in ( *tri[k]->third ) )
                    {
                        toAdd.push_back ( p1 ) ;
                    }

                    if ( count_2 == 0 && zonesVec[i].first[j]->in ( *tri[k]->second ) && zonesVec[i].first[j]->in ( *tri[k]->third ) )
                    {
                        toAdd.push_back ( p2 ) ;
                    }

// 				}

                }


                std::sort ( toAdd.begin(), toAdd.end() ) ;
                auto e = std::unique ( toAdd.begin(), toAdd.end() );
                toAdd.erase ( e, toAdd.end() ) ;

// 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;
                std::random_shuffle ( toAdd.begin(), toAdd.end() ) ;

                for ( size_t k = 0 ; k < toAdd.size() ; k++ )
                {
                    std::cerr << "\r refining feature " << i + 1 << "/" << zonesVec.size() << "...added " << points_added << " points" << std::flush ;
                    Point *p = new Point ( toAdd[k] ) ;
                    points_added++ ;
                    meshPoints.push_back ( std::pair<Point *, const Feature *> ( p, zonesVec[i].second->getFather() ) ) ;
                    dtree->insert ( p ) ;
                }
            }
        }

        std::cerr << "...done" << std::endl ;
    }
}

Form * FeatureTree::getElementBehaviour ( Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator & t, int layer,  bool onlyUpdate ) const
{
    if(structuredMesh)
        return nullptr;

    int root_box = 0 ;

    if ( !inRoot ( t->getCenter() ) )
    {
        return new VoidForm();
    }

    if ( t->getBoundingPoints().size() % 3 != 0 )
    {
        return new VoidForm() ;
    }

    std::vector<Feature *> targets ;

    if ( tree.size() > 32 )
    {
        std::vector<const Geometry *> targetstmp = grid->coOccur ( t->getPrimitive() ) ;

        for ( size_t i = 0 ; i < targetstmp.size() ; i++ )
        {
            const Feature * tmp = dynamic_cast<const Feature *> ( targetstmp[i] ) ;
            if ( tmp->getLayer() == layer && ( tmp->getBoundingPoints().size() || tmp->isVirtualFeature || samplingRestriction == 0 ) )
            {
                targets.push_back ( const_cast<Feature *> ( tmp ) ) ;
            }
        }
    }
    else
    {
        for ( size_t i = 0 ; i < tree.size() ; i++ )
        {
            if ( tree[i]->getLayer() == layer && ( tree[i]->getBoundingPoints().size() || tree[i]->isVirtualFeature || samplingRestriction == 0 ) )
            {
                targets.push_back ( tree[i] ) ;
            }
        }
    }

    Form * found = nullptr ;


    if ( !targets.empty() )
    {

        for ( int i = targets.size() - 1 ; i >= 0  ; i-- )
        {
            if ( !targets[i]->isEnrichmentFeature && targets[i]->in( t->getCenter() ) && targets[i]->inMask( t->getCenter(), 0. ) && ( !onlyUpdate || (onlyUpdate && targets[i]->isUpdated) ) )
            {
                bool notInChildren  = true ;

                std::vector<Feature *> descendants = targets[i]->getDescendants() ;

                for ( size_t j = 0 ; j < descendants.size() ; j++ )
                {
                    if ( descendants[j]->getLayer() == layer && !descendants[j]->isEnrichmentFeature && descendants[j]->in ( t->getCenter() ) && descendants[j]->inMask ( t->getCenter(), 0. ) )
                    {
                        notInChildren = false ;
                        break ;
                    }
                }

                for(size_t k = 0 ; k < targets.size() ; k++)
                {
                    if( ((int) k != i) && targets[k]->getLayer() == layer && !targets[k]->isEnrichmentFeature && targets[k]->isMaskedBy( targets[i] ) && targets[k]->in( t->getCenter() ) && targets[k]->inMask( t->getCenter(), 0. )  )
                    {
                        notInChildren = false ;
                        break ;
                    }
                }

                if ( notInChildren )
                {
                    if ( targets[i]->getBehaviour ( t->getCenter() )->timeDependent() )
                    {
                        if ( !targets[i]->getBehaviour ( t->getCenter() )->spaceDependent() )
                        {
                            Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                            if ( targets[i]->getBehaviourSource() )
                            {
                                b->setSource ( targets[i]->getBehaviourSource() );
                            }
                            else
                            {
                                b->setSource (targets[i]);
                            }
                            found = b ;
                        }
                        else
                        {
                            Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                            b->transform ( t ) ;
                            if ( targets[i]->getBehaviourSource() )
                            {
                                b->setSource ( targets[i]->getBehaviourSource() );
                            }
                            else
                            {
                                b->setSource (  targets[i] );
                            }
                            found = b ;
                        }
                    }
                    else if ( !targets[i]->getBehaviour ( t->getCenter() )->spaceDependent() )
                    {
                        Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                        if ( targets[i]->getBehaviourSource() )
                        {
                            b->setSource ( targets[i]->getBehaviourSource() );
                        }
                        else
                        {
                            b->setSource ( targets[i] );
                        }
                        found = b ;
                    }
                    else
                    {
                        Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                        b->transform ( t ) ;
                        if ( targets[i]->getBehaviourSource() )
                        {
                            b->setSource ( targets[i]->getBehaviourSource() );
                        }
                        else
                        {
                            b->setSource ( targets[i] );
                        }
                        found = b ;
                    }
                }
            }
        }
    }

    if ( found )
    {
        return found ;
    }

    if(!tree[root_box]->getBehaviour())
        return new VoidForm() ;
    if ( !onlyUpdate && tree[root_box]->getLayer() == layer )
    {
        if ( tree[root_box]->getBehaviour ( t->getCenter() )->timeDependent() )
        {
            if ( !tree[root_box]->getBehaviour ( t->getCenter() )->spaceDependent() )
            {
                Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
                if ( tree[root_box]->getBehaviourSource() )
                {
                    b->setSource ( tree[root_box]->getBehaviourSource() );
                }
                else
                {
                    b->setSource ( tree[root_box] );
                }
                return b ;
            }
            else
            {
                Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
                if ( tree[root_box]->getBehaviourSource() )
                {
                    b->setSource ( tree[root_box]->getBehaviourSource() );
                }
                else
                {
                    b->setSource ( tree[root_box] );
                }
                b->transform ( t ) ;
                return b ;
            }
        }
        else if ( !tree[root_box]->getBehaviour ( t->getCenter() )->spaceDependent() )
        {
            Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
            if ( tree[root_box]->getBehaviourSource() )
            {
                b->setSource ( tree[root_box]->getBehaviourSource() );
            }
            else
            {
                b->setSource ( tree[root_box] );
            }
            return b ;
        }
        else
        {
            Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
            if ( tree[root_box]->getBehaviourSource() )
            {
                b->setSource ( tree[root_box]->getBehaviourSource() );
            }
            else
            {
                b->setSource ( tree[root_box] );
            }
            b->transform ( t ) ;
            return b ;
        }
        Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
        if ( tree[root_box]->getBehaviourSource() )
        {
            b->setSource ( tree[root_box]->getBehaviourSource() );
        }
        else
        {
            b->setSource ( tree[root_box] );
        }
        return b ;
    }
    else if ( !onlyUpdate && tree[root_box]->getLayer() == layer )
    {
        return new VoidForm() ;
    }

    return new VoidForm() ;

}

Form * FeatureTree::getElementBehaviour ( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>::iterator & t, int layer,  bool onlyUpdate ) const
{
    if(structuredMesh)
        return nullptr;

    int root_box = 0 ;

    if ( !inRoot ( t->getCenter() ) )
    {
        return new VoidForm() ;
    }

// 	std::vector<Geometry *> targetstmp = grid3d->coOccur(t->getPrimitive()) ;
    std::vector<Feature *> targets ;

    if ( tree.size() > 32 )
    {
        std::vector<const Geometry *> targetstmp = grid3d->coOccur ( t->getPrimitive() ) ;

        for ( size_t i = 0 ; i < targetstmp.size() ; i++ )
        {
            const Feature * tmp = dynamic_cast<const Feature *> ( targetstmp[i] ) ;
            if ( tmp->getLayer() == layer )
            {
                targets.push_back ( const_cast<Feature *> ( tmp ) ) ;
            }
        }
    }
    else
    {
        std::vector<Feature *> targetstmp = tree ;
        for ( size_t i = 0 ; i < targetstmp.size() ; i++ )
        {
            if ( targetstmp[i]->getLayer() == layer )
            {
                targets.push_back ( targetstmp[i] ) ;
            }
        }
    }

    Form * found = nullptr ;


    if ( !targets.empty() )
    {

        for ( int i = targets.size() - 1 ; i >= 0  ; i-- )
        {

            if ( !targets[i]->isEnrichmentFeature && targets[i]->in ( t->getCenter() ) && ( !onlyUpdate || (onlyUpdate && targets[i]->isUpdated) ) )
            {
                bool notInChildren  = true ;

                std::vector<Feature *> descendants = targets[i]->getDescendants() ;

                for ( size_t j = 0 ; j < descendants.size() ; j++ )
                {
                    if ( descendants[j]->getLayer() == layer && !descendants[j]->isEnrichmentFeature && descendants[j]->in ( t->getCenter() ) )
                    {
                        notInChildren = false ;
                        break ;
                    }
                }

                if ( notInChildren )
                {
                    if ( targets[i]->getBehaviour ( t->getCenter() )->timeDependent() )
                    {
                        if ( !targets[i]->getBehaviour ( t->getCenter() )->spaceDependent() )
                        {
                            Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                            if ( targets[i]->getBehaviourSource() )
                            {
                                b->setSource ( targets[i]->getBehaviourSource() );
                            }
                            else
                            {
                                b->setSource ( targets[i] );
                            }
                            found = b ;
                        }
                        else
                        {
                            Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                            b->transform ( t ) ;
                            if ( targets[i]->getBehaviourSource() )
                            {
                                b->setSource ( targets[i]->getBehaviourSource() );
                            }
                            else
                            {
                                b->setSource ( targets[i] );
                            }
                            found = b ;
                        }
                    }
                    else if ( !targets[i]->getBehaviour ( t->getCenter() )->spaceDependent() )
                    {
                        Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                        if ( targets[i]->getBehaviourSource() )
                        {
                            b->setSource ( targets[i]->getBehaviourSource() );
                        }
                        else
                        {
                            b->setSource ( targets[i] );
                        }
                        found = b ;
                    }
                    else
                    {
                        Form *b = targets[i]->getBehaviour ( t->getCenter() )->getCopy() ;
                        b->transform ( t ) ;
                        if ( targets[i]->getBehaviourSource() )
                        {
                            b->setSource ( targets[i]->getBehaviourSource() );
                        }
                        else
                        {
                            b->setSource ( targets[i] );
                        }
                        found = b ;
                    }
                }
            }
        }
    }

    if ( found )
    {
        return found ;
    }

    if ( !onlyUpdate && tree[root_box]->getLayer() == layer )
    {
        if ( tree[root_box]->getBehaviour ( t->getCenter() )->timeDependent() )
        {
            if ( !tree[root_box]->getBehaviour ( t->getCenter() )->spaceDependent() )
            {
                Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
                if ( tree[root_box]->getBehaviourSource() )
                {
                    b->setSource ( tree[root_box]->getBehaviourSource() );
                }
                else
                {
                    b->setSource ( tree[root_box] );
                }
                return b ;
            }
            else
            {
                Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
                if ( tree[root_box]->getBehaviourSource() )
                {
                    b->setSource ( tree[root_box]->getBehaviourSource() );
                }
                else
                {
                    b->setSource ( tree[root_box] );
                }
                b->transform ( t ) ;
                return b ;
            }
        }
        else if ( !tree[root_box]->getBehaviour ( t->getCenter() )->spaceDependent() )
        {
            Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
            if ( tree[root_box]->getBehaviourSource() )
            {
                b->setSource ( tree[root_box]->getBehaviourSource() );
            }
            else
            {
                b->setSource ( tree[root_box] );
            }
            return b ;
        }
        else
        {
            Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
            if ( tree[root_box]->getBehaviourSource() )
            {
                b->setSource ( tree[root_box]->getBehaviourSource() );
            }
            else
            {
                b->setSource ( tree[root_box] );
            }
            b->transform ( t ) ;
            return b ;
        }
        Form *b = tree[root_box]->getBehaviour ( t->getCenter() )->getCopy() ;
        if ( tree[root_box]->getBehaviourSource() )
        {
            b->setSource ( tree[root_box]->getBehaviourSource() );
        }
        else
        {
            b->setSource ( tree[root_box] );
        }
        return b ;
    }
    else if ( !onlyUpdate && tree[root_box]->getLayer() == layer )
    {
        return new VoidForm() ;
    }

    return new VoidForm() ;

}

Point *FeatureTree::checkElement ( const DelaunayTetrahedron *t ) const
{
    if(structuredMesh)
        return nullptr;

    double pointDensity = 0 ;

    if ( is2D() )
    {
        pointDensity = .6 * sqrt ( tree[0]->area() ) / meshPoints.size() ;
    }
    else
    {
        pointDensity = .6 * pow ( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;
    }

    if ( !inRoot ( t->getCenter() ) )
    {
        return nullptr;
    }

    for ( int i = tree.size() - 1 ; i >= 0 ; i-- )
    {
        int inCount = tree[i]->in ( *t->first )
                      + tree[i]->in ( *t->second )
                      + tree[i]->in ( *t->third )
                      + tree[i]->in ( *t->fourth )
                      + tree[i]->in ( t->getCenter() );

        if ( inCount > 2 && inRoot ( t->getCenter() ) )
        {
            bool inChild = false ;

            std::vector<Feature *> tocheck = tree[i]->getDescendants();

            for ( size_t j = 0 ; j < tocheck.size() ; j++ )
            {
                if ( tocheck[j]->in ( t->getCenter() ) )
                {
                    inChild = true ;
                    break ;
                }
            }


            if ( !inChild )
            {

                size_t count_in = 0 ;

                count_in += tree[i]->inBoundary ( *t->first, pointDensity ) ;
                count_in += tree[i]->inBoundary ( *t->second, pointDensity ) ;
                count_in += tree[i]->inBoundary ( *t->third, pointDensity ) ;
                count_in += tree[i]->inBoundary ( *t->fourth, pointDensity ) ;

                if ( count_in == 4 && tree[i]->in ( t->getCenter() ) )
                {
                    return nullptr;
                }

                else
                {
                    Point p1 ( t->getCenter() ) ;
                    Point p0 ( t->getCircumCenter() ) ;
                    tree[i]->project ( &p1 );
                    tree[i]->project ( &p0 );
                    double d0 = std::min ( dist ( &p0, t->first ), std::min ( dist ( &p0, t->second ), dist ( &p0, t->third ) ) ) ;
                    double d1 = std::min ( dist ( &p1, t->first ), std::min ( dist ( &p1, t->second ), dist ( &p1, t->third ) ) ) ;

                    if ( t->inCircumSphere ( p0 )
                            && d0 > 1e-8
                            && d0 > d1
                       )
                    {
                        return new Point ( p0 );
                    }
                    else if ( t->inCircumSphere ( p1 )
                              && d1 > 1e-8
                              && d1 > d0
                            )
                    {
                        return new Point ( p1 );
                    }
                    else
                    {
                        return nullptr ;
                    }
                }

            }
        }
    }

    return nullptr;
}

Point *FeatureTree::checkElement ( const DelaunayTriangle *t ) const
{
    if(structuredMesh)
        return nullptr;

    double pointDensity = 0 ;

    if ( is2D() )
    {
        pointDensity = .6 * sqrt ( tree[0]->area() ) / meshPoints.size() ;
    }
    else
    {
        pointDensity = .6 * pow ( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;
    }

    if ( !inRoot ( t->getCenter() ) )
    {
        return nullptr;
    }

    for ( int i = tree.size() - 1 ; i >= 0 ; i-- )
    {
        int inCount = tree[i]->in ( *t->first )
                      + tree[i]->in ( *t->second )
                      + tree[i]->in ( *t->third ) + tree[i]->in ( t->getCenter() );

        if ( inCount > 1 && inRoot ( t->getCenter() ) )
        {
            bool inChild = false ;

            std::vector<Feature *> tocheck = tree[i]->getDescendants();

            for ( size_t j = 0 ; j < tocheck.size() ; j++ )
            {
                if ( tocheck[j]->in ( t->getCenter() ) )
                {
                    inChild = true ;
                    break ;
                }
            }


            if ( !inChild )
            {

                size_t count_in = 0 ;

                count_in += tree[i]->inBoundary ( *t->first, pointDensity ) ;
                count_in += tree[i]->inBoundary ( *t->second, pointDensity ) ;
                count_in += tree[i]->inBoundary ( *t->third, pointDensity ) ;

                if ( count_in == 3 && tree[i]->in ( t->getCenter() ) )
                {
                    return nullptr;
                }

                else
                {
                    Point p1 ( t->getCenter() ) ;
                    Point p0 ( t->getCircumCenter() ) ;
                    tree[i]->project ( &p1 );
                    tree[i]->project ( &p0 );
                    double d0 = std::min ( dist ( &p0, t->first ), std::min ( dist ( &p0, t->second ), dist ( &p0, t->third ) ) ) ;
                    double d1 = std::min ( dist ( &p1, t->first ), std::min ( dist ( &p1, t->second ), dist ( &p1, t->third ) ) ) ;

                    if ( t->inCircumCircle ( p0 )
                            && d0 > 1e-4
                            && d0 > d1
                       )
                    {
                        return new Point ( p0 );
                    }
                    else if ( t->inCircumCircle ( p1 )
                              && d1 > 1e-4
                              && d1 > d0
                            )
                    {
                        return new Point ( p1 );
                    }
                    else
                    {
                        return nullptr ;
                    }
                }

            }
        }
    }

    return nullptr;
}

Feature *FeatureTree::getFeatForTetra ( const DelaunayTetrahedron *t ) const
{
    if(structuredMesh)
        return nullptr;

    if ( !inRoot ( t->getCenter() ) )
    {
        return nullptr;
    }

    for ( int i = tree.size() - 1 ; i >= 0 ; i-- )
    {
        if ( tree[i]->in ( t->getCenter() ) && ( inRoot ( t->getCenter() ) ) )
        {
            std::vector<Feature *> tocheck = tree[i]->getChildren();
            std::vector<Feature *> tocheckNew =  tocheck;

            while ( !tocheckNew.empty() )
            {
                std::vector<Feature *> tocheckTemp ;

                for ( size_t k = 0 ; k < tocheckNew.size() ; k++ )
                {
                    tocheckTemp.insert (
                        tocheckTemp.end(),
                        tocheckNew[k]->getChildren().begin(),
                        tocheckNew[k]->getChildren().end()
                    ) ;
                }

                tocheck.insert ( tocheck.end(), tocheckTemp.begin(), tocheckTemp.end() ) ;
                tocheckNew = tocheckTemp ;
            }

            for ( size_t j = 0 ; j < tocheck.size() ; j++ )
            {
                if ( tocheck[j]->in ( t->getCenter() ) )
                {
                    break ;
                }
            }

            return tree[i];
        }
    }

    return nullptr;
}

void FeatureTree::setElementBehaviours()
{

    if ( !father3D )
    {
        father3D = new TetrahedralElement ( elemOrder ) ;
    }

//     father3D->compileAndPrecalculate() ;

    if ( !father2D )
    {
        father2D = new TriElement ( elemOrder ) ;
    }

//     father2D->compileAndPrecalculate() ;

    if(structuredMesh)
    {
        std::cerr << " setting behaviours..." << std::flush ;
        int setcount = 0 ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( setcount % 1000 == 0 )
            {
                std::cerr << "\r setting behaviours : base layer : tet " << setcount << "/" << i.size() << std::flush ;
            }

            i->refresh ( father3D ) ;

            if ( !i->getBehaviour() )
            {
                i->setBehaviour ( dtree3D, new VoidForm() ) ;
            }
            setcount++ ;
        }


        std::cerr << " ...done" << std::endl ;
        return ;
    }

    if ( is2D() )
    {
        double remainder = 0 ;
        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            if ( i->first != -1 )
            {
                remainder += scalingFactors[i->first] ;
            }
        }
        scalingFactors[-1] = 1.-remainder ;
        layer2d[-1] = dtree ;



        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            int setcount = 0 ;
            std::cerr << "\r setting behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount++ << "/" << i->second->begin().size() << "    " << std::flush ;
            for ( auto j = i->second->begin() ; j != i->second->end() ; j++ )
            {
                if ( setcount++ % 1000 == 0 )
                {
                    std::cerr << "\r setting behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount << "/" << j.size() << "    " << std::flush ;
                }

                j->refresh ( father2D ) ;
                Form * bf =  getElementBehaviour ( j, i->first );
                if(bf && bf->getDamageModel() && bf->getDamageModel()->alternating)
                    setAlternating(true) ;

                j->setBehaviour ( i->second, bf ) ;


                if ( !j->getBehaviour() )
                {
                    j->setBehaviour ( i->second, new VoidForm() ) ;
                }

            }
            std::cerr << " ...done" << std::endl ;

        }


    }
    else
    {

        std::cerr << " setting behaviours..." << std::flush ;
        int setcount = 0 ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
//             if ( setcount % 1000 == 0 )
//             {
//                 std::cerr << "\r setting behaviours : base layer : tet " << setcount << "/" << i.size() << std::flush ;
//             }
            i->refresh ( father3D ) ;
            Form * bf =  getElementBehaviour ( i );
            if(bf && bf->getDamageModel() && bf->getDamageModel()->alternating)
                setAlternating(true) ;

            i->setBehaviour ( dtree3D, bf ) ;

            if ( !i->getBehaviour() )
            {
                i->setBehaviour ( dtree3D, new VoidForm() ) ;
            }
            setcount++ ;
        }


        std::cerr << " ...done" << std::endl ;


    }
    behaviourSet = true ;
}

void FeatureTree::updateElementBehaviours()
{
    if(structuredMesh)
        return ;
    double n_void ;

    if ( is2D() )
    {
        double remainder = 0 ;
        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            if ( i->first != -1 )
            {
                remainder += scalingFactors[i->first] ;
            }
        }
        scalingFactors[-1] = 1.-remainder ;
        layer2d[-1] = dtree ;


        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            std::vector<double> scales ;
            if ( layer2d.size() > 1 )
            {
                for ( auto j = ++layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    if ( i != j )
                    {
                        scales.push_back ( scalingFactors[j->first]/scalingFactors[i->first] );
                    }
                }
            }

            int setcount = 0 ;
            std::cerr << "\r updating behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount++ << "/" << i->second->begin().size() << "    " << std::flush ;
            for ( auto j = i->second->begin() ; j != i->second->end() ; j++ )
            {

                if ( setcount++ % 1000 == 0 )
                {
                    std::cerr << "\r updating behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount << "/" << j.size() << "    " << std::flush ;
                }

                j->refresh ( father2D ) ;
                Form * bf =  getElementBehaviour ( j, i->first, true );

                if ( bf )
                {
                    j->setBehaviour ( dtree, bf ) ;
                }
            }
            std::cerr << " ...done" << std::endl ;

        }

        setBehaviours = true ;
    }
    else
    {
        std::cerr << " setting behaviours..." << std::flush ;
        int setcount = 0 ;

        for ( auto i = dtree3D->begin() ; i < dtree3D->end() ; i++ )
        {
            i->refresh ( father3D ) ;
        }



        #pragma omp parallel for schedule(runtime)
        for ( auto i = dtree3D->begin() ; i < dtree3D->end() ; i++ )
        {
            if ( setcount % 1000 == 0 )
            {
                std::cerr << "\r setting behaviours... tet " << setcount << "/" << i.size() << std::flush ;
            }

            Form *b = getElementBehaviour ( i, true ) ;

            if ( b )
            {
                i->setBehaviour ( dtree3D, b ) ;
            }

            if ( !i->getBehaviour() || i->getBehaviour()->type == VOID_BEHAVIOUR )
            {
                i->setBehaviour ( dtree3D, getElementBehaviour ( i ) ) ;
            }

            n_void++ ;
            setcount++ ;
        }

        std::cerr << " ...done" << std::endl ;

        setBehaviours = true ;
    }


    for ( size_t i = 0 ; i < tree.size() ; i++ )
    {
        if ( !tree[i]->isEnrichmentFeature )
        {
            tree[i]->isUpdated = false ;
        }
    }
}

void FeatureTree::enrich()
{

    enrichmentChange = false ;
    lastEnrichmentId = getNodes().size() ;

    std::cerr << "\r enriching... feature " << 0 << "/" << tree.size() << std::flush ;

    for ( size_t i = 1 ; i < tree.size() ; i++ )
    {
        if ( is3D() )
        {
            std::vector<Mesh <DelaunayTetrahedron, DelaunayTreeItem3D > *> extra3dMeshes ;

            if ( tree[i]->isEnrichmentFeature && ( dynamic_cast<EnrichmentFeature *> ( tree[i] )->moved() || !state.enriched ) )
            {
                if ( !state.enriched )
                {
                    dynamic_cast<EnrichmentFeature *> ( tree[i] )->update ( dtree3D ) ;
                }

                dynamic_cast<EnrichmentFeature *> ( tree[i] )->enrich ( lastEnrichmentId, dtree3D ) ;

                enrichmentChange = true ;
                reuseDisplacements = false ;

            }

            if ( i % 10 == 0 )
            {
                std::cerr << "\r enriching... feature " << i + 1 << "/" << tree.size() << std::flush ;
            }
        }
        else
        {
            std::vector<Mesh <DelaunayTriangle, DelaunayTreeItem > *> extra2dMeshes ;
            if ( layer2d.size() > 1 )
            {
                for ( auto i = ++layer2d.begin() ; i != layer2d.end() ; i++ )
                {
                    extra2dMeshes.push_back ( i->second ) ;
                }
            }

            if ( tree[i]->isEnrichmentFeature && ( dynamic_cast<EnrichmentFeature *> ( tree[i] )->moved() || !state.enriched ) )
            {
                if ( !state.enriched )
                {
                    dynamic_cast<EnrichmentFeature *> ( tree[i] )->update ( dtree ) ;
                }
                dynamic_cast<EnrichmentFeature *> ( tree[i] )->enrich ( lastEnrichmentId, dtree ) ;

                enrichmentChange = true ;
                reuseDisplacements = false ;

            }

            if ( i % 10 == 0 )
            {
                std::cerr << "\r enriching... feature " << i + 1 << "/" << tree.size() << std::flush ;
            }
        }
    }

    if(enrichmentChange)
        K->clear();

    std::cerr << " ...done" << std::endl ;
}

void FeatureTree::assemble()
{
    K->clearElements();

    if(father3D && father3D->getOrder() != elemOrder)
        delete father3D ;

    if ( !father3D )
    {
        father3D = new TetrahedralElement ( elemOrder ) ;
    }

//     father3D->compileAndPrecalculate() ;

    if(father2D && father2D->getOrder() != elemOrder)
        delete father2D ;

    if ( !father2D )
    {
        father2D = new TriElement ( elemOrder ) ;
    }


    if ( is2D() )
    {
        numdofs = dtree->getLastNodeId() ;
        for ( auto i = layer2d.begin() ; i != layer2d.end() ; i++ )
        {
            for ( auto j = i->second->begin() ; j != i->second->end() ; j++ )
            {
                if ( j.getPosition() % 1000 == 0 )
                {
                    std::cerr << "\r assembling stiffness matrix, layer "<< i->first << " ... triangle " << j.getPosition() + 1 << "/" << j.size() << std::flush ;
                }

                if (j->getBehaviour() && j->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    j->refresh ( father2D ) ;
                    j->getBehaviour()->preProcess ( deltaTime, j->getState() ) ;
// 					if(!tris[j]->getBehaviour()->fractured())
                    K->add ( j, scalingFactors[i->first] ) ;
                }

            }
            K->setMaxDof ( std::max ( getNodes().size(),lastEnrichmentId ) ) ;
            std::cerr << " ...done." << std::endl ;
        }

    }
    else
    {
        numdofs = dtree3D->getLastNodeId() ;

        for ( auto j = dtree3D->begin() ; j != dtree3D->end() ; j++ )
        {
            if ( j->getBehaviour() && j->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                if ( j.getPosition() % 1000 == 0 )
                {
                    std::cerr << "\r assembling stiffness matrix... tetrahedron " << j.getPosition() + 1 << "/" << j.size() << std::flush ;
                }

                j->refresh ( father3D ) ;
                j->getBehaviour()->preProcess ( deltaTime, j->getState() ) ;
                K->add ( j ) ;
            }
        }
        K->setMaxDof ( std::max ( getNodes().size(),lastEnrichmentId ) ) ;
        std::cerr << " ...done." << std::endl ;

    }
}

std::vector<DelaunayTriangle> FeatureTree::getSnapshot2D() const
{
    std::vector<DelaunayTriangle> copy ;

    for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
    {
        copy.push_back ( *i ) ;
        copy.back().setBehaviour ( dtree, i->getBehaviour()->getCopy() ) ;
        copy.back().getState().initialize ( dtree ) ;
    }

    return copy ;
}

Vector FeatureTree::stressFromDisplacements()
{
    state.setStateTo ( XFEM_STEPPED, false ) ;

    VirtualMachine vm ;
    if ( dtree )
    {
        Vector stress ( 0.f, 3 * 3 * dtree->begin().size() ) ;

        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                std::valarray<Point *> pts ( 3 ) ;
                pts[0] =  i->first ;
                pts[1] =  i->second ;
                pts[2] =  i->third ;

                Vector str ( 0., 3*3 ) ;
                i->getState().getField ( REAL_STRESS_FIELD, pts, str, false, &vm ) ;

                for ( size_t j = 0 ; j < 9 ; j++ )
                {
                    stress[i.getPosition() * 3 * 3 + j] = str[j] ;
                }

                std::cerr << "\r computing stress... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
            }

        }

        std::cerr << " ...done." << std::endl ;
        return stress ;
    }
    else
    {
        Vector stress ( 0., 4 * 6 * dtree3D->begin().size() ) ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            std::valarray<Point *> pts ( 4 ) ;
            pts[0] =  i->first ;
            pts[1] =  i->second ;
            pts[2] =  i->third ;
            pts[2] =  i->fourth ;

            Vector str ( 0., 24 ) ;
            i->getState().getField ( REAL_STRESS_FIELD, pts, str, false, &vm ) ;

            for ( size_t j = 0 ; j < 24 ; j++ )
            {
                stress[i.getPosition() * 4 * 6 + j] = str[j] ;
            }

            std::cerr << "\r computing stress... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;

        }

        std::cerr << " ...done." << std::endl ;
        return stress ;
    }
}

const Vector &FeatureTree::getDisplacements ( int g, bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    return  K->getDisplacements() ;
}

Vector FeatureTree::getDisplacements ( Point * pt, int g , bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    if ( is2D() )
    {
        Vector ret ;
        std::vector<DelaunayTriangle *> elements = dtree->getConflictingElements ( pt ) ;
        for ( size_t i = 0 ; i < elements.size() ; i++ )
        {

            if ( elements[i]->in ( *pt ) )
            {
                ret.resize ( elements[i]->getBehaviour()->getNumberOfDegreesOfFreedom() , 0. ) ;
                elements[i]->getState().getField ( DISPLACEMENT_FIELD, *pt, ret, false ) ;
                return ret ;
            }
        }
        ret.resize ( 2 , 0. ) ;
        return ret ;
    }
    else
    {
        Vector ret ;
        std::vector<DelaunayTetrahedron *> elements = dtree3D->getConflictingElements ( pt ) ;
        for ( size_t i = 0 ; i < elements.size() ; i++ )
        {
            if ( elements[i]->in ( *pt ) )
            {
                ret.resize ( elements[i]->getBehaviour()->getNumberOfDegreesOfFreedom() , 0. ) ;
                elements[i]->getState().getField ( DISPLACEMENT_FIELD, *pt, ret, false ) ;
                return ret ;
            }
        }
        ret.resize ( 3 , 0. ) ;
        return ret ;
    }
}

Vector FeatureTree::getField ( FieldType f, Point * pt, int g , bool stepTree , bool localTime)
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    if ( is2D() )
    {
        Vector ret ;
        std::vector<DelaunayTriangle *> elements = dtree->getConflictingElements ( pt ) ;
        for ( size_t i = 0 ; i < elements.size() ; i++ )
        {

            if ( elements[i]->in ( *pt ) )
            {
                int b = 1 ;
                if(dynamic_cast<Viscoelasticity *>(elements[i]->getBehaviour()))
                    b = dynamic_cast<Viscoelasticity *>(elements[i]->getBehaviour())->blocks ;
                ret.resize ( fieldTypeElementarySize( f, SPACE_TWO_DIMENSIONAL, b ) , 0. ) ;
                double prevT = pt->getT() ;
                if(localTime)
                {
                    double t0 = elements[i]->getBoundingPoint(0).getT() ;
                    double t1 = elements[i]->getBoundingPoint( elements[i]->getBoundingPoints().size() -1 ).getT() ;

                    pt->setT( t0 + (t1-t0)*(prevT+1)*0.5 ) ;
                }
                elements[i]->getState().getField ( f, *pt, ret, false ) ;
                pt->setT( prevT ) ;
                return ret ;
            }
        }
        ret.resize ( fieldTypeElementarySize( f, SPACE_TWO_DIMENSIONAL, 1 ) , 0. ) ;
        return ret ;
    }
    else
    {
        Vector ret ;
        std::vector<DelaunayTetrahedron *> elements = dtree3D->getConflictingElements ( pt ) ;
        for ( size_t i = 0 ; i < elements.size() ; i++ )
        {
            if ( elements[i]->in ( *pt ) )
            {
                int b = 1 ;
                if(dynamic_cast<Viscoelasticity *>(elements[i]->getBehaviour()))
                    b = dynamic_cast<Viscoelasticity *>(elements[i]->getBehaviour())->blocks ;
                ret.resize ( fieldTypeElementarySize( f, SPACE_THREE_DIMENSIONAL, b ) , 0. ) ;
                elements[i]->getState().getField ( f, *pt, ret, false ) ;
                return ret ;
            }
        }
        ret.resize ( fieldTypeElementarySize( f, SPACE_THREE_DIMENSIONAL, 1 ) , 0. ) ;
        return ret ;
    }
}

double FeatureTree::getField ( std::string f, Point * pt, int g , bool stepTree , bool localTime)
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    if ( is2D() )
    {
        std::vector<DelaunayTriangle *> elements = dtree->getConflictingElements ( pt ) ;
        for ( size_t i = 0 ; i < elements.size() ; i++ )
        {

            if ( elements[i]->in ( *pt ) )
            {
                std::map<std::string, double > dummy ;
                if( dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables *>( &(elements[i]->getState()) ) )
                    return dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables *>( &(elements[i]->getState()) )->get( f, dummy ) ;
                else
                    return 0. ;
            }
        }
    }
    return 0. ;
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain (bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    VirtualMachine vm ;

    if ( dtree )
    {
        std::pair<Vector , Vector > stress_strain ( Vector ( 0., dtree->begin()->getBoundingPoints().size() * 3 * dtree->begin().size() ), Vector ( 0., dtree->begin()->getBoundingPoints().size() * 3 * dtree->begin().size() ) ) ;
        for(auto l = layer2d.begin() ; l !=layer2d.end() ; l++)
        {


            int donecomputed = 0 ;
//
//         #pragma omp parallel for shared(donecomputed) schedule(runtime)
            #pragma omp parallel
            {
                #pragma omp single
                for ( auto i = l->second->begin() ; i != l->second->end() ; i++ )
                {
                    #pragma omp task firstprivate(i)
                    if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                    {
                        // 				std::valarray<Point *> pts(3) ;
                        // 				pts[0] =  elements[i]->first ;
                        // 				pts[1] =  elements[i]->second ;
                        // 				pts[2] =  elements[i]->third ;

                        Vector strain ( 0., 3*i->getBoundingPoints().size() ) ;
                        Vector stress ( 0., 3*i->getBoundingPoints().size() ) ;
                        i->getState().getField ( STRAIN_FIELD, REAL_STRESS_FIELD, i->getBoundingPoints(), strain, stress, false, &vm ) ;

                        for ( size_t j = 0 ; j < i->getBoundingPoints().size() * 3 ; j++ )
                        {
                            stress_strain.first[i.getPosition() * i->getBoundingPoints().size() * 3 + j] = stress[j] ;
                            stress_strain.second[i.getPosition() * i->getBoundingPoints().size() * 3 + j] = strain[j] ;
                        }

                        if ( donecomputed % 10000 == 0 )
                        {
                            std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << i.size() << std::flush ;
                        }
                    }

                    donecomputed++ ;
                }
            }

            std::cerr << " ...done." << std::endl ;
        }
        return stress_strain ;

    }
    else
    {

        std::pair<Vector , Vector > stress_strain ( Vector ( 0.f, 4 * 6 * dtree3D->begin().size() ), Vector ( 0.f, 4 * 6 * dtree3D->begin().size() ) ) ;
        int donecomputed = 0 ;

        #pragma omp parallel
        {
            #pragma omp single
            for ( auto i  = dtree3D->begin() ; i != dtree3D->end() ; i++ )
            {
                #pragma omp task firstprivate(i)
                if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    std::valarray<Point *> pts ( 4 ) ;
                    pts[0] =  i->first ;
                    pts[1] =  i->second ;
                    pts[2] =  i->third ;
                    pts[3] =  i->fourth ;

                    Vector strain ( 0., 24 ) ;
                    Vector stress ( 0., 24 ) ;
                    i->getState().getField ( STRAIN_FIELD, REAL_STRESS_FIELD, pts, strain, stress, false, &vm ) ;

                    for ( size_t j = 0 ; j < 24 ; j++ )
                    {
                        stress_strain.first[i.getPosition() * 4 * 6 + j] = strain[j] ;
                        stress_strain.second[i.getPosition() * 4 * 6 + j] = stress[j] ;
                    }
                }

                if ( donecomputed % 1000 == 0 )
                {
                    std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << i.size() << std::flush ;
                }

                donecomputed++ ;
            }
        }

        std::cerr << " ...done." << std::endl ;
        return stress_strain ;
    }
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrainInLayer ( int g, bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    VirtualMachine vm ;

    if ( dtree )
    {
        Mesh<DelaunayTriangle, DelaunayTreeItem> * msh = dtree ;
        if ( g != -1 && layer2d.find ( g ) != layer2d.end() )
        {
            msh = layer2d[g] ;
        }

        std::pair<Vector , Vector > stress_strain ( Vector ( 0., msh->begin()->getBoundingPoints().size() * 3 * msh->begin().size() ), Vector ( 0., msh->begin()->getBoundingPoints().size() * 3 * msh->begin().size() ) ) ;
        int donecomputed = 0 ;

        #pragma omp parallel
        {
            #pragma omp single
            for ( auto i  = msh->begin() ; i != msh->end() ; i++ )
            {
                #pragma omp task firstprivate(i)
                if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    // 				std::valarray<Point *> pts(3) ;
                    // 				pts[0] =  elements[i]->first ;
                    // 				pts[1] =  elements[i]->second ;
                    // 				pts[2] =  elements[i]->third ;

                    Vector strain ( 0., 3*i->getBoundingPoints().size() ) ;
                    Vector stress ( 0., 3*i->getBoundingPoints().size() ) ;
                    i->getState().getField ( STRAIN_FIELD, REAL_STRESS_FIELD, i->getBoundingPoints(), strain, stress, false, &vm ) ;

                    for ( size_t j = 0 ; j < i->getBoundingPoints().size() * 3 ; j++ )
                    {
                        stress_strain.first[i.getPosition() * i->getBoundingPoints().size() * 3 + j] = stress[j] ;
                        stress_strain.second[i.getPosition() * i->getBoundingPoints().size() * 3 + j] = strain[j] ;
                    }

                    if ( donecomputed % 10000 == 0 )
                    {
                        std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << i.size() << std::flush ;
                    }
                }

                donecomputed++ ;
            }
        }

        std::cerr << " ...done." << std::endl ;
        return stress_strain ;
    }
    else
    {
//         std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

//         if ( g != -1 && layer3d.find ( g ) != layer3d.end() )
//         {
//             tets = layer3d[g]->getElements() ;
//         }

        std::pair<Vector , Vector > stress_strain ( Vector ( 0.f, 4 * 6 * dtree3D->begin().size() ), Vector ( 0.f, 4 * 6 * dtree3D->begin().size() ) ) ;
        int donecomputed = 0 ;

        #pragma omp parallel
        {
            #pragma omp single
            for ( auto i  = dtree3D->begin() ; i != dtree3D->end() ; i++ )
            {
                #pragma omp task firstprivate(i)
                if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    std::valarray<Point *> pts ( 4 ) ;
                    pts[0] =  i->first ;
                    pts[1] =  i->second ;
                    pts[2] =  i->third ;
                    pts[3] =  i->fourth ;

                    Vector strain ( 0., 24 ) ;
                    Vector stress ( 0., 24 ) ;
                    i->getState().getField ( STRAIN_FIELD, REAL_STRESS_FIELD, pts, strain, stress, false, &vm ) ;

                    for ( size_t j = 0 ; j < 24 ; j++ )
                    {
                        stress_strain.first[i.getPosition() * 4 * 6 + j] = stress[j] ;
                        stress_strain.second[i.getPosition() * 4 * 6 + j] = strain[j] ;
                    }
                }

                if ( donecomputed % 1000 == 0 )
                {
                    std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << i.size() << std::flush ;
                }

                donecomputed++ ;
            }
        }

        std::cerr << " ...done." << std::endl ;
        return stress_strain ;
    }
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFlux ( bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    if ( dtree != nullptr )
    {

        std::pair<Vector , Vector > grad_flux ( Vector ( 0., dtree->begin()->getBoundingPoints().size() * 2 * dtree->begin().size() ), Vector ( 0., dtree->begin()->getBoundingPoints().size() * 2 * dtree->begin().size() ) ) ;

        for ( auto i  = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

                Vector gradient ( 0., 2*i->getBoundingPoints().size() ) ;
                Vector flux ( 0., 2*i->getBoundingPoints().size() ) ;
                i->getState().getField ( GRADIENT_FIELD, FLUX_FIELD, i->getBoundingPoints(), gradient, flux, false ) ;

                for ( size_t j = 0 ; j < i->getBoundingPoints().size() * 2 ; j++ )
                {
                    grad_flux.first[i.getPosition() * i->getBoundingPoints().size() * 2 + j] = gradient[j] ;
                    grad_flux.second[i.getPosition() * i->getBoundingPoints().size() * 2 + j] = flux[j] ;
                }

                if ( i.getPosition() % 1000 == 0 )
                {
                    std::cerr << "\r computing gradient+flux... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
                }
            }
        }

        std::cerr << " ...done." << std::endl ;
        return grad_flux ;
    }
    else
    {

        size_t npoints = dtree3D->begin()->getBoundingPoints().size() ;
        std::pair<Vector , Vector > grad_flux ( Vector ( 0.f, npoints * 3 * dtree3D->begin().size() ), Vector ( 0.f, npoints * 3 * dtree3D->begin().size() ) ) ;

        for ( auto i  = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {

            Vector gradient ( 0., 3*i->getBoundingPoints().size() ) ;
            Vector flux ( 0., 3*i->getBoundingPoints().size() ) ;
            i->getState().getField ( GRADIENT_FIELD, FLUX_FIELD, i->getBoundingPoints(), gradient, flux, false ) ;

            for ( size_t j = 0 ; j < npoints * 3 ; j++ )
            {
                grad_flux.first[i.getPosition() * npoints * 3 + j] = gradient[j] ;
                grad_flux.second[i.getPosition() * npoints * 3 + j] = flux[j] ;
            }

            if ( i.getPosition() % 1000 == 0 )
            {
                std::cerr << "\r computing gradient+flux... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
            }

        }

        std::cerr << " ...done." << std::endl ;
        return grad_flux ;
    }
}

std::vector<int>FeatureTree:: listLayers() const
{
    std::vector<int> ret ;
    if ( is2D() )
    {
        for ( auto i = layer2d.begin() ; i!= layer2d.end() ; ++i )
        {
            ret.push_back ( i->first );
        }
    }
    else
    {

    }

    return ret ;
}

void FeatureTree::setDiscretizationParameters ( ConfigTreeItem * config, ConfigTreeItem * def )
{
    if ( def == nullptr )
    {
        def = new ConfigTreeItem() ;
        def->addChild ( new ConfigTreeItem ( def, "sampling_number", 4 ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "order", "LINEAR" ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "sampling_restriction", 0 ) ) ;
    }
    double sampling = config->getData ( "sampling_number", def->getData ( "sampling_number" ) ) ;
    std::string order = config->getStringData ( "order", def->getStringData ( "order" ) ) ;
    int restriction = config->getData ( "sampling_restriction", def->getData ( "sampling_restriction" ) ) ;
    setSamplingNumber ( sampling ) ;
    setOrder ( Enum::getOrder ( order ) ) ;
    setSamplingRestriction ( restriction ) ;

    std::vector<ConfigTreeItem *> zones = config->getAllChildren ( "refinement_zone" ) ;
    for ( size_t i = 0 ; i < zones.size() ; i++ )
    {
        Sample * s = zones[i]->getSample() ;
        addRefinementZone ( dynamic_cast<Rectangle *> ( s ) ) ;
    }

}

Vector FeatureTree::setSteppingParameters ( ConfigTreeItem * config, ConfigTreeItem * def )
{
    if ( def == nullptr )
    {
        def = new ConfigTreeItem() ;
        def->addChild ( new ConfigTreeItem ( def, "time_step", 1. ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "minimum_time_step", 1e-9 ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "maximum_iterations_per_step", 256 ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "number_of_time_steps", 1 ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "first_time_step", 0.1 ) ) ;
        def->addChild ( new ConfigTreeItem ( def, "logarithmic", "FALSE" ) ) ;
    }
    double deltaTime = config->getData ( "time_step", def->getData ( "time_step" ) ) ;
    double minDeltaTime = config->getData ( "minimum_time_step", def->getData ( "minimum_time_step" ) ) ;
    int maxIter = config->getData ( "maximum_iterations_per_step", def->getData ( "maximum_iterations_per_step" ) ) ;
    int nSteps = config->getData ( "number_of_time_steps", def->getData ( "number_of_time_steps" ) ) ;
    setDeltaTime ( deltaTime ) ;
    setMinDeltaTime ( minDeltaTime ) ;
    setMaxIterationsPerStep ( maxIter ) ;
    Vector cinstants ( nSteps+1 ) ;
    if( config->hasChild( "solver_precision" ) )
        K->setEpsilon( config->getData ( "solver_precision", POINT_TOLERANCE ) ) ;
    if ( config->hasChild ( "list_of_time_steps" ) )
    {
        if(config->getStringData("list_of_time_steps").find(',') != std::string::npos)
            cinstants = ConfigTreeItem::readLineAsVector( config->getStringData("list_of_time_steps") ) ;
        else
            cinstants = config->getChild ( "list_of_time_steps" )->readVectorFromFile() ;
    }
    else
    {
        if ( config->hasChild ( "next_time_step" ) )
        {
            Function f = config->getChild ( "next_time_step" )->getFunction() ;
            double t = 0. ;
            double x = deltaTime ;
            VirtualMachine vm ;
            cinstants[0] = 0. ;
            cinstants[1] = deltaTime ;
            for ( size_t i = 2 ; i < cinstants.size() ; i++ )
            {
                x = vm.eval ( f, x, 0.,0., t ) ;
                cinstants[i] = cinstants[i-1]+x ;
                t = t+x ;
            }
        }
        else
        {
            if( config->getStringData("logarithmic", def->getStringData("logarithmic")) == "TRUE" )
            {
                cinstants[0] = 0. ;
                cinstants[1] = config->getData ( "first_time_step", def->getData ( "first_time_step" ) ) ;
                for(size_t i = 2 ; i < cinstants.size() ; i++)
                {
                    double logprev = std::log10(cinstants[i-1]) ;
                    double lognext = logprev+deltaTime ;
                    cinstants[i] = std::pow(10., lognext) ;
                }
            }
            else
            {
                for ( size_t i = 0 ; i < cinstants.size() ; i++ )
                {
                    cinstants[i] = deltaTime*i ;
                }
            }
        }
    }
    return cinstants ;
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFluxInLayer ( int g, bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }

    if ( dtree != nullptr )
    {

        std::pair<Vector , Vector > grad_flux ( Vector ( 0., layer2d[g]->begin()->getBoundingPoints().size() * 2 * layer2d[g]->begin().size() ), Vector ( 0., layer2d[g]->begin()->getBoundingPoints().size() * 2 * layer2d[g]->begin().size() ) ) ;

        for ( auto i  = layer2d[g]->begin() ; i < layer2d[g]->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

                Vector gradient ( 0., 2*i->getBoundingPoints().size() ) ;
                Vector flux ( 0., 2*i->getBoundingPoints().size() ) ;
                i->getState().getField ( GRADIENT_FIELD, FLUX_FIELD, i->getBoundingPoints(), gradient, flux, false ) ;

                for ( size_t j = 0 ; j < layer2d[g]->begin()->getBoundingPoints().size() * 2 ; j++ )
                {
                    grad_flux.first[i.getPosition() * layer2d[g]->begin()->getBoundingPoints().size() * 2 + j] = gradient[j] ;
                    grad_flux.second[i.getPosition() * layer2d[g]->begin()->getBoundingPoints().size() * 2 + j] = flux[j] ;
                }

                if ( i.getPosition() % 1000 == 0 )
                {
                    std::cerr << "\r computing gradient+flux... element " << i.getPosition() + 1 << "/" << layer2d[g]->begin().size() << std::flush ;
                }
            }
        }

        std::cerr << " ...done." << std::endl ;
        return grad_flux ;
    }
    else
    {

//         tets = layer3d[g]->getElements() ;


        size_t npoints = dtree3D->begin()->getBoundingPoints().size() ;
        std::pair<Vector , Vector > grad_flux ( Vector ( 0.f, npoints * 3 * dtree3D->begin().size() ), Vector ( 0.f, npoints * 3 * dtree3D->begin().size() ) ) ;

        for ( auto i  = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {

            Vector gradient ( 0., 3*i->getBoundingPoints().size() ) ;
            Vector flux ( 0., 3*i->getBoundingPoints().size() ) ;
            i->getState().getField ( GRADIENT_FIELD, FLUX_FIELD, i->getBoundingPoints(), gradient, flux, false ) ;

            for ( size_t j = 0 ; j < npoints * 3 ; j++ )
            {
                grad_flux.first[i.getPosition() * npoints * 3 + j] = gradient[j] ;
                grad_flux.second[i.getPosition() * npoints * 3 + j] = flux[j] ;
            }

            if ( i.getPosition() % 1000 == 0 )
            {
                std::cerr << "\r computing gradient+flux... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
            }

        }

        std::cerr << " ...done." << std::endl ;
        return grad_flux ;
    }
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFlux ( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>::iterator begin, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>::iterator end , bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }
    std::pair<Vector , Vector > stress_strain ( Vector ( 4 * 3 * begin.size() ), Vector ( 4 * 3 * begin.size() ) ) ;

    for ( auto i  = begin ; i != end ; i++ )
    {
        std::valarray<Point *> pts ( 4 ) ;
        pts[0] =  i->first ;
        pts[1] =  i->second ;
        pts[2] =  i->third ;
        pts[3] =  i->fourth ;

        Vector gradient ( 0., 12 ) ;
        Vector flux ( 0., 12 ) ;
        i->getState().getField ( GRADIENT_FIELD, FLUX_FIELD, i->getBoundingPoints(), gradient, flux, false ) ;

        for ( size_t j = 0 ; j < 4 ; j++ )
        {
            for ( size_t k = 0 ; k < 3 ; k++ )
            {
                stress_strain.first[i.getPosition() * 4 * 3 + j * 3 + k] = gradient[j * 3 + k] ;
                stress_strain.second[i.getPosition() * 4 * 3 + j * 3 + k] = flux[j * 3 + k] ;
            }
        }

//			if(i%1000 == 0)
//				std::cerr << "\r computing gradient+flux... element " << i+1 << "/" << tets.size() << std::flush ;
    }

    std::cerr << " ...done." << std::endl ;
    return stress_strain ;
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain ( Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >::iterator begin, Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >::iterator end, bool stepTree )
{
    if ( stepTree )
    {
        state.setStateTo ( XFEM_STEPPED, false ) ;
    }
    std::pair<Vector , Vector > stress_strain ( Vector ( begin->getBoundingPoints().size() * 6 * begin.size() ), Vector ( begin->getBoundingPoints().size() * 6 * begin.size() ) ) ;

    for ( auto i = begin ; i != end ; i++ )
    {
        std::valarray<Point *> pts ( 4 ) ;
        pts[0] =  i->first ;
        pts[1] =  i->second ;
        pts[2] =  i->third ;
        pts[3] =  i->fourth ;

        Vector strain ( 0., i->getBoundingPoints().size() * 6 ) ;
        Vector stress ( 0., i->getBoundingPoints().size() * 6 ) ;
        i->getState().getField ( STRAIN_FIELD, REAL_STRESS_FIELD, i->getBoundingPoints(), strain, stress, false ) ;

        for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
        {
            for ( size_t k = 0 ; k < 6 ; k++ )
            {
                stress_strain.first[i.getPosition() * i->getBoundingPoints().size() * 6 + j * 6 + k] = stress[j * 6 + k] ;
                stress_strain.second[i.getPosition() * i->getBoundingPoints().size() * 6 + j * 6 + k] = strain[j * 6 + k] ;
            }
        }

        if ( i.getPosition() % 1000 == 0 )
        {
            std::cerr << "\r computing strain+stress... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
        }
    }

    std::cerr << " ...done." << std::endl ;
    return stress_strain ;
}

Vector FeatureTree::strainFromDisplacements()
{
    state.setStateTo ( XFEM_STEPPED, false ) ;

    if ( dtree )
    {
        Vector strain ( 0.f, 3 * 3 * dtree->begin().size() ) ;

        for ( auto i  = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                std::valarray<Point *> pts ( 3 ) ;
                pts[0] =  i->first ;
                pts[1] =  i->second ;
                pts[2] =  i->third ;

                Vector str ( 0., 9 ) ;
                i->getState().getField ( STRAIN_FIELD, pts, str, false ) ;

                for ( size_t j = 0 ; j < 9 ; j++ )
                {
                    strain[i.getPosition() * 3 * 3 + j] = str[j] ;
                }

                std::cerr << "\r computing strain... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
            }
        }

        std::cerr << " ...done." << std::endl ;
        return strain ;
    }
    else
    {
//         std::vector<DelaunayTetrahedron *> elements3D = dtree3D->getElements() ;
        Vector strain ( 0., 4 * 6 * dtree3D->begin().size() ) ;

        for ( auto i  = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            std::valarray<Point *>  pts ( 4 ) ;
            pts[0] =  i->first ;
            pts[1] =  i->second ;
            pts[2] =  i->third ;
            pts[3] =  i->fourth ;

            Vector str ( 0., 24 ) ;
            i->getState().getField ( STRAIN_FIELD, pts, str, false ) ;

            for ( size_t j = 0 ; j < 24 ; j++ )
            {
                strain[i.getPosition() * 4 * 6 + j] = str[j] ;
            }

            std::cerr << "\r computing strain... element " << i.getPosition() + 1 << "/" << i.size() << std::flush ;
        }

        std::cerr << " ...done." << std::endl ;
        return strain ;
    }

}

Assembly *FeatureTree::getAssembly ( bool forceReassembly )
{
    if ( forceReassembly )
    {
        state.setStateTo ( ASSEMBLED, false ) ;
    }
    return K ;
}

void FeatureTree::insert ( Point *p )
{
    if(structuredMesh)
        return ;
    double pointDensity = 0 ;

    if ( is2D() )
    {
        pointDensity = .6 * sqrt ( tree[0]->area() ) / meshPoints.size() ;
    }
    else
    {
        pointDensity = .6 * pow ( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;
    }

    Feature *mother = nullptr;

    for ( size_t i  = 0 ; i < this->tree.size() ; i ++ )
    {
        if ( this->tree[i]->in ( ( *p ) ) )
        {
            bool yes = true ;

            for ( size_t k  =  0 ; k < this->tree[i]->getChildren().size() && yes; k++ )
            {
                if ( this->tree[i]->getChild ( k )->in ( ( *p ) ) )
                {
                    yes = false ;
                }
            }

            if ( yes )
            {
                mother = this->tree[i] ;
                break ;
            }
        }
    }

    bool yes = true ;

    for ( size_t k  =  0 ; k <  mother->getChildren().size() && yes; k++ )
    {
        if ( mother->getChild ( k )->inBoundary ( *p, pointDensity ) )
        {
            yes = false ;
        }
    }

    if ( yes )
    {
        this->meshPoints.push_back ( std::pair<Point *, Feature *> ( p, mother ) ) ;

        if ( dtree != nullptr )
        {
            this->dtree->insert ( p ) ;
        }
        else
        {
            this->dtree3D->insert ( p ) ;
        }
    }
}

bool FeatureTree::solverConverged() const
{
    return solverConvergence ;
}

bool FeatureTree::behaviourChanged() const
{
    return behaviourChange ;
}

bool FeatureTree::enrichmentChanged() const
{
    return enrichmentChange ;
}

void FeatureTree::forceEnrichmentChange()
{
    enrichmentChange = true ;
}

void FeatureTree::elasticStep()
{
    Vector lastx ( K->getDisplacements() ) ;
    K->clear() ;
    assemble() ;
    solve() ;
    bool prevElastic = elastic ;
    elastic = true ;
    stepElements() ;
    elastic = prevElastic ;
}

void FeatureTree::solve()
{
    Vector lastx ( K->getDisplacements() ) ;

    if ( enrichmentChange || needMeshing )
    {
        K->clear() ;
    }

    timeval time0, time1 ;
    gettimeofday ( &time0, nullptr );

    if ( dtree )
    {
        std::cerr << "finding nodes for boundary conditions... " << std::flush ;
        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
        {
            for ( auto i = j->second->begin() ; i != j->second->end() ; i++ )
            {
                i->applyBoundaryCondition ( K ) ;
            }
        }

    }
    else
    {
        std::cerr << "finding nodes for boundary conditions... " << std::flush ;
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            i->applyBoundaryCondition ( K ) ;
        }
    }


    gettimeofday ( &time1, nullptr );
    double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << "...done. Time (s) " << delta / 1e6 << std::endl ;

    gettimeofday ( &time0, nullptr );

    std::cerr << "Applying boundary conditions... " << std::flush ;
    for ( size_t i = 0 ; i < boundaryCondition.size() ; ++i )
    {
//         if ( i%20 == 0 )
//         {
        std::cerr << "\rApplying boundary conditions... " << i+1 << "/" << boundaryCondition.size() << std::flush ;
//         }

        boundaryCondition[i]->step( now, realDeltaTime ) ;

        if ( dtree )
        {
            boundaryCondition[i]->apply ( K, dtree ) ;
        }

        if ( dtree3D )
        {
            boundaryCondition[i]->apply ( K, dtree3D ) ;
        }
    }
    gettimeofday ( &time1, nullptr );
    delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << "\rApplying boundary conditions... " << boundaryCondition.size() << "/" << boundaryCondition.size() << "...done  Time (s) " << delta / 1e6 << std::endl  ;

    if ( solverConvergence || reuseDisplacements )
    {

        solverConvergence = K->cgsolve ( lastx ) ;

        Vector r = K->getMatrix() * K->getDisplacements() - K->getForces() ;
        double perror = residualError ;
        residualError = sqrt ( parallel_inner_product ( &r[0], &r[0], r.size() ) ) ;

        if ( perror > residualError || solverConvergence )
        {
            reuseDisplacements = true;
        }
    }
    else
    {
        lastx = 0 ;


        solverConvergence = K->cgsolve() ;

// 		dtree->project(coarseTrees[3], K->getDisplacements(), coarseAssemblies[3]->getDisplacements(), false) ;
        Vector r = K->getMatrix() * K->getDisplacements() - K->getForces() ;
        double perror = residualError ;
        residualError = sqrt ( parallel_inner_product ( &r[0], &r[0], r.size() ) ) ;

        if ( perror > residualError || solverConvergence )
        {
            reuseDisplacements = true;
        }
    }

}

void FeatureTree::stepXfem()
{
    enrichmentChange = false ;

    if ( solverConvergence )
    {
        std::vector<EnrichmentFeature *> featuresToStep ;
        for ( size_t i = 0 ; i < tree.size() ; i++ )
        {
            if ( tree[i]->isEnrichmentFeature )
            {

                featuresToStep.push_back(dynamic_cast<EnrichmentFeature *> ( tree[i] ));
            }
        }
        if ( is2D() )
        {

// 			#pragma omp parallel for schedule(runtime)

            for ( size_t i = 0 ; i < featuresToStep.size() ; i++ )
            {
                featuresToStep[i]->step( deltaTime, &K->getForces(), dtree ) ;
                bool moved = featuresToStep[i]->moved() ;
                enrichmentChange = enrichmentChange || moved;

                if ( moved )
                {
                    reuseDisplacements = false ;
                }
            }
            for ( size_t i = 0 ; i < enrichmentManagers.size() ; i++ )
            {
                bool moved = enrichmentManagers[i]->step(deltaTime, &K->getForces(), dtree)  ;
                if ( moved )
                {
                    reuseDisplacements = false ;
                }
            }

            for ( size_t i = 0 ; i < tree.size() ; i++ )
            {
                if ( !tree[i]->isEnrichmentFeature && tree[i]->isUpdated )
                {
                    needMeshing = true ;
                    reuseDisplacements = false ;
                }
            }

        }
        else if ( is3D() )
        {

            for ( size_t i = 0 ; i < featuresToStep.size() ; i++ )
            {
                featuresToStep[i]->step( deltaTime, &K->getForces(), dtree3D ) ;
                bool moved = featuresToStep[i]->moved() ;
                enrichmentChange = enrichmentChange || moved;

                if ( moved )
                {
                    reuseDisplacements = false ;
                }
            }
            for ( size_t i = 0 ; i < enrichmentManagers.size() ; i++ )
            {
                bool moved = enrichmentManagers[i]->step(deltaTime, &K->getForces(), dtree3D) ;
                if ( moved )
                {
                    reuseDisplacements = false ;
                }
            }
            for ( size_t i = 0 ; i < tree.size() ; i++ )
            {
                if ( !tree[i]->isEnrichmentFeature && tree[i]->isUpdated )
                {
                    needMeshing = true ;
                    reuseDisplacements = false ;
                }
            }
        }
    }

    if ( enrichmentChange )
    {
        residualError = 1e9 ;
    }
}

bool sortByScore ( DelaunayTriangle * tri1, DelaunayTriangle * tri2 )
{
    if ( tri1->getBehaviour()->getFractureCriterion() && tri2->getBehaviour()->getFractureCriterion() )
    {
        return tri1->getBehaviour()->getFractureCriterion()->getScoreAtState() > tri2->getBehaviour()->getFractureCriterion()->getScoreAtState();
    }
    return false ;
}

void FeatureTree::stepMesh()
{

    if ( is2D() )
    {
        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            j->second->step(deltaTime) ;
    }
    else if(is3D())
    {
        if(dtree3D)
            dtree3D->step(deltaTime);
    }
}

bool FeatureTree::stepElements()
{
    behaviourChange = false ;
    stateConverged = false ;
    maxScore = -1 ;
    double maxTolerance = 1e-6 ;
    foundCheckPoint = true ;
    if ( resetcalcul )
    {
        deltaTime = 0 ;
        now = 0 ;
        K->getDisplacements() = 0 ;
        K->getPreviousDisplacements() = 0 ;
        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
        {
            #pragma omp parallel
            {
                #pragma omp single
                {
                    for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                    {
                        #pragma omp task firstprivate(i)
                        {
                            if ( i->getBehaviour()->getDamageModel() )
                            {
                                i->getBehaviour()->getDamageModel()->getState(true) = 0 ;
                            }
                        }
                    }
                }
            }
        }
        resetcalcul = false ;
    }

    if ( solverConvergence )
    {
        if ( is2D() )
        {
            if(cachedVolumes.empty())
            {
                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                    cachedVolumes.push_back(std::vector<double>());
            }
            int lcounter = 0 ;
            double volume = 0 ;

            for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            {
                if ( cachedVolumes[lcounter].empty() )
                {
                    for ( auto i = j->second->begin() ; i != j->second->end() ; i++ )
                    {
                        if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                        {
                            cachedVolumes[lcounter].push_back ( i->area() ) ;
                        }
                        else
                        {
                            cachedVolumes[lcounter].push_back ( 0. ) ;
                        }
                    }
                }
                volume += std::accumulate ( cachedVolumes[lcounter].begin(), cachedVolumes[lcounter].end(), double ( 0 ) ) ;
                lcounter++ ;
            }

            double adamage = 0 ;
            if ( !elastic )
            {
                crackedVolume = 0 ;
                damagedVolume = 0 ;
                averageDamage = 0. ;
            }
            //this will update the state of all elements. This is necessary as
            //the behaviour updates might depend on the global state of the
            //simulation.

            for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            {
                std::cerr << " stepping through elements... " << std::flush ;
                #pragma omp parallel
                {
                    #pragma omp single
                    {
                        for (  auto i = j->second->begin() ; i != j->second->end() ; i++ )
                        {
                            #pragma omp task firstprivate(i)
                            {
                                if ( i.getPosition() % 1000 == 0 )
                                {
                                    std::cerr << "\r stepping through elements... " << i.getPosition() << "/" << i.size() << std::flush ;
                                }
                                i->step ( deltaTime, &K->getDisplacements() ) ;


                            }
                        }
                    }
                }
                std::cerr << " ...done" << std::endl ;
            }



            int fracturedCount = 0 ;
            int ccount = 0 ;


            if ( !elastic )
            {
                double maxScoreInit = -1;
                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                            {
                                #pragma omp task firstprivate(i)
                                {
                                    if ( i.getPosition() % 200 == 0 )
                                        std::cerr << "\r checking for fractures (1)... " << i.getPosition() << "/" << i.size() << std::flush ;
                                    if ( i->getBehaviour()->getFractureCriterion() )
                                    {
                                        i->getBehaviour()->getFractureCriterion()->step ( i->getState() ) ;
                                        #pragma omp flush(maxScoreInit)
                                        double tmpmax = std::max ( i->getBehaviour()->getFractureCriterion()->getScoreAtState(), maxScoreInit ) ;
                                        #pragma omp atomic write
                                        maxScoreInit = tmpmax ;
                                    }
                                }
                            }
                            std::cerr << ". Maxscore = " << maxScoreInit <<" ...done. " << std::endl ;
                        }
                    }
                }



// 				std::stable_sort(elements.begin(), elements.end(), sortByScore) ;

                lcounter = 0 ;
                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                            {
                                #pragma omp task firstprivate(i)
                                {
                                    double are = cachedVolumes[lcounter][i.getPosition()] ;

                                    if ( i.getPosition() % 10000 == 0 )
                                    {
                                        std::cerr << "\r checking for fractures (2)... " << i.getPosition() << "/" << i.size() << std::flush ;
                                    }

                                    if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
                                    {
                                        DamageModel * dmodel = i->getBehaviour()->getDamageModel() ;
                                        bool wasFractured = i->getBehaviour()->fractured() ;
                                        i->getBehaviour()->step ( deltaTime, i->getState(), maxScoreInit ) ;

                                        if ( dmodel )
                                        {
                                            if ( !i->getBehaviour()->fractured() )
                                            {
                                                #pragma omp atomic update
                                                adamage += are  * dmodel->getState().max() ;
                                            }
                                            else
                                            {
                                                #pragma omp atomic update
                                                adamage += are ;    // * dmodel->getState().max() ;
                                            }

                                        }

                                        if ( i->getBehaviour()->changed() )
                                        {
                                            #pragma omp atomic write
                                            behaviourChange = true ;
                                            #pragma omp atomic update
                                            ccount++ ;
                                        }

                                        if ( i->getBehaviour()->fractured() )
                                        {
                                            #pragma omp atomic update
                                            fracturedCount++ ;
                                            #pragma omp atomic update
                                            crackedVolume += are ;

                                            if ( !wasFractured )
                                            {
                                                #pragma omp atomic write
                                                behaviourChange = true ;
                                            }
                                        }
                                        else if ( dmodel && dmodel->getState().max() > POINT_TOLERANCE )
                                        {
                                            #pragma omp atomic update
                                            damagedVolume += are ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    lcounter++ ;
                }
                averageDamage = adamage/volume ;

                std::cerr << ". Average damage = " << averageDamage << " ...done. " << ccount << " elements changed."  << std::endl ;

                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                            {
                                #pragma omp task firstprivate(i)
                                {
                                    if ( i.getPosition() % 1000 == 0 )
                                    {
                                        std::cerr << "\r checking for fractures (3)... " << i.getPosition() << "/" << i.size() << std::flush ;
                                    }

                                    if ( i->getBehaviour()->getDamageModel() )
                                    {
                                        i->getBehaviour()->getDamageModel()->postProcess() ;
                                        if ( i->getBehaviour()->changed() )
                                        {
                                            #pragma omp atomic write
                                            behaviourChange = true ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                std::cerr << ". Average damage = " << averageDamage << " ...done. " << ccount << " elements changed." << std::endl ;

                foundCheckPoint = true ;
                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                            {
                                bool done = false ;
                                #pragma omp task firstprivate(i)
                                if (i->getBehaviour()->getDamageModel() && !i->getBehaviour()->getDamageModel()->converged )
                                {
                                    #pragma omp critical
                                    {
                                        foundCheckPoint = false ;
                                        done = true ;
                                    }
                                }
                                if(done)
                                {
                                    break ;
                                }
                            }
                        }
                    }
                }

                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    #pragma omp parallel
                    {
                        #pragma omp single
                        {
                            if ( !behaviourChange )
                            {
                                for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                                {
                                    bool done = false ;
                                    #pragma omp task firstprivate(i)
                                    if ( i->getBehaviour()->getFractureCriterion() && i->getBehaviour()->getFractureCriterion()->met() )
                                    {
                                        #pragma omp critical
                                        {
                                            behaviourChange = true ;
                                            done = true ;
                                        }
                                    }
                                    if(done)
                                    {
                                        break ;
                                    }
                                }
                            }
                        }
                    }
                }

                maxScore = -1. ;

                if (foundCheckPoint )
                {
                    std::cerr << "[" << averageDamage << " ; " << ccount << " ; " <<  std::flush ;

                    maxTolerance = 1 ;
                    for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                    {
                        #pragma omp parallel
                        {
                            #pragma omp single
                            {
                                for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                                {
                                    #pragma omp task firstprivate(i)
                                    if ( i->getBehaviour()->getFractureCriterion() )
                                    {
                                        #pragma omp critical
                                        {
                                            maxScore = std::max (i->getBehaviour()->getFractureCriterion()->getScoreAtState(), maxScore ) ;
                                            maxTolerance = std::max (i->getBehaviour()->getFractureCriterion()->getScoreTolerance(), maxTolerance ) ;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    std::cerr << maxScore << "]" << std::flush ;
                    for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                    {
                        if ( j->second->begin()->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0 && maxScore < 1.-POINT_TOLERANCE )
                        {
                            std::cerr << "adjusting time step..." << std::endl ;
                            double begin = j->second->cbegin()->getBoundingPoint ( 0 ).getT() ;
                            double end = j->second->cbegin()->getBoundingPoint ( j->second->cbegin()->getBoundingPoints().size() -1 ).getT() ;
                            if ( maxScore* ( end-begin ) > minDeltaTime )
                            {
                                moveFirstTimePlanes ( ( 1.-maxScore ) * ( end-begin ) , j->second->begin(), j->second->end() ) ;
                            }
                            else if ( end - begin > minDeltaTime )
                            {
                                moveFirstTimePlanes ( end-begin-minDeltaTime , j->second->begin(), j->second->end() ) ;
                            }
                            else
                            {
                                std::cerr << "time step too small: setting to " << minDeltaTime << std::endl ;
                                moveFirstTimePlanes ( 0., j->second->begin(), j->second->end() ) ;
                            }
                        }

                        if ( maxScore > 1.-POINT_TOLERANCE || maxScore < POINT_TOLERANCE )
                        {
                            moveFirstTimePlanes ( 0., j->second->begin(), j->second->end() ) ;
                        }
                    }
                }

                std::cerr << " ...done. " << std::endl ;
            }
        }
        else if ( is3D() )
        {

            if ( cachedVolumes.empty() )
            {
                for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
                {
                    cachedVolumes.push_back(std::vector<double>());
                    if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                    {
                        cachedVolumes[0].push_back ( i->volume() ) ;
                    }
                    else
                    {
                        cachedVolumes[0].push_back ( 0. ) ;
                    }
                }
            }
            double volume = std::accumulate ( cachedVolumes[0].begin(), cachedVolumes[0].end(), double ( 0 ) ) ;
            double adamage = 0 ;
            if ( !elastic )
            {
                crackedVolume = 0 ;
                damagedVolume = 0 ;
                averageDamage = 0. ;
            }
            //this will update the state of all elements. This is necessary as
            //the behaviour updates might depend on the global state of the
            //simulation.
            std::cerr << " stepping through elements... " << std::flush ;
            #pragma omp parallel
            {
                #pragma omp single
                for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
                {
                    #pragma omp task firstprivate(i)
                    {
                        if ( i.getPosition() % 1000 == 0 )
                        {
                            std::cerr << "\r stepping through elements... " << i.getPosition() << "/" << i.size() << std::flush ;
                        }

                        i->step ( deltaTime, &K->getDisplacements() ) ;
                    }
                }
            }

            std::cerr << " ...done" << std::endl ;

            int fracturedCount = 0 ;
            int ccount = 0 ;

            if ( !elastic )
            {
                double maxScoreInit = -1;
                #pragma omp parallel
                {
                    #pragma omp single
                    for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
                    {
                        #pragma omp task firstprivate(i)
                        {
                            if ( i.getPosition() % 200 == 0 )
                            {
                                std::cerr << "\r checking for fractures (1)... " << i.getPosition() << "/" << i.size() << std::flush ;
                            }

                            if ( i->getBehaviour()->getFractureCriterion() )
                            {
                                i->getBehaviour()->getFractureCriterion()->step ( i->getState() ) ;
                                #pragma omp flush(maxScoreInit)
                                double tmpmax = std::max ( i->getBehaviour()->getFractureCriterion()->getScoreAtState(), maxScoreInit ) ;
                                #pragma omp atomic write
                                maxScoreInit = tmpmax ;

                            }
                        }
                    }
                }

                std::cerr << " ...done. " << std::endl ;

// 				std::stable_sort(elements.begin(), elements.end(), sortByScore) ;

                #pragma omp parallel
                {
                    #pragma omp single
                    for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
                    {
                        #pragma omp task firstprivate(i)
                        {

                            double are = cachedVolumes[0][i.getPosition()] ;

                            if ( i.getPosition() % 10000 == 0 )
                            {
                                std::cerr << "\r checking for fractures (2)... " << i.getPosition() << "/" << i.size() << std::flush ;
                            }

                            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
                            {
                                DamageModel * dmodel = i->getBehaviour()->getDamageModel() ;
                                bool wasFractured = i->getBehaviour()->fractured() ;

                                i->getBehaviour()->step ( deltaTime, i->getState(), maxScoreInit ) ;
                                if ( dmodel )
                                {
                                    if ( !i->getBehaviour()->fractured() )
                                    {
                                        adamage += are * dmodel->getState().max() ;
                                    }
                                }
                                #pragma omp critical
                                if ( i->getBehaviour()->changed() )
                                {
                                    behaviourChange = true ;
                                    ccount++ ;
                                }
                                #pragma omp critical
                                if ( i->getBehaviour()->fractured() )
                                {
                                    fracturedCount++ ;
                                    crackedVolume += are ;

                                    if ( !wasFractured )
                                    {
                                        behaviourChange = true ;
                                    }
                                }
                                else if ( dmodel && dmodel->getState().max() > POINT_TOLERANCE )
                                {
                                    damagedVolume += are ;
                                }
                            }
                        }
                    }
                }
                averageDamage = adamage/volume ;

                std::cerr << " ...done. " << ccount << " elements changed." << std::endl ;

                #pragma omp parallel
                {
                    #pragma omp single
                    for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
                    {
                        #pragma omp task firstprivate(i)
                        {
                            if ( i.getPosition() % 1000 == 0 )
                            {
                                std::cerr << "\r checking for fractures (3)... " << i.getPosition() << "/" << i.size() << std::flush ;
                            }

                            if ( i->getBehaviour()->getDamageModel() )
                            {
                                i->getBehaviour()->getDamageModel()->postProcess() ;
                                if ( i->getBehaviour()->changed() )
                                {
                                    #pragma omp critical
                                    {
                                        behaviourChange = true ;
                                    }
                                }
                            }
                        }
                    }
                }
                foundCheckPoint = true ;
                for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
                {
                    if ( i->getBehaviour()->getDamageModel() && !i->getBehaviour()->getDamageModel()->converged )
                    {
                        foundCheckPoint = false ;
                        maxScore = maxScoreInit ;
                        break ;
                    }
                }

                if ( !behaviourChange )
                {
                    for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
                    {
                        if ( i->getBehaviour()->getFractureCriterion() && i->getBehaviour()->getFractureCriterion()->met() )
                        {
                            behaviourChange = true ;
                            break ;
                        }
                    }
                }
            }

            maxScore = -1. ;

            if ( !elastic && foundCheckPoint )
            {
                std::cerr << "[" << averageDamage << " ; " << std::flush ;

                maxTolerance = 1 ;
// 				double maxs = -1 ;
// 				double maxtol = 1 ;
// 				#pragma omp parallel lastshared(maxs,maxtol)
// 				{
//
// 				  #pragma omp for nowait
                for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
                {
                    if (i->getBehaviour()->getFractureCriterion() )
                    {
//                         i->getBehaviour()->getFractureCriterion()->setCheckpoint ( true ) ;
                        maxScore = std::max ( i->getBehaviour()->getFractureCriterion()->getScoreAtState(), maxScore ) ;
                        maxTolerance = std::max ( i->getBehaviour()->getFractureCriterion()->getScoreTolerance(), maxTolerance ) ;

                    }
                }
// 				  #pragma omp critical
// 				  {
// 				    maxScore = std::max(maxScore, maxs) ;
// 				    maxTolerance = std::max(maxTolerance, maxtol) ;
// 				  }
// 				}

                std::cerr << maxScore << "]" << std::flush ;
                if ( dtree3D->begin()->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0. && maxScore < 1. )
                {
                    std::cerr << "adjusting time step..." << std::endl ;
                    double begin = dtree3D->begin()->getBoundingPoint ( 0 ).getT() ;
                    double end = dtree3D->begin()->getBoundingPoint ( dtree3D->begin()->getBoundingPoints().size() -1 ).getT() ;
                    if ( maxScore * ( end-begin ) > minDeltaTime )
                    {
                        moveFirstTimePlanes ( ( 1.-maxScore ) * ( end-begin ) , dtree3D->begin(), dtree3D->end() ) ;
                    }
                    else if ( end - begin > minDeltaTime )
                    {
                        moveFirstTimePlanes ( end-begin-minDeltaTime , dtree3D->begin(), dtree3D->end() ) ;
                    }
                    else
                    {
                        std::cerr << "time step too small: setting to " << minDeltaTime << std::endl ;
                        this->moveFirstTimePlanes ( 0., dtree3D->begin(), dtree3D->end() ) ;
                    }
                }


            }
            else if ( !elastic )
            {

                if ( dtree3D->begin()->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0 )
                {
                    moveFirstTimePlanes ( 0. , dtree3D->begin(), dtree3D->end() ) ;
                }


//                 #pragma omp parallel for schedule(runtime)
//                 for (  auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
//                 {
//                     if ( i->getBehaviour()->getFractureCriterion() )
//                     {
//                         i->getBehaviour()->getFractureCriterion()->setCheckpoint ( false ) ;
//                     }
//                 }



            }

            std::cerr << " ...done. " << std::endl ;

        }
    }
    else
    {
        Vector dummyx ( 0., K->getDisplacements().size() ) ;

        if ( is2D() )
        {   for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
            {
                crackedVolume = 0 ;
                damagedVolume = 0 ;
                //this will update the state of all elements. This is necessary as
                //the behaviour updates might depend on the global state of the
                //simulation.
                std::cerr << " stepping through elements... " << std::flush ;

                for (auto i = j->second->begin() ; i != j->second->end() ; i++)
                {
                    if ( i.getPosition() % 1000 == 0 )
                    {
                        std::cerr << "\r stepping through elements... " << i.getPosition() << "/" << i.size() << std::flush ;
                    }

                    i->step ( 0., &dummyx ) ;
                }

                std::cerr << " ...done" << std::endl ;

            }
        }
        else if ( is3D() )
        {

            //this will update the state of all elements. This is necessary as
            //the behaviour updates might depend on the global state of the
            //simulation.

            //this will update the state of all elements. This is necessary as
            //the behaviour updates might depend on the global state of the
            //simulation.
            std::cerr << " stepping through elements... " << std::flush ;

            for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
            {
                if ( i.getPosition() % 1000 == 0 )
                {
                    std::cerr << "\r stepping through elements... " << i.getPosition() << "/" << i.size() << std::flush ;
                }

                i->step ( 0., &dummyx ) ;
            }

            std::cerr << " ...done" << std::endl ;

        }
    }

    stateConverged = foundCheckPoint && maxScore < maxTolerance ;

    if ( behaviourChange )
    {
        residualError = 1e9 ;
    }
    return stateConverged;
}


void FeatureTree::State::setStateTo ( StateType s, bool stepChanged )
{
    if( !spaceTimeSwitched )
    {
        ft->checkSpaceTimeConsistency() ;
        spaceTimeSwitched = true ;
    }

    bool behaviourChanged = ft->behaviourChanged() ;
    bool xfemChanged = ft->enrichmentChanged() ;

    bool samplingChanged = ft->needMeshing ;

    if ( samplingChanged )
    {
        sampled = false ;
        meshed = false ;
        behaviourSet = behaviourSet ;
        behaviourUpdated = false;
        stitched = false ;
        renumbered = false ;
        initialised = false ;
        enriched = false ;
        assembled = false ;
        solved = false ;
        behaviourStepped = false;
        xfemStepped = false ;
        featureStepped = false;
    }

    if ( stepChanged )
    {
        if ( ft->deltaTime > POINT_TOLERANCE )
        {
            enriched = false ;
        }

        assembled = false ;
        solved = false ;
        behaviourStepped = false;
        xfemStepped = false ;
        featureStepped = false;
    }

    if ( xfemChanged )
    {
        initialised = false ;
        enriched = false ;
        assembled = false ;
        solved = false ;
        behaviourStepped = false;
        xfemStepped = false ;
        featureStepped = false;
    }

    if ( behaviourChanged )
    {
        assembled = false ;
        solved = false ;
        behaviourStepped = false;
        xfemStepped = false ;
        featureStepped = false;
    }
    ft->stepMesh();

#ifdef HAVE_OMP
    if ( !sampled )
    {
        double t0 = omp_get_wtime() ;
        ft->sample();
        sampled = true ;
        std::cerr << "... time to sample (s) " << omp_get_wtime()-t0 << std::endl ;
    }
#else
    if ( !sampled )
    {
        ft->sample();
        sampled = true ;
    }
#endif

    if ( s == SAMPLED )
    {
        return ;
    }

#ifdef HAVE_OMP
    if ( !meshed )
    {
        double t0 = omp_get_wtime() ;
        ft->generateElements();
        meshed = true ;
        std::cerr << "... time to mesh (s) " << omp_get_wtime()-t0 << std::endl ;
    }
#else
    if ( !meshed )
    {
        ft->generateElements();
        meshed = true ;
    }
#endif

    if ( s == MESHED )
    {
        return ;
    }

    if ( !stitched )
    {
        ft->stitch();
        stitched = true ;
    }
    if ( s == STITCHED )
    {
        return ;
    }

    if ( !behaviourSet )
    {
        ft->setElementBehaviours() ;
        behaviourSet = true;
        behaviourUpdated = true;
        initialised = false ;
    }
    else if ( !behaviourUpdated && behaviourSet )
    {
        ft->updateElementBehaviours();
        behaviourUpdated = true;
    }
    if ( s == BEHAVIOUR_SET )
    {
        return ;
    }

    if ( !renumbered )
    {
        ft->renumber();
        renumbered = true ;
    }
    if ( s == RENUMBERED )
    {
        return ;
    }

    if ( !initialised )
    {
        ft->initializeElements( );
        initialised = true ;
    }
    if ( s == INITIALISED )
    {
        return ;
    }

    if ( !enriched )
    {
        ft->enrich() ;
        enriched = true ;
    }
    if ( s == ENRICHED )
    {
        return ;
    }

    if ( !assembled )
    {
        ft->assemble();
        assembled = true ;
    }
    if ( s == ASSEMBLED )
    {
        return ;
    }

    if ( !solved )
    {
        ft->solve() ;
        solved = true ;
    }
    if ( s == SOLVED )
    {
        return ;
    }

    if ( !behaviourStepped )
    {
        ft->stepElements();
        behaviourStepped = true ;
    }
    if ( s == BEHAVIOUR_STEPPED )
    {
        return ;
    }
    if ( !xfemStepped )
    {
        ft->stepXfem();
        xfemStepped = true ;
    }
}

void FeatureTree::resetBoundaryConditions()
{
    boundaryCondition.clear() ;
    if ( K )
    {
        K->clear() ;
    }
}

void FeatureTree::checkSpaceTimeConsistency()
{
    bool spaceTimeElemOrder = (elemOrder >= CONSTANT_TIME_LINEAR && elemOrder < QUADTREE_REFINED);
    bool spaceTimeElemBehaviour = false ;
    bool isDiffusion = false ;
    int maxBlocks = 1 ;
    for(size_t i = 0 ; i < tree.size() ; i++)
    {
        Form * behaviour = tree[i]->getBehaviour() ;
        if(dynamic_cast<Viscoelasticity*>(behaviour))
        {
            spaceTimeElemBehaviour = true ;
            if( dynamic_cast<Viscoelasticity*>(behaviour)->blocks > maxBlocks )
                maxBlocks = dynamic_cast<Viscoelasticity*>(behaviour)->blocks ;
        }
        if(dynamic_cast<FiniteDifferenceViscoelasticity*>(behaviour))
        {
            ViscoelasticModel model = dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->model ;
            int tensors = dynamic_cast<FiniteDifferenceViscoelasticity *>(behaviour)->tensors.size() ;
            int blocks = 1 ;
            switch(model)
            {
            case PURE_ELASTICITY:
            case PURE_VISCOSITY:
            case KELVIN_VOIGT:
                blocks = 1 ;
                break ;
            case MAXWELL:
                blocks = 2 ;
                break ;
            case BURGER:
                blocks = 3 ;
                break ;
            case GENERALIZED_KELVIN_VOIGT:
            case GENERALIZED_MAXWELL:
            case GENERAL_VISCOELASTICITY:
                blocks = (tensors-1)/2 ;
                break ;
            }
            if(blocks > maxBlocks)
                maxBlocks = blocks ;
        }
        if ( (dynamic_cast<AggregateBehaviour *>(behaviour) || dynamic_cast<PasteBehaviour *>(behaviour) || dynamic_cast<GelBehaviour *>(behaviour)) && maxBlocks < 3)
        {
            maxBlocks = 3 ;
        }
        if( dynamic_cast<ViscoElasticOnlyAggregateBehaviour*>(behaviour) || dynamic_cast<ViscoDamageAggregateBehaviour*>(behaviour) ||
                dynamic_cast<ViscoElasticOnlyPasteBehaviour*>(behaviour) || dynamic_cast<ViscoDamagePasteBehaviour*>(behaviour) ||
                dynamic_cast<ViscoElasticOnlyGelBehaviour*>(behaviour))
        {
            spaceTimeElemBehaviour = true ;
            int fblocks = 3 ;
            if( dynamic_cast<ViscoElasticOnlyAggregateBehaviour*>(behaviour) ) {
                fblocks += dynamic_cast<ViscoElasticOnlyAggregateBehaviour*>(behaviour)->freeblocks ;
            }
            else if( dynamic_cast<ViscoDamageAggregateBehaviour*>(behaviour) ) {
                fblocks += dynamic_cast<ViscoDamageAggregateBehaviour*>(behaviour)->freeblocks ;
            }
            else if( dynamic_cast<ViscoElasticOnlyPasteBehaviour*>(behaviour) ) {
                fblocks += dynamic_cast<ViscoElasticOnlyPasteBehaviour*>(behaviour)->freeblocks ;
            }
            else if( dynamic_cast<ViscoDamagePasteBehaviour*>(behaviour) ) {
                fblocks += dynamic_cast<ViscoDamagePasteBehaviour*>(behaviour)->freeblocks ;
            }
            else if( dynamic_cast<ViscoElasticOnlyGelBehaviour*>(behaviour) ) {
                fblocks += dynamic_cast<ViscoElasticOnlyGelBehaviour*>(behaviour)->freeblocks ;
            }
            if(fblocks > maxBlocks)
                maxBlocks = fblocks ;

        }
        Diffusion * dbehaviour = dynamic_cast<Diffusion *>(tree[i]->getBehaviour()) ;
        if(dbehaviour)
        {
            spaceTimeElemBehaviour = true ;
            isDiffusion = true ;
        }
        if(isDiffusion && !dbehaviour)
        {
            std::cout << "unable to do transport-mechanical analysis" << std::endl ;
            exit(0) ;
        }

    }

    if( spaceTimeElemBehaviour && !spaceTimeElemOrder )
    {
        switch(elemOrder)
        {
        case CONSTANT:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to CONSTANT_TIME_LINEAR element order" << std::endl ;
            setOrder(CONSTANT_TIME_LINEAR) ;
            break ;
        case LINEAR:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to LINEAR_TIME_LINEAR element order" << std::endl ;
            setOrder(LINEAR_TIME_LINEAR) ;
            break ;
        case QUADRATIC:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to QUADRATIC_TIME_LINEAR element order" << std::endl ;
            setOrder(QUADRATIC_TIME_LINEAR) ;
            break ;
        case CUBIC:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to CUBIC_TIME_LINEAR element order" << std::endl ;
            setOrder(CUBIC_TIME_LINEAR) ;
            break ;
        case QUADRIC:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to QUADRIC_TIME_LINEAR element order" << std::endl ;
            setOrder(QUADRIC_TIME_LINEAR) ;
            break ;
        case QUINTIC:
            std::cerr << "element type inappropriate for space-time finite element analysis, switching to QUINTIC_TIME_LINEAR element order" << std::endl ;
            setOrder(QUINTIC_TIME_LINEAR) ;
            break ;
        case QUADTREE_REFINED:
        case REGULAR_GRID:
            std::cerr << "unable to convert " << (elemOrder == QUADTREE_REFINED ? "QUADTREE_REFINED" : "REGULAR_GRID") << " to space-time finite elements" << std::endl ;
            exit(0) ;
            break ;
        default:
            break ;
        }
        spaceTimeElemOrder = true ;
    }

    if(spaceTimeElemOrder && !isDiffusion)
    {
        for(size_t i = 0 ; i < tree.size() ; i++)
        {
            std::cerr << "\rchecking space-time finite element behaviour (feature " << i << "/" << tree.size() << ")" << std::flush ;
            if(! BehaviourConverter::toSpaceTimeBehaviour( tree[i], maxBlocks ) )
            {
                std::cerr << "conversion of behaviour to space-time failed for feature " << i << "/" << tree.size() << std::endl ;
                exit( 0 ) ;
            }

            if(tree[i]->isEnrichmentFeature && !dynamic_cast<TimeDependentEnrichmentInclusion*>(tree[i]) )
            {
                std::cerr << "enrichment inconsistent with space-time finite elements for feature " << i << "/" << tree.size() << std::endl ;
                exit( 0 ) ;
            }


        }
        std::cerr << "\rchecking space-time finite element behaviour (feature " << tree.size() << "/" << tree.size() << ") ...done" << std::endl ;
    }


}

bool FeatureTree::stepInternal(bool guided, bool xfemIteration)
{
    StateType lastStep = XFEM_STEPPED ;

    bool ret = true ;
    size_t totit = 0 ;
    int notConvergedCounts = 0 ;

    if ( stateConverged && state.meshed && solverConverged() && !behaviourChanged() )
    {
        deltaTime = realDeltaTime ;
        now += deltaTime ;
        setDeltaTime( deltaTime, true ) ;
        for ( size_t i = 0 ; i < nodes.size() ; i++ )
        {
            nodes[i]->getT() += deltaTime ;
        }
    }
    else
    {
        deltaTime = 0 ;
    }

    for ( size_t i = 0 ; i < boundaryCondition.size() ; i++ )
    {
        TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *> ( boundaryCondition[i] ) ;
        if ( timec != nullptr )
        {
            timec->goToNext = stateConverged && !behaviourChanged();
            if(!behaviourChanged())
                timec->instant = 1. ;
            else if(!damageConverged)
                timec->instant = 0. ;
        }
    }

    size_t it = 0 ;
    bool needexit = false ;
    do
    {
        if ( it == 1 )
        {
            for ( size_t k = 0 ; k < boundaryCondition.size() ; k++ )
            {
                TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *> ( boundaryCondition[k] ) ;
                if ( timec )
                {
                    timec->goToNext = stateConverged && !behaviourChanged() ;
                    if(!behaviourChanged())
                        timec->instant = 1. ;
                    else if(!damageConverged)
                        timec->instant = 0. ;
                }
            }
        }
        state.setStateTo ( lastStep, true ) ;
        ++it ;
        deltaTime = 0 ;
        if ( solverConverged() )
        {
            std::cout << "." << std::flush ;
            notConvergedCounts = 0 ;

            if(foundCheckPoint  && !(enrichmentChange || behaviourChanged() ) )
                needexit = true ;
        }
        else
        {
            notConvergedCounts++ ;
            if(notConvergedCounts > 5)
                needexit = true ;
            std::cerr << "+" << std::flush ;
        }

        if ( enrichmentChange || needMeshing )
        {
            K->clear() ;
        }


        if ( it > maxitPerStep && foundCheckPoint )
        {
            ret = false ;
            needexit = true ;
        }

/*        std::cout << foundCheckPoint << std::endl ;
        std::cout << enrichmentChange << std::endl ;
        std::cout << behaviourChanged() << std::endl ;
        std::cout << needexit << std::endl ;*/

    }
    while ( !needexit ) ;

    totit += it ;

    if(is2D())
    {
        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
        {
            #pragma omp parallel
            {
                #pragma omp single
                {
                    for (  auto i = j->second->begin() ; i != j->second->end() ; i++  )
                    {
                        #pragma omp task firstprivate(i)
                        if ( i->getBehaviour()->getDamageModel() )
                        {
                            i->getBehaviour()->getFractureCriterion()->setCheckpoint ( true ) ;
                            i->getBehaviour()->getDamageModel()->alternateCheckpoint = false  ;
                            i->getBehaviour()->getDamageModel()->converged = true ;
                            i->getBehaviour()->getDamageModel()->alternate = true ;
                        }
                    }
                }
            }
        }
    }
    else
    {
        #pragma omp parallel
        {
            #pragma omp single
            {
                for (  auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
                {
                    #pragma omp task firstprivate(i)
                    if ( i->getBehaviour()->getDamageModel() )
                    {
                        i->getBehaviour()->getFractureCriterion()->setCheckpoint ( true ) ;
                        i->getBehaviour()->getDamageModel()->alternateCheckpoint = false  ;
                        i->getBehaviour()->getDamageModel()->converged = true ;
                        i->getBehaviour()->getDamageModel()->alternate = true ;
                    }
                }
            }
        }
    }

    if ( notConvergedCounts >= 20 )
    {
        ret = false ;
    }
    std::cerr << std::endl ;

    if ( ret && !guided)
    {
        setDeltaTime ( realDeltaTime, false ) ;
    }
//     std::cout << totit << std::endl ;
    std::cerr << totit << "/" << maxitPerStep << "." << std::flush ;
    damageConverged = solverConverged() && !behaviourChanged() /*stateConverged*/ && ret && ( totit < maxitPerStep ) ;

    if(damageConverged)
        K->setPreviousDisplacements() ;

    return solverConverged() && stateConverged && maxScore < 0 ;

}

bool FeatureTree::step(bool guided)
{
    return stepInternal(guided, true) ;
}

bool FeatureTree::stepToCheckPoint( int iterations, double precision)
{
//     setDeltaTime ( realDeltaTime, false ) ;
    double initialscale = 1. ;
    for ( const auto & bc : boundaryCondition)
    {
        initialscale  =std::min(bc->getScale(), initialscale) ;
    }

    int prevmaxit = maxitPerStep ;
    maxitPerStep = 2 ;
    scaleBoundaryConditions ( 1. );

    for(int iter = 0 ; iter < iterations ; iter++)
    {
//         setDeltaTime ( realDeltaTime, false ) ;
        std::cout <<"["<<iter+1<<"/"<<iterations<< " : "<< maxScore << "]" << std::flush ;
        step(true) ;
        if(maxScore < 0)
            break ;
    }
    maxitPerStep = prevmaxit ;
    //at this point, we should have found a checkpoint.
    if(maxScore > 0)
    {
        for ( size_t k = 0 ; k < boundaryCondition.size() ; k++ )
        {
            TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *> ( boundaryCondition[k] ) ;
            if ( timec )
            {
                timec->goToNext = false ;
                timec->instant = 0. ;
            }
        }

        elastic = true ;
        double currentScale = 0.8 ;
        double highscale = 1. ;
        double bottomscale= 0. ;

        while(highscale-bottomscale > precision*initialscale)
        {
            currentScale = highscale*.8+bottomscale*.2 ;
            scaleBoundaryConditions(currentScale) ;

            bool met = false ;
            state.setStateTo ( XFEM_STEPPED, true ) ;

            if(is2D())
            {
                for ( auto j = layer2d.begin() ; j != layer2d.end() && !met; j++ )
                {

                    for (  auto i = j->second->begin() ; i != j->second->end()  && !met; i++  )
                    {

                        if ( i->getBehaviour()->getFractureCriterion() )
                        {
                            i->getBehaviour()->getFractureCriterion()->step ( i->getState() ) ;
                            if(i->getBehaviour()->getFractureCriterion()->met())
                            {
                                met = true ;
                            }
                        }
                    }
                }
            }
            else
            {

                for (  auto i = dtree3D->begin() ; i != dtree3D->end() && !met; i++  )
                {

                    if ( i->getBehaviour()->getFractureCriterion() )
                    {
                        i->getBehaviour()->getFractureCriterion()->step ( i->getState() ) ;
                        if(i->getBehaviour()->getFractureCriterion()->met())
                        {
                            met = true ;
                        }
                    }
                }
            }

            if(met)
                highscale = currentScale ;
            else
                bottomscale = currentScale ;
        }
        std::cout << "Scale = " << bottomscale << std::endl ;
        scaleBoundaryConditions(bottomscale) ;
        elastic = false ;

        state.setStateTo ( XFEM_STEPPED, true ) ;
        if(damageConverged)
            K->setPreviousDisplacements() ;
    }
    else
    {
        if(damageConverged)
            K->setPreviousDisplacements() ;
    }

    return true ;
}

bool orderPointsByID ( Point * p1, Point * p2 )
{
    return p1->getId() < p2->getId() ;
}

std::vector<Point *> FeatureTree::getNodes ()
{
    if ( nodes.size() > 0 )
    {
        return nodes ;
    }

    std::vector<Point *> pts ;
    if ( is2D() )
    {
        std::valarray<bool> done ( dtree->begin().size() *dtree->begin()->getBoundingPoints().size() ) ;
        done = false ;
        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                if ( !done[ i->getBoundingPoint ( j ).getId()] )
                {
                    done[ i->getBoundingPoint ( j ).getId()] = true ;
                    pts.push_back ( &i->getBoundingPoint ( j ) ) ;
                }
            }
        }
    }
    if ( is3D() )
    {
        std::valarray<bool> done ( dtree3D->begin().size() *dtree3D->begin()->getBoundingPoints().size() ) ;
        done = false ;
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
            {
                if ( !done[ i->getBoundingPoint ( j ).getId()] )
                {
                    done[ i->getBoundingPoint ( j ).getId()] = true ;
                    pts.push_back ( &i->getBoundingPoint ( j ) ) ;
                }
            }
        }
    }
    std::sort ( pts.begin(), pts.end(), orderPointsByID ) ;
    return pts ;
}

Vector FeatureTree::getAverageField ( FieldType f, int grid , double t )
{
    if ( is2D() )
    {
        Vector ret = get2DMesh()->getField ( f, -1, t ) ;
        for ( auto layer = layer2d.begin() ; layer!=layer2d.end() ; layer++ )
        {
            if(layer->second != dtree)
                ret += layer->second->getField ( f, -1, t ) ;
        }
        return ret ;
    }


    return get3DMesh()->getField ( f, -1, t ) ;

}

Vector FeatureTree::getAverageFieldOnBoundary ( BoundingBoxPosition edge, FieldType f, int dummy, double t )
{
    if(boundaryCache.find(edge) == boundaryCache.end())
    {
        if(is2D())
        {
            Segment * pos = nullptr ;
            switch(edge)
            {
            case TOP:
                pos = new Segment( tree[0]->getBoundingBox()[0], tree[0]->getBoundingBox()[1] ) ;
                break ;
            case BOTTOM:
                pos = new Segment( tree[0]->getBoundingBox()[2], tree[0]->getBoundingBox()[3] ) ;
                break ;
            case LEFT:
                pos = new Segment( tree[0]->getBoundingBox()[3], tree[0]->getBoundingBox()[0] ) ;
                break ;
            case RIGHT:
                pos = new Segment( tree[0]->getBoundingBox()[1], tree[0]->getBoundingBox()[2] ) ;
                break ;
            default:
                std::cout << "cannot calculate field on selected boundary" << std::endl ;
                return get2DMesh()->getField ( f, -1, t ) ;
            }
            unsigned int id = get2DMesh()->generateCache( pos ) ;
            boundaryCache[ edge ] = std::make_pair( pos, id ) ;
        }
        else
        {
            std::cout << "cannot calculate field on selected boundary" << std::endl ;
            return get3DMesh()->getField ( f, -1, t ) ;
        }
    }


    if ( is2D() )
    {
        return get2DMesh()->getField ( f, boundaryCache[edge].first, boundaryCache[edge].second, dummy, t ) ;
    }

    return get3DMesh()->getField ( f, boundaryCache[edge].first, boundaryCache[edge].second, dummy, t ) ;

}

double FeatureTree::getAverageFieldOnBoundary ( BoundingBoxPosition edge, std::string f, int dummy, double t )
{
    if(boundaryCache.find(edge) == boundaryCache.end())
    {
        if(is2D())
        {
            Segment * pos = nullptr ;
            switch(edge)
            {
            case TOP:
                pos = new Segment( tree[0]->getBoundingBox()[0], tree[0]->getBoundingBox()[1] ) ;
                break ;
            case BOTTOM:
                pos = new Segment( tree[0]->getBoundingBox()[2], tree[0]->getBoundingBox()[3] ) ;
                break ;
            case LEFT:
                pos = new Segment( tree[0]->getBoundingBox()[3], tree[0]->getBoundingBox()[0] ) ;
                break ;
            case RIGHT:
                pos = new Segment( tree[0]->getBoundingBox()[1], tree[0]->getBoundingBox()[2] ) ;
                break ;
            default:
                std::cout << "cannot calculate field on selected boundary" << std::endl ;
                return get2DMesh()->getField ( f, -1 ) ;
            }
            unsigned int id = get2DMesh()->generateCache( pos ) ;
            boundaryCache[ edge ] = std::make_pair( pos, id ) ;
        }
        else
        {
            std::cout << "cannot calculate field on selected boundary" << std::endl ;
            return get3DMesh()->getField ( f, -1 ) ;
        }
    }


    if ( is2D() )
    {
        return get2DMesh()->getField ( f, boundaryCache[edge].first, boundaryCache[edge].second ) ;
    }

    return get3DMesh()->getField ( f, boundaryCache[edge].first, boundaryCache[edge].second ) ;

}

std::vector<double>  FeatureTree::getMacroscopicStrain ( const Geometry * base, double tol )
{
    Vector disps = getDisplacements() ;
    std::vector<double> ret ;
    if ( !disps.size() )
    {
        return ret ;
    }
    if ( is2D() )
    {
        double size_x = dynamic_cast<const Rectangle *> ( base )->width() ;
        double size_y = dynamic_cast<const Rectangle *> ( base )->height() ;

        double dxp = 0;
        double dxpc = 0 ;
        double dyp = 0;
        double dypc = 0 ;
        double dxm = 0;
        double dxmc = 0 ;
        double dym = 0;
        double dymc = 0 ;
        std::set<int> doneIds ;
        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {

                    Point test ( i->getBoundingPoint ( j ) ) ;

                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX() +size_x*.5 ) ) < tol && doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp += disps[2*i->getBoundingPoint ( j ).getId()] ;
                        dxpc++ ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY() +size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dyp += disps[2*i->getBoundingPoint ( j ).getId() +1] ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                        dypc++ ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm += disps[2*i->getBoundingPoint ( j ).getId()] ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                        dxmc++ ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dym += disps[2*i->getBoundingPoint ( j ).getId() +1] ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                        dymc++ ;
                    }
                }
            }
        }

        ret.push_back ( ( dxp/dxpc-dxm/dxmc ) /size_x );
        ret.push_back ( ( dyp/dypc-dym/dymc ) /size_y );
    }
    else
    {
        double size_x = dynamic_cast<const Hexahedron *> ( base )->getXSize() ;
        double size_y = dynamic_cast<const Hexahedron *> ( base )->getYSize() ;
        double size_z = dynamic_cast<const Hexahedron *> ( base )->getZSize() ;

        double dxp = 0;
        double dxpc = 0 ;
        double dyp = 0;
        double dypc = 0 ;
        double dzp = 0;
        double dzpc = 0 ;
        double dxm = 0;
        double dxmc = 0 ;
        double dym = 0;
        double dymc = 0 ;
        double dzm = 0;
        double dzmc = 0 ;
        std::set<int> doneIds ;
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {
                    Point test ( i->getBoundingPoint ( j ) ) ;
                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp = disps[3*i->getBoundingPoint ( j ).getId()] ;
                        dxpc++ ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol&& doneIds.find ( 3*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dyp = disps[3*i->getBoundingPoint ( j ).getId() +1] ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() +1 );
                        dypc++ ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getZ() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol&& doneIds.find ( 3*i->getBoundingPoint ( j ).getId() +2 ) == doneIds.end() )
                    {
                        dzp = disps[3*i->getBoundingPoint ( j ).getId() +2] ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() +2 );
                        dzpc++ ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol&& doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm = disps[3*i->getBoundingPoint ( j ).getId()] ;
                        dxmc++ ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol&& doneIds.find ( 3*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dym = disps[3*i->getBoundingPoint ( j ).getId() +1] ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() +1 );
                        dymc++ ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol&& doneIds.find ( 3*i->getBoundingPoint ( j ).getId() +2 ) == doneIds.end() )
                    {
                        dzm = disps[3*i->getBoundingPoint ( j ).getId() +2] ;
                        doneIds.insert ( 3*i->getBoundingPoint ( j ).getId() +2 );
                        dzmc++ ;
                    }
                }
            }
        }


        ret.push_back ( ( dxp/dxpc-dxm/dxmc ) /size_x );
        ret.push_back ( ( dyp/dypc-dym/dymc ) /size_y );
        ret.push_back ( ( dzp/dzpc-dzm/dzmc ) /size_z );
    }

    return ret ;
}

std::vector<double> FeatureTree::getCutMacroscopicStrain ( const Geometry * base, double tol, double cut )
{
    Vector disps = getDisplacements() ;
    std::vector<double> ret ;
    if ( !disps.size() )
    {
        return ret ;
    }
    if ( is2D() )
    {
        double size_x = dynamic_cast<const Rectangle *> ( base )->width() ;
        double size_y = dynamic_cast<const Rectangle *> ( base )->height() ;
        tol = .05*std::min ( size_x, size_y ) ;
        std::vector<double> dxp ;
        std::vector<double> dyp ;
        std::vector<double> dxm ;
        std::vector<double> dym ;
        std::set<int> doneIds ;
        for (auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {

                    Point test ( i->getBoundingPoint ( j ) ) ;

                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX() +size_x*.5 ) ) < tol && doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp.push_back ( disps[2*i->getBoundingPoint ( j ).getId()] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY() +size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dyp.push_back ( disps[2*i->getBoundingPoint ( j ).getId() +1] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm.push_back ( disps[2*i->getBoundingPoint ( j ).getId()] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dym.push_back ( disps[2*i->getBoundingPoint ( j ).getId() +1] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                    }
                }
            }
        }
        std::sort ( dxp.begin(), dxp.end() ) ;
        std::sort ( dxm.begin(), dxm.end() ) ;
        std::sort ( dyp.begin(), dyp.end() ) ;
        std::sort ( dym.begin(), dym.end() ) ;
        double xpmax = dxp[dxp.size() *cut];
        double xpmin = dxp[dxp.size() * ( 1.-cut )];
        double xmmax = dxm[dxm.size() *cut];
        double xmmin = dxm[dxm.size() * ( 1.-cut )];
        double ypmax = dyp[dyp.size() *cut];
        double ypmin = dyp[dyp.size() * ( 1.-cut )];
        double ymmax = dym[dym.size() *cut];
        double ymmin = dym[dym.size() * ( 1.-cut )];

        double dxpv = 0;
        double dxpc = 0;
        double dxmv = 0;
        double dxmc = 0;
        double dypv = 0;
        double dypc = 0;
        double dymv = 0;
        double dymc = 0;

        for ( size_t i = 0 ; i < dxp.size() ; i++ )
        {
            if ( dxp[i] > xpmin && dxp[i] < xpmax )
            {
                dxpv += dxp[i] ;
                dxpc++ ;
            }
        }
        for ( size_t i = 0 ; i < dxm.size() ; i++ )
        {
            if ( dxm[i] > xmmin && dxm[i] < xmmax )
            {
                dxmv += dxm[i] ;
                dxmc++ ;
            }
        }

        for ( size_t i = 0 ; i < dyp.size() ; i++ )
        {
            if ( dyp[i] > ypmin && dyp[i] < ypmax )
            {
                dypv += dyp[i] ;
                dypc++ ;
            }
        }
        for ( size_t i = 0 ; i < dym.size() ; i++ )
        {
            if ( dym[i] > ymmin && dym[i] < ymmax )
            {
                dymv += dym[i] ;
                dymc++ ;
            }
        }

        ret.push_back ( ( dxpv/dxpc-dxmv/dxmc ) /size_x );
        ret.push_back ( ( dypv/dypc-dymv/dymc ) /size_y );
    }
    else
    {
        double size_x = dynamic_cast<const Hexahedron *> ( base )->getXSize() ;
        double size_y = dynamic_cast<const Hexahedron *> ( base )->getYSize() ;
        double size_z = dynamic_cast<const Hexahedron *> ( base )->getZSize() ;
        std::vector<double> dxp ;
        std::vector<double> dyp ;
        std::vector<double> dzp ;
        std::vector<double> dxm ;
        std::vector<double> dym ;
        std::vector<double> dzm ;
        std::set<int> doneIds ;
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {
                    Point test ( i->getBoundingPoint ( j ) ) ;
                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp.push_back ( disps[3*i->getBoundingPoint ( j ).getId()] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dyp.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +1] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getZ() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dzp.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +2] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm.push_back ( disps[3*i->getBoundingPoint ( j ).getId()] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dym.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +1] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dzm.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +2] ) ;
                    }
                }
            }
        }

        std::sort ( dxp.begin(), dxp.end() ) ;
        std::sort ( dxm.begin(), dxm.end() ) ;
        std::sort ( dyp.begin(), dyp.end() ) ;
        std::sort ( dym.begin(), dym.end() ) ;
        std::sort ( dzp.begin(), dzp.end() ) ;
        std::sort ( dzm.begin(), dzm.end() ) ;

        double xpmax = dxp[dxp.size() *cut];
        double xpmin = dxp[dxp.size() * ( 1.-cut )];
        double xmmax = dxm[dxm.size() *cut];
        double xmmin = dxm[dxm.size() * ( 1.-cut )];
        double ypmax = dyp[dyp.size() *cut];
        double ypmin = dyp[dyp.size() * ( 1.-cut )];
        double ymmax = dym[dym.size() *cut];
        double ymmin = dym[dym.size() * ( 1.-cut )];
        double zpmax = dzp[dzp.size() *cut];
        double zpmin = dzp[dzp.size() * ( 1.-cut )];
        double zmmax = dzm[dzm.size() *cut];
        double zmmin = dzm[dzm.size() * ( 1.-cut )];

        double dxpv = 0;
        double dxpc = 0;
        double dxmv = 0;
        double dxmc = 0;
        double dypv = 0;
        double dypc = 0;
        double dymv = 0;
        double dymc = 0;
        double dzpv = 0;
        double dzpc = 0;
        double dzmv = 0;
        double dzmc = 0;

        for ( size_t i = 0 ; i < dxp.size() ; i++ )
        {
            if ( dxp[i] > xpmin && dxp[i] < xpmax )
            {
                dxpv += dxp[i] ;
                dxpc++ ;
            }
        }
        for ( size_t i = 0 ; i < dxm.size() ; i++ )
        {
            if ( dxm[i] > xmmin && dxm[i] < xmmax )
            {
                dxmv += dxm[i] ;
                dxmc++ ;
            }
        }

        for ( size_t i = 0 ; i < dyp.size() ; i++ )
        {
            if ( dyp[i] > ypmin && dyp[i] < ypmax )
            {
                dypv += dyp[i] ;
                dypc++ ;
            }
        }
        for ( size_t i = 0 ; i < dym.size() ; i++ )
        {
            if ( dym[i] > ymmin && dym[i] < ymmax )
            {
                dymv += dym[i] ;
                dymc++ ;
            }
        }

        for ( size_t i = 0 ; i < dzp.size() ; i++ )
        {
            if ( dzp[i] > zpmin && dzp[i] < zpmax )
            {
                dzpv += dzp[i] ;
                dzpc++ ;
            }
        }
        for ( size_t i = 0 ; i < dym.size() ; i++ )
        {
            if ( dzm[i] > zmmin && dym[i] < zmmax )
            {
                dzmv += dzm[i] ;
                dzmc++ ;
            }
        }


        ret.push_back ( ( dxpv/dxpc-dxmv/dxmc ) /size_x );
        ret.push_back ( ( dypv/dypc-dymv/dymc ) /size_y );
        ret.push_back ( ( dzpv/dzpc-dzmv/dzmc ) /size_z );
    }

    return ret ;
}

std::vector<double> FeatureTree::getMedianMacroscopicStrain ( const Geometry * base, double tol )
{
    Vector disps = getDisplacements() ;
    std::vector<double> ret ;
    if ( !disps.size() )
    {
        return ret ;
    }
    if ( is2D() )
    {
        double size_x = dynamic_cast<const Rectangle *> ( base )->width() ;
        double size_y = dynamic_cast<const Rectangle *> ( base )->height() ;

        tol = .05*std::min ( size_x, size_y ) ;
        std::vector<double> dxp ;
        std::vector<double> dyp ;
        std::vector<double> dxm ;
        std::vector<double> dym ;
        std::set<int> doneIds ;
        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {

                    Point test ( i->getBoundingPoint ( j ) ) ;

                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX() +size_x*.5 ) ) < tol && doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp.push_back ( disps[2*i->getBoundingPoint ( j ).getId()] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY() +size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dyp.push_back ( disps[2*i->getBoundingPoint ( j ).getId() +1] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm.push_back ( disps[2*i->getBoundingPoint ( j ).getId()] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() );
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol&& doneIds.find ( 2*i->getBoundingPoint ( j ).getId() +1 ) == doneIds.end() )
                    {
                        dym.push_back ( disps[2*i->getBoundingPoint ( j ).getId() +1] ) ;
                        doneIds.insert ( 2*i->getBoundingPoint ( j ).getId() +1 );
                    }
                }
            }
        }
        std::sort ( dxp.begin(), dxp.end() ) ;
        std::sort ( dxm.begin(), dxm.end() ) ;
        std::sort ( dyp.begin(), dyp.end() ) ;
        std::sort ( dym.begin(), dym.end() ) ;

        double dxpv = dxp[dxp.size() /2];
        double dxpc = 1;
        double dxmv = dxm[dxm.size() /2];
        double dxmc = 1;
        double dypv = dyp[dyp.size() /2];
        double dypc = 1;
        double dymv = dym[dym.size() /2] ;
        double dymc = 1;


        ret.push_back ( ( dxpv/dxpc-dxmv/dxmc ) /size_x );
        ret.push_back ( ( dypv/dypc-dymv/dymc ) /size_y );
    }
    else
    {
        double size_x = dynamic_cast<const Hexahedron *> ( base )->getXSize() ;
        double size_y = dynamic_cast<const Hexahedron *> ( base )->getYSize() ;
        double size_z = dynamic_cast<const Hexahedron *> ( base )->getZSize() ;

        std::vector<double> dxp ;
        std::vector<double> dyp ;
        std::vector<double> dzp ;
        std::vector<double> dxm ;
        std::vector<double> dym ;
        std::vector<double> dzm ;
        std::set<int> doneIds ;
        for ( auto i = dtree3D->begin() ; i  != dtree3D->end() ; i++ )
        {
            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                for ( size_t j = 0 ; j < i->getBoundingPoints().size() ; j++ )
                {
                    Point test ( i->getBoundingPoint ( j ) ) ;
                    base->project ( &test );
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxp.push_back ( disps[3*i->getBoundingPoint ( j ).getId()] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dyp.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +1] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getZ() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dzp.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +2] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getX() - ( base->getCenter().getX()-size_x*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dxm.push_back ( disps[3*i->getBoundingPoint ( j ).getId()] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getY()-size_y*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dym.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +1] ) ;
                    }
                    if ( dist ( test, i->getBoundingPoint ( j ) ) < tol && std::abs ( i->getBoundingPoint ( j ).getY() - ( base->getCenter().getZ()-size_z*.5 ) ) < tol && doneIds.find ( 3*i->getBoundingPoint ( j ).getId() ) == doneIds.end() )
                    {
                        dzm.push_back ( disps[3*i->getBoundingPoint ( j ).getId() +2] ) ;
                    }
                }
            }
        }

        std::sort ( dxp.begin(), dxp.end() ) ;
        std::sort ( dxm.begin(), dxm.end() ) ;
        std::sort ( dyp.begin(), dyp.end() ) ;
        std::sort ( dym.begin(), dym.end() ) ;
        std::sort ( dzp.begin(), dzp.end() ) ;
        std::sort ( dzm.begin(), dzm.end() ) ;


        double dxpv = dxp[dxp.size() /2];
        double dxpc = 1;
        double dxmv = dxm[dxm.size() /2];
        double dxmc = 1;
        double dypv = dyp[dyp.size() /2];
        double dypc = 1;
        double dymv = dym[dym.size() /2] ;
        double dymc = 1;
        double dzpv = dzp[dzp.size() /2];
        double dzpc = 1;
        double dzmv = dzm[dzm.size() /2] ;
        double dzmc = 1;

        ret.push_back ( ( dxpv/dxpc-dxmv/dxmc ) /size_x );
        ret.push_back ( ( dypv/dypc-dymv/dymc ) /size_y );
        ret.push_back ( ( dzpv/dzpc-dzmv/dzmc ) /size_z );
    }

    return ret ;
}


// Vector FeatureTree::getAverageField( FieldType f, const std::vector<DelaunayTriangle *> & tri )
// {
//     Vector avg ;
//     Vector buffer ;
//     double volume = 0 ;
//     avg.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ;
//     buffer.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ;
//     avg = 0 ;
//     buffer = 0 ;
//     for(size_t i = 0 ; i < tri.size() ; i++)
//     {
//         tri[i]->getState().getAverageField( f, buffer ) ;
//         avg += buffer * tri[i]->area() ;
//         volume += tri[i]->area() ;
//     }
//     return avg/volume ;
// }

// Vector FeatureTree::getAverageField( FieldType f, const std::vector<DelaunayTetrahedron *> & tet )
// {
//     Vector avg ;
//     Vector buffer ;
//     double volume = 0 ;
//     avg.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ;
//     buffer.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ;
//     avg = 0 ;
//     buffer = 0 ;
//     for(size_t i = 0 ; i < tet.size() ; i++)
//     {
//         tet[i]->getState().getAverageField( f, buffer ) ;
//         avg += buffer * tet[i]->volume() ;
//         volume += tet[i]->volume() ;
//     }
//     return avg/volume ;
// }

std::pair<Vector, Vector> FeatureTree::getFieldMinMax ( FieldType f, int grid, double t )
{
    Vector min ;
    Vector max ;
    Vector buffer ;
    VirtualMachine vm ;
    if ( is2D() )
    {
        for ( auto layer = layer2d.begin() ; layer!=layer2d.end() ; layer++ )
        {
            size_t blocks = layer->second->begin()->getBehaviour()->getNumberOfDegreesOfFreedom() /2 ;
            min.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL, blocks ) ) ;
            max.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL, blocks ) ) ;
            buffer.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL, blocks ) ) ;
            buffer = 0 ;
            auto start = layer->second->begin() ;
            while ( start->getBehaviour()->type == VOID_BEHAVIOUR && start != layer->second->end() )
            {
                start++ ;
            }
            start->getState().getAverageField ( f, buffer ) ;
            min = buffer ;
            max = buffer ;
            start++ ;
            for (  ; start != layer->second->end() ; start++ )
            {
                if ( start->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    start->getState().getAverageField ( f, buffer,&vm, -1, t ) ;
                    for ( size_t j = 0 ; j < min.size() ; j++ )
                    {
                        min[j] = std::min ( min[j], buffer[j] ) ;
                        max[j] = std::max ( max[j], buffer[j] ) ;
                    }
                }
            }
        }
    }
    else
    {
        size_t blocks = dtree3D->begin()->getBehaviour()->getNumberOfDegreesOfFreedom() /3 ;
        min.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL, blocks ) ) ;
        max.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL, blocks ) ) ;
        buffer.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL, blocks ) ) ;
        buffer = 0 ;
        auto start = dtree3D->begin() ;
        while ( start->getBehaviour()->type == VOID_BEHAVIOUR && start != dtree3D->end() )
        {
            start++ ;
        }
        start->getState().getAverageField ( f, buffer ) ;
        min = buffer ;
        max = buffer ;
        for ( ; start != dtree3D->end() ; start++ )
        {
            if ( start->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                start->getState().getAverageField ( f, buffer,&vm, -1, t ) ;
                for ( size_t j = 0 ; j < min.size() ; j++ )
                {
                    min[j] = std::min ( min[j], buffer[j] ) ;
                    max[j] = std::max ( max[j], buffer[j] ) ;
                }
            }
        }
    }
    return std::make_pair ( min, max ) ;
}


std::pair<Vector, Vector> FeatureTree::getFieldMinMax ( FieldType f, const std::vector<DelaunayTriangle *> & tri )
{
    Vector min ;
    Vector max ;
    Vector buffer ;
    min.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL ) ) ;
    max.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL ) ) ;
    buffer.resize ( fieldTypeElementarySize ( f, SPACE_TWO_DIMENSIONAL ) ) ;
    tri[0]->getState().getAverageField ( f, buffer ) ;
    min = buffer ;
    max = buffer ;
    for ( size_t i = 1 ; i < tri.size() ; i++ )
    {
        tri[i]->getState().getAverageField ( f, buffer ) ;
        for ( size_t j = 0 ; j < min.size() ; j++ )
        {
            min[j] = std::min ( min[j], buffer[j] ) ;
            max[j] = std::max ( max[j], buffer[j] ) ;
        }
    }
    return std::make_pair ( min, max ) ;
}

std::pair<Vector, Vector> FeatureTree::getFieldMinMax ( FieldType f, const std::vector<DelaunayTetrahedron *> & tet )
{
    Vector min ;
    Vector max ;
    Vector buffer ;
    min.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL ) ) ;
    max.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL ) ) ;
    buffer.resize ( fieldTypeElementarySize ( f, SPACE_THREE_DIMENSIONAL ) ) ;
    tet[0]->getState().getAverageField ( f, buffer ) ;
    min = buffer ;
    max = buffer ;
    for ( size_t i = 1 ; i < tet.size() ; i++ )
    {
        tet[i]->getState().getAverageField ( f, buffer ) ;
        for ( size_t j = 0 ; j < min.size() ; j++ )
        {
            min[j] = std::min ( min[j], buffer[j] ) ;
            max[j] = std::max ( max[j], buffer[j] ) ;
        }
    }
    return std::make_pair ( min, max ) ;
}

bool FeatureTree::isStable()
{
    bool meshChangeinit = behaviourChange ;
    bool enrichmentChangeinit = enrichmentChange ;
    double crackedVolumeinit = crackedVolume ;
    double damagedVolumeinit = damagedVolume ;
    setMaxIterationsPerStep ( 0 ) ;
    elasticStep();
    bool stable = true ;


    if ( is2D() )
    {
        for ( auto j = layer2d.begin() ; j!=layer2d.end() ; ++j )
        {

            for ( auto i = j->second->begin() ; i != j->second->end() && stable; i++ )
            {
                if ( i->getBehaviour() && i->getBehaviour()->getFractureCriterion() && !i->getBehaviour()->fractured() )
                {
                    if ( i->getBehaviour()->getFractureCriterion()->grade ( i->getState() ) > 0 )
                    {
                        stable = false ;
                        break ;
                    }
                }
            }

        }
    }
    else
    {
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() && stable; i++ )
        {
            if ( i->getBehaviour() && i->getBehaviour()->getFractureCriterion() )
            {
                if ( i->getBehaviour()->getFractureCriterion()->grade ( i->getState() ) > 0 )
                {
                    stable = false ;
                    break ;
                }
            }
        }
    }

    setMaxIterationsPerStep ( maxitPerStep ) ;
    behaviourChange = meshChangeinit ;
    enrichmentChange = enrichmentChangeinit ;
    crackedVolume = crackedVolumeinit ;
    damagedVolume = damagedVolumeinit ;

    return stable ;
}

double FeatureTree::getMaximumDisplacement()
{
    state.setStateTo ( RENUMBERED, false ) ;

    if ( is2D() )
    {

        double max = 0 ;

        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( ! i->getBehaviour() )
            {
                continue ;
            }

            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                max = std::max ( max, i->getState().getDisplacements().max() ) ;
            }
        }

        return max ;
    }
    else if ( is3D() )
    {

        double max = 0 ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( ! i->getBehaviour() )
            {
                continue ;
            }

            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                max = std::max ( max, i->getState().getDisplacements().max() ) ;
            }
        }

        return max ;
    }

    return 0 ;
}

double FeatureTree::getMinimumDisplacement()
{
    state.setStateTo ( RENUMBERED, false ) ;

    if ( is2D() )
    {

        double max = 0 ;

        for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
        {
            if ( ! i->getBehaviour() )
            {
                continue ;
            }
            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                max = std::min ( max, i->getState().getDisplacements().min() ) ;
            }
        }

        return max ;
    }
    else if ( is3D() )
    {

        double max = 0 ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( ! i->getBehaviour() )
            {
                continue ;
            }
            if ( i->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                max = std::min ( max, i->getState().getDisplacements().min() ) ;
            }
        }

        return max ;
    }

    return 0 ;
}

size_t FeatureTree::numPoints() const
{
    return this->numdofs ;
}

void FeatureTree::print() const
{
    printForFeature ( tree[0] );

}

void FeatureTree::printForFeature ( const Feature *f ) const
{
    f->print();
    std::vector<Feature *> children = f->getChildren();

    for ( size_t i = 0; i != children.size(); ++i )
    {
// 		if ( !(*children)[i]->getChildren().empty())
        printForFeature ( children[i] );
    }

}

void FeatureTree::reMesh()
{
    needMeshing = true ;
}

bool FeatureTree::is3D() const
{
    if(structuredMesh)
        return true ;
    return tree[0]->spaceDimensions() == SPACE_THREE_DIMENSIONAL ;
}

bool FeatureTree::is2D() const
{
    if(structuredMesh)
        return false ;
    return tree[0]->spaceDimensions() == SPACE_TWO_DIMENSIONAL ;
}

void FeatureTree::initializeElements( )
{

    if(father3D && father3D->getOrder() != elemOrder)
        delete father3D ;

    if ( !father3D )
    {
        father3D = new TetrahedralElement ( elemOrder ) ;
    }

//     father3D->compileAndPrecalculate() ;

    if(father2D && father2D->getOrder() != elemOrder)
        delete father2D ;

    if ( !father2D )
    {
        father2D = new TriElement ( elemOrder ) ;
    }

//     father2D->compileAndPrecalculate() ;

    timeval time0, time1 ;
    gettimeofday ( &time0, nullptr );


    if ( is2D() )
    {

        std::cerr << " initialising..." << std::flush;

        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
        {

            int ecounter = 0 ;
            for ( auto i = j->second->begin() ; i != j->second->end() ; i++ )
            {
                if ( !i->getBehaviour() )
                {
                    std::cerr << "ouch" << std::endl ;
                }
                else
                {
                    i->refresh ( father2D );
                    i->getState().initialize ( j->second ) ;
                    ecounter++ ;
                    if ( ecounter % 100 == 0 )
                    {
                        std::cerr << "\r initialising... element " << ecounter << "/" << i.size() << std::flush ;
                    }
                }
            }
        }

        gettimeofday ( &time1, nullptr );
        int numtris = layer2d.begin()->second->begin().size() ;
        double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
        std::cerr << "\r initialising... element " << numtris << "/" << numtris << ". Time to initialise (s) " << delta / 1e6 << std::endl ;


    }

    if ( is3D() )
    {

        std::cerr << " initialising..." ;

//         #pragma omp parallel for schedule(runtime)
        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++  )
        {
            if ( !i->getBehaviour() )
            {
                continue ;
            }
            i->refresh ( father3D );
            i->getState().initialize ( dtree3D ) ;
        }

        gettimeofday ( &time1, nullptr );
        double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
        std::cerr << "\r initialising... element " << 0 << "/" << dtree3D->begin().size() << ". Time to initialise (s) " << delta / 1e6 << std::endl ;


        int ecounter = 0 ;

        for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
        {
            if ( !i->getBehaviour() )
            {
                std::cerr << "ouch" << std::endl ;
            }
            i->refresh ( father3D );
            i->getState().initialize ( dtree3D ) ;
            #pragma omp critical
            {
                ecounter++ ;
                if ( ecounter % 100 == 0 )
                {
                    std::cerr << "\r initialising... element " << ecounter << "/" << i.size() << std::flush ;
                }
            }
        }
    }

}

double FeatureTree::getDeltaTime () const
{
    return deltaTime ;
}

void FeatureTree::setMinDeltaTime ( double d )
{
    minDeltaTime = d ;
}

void FeatureTree::setMaxIterationsPerStep ( size_t its )
{
    maxitPerStep = its ;
}

void FeatureTree::setDeltaTime ( double d, bool isreal )
{
    previousDeltaTime = deltaTime ;
    deltaTime = d ;
    if(isreal)
        realDeltaTime = d ;
    if(needMeshing)
        return ;

    if ( dtree )
    {
        for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
        {
            previousDeltaTime = j->second->begin()->getBoundingPoint ( j->second->begin()->getBoundingPoints().size() -1 ).getT() - j->second->begin()->getBoundingPoint ( 0 ).getT() ;
            double end = j->second->begin()->getBoundingPoint ( j->second->begin()->getBoundingPoints().size() -1 ).getT() ;
            if ( j->second->begin().size() && j->second->begin()->timePlanes() > 1 )
            {
                for ( auto i = j->second->begin() ; i != j->second->end() ; i++ )
                {
                    size_t k0 = i->getBoundingPoints().size() /i->timePlanes() ;
                    for ( size_t t = 0 ; t < i->timePlanes() -1 ; t++ )
                    {
                        for ( size_t k = 0 ; k < k0 ; k++ )
                        {
                            i->getBoundingPoint ( k+k0*t ).getT() = end - d + d*t/ ( i->timePlanes()-1 ) ;
                        }
                    }

                    if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR)
                    {
                        i->adjustElementaryMatrix ( previousDeltaTime, d ) ;

//                        i->moved = true ;
                    }
                }
            }
        }
    }
    else
    {

        previousDeltaTime = dtree3D->begin()->getBoundingPoint ( dtree3D->begin()->getBoundingPoints().size() -1 ).getT() - dtree3D->begin()->getBoundingPoint ( 0 ).getT() ;
        double end = dtree3D->begin()->getBoundingPoint ( dtree3D->begin()->getBoundingPoints().size() -1 ).getT() ;
        if ( dtree3D->begin().size() && dtree3D->begin()->timePlanes() > 1 )
        {
            for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
            {
                size_t k0 = i->getBoundingPoints().size() /i->timePlanes() ;
                for ( size_t t = 0 ; t < i->timePlanes() -1 ; t++ )
                {
                    for ( size_t k = 0 ; k < k0 ; k++ )
                    {
                        i->getBoundingPoint ( k+k0*t ).getT() = end - d + d*t/ ( i->timePlanes()-1 ) ;
                    }
                }

                if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    i->adjustElementaryMatrix ( previousDeltaTime, d ) ;
                }
            }
        }
    }
}

void FeatureTree::moveFirstTimePlanes ( double d, const Mesh<DelaunayTriangle,  DelaunayTreeItem >::iterator & begin,  const Mesh<DelaunayTriangle, DelaunayTreeItem >::iterator & end)
{
    Mesh<DelaunayTriangle,  DelaunayTreeItem >::iterator i(begin) ;
    double prev = i->getBoundingPoint ( i->getBoundingPoints().size() -1 ).getT() - i->getBoundingPoint ( 0 ).getT() ;
    size_t ndof = 0 ;
    while ( ndof == 0 && i != end )
    {
        ndof = i->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        i++ ;
    }
    Vector buff ( 0.,ndof ) ;

    if ( dtree )
    {


        if ( std::abs ( d ) > POINT_TOLERANCE )
        {
            setDeltaTime ( prev - d, false ) ;
        }

        for(size_t i = 0 ; i < boundaryCondition.size() ; i++)
        {
            TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *>(boundaryCondition[i]) ;
            if(timec)
            {
                timec->instant = d/prev ;
            }
        }

    }
}

void FeatureTree::moveFirstTimePlanes ( double d, const Mesh<DelaunayTetrahedron,DelaunayTreeItem3D >::iterator & begin,  const Mesh<DelaunayTetrahedron, DelaunayTreeItem3D >::iterator & end)
{
    double prev = 0.;
    auto i = begin ;
    prev = i->getBoundingPoint ( i->getBoundingPoints().size() -1 ).getT() - i->getBoundingPoint ( 0 ).getT() ;
    size_t ndof = 0 ;
    while ( ndof == 0 && i != end )
    {
        ndof = i->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        i++ ;
    }
    Vector buff ( 0.,ndof ) ;

    if ( dtree3D )
    {


        if ( std::abs ( d ) > POINT_TOLERANCE )
        {
            setDeltaTime ( prev - d, false ) ;
        }

        for(size_t i = 0 ; i < boundaryCondition.size() ; i++)
        {
            TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *>(boundaryCondition[i]) ;
            if(timec)
            {
                timec->instant = d/prev ;
            }
        }

    }
}

void FeatureTree::generateElements()
{

    for ( size_t i = 0 ; i < boundaryCondition.size() ; i++ )
    {
        boundaryCondition[i]->clearCache() ;
    }

    if ( dtree || dtree3D )
    {
        if ( K )
        {
            K->clear() ;
        }

    }

    needMeshing = false ;

    if(structuredMesh)
        return ;

    double pointDensity = 0 ;

    if ( is2D() )
    {
        pointDensity = sqrt ( tree[0]->area() / ( tree[0]->getBoundingPoints().size() + tree[0]->getInPoints().size() ) ) ;
    }
    else
    {
        pointDensity = pow ( tree[0]->volume() / ( tree[0]->getBoundingPoints().size() + tree[0]->getInPoints().size() ), .33333333333 ) ;
    }

    std::cerr << "space meshed with " << pointDensity << " points per unit length" << std::endl ;

    std::vector<Point> bbox = tree[0]->getBoundingBox() ;
    double min_x = 0, min_y = 0, max_x = 0, max_y = 0, max_z = 0, min_z = 0;

    for ( size_t j  =  0 ; j <  bbox.size() ; j++ )
    {
        if ( bbox[j].getY() < min_y )
        {
            min_y = bbox[j].getY() ;
        }

        if ( bbox[j].getY() > max_y )
        {
            max_y = bbox[j].getY() ;
        }

        if ( bbox[j].getX() < min_x )
        {
            min_x = bbox[j].getX() ;
        }

        if ( bbox[j].getX() > max_x )
        {
            max_x = bbox[j].getX() ;
        }

        if ( bbox[j].getZ() < min_z )
        {
            min_z = bbox[j].getZ() ;
        }

        if ( bbox[j].getZ() > max_z )
        {
            max_z = bbox[j].getZ() ;
        }
    }

    if ( is3D() )
    {
        bbox[0] = Point ( min_x, min_y, min_z ) ;
        bbox[1] = Point ( min_x, min_y, max_z ) ;
        bbox[2] = Point ( min_x, max_y, min_z ) ;
        bbox[3] = Point ( min_x, max_y, max_z ) ;
        bbox[4] = Point ( max_x, min_y, min_z ) ;
        bbox[5] = Point ( max_x, min_y, max_z ) ;
        bbox[6] = Point ( max_x, max_y, min_z ) ;
        bbox[7] = Point ( max_x, max_y, max_z ) ;
    }
    else
    {
        bbox.push_back ( Point() );
        bbox.push_back ( Point() );
        bbox.push_back ( Point() );
        bbox.push_back ( Point() );
        bbox[0] = Point ( min_x, min_y, min_z ) ;
        bbox[1] = Point ( min_x, min_y, max_z ) ;
        bbox[2] = Point ( min_x, max_y, min_z ) ;
        bbox[3] = Point ( min_x, max_y, max_z ) ;
        bbox[4] = Point ( max_x, min_y, min_z ) ;
        bbox[5] = Point ( max_x, min_y, max_z ) ;
        bbox[6] = Point ( max_x, max_y, min_z ) ;
        bbox[7] = Point ( max_x, max_y, max_z ) ;
    }

    size_t basepoints = 0 ;
    std::cerr << " getting mesh points..." << std::flush ;
    std::vector<Feature *> nullFatherFeatures ;

    for ( size_t i  = 1 ; i < tree.size() ; i++ )
    {
        if ( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature )
        {
            if ( !tree[i]->getFather() )
            {
                nullFatherFeatures.push_back ( tree[i] );
            }
        }
    }


    for ( size_t i  = 0 ; i < tree.size() ; i++ )
    {
        std::cerr << "\r getting mesh points... feature " << i << "/" << tree.size() << std::flush ;

        if ( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature )
        {
            std::vector<Feature *> descendants = tree[i]->getDescendants() ;
            std::stable_sort ( descendants.begin(), descendants.end() ) ;

            for ( size_t j  =  0 ; j <  tree[i]->getBoundingPoints().size() ; j++ )
            {
                bool isIn = false ;

                std::vector<const Geometry *> potentialFeaturestmp  ;
                std::vector<Feature *> potentialFeatures ;

                if ( is2D() )
                {
                    potentialFeaturestmp = grid->coOccur ( tree[i]->getBoundingPoint ( j ) ) ;
                }
                else
                {
                    potentialFeaturestmp = grid3d->coOccur ( tree[i]->getBoundingPoint ( j ) ) ;
                }

                for ( size_t l = 0 ; l < potentialFeaturestmp.size() ; l++ )
                {
                    potentialFeatures.push_back ( const_cast<Feature *> ( dynamic_cast<const Feature *> ( potentialFeaturestmp[l] ) ) ) ;
                }

                std::vector<Feature *> potentialChildren ;

                for ( size_t l = 0 ; l < potentialFeatures.size() ; l++ )
                {
                    if ( !potentialFeatures[l]->isEnrichmentFeature
                            && std::binary_search ( descendants.begin(), descendants.end(), potentialFeatures[l] ) )
                    {
                        potentialChildren.push_back ( potentialFeatures[l] ) ;
                    }
                }

                if ( tree.size() < 128 )
                {
                    potentialChildren = descendants ;
                }

                for ( size_t k  =  0 ; k <  potentialFeatures.size() ; k++ )
                {
                    if ( i && tree[i] != potentialFeatures[k]
                            && !potentialFeatures[k]->isVirtualFeature
                            && potentialFeatures[k]->onBoundary ( tree[i]->getBoundingPoint ( j ), POINT_TOLERANCE )
                            && potentialFeatures[k]->getBoundingPoints().size() )
                    {
                        isIn = true ;
                        break ;
                    }
                }
                for ( size_t k  =  0 ; k <  potentialChildren.size() ; k++ )
                {
                    if ( ( !potentialChildren[k]->isVirtualFeature
                            && (potentialChildren[k]->inBoundary ( tree[i]->getBoundingPoint ( j ), pointDensity*(i == 0 ? 0.1 : 0.33) ))
                            && potentialChildren[k]->inMask( tree[i]->getBoundingPoint ( j ), pointDensity*(i == 0 ? 0.1 : 0.33)) )
                            || ( potentialChildren[k]->isVirtualFeature
                                 && tree[i]->isVirtualFeature
                                 && ( dynamic_cast<VirtualFeature *> ( potentialChildren[k] )->getSource()
                                      != dynamic_cast<VirtualFeature *> ( tree[i] )->getSource() )
                                 && potentialChildren[k]->inBoundary ( tree[i]->getBoundingPoint ( j ), pointDensity*0.33 )
                               )
                       )
                    {
                        if ( potentialChildren[k]->getBoundingPoints().size() )
                        {
                            isIn = true ;
                            break ;
                        }
                    }
                }

                if( i && ((!tree[i]->inMask( tree[i]->getBoundingPoint ( j ), pointDensity*0.1 )) || tree[i]->onMaskBoundary( tree[i]->getBoundingPoint ( j ), pointDensity*0.1 )))
                {
                    isIn = true ;
                }

                if ( i && tree[i]->getFather() && !inRoot ( tree[i]->getBoundingPoint ( j ) ) )
                {
                    isIn = true ;
                }

                if ( !isIn && tree[i]->isVirtualFeature && !tree[i]->in ( tree[i]->getBoundingPoint ( j ) ) )
                {
                    isIn = true ;
                }

//                if ( tree[i]->getFather() && tree[i]->getFather()->onBoundary ( tree[i]->getBoundingPoint ( j ), pointDensity*0.33 ) )
//                {
//                    isIn = true ;
//                }

                if ( tree[i]->getFather() && !isIn && i && tree[0]->onBoundary ( tree[i]->getBoundingPoint ( j ), pointDensity*0.33 ) )
                {
                    isIn = true ;
                }

                if ( !tree[i]->getFather() && i )
                {
                    isIn = false ;
                    for ( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
                    {
                        if ( nullFatherFeatures[k] == tree[i] )
                        {
                            break ;
                        }

                        if ( nullFatherFeatures[k]->inBoundary ( tree[i]->getBoundingPoint ( j ), 2.*POINT_TOLERANCE ) )
                        {
                            isIn = true ;
                            break ;
                        }
                    }

                    Point proj ( tree[i]->getBoundingPoint ( j ) ) ;
                    tree[0]->project ( &proj ) ;
                    if ( dist ( proj, tree[i]->getBoundingPoint ( j ) ) < 2.*POINT_TOLERANCE )
                    {
                        isIn = false ;
                    }
                }

                //the border is always defined by non-root features of nullptr father
                if ( !i && !isIn )
                {
                    for ( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
                    {
                        Point proj ( tree[i]->getBoundingPoint ( j ) ) ;
                        nullFatherFeatures[k]->project ( &proj ) ;
                        if ( dist ( proj, tree[i]->getBoundingPoint ( j ) ) < 2.*POINT_TOLERANCE )
                        {
                            isIn = true ;
                            break ;
                        }
                    }
                }


                if ( !isIn )
                {
                    meshPoints.push_back ( std::pair<Point *, Feature *> ( &tree[i]->getBoundingPoint ( j ), this->tree[i] ) ) ;

                    if ( i == 0 )
                    {
                        basepoints++ ;
                    }
                }
            }

            for ( size_t j  =  0 ; j <  tree[i]->getInPoints().size() ; j++ )
            {
                bool isIn = false ;
                std::vector<const Geometry *> potentialFeaturestmp  ;

                if ( is2D() )
                {
                    potentialFeaturestmp = grid->coOccur ( tree[i]->getInPoint ( j ) ) ;
                }
                else
                {
                    potentialFeaturestmp = grid3d->coOccur ( tree[i]->getInPoint ( j ) ) ;
                }

                std::vector<Feature *> potentialFeatures ;

                for ( size_t k = 0 ; k < potentialFeaturestmp.size() ; ++k )
                {
                    potentialFeatures.push_back ( const_cast<Feature *> ( dynamic_cast<const Feature *> ( potentialFeaturestmp[k] ) ) ) ;
                }

                std::vector<Feature *> potentialChildren ;

                for ( size_t l = 0 ; l < potentialFeatures.size() ; l++ )
                {
                    if ( !potentialFeatures[l]->isVirtualFeature
                            && !potentialFeatures[l]->isEnrichmentFeature
                            && std::binary_search ( descendants.begin(), descendants.end(), potentialFeatures[l] ) )
                    {
                        potentialChildren.push_back ( potentialFeatures[l] ) ;
                    }
                }

                if ( tree.size() < 128 )
                {
                    potentialChildren = descendants ;
                }

                for ( size_t l = 0 ; l < potentialFeatures.size() ; l++ )
                {
                    if ( !potentialFeatures[l]->isVirtualFeature
                            && !potentialFeatures[l]->isEnrichmentFeature
                            && potentialFeatures[l]->isMaskedBy( tree[i] ) )
                    {
                        potentialChildren.push_back ( potentialFeatures[l] ) ;
                    }
                }

                for ( size_t k  =  0 ; k <  potentialChildren.size() ; k++ )
                {
                    if (
                        (
                            !potentialChildren[k]->isVirtualFeature
                            && potentialChildren[k]->inBoundary ( tree[i]->getInPoint ( j ), pointDensity*0.5 )
                            && potentialChildren[k]->inMask ( tree[i]->getInPoint ( j ), pointDensity*0.1 )
                        )
                        ||
                        (
                            potentialChildren[k]->isVirtualFeature
                            && tree[i]->isVirtualFeature
                            && (
                                dynamic_cast<VirtualFeature *> ( potentialChildren[k] )->getSource()
                                != dynamic_cast<VirtualFeature *> ( tree[i] )->getSource()
                            )
                            && potentialChildren[k]->inBoundary ( tree[i]->getInPoint ( j ), pointDensity*0.33 )
                        )
                    )
                    {
                        if ( potentialChildren[k]->getBoundingPoints().size() )
                        {
                            isIn = true ;
                            break ;
                        }
                    }
                }

                if ( i && !inRoot ( tree[i]->getInPoint ( j ) ) )
                {
                    isIn = true ;
                }

                if( i && ( (!tree[i]->inMask( tree[i]->getInPoint ( j ) ))  || tree[i]->onMaskBoundary( tree[i]->getInPoint ( j ), pointDensity*0.5 )))
                {
                    isIn = true ;
                }

                if ( tree[i]->getFather() && tree[i]->getFather()->onBoundary ( tree[i]->getInPoint ( j ), pointDensity*0.33 ) && ! tree[i]->isMaskedBy( tree[i]->getFather() ) )
                {
                    isIn = true ;
                }

                if ( tree[i]->isVirtualFeature && !tree[i]->in ( tree[i]->getInPoint ( j ) ) )
                {
                    isIn = true ;
                }

                if ( tree[i]->getFather() && i && tree[0]->onBoundary ( tree[i]->getInPoint ( j ), pointDensity*0.33 ) )
                {
                    isIn = true ;
                }

                if ( !tree[i]->getFather() && i )
                {
                    isIn = false ;
                }


                for ( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
                {
                    if ( tree[i] == nullFatherFeatures[k] )
                    {
                        break ;
                    }

                    if ( nullFatherFeatures[k]->inBoundary ( tree[i]->getInPoint ( j ), 2.*POINT_TOLERANCE ) )
                    {
                        isIn = true ;
                        break ;
                    }
                }

                if ( !isIn )
                {
                    meshPoints.push_back ( std::pair<Point *, Feature *> ( &tree[i]->getInPoint ( j ), tree[i] ) ) ;

                    if ( i == 0 )
                    {
                        basepoints++ ;
                    }
                }
            }

        }
    }

    std::cerr << "...done" << std::endl ;

    size_t count  = 0 ;

    if ( computeIntersections )
    {
        std::set< std::pair<const Feature*, const Feature*> > done ;
        for ( size_t i = 1 ;  i < tree.size() ; i++ )
        {
            if ( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getFather() != nullptr )
            {
                std::vector<const Geometry *> coOccuringFeaturestmp ;
                std::vector<Feature *> descendants = tree[i]->getDescendants() ;

                if ( is3D() )
                {
                    coOccuringFeaturestmp = grid3d->coOccur ( tree[i] ) ;
                }
                else
                {
                    coOccuringFeaturestmp = grid->coOccur ( tree[i] ) ;
                }

                std::vector<const Feature *> coOccuringFeatures ;

                for ( size_t k = 0 ; k < coOccuringFeaturestmp.size() ; k++ )
                {
                    coOccuringFeatures.push_back ( dynamic_cast<const Feature *> ( coOccuringFeaturestmp[k] ) ) ;
                }

                /*                if(tree[i]->getGeometryType() == ELLIPSE)
                		{
                			for(size_t j = 0 ; j < coOccuringFeatures.size() ; j++ )
                			{
                				if(coOccuringFeatures[j]->getGeometryType() != ELLIPSE)
                				{
                					std::cout << i << " " << j << " " << std::string( tree[i]->intersects(coOccuringFeatures[j]) ? "yes" : "no" ) << " " << tree[i]->intersection(coOccuringFeatures[j]).size() << std::endl ;
                				}
                			}
                		}*/


                for ( size_t j  = 0 ; j < coOccuringFeatures.size() ; j++ )
                {
                    if ( !coOccuringFeatures[j]->isEnrichmentFeature
                            && !coOccuringFeatures[j]->isVirtualFeature
                            && tree[i] != coOccuringFeatures[j]
                            && tree[i]->intersects ( coOccuringFeatures[j] ) )
                    {
                        std::pair<const Feature*, const Feature*> testa = std::make_pair( tree[i], coOccuringFeatures[j] ) ;
                        std::pair<const Feature*, const Feature*> testb = std::make_pair( coOccuringFeatures[j], tree[i] ) ;
                        if(done.find( testa ) != done.end() || done.find( testb ) != done.end() )
                        {
                            continue ;
                        }
                        done.insert(testa) ;
                        done.insert(testb) ;

                        std::vector<Point> inter = tree[i]->intersection ( coOccuringFeatures[j] ) ;

                        int base = inter.size() ;
                        for( size_t k = 0 ; k < tree[i]->getBoundingPoints().size() ; k++)
                        {
                            if( coOccuringFeatures[j]->inBoundary(tree[i]->getBoundingPoint(k), pointDensity*0.1) && coOccuringFeatures[j]->inMask(tree[i]->getBoundingPoint(k), pointDensity*0.1))
                            {
                                bool alone = true ;
                                for(int l = 0 ; alone && l < base ; l++)
                                    alone = ( dist( inter[l], tree[i]->getBoundingPoint(k) ) > pointDensity*0.33 ) ;
                                if(alone)
                                    inter.push_back( tree[i]->getBoundingPoint(k) ) ;
                            }
                        }
//			std::cout << i << " NOPE " << tree[i]->getGeometryType()  << " "<< coOccuringFeatures[j]->getGeometryType() << " " << inter.size() << " " << tree[i]->getBoundingPoints().size() << std::endl ;
                        if(inter.size() < 3)
                        {
                            continue ;
                        }
                        for ( size_t k = 0 ;  k < inter.size() ; k++ )
                        {

                            bool indescendants = false ;

                            for ( size_t l = 0 ; l < descendants.size() ; l++ )
                            {
                                if ( !descendants[l]->isVirtualFeature && descendants[l]->inBoundary ( inter[k], pointDensity*0.33 ) && !descendants[l]->isMaskedBy( tree[i] ))
                                {
                                    indescendants = true ;
                                    break ;
                                }
                            }

//
                            if ( !indescendants )
                            {
                                if ( inRoot ( inter[k] ) && tree[i]->inMask( inter[k], pointDensity*0.1 ) && (!tree[i]->onMaskBoundary( inter[k], pointDensity*0.33) || tree[i]->isMaskedBy( coOccuringFeatures[j])) )
                                {
                                    Point *p = new Point ( inter[k] ) ;
                                    additionalPoints.push_back ( p ) ;
                                    ++count ;
                                    meshPoints.push_back ( std::make_pair ( p, tree[i] ) ) ;
                                }
                            }

                            if ( count % 100 == 0 )
                            {
                                std::cerr << "\r adding intersection points... " << count << std::flush ;
                            }
                        }
                    }
                }
            }
        }

        for ( auto & feature : tree )
        {
            if(feature == tree[0])
                continue ;

            if ( !feature->isEnrichmentFeature && feature->getBoundingPoints().size() && !feature->isVirtualFeature && tree[0]->intersects ( feature ) && feature->getFather() != nullptr )
            {
                std::vector<Point> inter = tree[0]->intersection ( feature ) ;
                if(inter.size() == 0)
                    inter = feature->intersection( tree[0] ) ;
                std::vector<Feature *> descendants = feature->getDescendants() ;
                std::vector<Feature *> fatherdescendants = tree[0]->getDescendants() ;

                size_t base = inter.size() ;
                for( size_t k = 0 ; k < tree[0]->getBoundingPoints().size() ; k++)
                {
                    if( feature->inBoundary(tree[0]->getBoundingPoint(k), pointDensity*0.1) && feature->inMask(tree[0]->getBoundingPoint(k), pointDensity*0.1))
                    {
                        bool alone = true ;
                        for(size_t l = 0 ; alone && l < base ; l++)
                            alone = ( dist( inter[l], tree[0]->getBoundingPoint(k) ) > pointDensity*0.33 ) ;
                        if(alone)
                            inter.push_back( tree[0]->getBoundingPoint(k) ) ;
                    }
                }
//                size_t start = ( base == inter.size() ? 0 : base ) ;

                for ( size_t k = 0 ;  k < inter.size() ; k++ )
                {



                    bool indescendants = false ;

                    for ( size_t l = 0 ; l < descendants.size() ; l++ )
                    {
                        if ( descendants[l]->inBoundary ( inter[k], pointDensity*0.33 ) && descendants[l]->inMask ( inter[k], pointDensity*0.33 ) )
                        {
                            indescendants = true ;
                            break ;
                        }
                    }

                    for ( size_t l = 0 ; l < fatherdescendants.size() ; l++ )
                    {
                        if ( fatherdescendants[l] != feature && !fatherdescendants[l]->isVirtualFeature && fatherdescendants[l]->getBoundingPoints().size() && fatherdescendants[l]->inBoundary ( inter[k], pointDensity*0.33 ) && !feature->onBoundary ( inter[k], pointDensity*0.33 ) )
                        {
//                            indescendants = true ;
                            break ;
                        }
                    }

                    if ( is3D() )
                    {

                        Point proj ( inter[k] ) ;
                        tree[0]->project ( &proj ) ;
                        Point proj0 ( inter[k] + Point ( 2.*POINT_TOLERANCE, 0, 0 ) ) ;
                        Point proj1 ( inter[k] + Point ( -2.*POINT_TOLERANCE, 0, 0 ) ) ;
                        Point proj2 ( inter[k] + Point ( 0, 2.*POINT_TOLERANCE, 0 ) ) ;
                        Point proj3 ( inter[k] + Point ( 0, -2.*POINT_TOLERANCE, 0 ) ) ;
                        Point proj4 ( inter[k] + Point ( 0, 0, 2.*POINT_TOLERANCE ) ) ;
                        Point proj5 ( inter[k] + Point ( 0, 0, -2.*POINT_TOLERANCE ) ) ;

                        int position = tree[0]->in ( proj0 )
                                       + tree[0]->in ( proj1 )
                                       + tree[0]->in ( proj2 )
                                       + tree[0]->in ( proj3 )
                                       + tree[0]->in ( proj4 )
                                       + tree[0]->in ( proj5 ) ;

                        bool onSurface = ( position == 5 ) ;
                        bool onEdge = ( position == 4 ) ;
                        bool onVertex = ( position == 3 ) ;
                        proj0 = ( inter[k] + Point ( pointDensity, 0, 0 ) ) ;
                        proj1 = ( inter[k] + Point ( -pointDensity, 0, 0 ) ) ;
                        proj2 = ( inter[k] + Point ( 0, pointDensity, 0 ) ) ;
                        proj3 = ( inter[k] + Point ( 0, -pointDensity, 0 ) ) ;
                        proj4 = ( inter[k] + Point ( 0, 0, pointDensity ) ) ;
                        proj5 = ( inter[k] + Point ( 0, 0, -pointDensity ) ) ;
                        int tooClose =  tree[0]->in ( proj0 )
                                        + tree[0]->in ( proj1 )
                                        + tree[0]->in ( proj2 )
                                        + tree[0]->in ( proj3 )
                                        + tree[0]->in ( proj4 )
                                        + tree[0]->in ( proj5 ) ;

                        // no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
                        if ( !indescendants && squareDist3D ( proj, inter[k] ) < POINT_TOLERANCE * POINT_TOLERANCE && /*inRoot(inter[k]) && */ ( ( onSurface && tooClose == 5 ) || ( onEdge && tooClose == 4 ) || onVertex ) )
                        {
                            Point *p = new Point ( inter[k] ) ;
                            additionalPoints.push_back ( p ) ;
                            ++count ;
                            meshPoints.push_back ( std::make_pair ( p, feature ) ) ;
                        }
                    }
                    else
                    {
                        Point proj ( inter[k] ) ;
                        tree[0]->project ( &proj ) ;
                        Point proj0 ( inter[k] + Point ( 2.*POINT_TOLERANCE, 0, 0 ) ) ;
                        Point proj1 ( inter[k] + Point ( -2.*POINT_TOLERANCE, 0, 0 ) ) ;
                        Point proj2 ( inter[k] + Point ( 0, 2.*POINT_TOLERANCE, 0 ) ) ;
                        Point proj3 ( inter[k] + Point ( 0, -2.*POINT_TOLERANCE, 0 ) ) ;


                        int position = tree[0]->in ( proj0 )
                                       + tree[0]->in ( proj1 )
                                       + tree[0]->in ( proj2 )
                                       + tree[0]->in ( proj3 );

                        bool onEdge = ( position == 3 ) ;
                        bool onVertex = ( position == 2 ) ;
                        proj0 = ( inter[k] + Point ( pointDensity, 0, 0 ) ) ;
                        proj1 = ( inter[k] + Point ( -pointDensity, 0, 0 ) ) ;
                        proj2 = ( inter[k] + Point ( 0, pointDensity, 0 ) ) ;
                        proj3 = ( inter[k] + Point ( 0, -pointDensity, 0 ) ) ;

                        int tooClose =  tree[0]->in ( proj0 )
                                        + tree[0]->in ( proj1 )
                                        + tree[0]->in ( proj2 )
                                        + tree[0]->in ( proj3 );


                        // no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
                        if ( ( feature->inBoundary ( inter[k], pointDensity*0.33 ) && feature->inMask ( inter[k], pointDensity*0.1 ) ) || ( !indescendants && squareDist3D ( proj, inter[k] ) < POINT_TOLERANCE * POINT_TOLERANCE && inRoot ( inter[k] ) && ( ( onEdge && tooClose == 3 ) || onVertex ) ) )
                        {
                            Point *p = new Point ( inter[k] ) ;
                            additionalPoints.push_back ( p ) ;
                            ++count ;
                            meshPoints.push_back ( std::make_pair ( p, feature ) ) ;
                        }
                        else
                        {

                        }
                    }

                }

                if ( count % 100 == 0 )
                {
                    std::cerr << "\r adding intersection points... " << count << std::flush ;
                }
            }

        }
    }

    std::cerr << "\r adding intersection points... " << count << " ...done." << std::endl ;
    count = 0 ;

    //shuffle for efficiency
    shuffleMeshPoints() ;

    if ( is2D() )
    {
        additionalPoints.push_back ( new Point ( bbox[0] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[2] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[4] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[6] ) ) ;
        for ( auto i = extraPoints.begin() ; i!= extraPoints.end() ; ++i )
        {
            meshPoints.push_front ( std::make_pair ( *i, tree[0] ) ) ;
        }
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 4], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 3], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 2], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 1], tree[0] ) ) ;

        Mesh<DelaunayTriangle, DelaunayTreeItem> * oldDtree = dtree ;

#ifdef HAVE_OPENMP
        double t0 = omp_get_wtime() ;
#endif

        dtree = new /*Parallel*/DelaunayTree ( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first/*, domains*/ ) ;
        dtree->insert ( meshPoints[3].first ) ;
        layer2d[-1] = dtree ;

        for ( size_t i  = 0 ; i < tree.size() ; i++ )
        {
            if ( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getLayer() != -1 )
            {
                if ( layer2d.find ( tree[i]->getLayer() ) == layer2d.end() )
                {
                    layer2d[tree[i]->getLayer()] = new /*Parallel*/DelaunayTree ( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first/*, domains*/ ) ;
                    layer2d[tree[i]->getLayer()]->insert ( meshPoints[3].first ) ;
                }
            }
        }

        std::vector<size_t> iterators ( meshPoints.size()-4 ) ;
        for ( size_t i = 0 ; i < iterators.size() ; i++ )
        {
            iterators[i] = i+4 ;
        }

        std::random_shuffle ( iterators.begin(), iterators.end() );

        for ( size_t i = 0 ; i < iterators.size() ; i++ )
        {
            if ( ( i ) % 1000 == 0 )
            {
                std::cerr << "\r generating triangles... point " << count << "/" << meshPoints.size() << std::flush ;
            }

            ++count ;

            if ( *meshPoints[iterators[i]].first != bbox[0] &&
                    *meshPoints[iterators[i]].first != bbox[2] &&
                    *meshPoints[iterators[i]].first != bbox[4] &&
                    *meshPoints[iterators[i]].first != bbox[6] && ( inRoot ( *meshPoints[iterators[i]].first ) || meshPoints[iterators[i]].second->getFather() == nullptr )
               )
            {
                for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                {
                    j->second->insert ( meshPoints[iterators[i]].first ) ;
                }
            }
        }

        std::cerr << "\r generating triangles.... point " << meshPoints.size() - 3 << "/" << meshPoints.size() - 4 << std::flush ;


        bool correct = false ;
        int tries = correctionSteps ;

        while ( !correct && tries )
        {

            std::vector< Point *> to_insert ;

            for ( auto i = dtree->begin() ; i != dtree->end() ; i++ )
            {
                Point *test = checkElement ( i );

                if ( test )
                {
                    to_insert.push_back ( test );
                }
            }

            if ( to_insert.empty() )
            {
                correct = true ;
            }

            for ( size_t i = 0 ; i < to_insert.size() ; i++ )
            {
                std::cerr << "\r generating triangles.... point " << ++count << "/" << meshPoints.size() - 3 << std::flush ;

                if ( *to_insert[i] != bbox[0] &&
                        *to_insert[i] != bbox[2] &&
                        *to_insert[i] != bbox[4] &&
                        *to_insert[i] != bbox[6] &&
                        inRoot ( *to_insert[i] )
                   )
                {
                    for ( auto j = layer2d.begin() ; j != layer2d.end() ; j++ )
                    {
                        j->second->insert ( to_insert[i] ) ;
                    }
                }

                if ( to_insert[i]->getId() == -1 )
                {
                    delete to_insert[i] ;
                }
            }

            tries-- ;

        }

#ifdef HAVE_OPENMP
        std::cerr << " ...done. " << omp_get_wtime() - t0 <<  " (s)"<< std::endl ;
#else
        std::cerr << " ...done. " << std::endl ;
#endif

        for ( size_t i = 0 ; i < refinementZones.size() ; i++ )
        {
            quadTreeRefine ( refinementZones[i] ) ;
        }

        if ( oldDtree )
        {
            setElementBehavioursFromMesh<Mesh<DelaunayTriangle, DelaunayTreeItem>, DelaunayTriangle> ( oldDtree, dtree ) ;
            delete oldDtree ;
        }

    }
    else if ( is3D() )
    {
        additionalPoints.push_back ( new Point ( bbox[0] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[1] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[7] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[2] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[3] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[4] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[5] ) ) ;
        additionalPoints.push_back ( new Point ( bbox[6] ) ) ;
        for ( auto i = extraPoints.begin() ; i!= extraPoints.end() ; ++i )
        {
            meshPoints.push_front ( std::make_pair ( *i, tree[0] ) ) ;
        }
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 8], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 7], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 6], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 5], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 4], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 3], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 2], tree[0] ) ) ;
        meshPoints.push_front ( std::make_pair ( additionalPoints[additionalPoints.size() - 1], tree[0] ) ) ;


        Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * oldDtree = dtree3D ;

#ifdef HAVE_OPENMP
        double t0 = omp_get_wtime() ;
#endif
        dtree3D = new /*Parallel*/DelaunayTree3D ( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first, meshPoints[3].first/*, domains*/ ) ;
        dtree3D->insert ( meshPoints[4].first ) ;
        dtree3D->insert ( meshPoints[5].first ) ;
        dtree3D->insert ( meshPoints[6].first ) ;
        dtree3D->insert ( meshPoints[7].first ) ;


        std::pair<std::vector<int>, std::vector<int> > pb ;

        for ( size_t i = 0 ; i < meshPoints.size() ; i++ )
        {
            if ( meshPoints[i].first == nullptr )
            {
                pb.first.push_back ( i ) ;
            }

            if ( meshPoints[i].second == nullptr )
            {
                pb.second.push_back ( i ) ;
            }
        }

        std::vector<Point *> toInsert ;

        for ( auto i = meshPoints.begin() + 8 ; i != meshPoints.end(); ++i )
        {

            if ( ( i - meshPoints.begin() ) % 1000 == 0 )
            {
                std::cerr << "\r generating tetrahedrons... point " << ++count * 1000 << "/" << meshPoints.size() - 8 <<  std::flush ;
            }

            if ( *i->first != bbox[0] &&
                    *i->first != bbox[1] &&
                    *i->first != bbox[2] &&
                    *i->first != bbox[3] &&
                    *i->first != bbox[4] &&
                    *i->first != bbox[5] &&
                    *i->first != bbox[6] &&
                    *i->first != bbox[7]
               )
            {

                dtree3D->insert ( i->first ) ;


                if ( i->first->getId() == -1 )
                {
                    toInsert.push_back ( i->first ) ;
                }
            }
        }

        bool correct = false ;
        int tries = correctionSteps ;

        while ( !correct && tries )
        {

            std::vector< Point *> to_insert ;

            for ( auto i = dtree3D->begin() ; i != dtree3D->end() ; i++ )
            {
                Point *test = checkElement ( i );

                if ( test )
                {
                    to_insert.push_back ( test );
                }
            }

            if ( to_insert.empty() )
            {
                correct = true ;
            }

            for ( size_t i = 0 ; i < to_insert.size() ; i++ )
            {
                std::cerr << "\r generating tetrahedrons.... point " << ++count << "/" << meshPoints.size() - 3 << std::flush ;

                if ( *to_insert[i] != bbox[0] &&
                        *to_insert[i] != bbox[1] &&
                        *to_insert[i] != bbox[2] &&
                        *to_insert[i] != bbox[3] &&
                        *to_insert[i] != bbox[4] &&
                        *to_insert[i] != bbox[5] &&
                        *to_insert[i] != bbox[6] &&
                        *to_insert[i] != bbox[7] &&
                        inRoot ( *to_insert[i] )
                   )
                {
                    dtree3D->insert ( to_insert[i] ) ;
                }

                if ( to_insert[i]->getId() == -1 )
                {
                    delete to_insert[i] ;
                }
            }

            tries-- ;

        }

#ifdef HAVE_OPENMP
        std::cerr << " ...done. " << omp_get_wtime() - t0 <<  " (s)"<< std::endl ;
#else
        std::cerr << " ...done. " << std::endl ;
#endif
//		dtree3D->purge() ;

        if ( oldDtree )
        {
            setElementBehavioursFromMesh<Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>, DelaunayTetrahedron> ( oldDtree, dtree3D ) ;
            delete oldDtree ;
        }

    }
}

void FeatureTree::shuffleMeshPoints()
{
    std::random_shuffle ( meshPoints.begin(), meshPoints.end() ) ;
}

void FeatureTree::printReport(bool printHeader, bool vertical)
{
    if(is2D())
    {
        Vector x = getDisplacements() ;

        std::cout << "unknowns :" << x.size() << std::endl ;

        int npoints = get2DMesh()->begin()->getBoundingPoints().size() ;

        double volume = 0 ;

        double xavg = 0 ;

        for(auto k = get2DMesh()->begin() ; k != get2DMesh()->end() ; k++)
        {
            if(k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                double ar = k->area() ;
                volume += ar ;
                for(int l = 0 ; l < npoints ; l++)
                {
                    xavg += x[k->getBoundingPoint(l).getId()*2]*ar/npoints ;
                }
            }
        }

        xavg /= volume ;
        reportValues.push_back(Vector(23));
        std::pair<Vector, Vector> stempm = getFieldMinMax(REAL_STRESS_FIELD) ;
        std::pair<Vector, Vector> etempm = getFieldMinMax(STRAIN_FIELD) ;
        std::pair<Vector, Vector> vmm = getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
        Vector stemp = getAverageField(REAL_STRESS_FIELD) ;
        Vector etemp = getAverageField(STRAIN_FIELD) ;
        reportValues.back()[0] = x.max() ;
        reportValues.back()[1] = x.min() ;
        reportValues.back()[2] = xavg    ;
        reportValues.back()[3] = stempm.second[0] ;
        reportValues.back()[4] = stempm.first[0]  ;
        reportValues.back()[5] = stempm.second[1] ;
        reportValues.back()[6] = stempm.first[1]  ;
        reportValues.back()[7] = stempm.second[2] ;
        reportValues.back()[8] = stempm.first[2]   ;
        reportValues.back()[9] =  etempm.second[0] ;
        reportValues.back()[10] =  etempm.first[0]  ;
        reportValues.back()[11] =  etempm.second[1] ;
        reportValues.back()[12] =  etempm.first[1]  ;
        reportValues.back()[13] =  etempm.second[2]  ;
        reportValues.back()[14] =  etempm.first[2]  ;
        reportValues.back()[15] =  vmm.second[0]    ;
        reportValues.back()[16] =  vmm.first[0]     ;
        reportValues.back()[17] =   stemp[0]        ;
        reportValues.back()[18] =   stemp[1]        ;
        reportValues.back()[19] =   stemp[2]        ;
        reportValues.back()[20] =   etemp[0]        ;
        reportValues.back()[21] =   etemp[1]        ;
        reportValues.back()[22] =   etemp[2]      ;
        std::cout << std::endl ;
        if(vertical)
        {
            if(printHeader)
                std::cout << "        max value\t" <<  std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][0] << "\t" << std::flush ;
            std::cout <<std::endl ;
            if(printHeader)
                std::cout << "        min value\t" <<  std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][1] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "        avg value\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][2] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][3] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][4] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][5] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][6] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][7] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][8] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][9] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][10] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][11] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][12] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][13] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][14] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max von Mises\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][15] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min von Mises\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][16] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][17] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][18] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][19] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][20] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][21] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][22] << "\t" << std::flush ;
            std::cout << std::endl ;
        }
        else
        {
            if(printHeader)
                std::cout << "max value        \tmin value        \tavg value        \tmax sigma11      \tmin sigma11      \tmax sigma22      \tmin sigma22      \tmax sigma12      \tmin sigma12      \tmax epsilon11    \tmin epsilon11    \tmax epsilon22    \tmin epsilon22    \tmax epsilon12    \tmin epsilon12    \tmax von Mises    \tmin von Mises    \taverage sigma11  \taverage sigma22  \taverage sigma12  \taverage epsilon11\taverage epsilon22\taverage epsilon12\t" << std::endl ;
            for(size_t j = 0 ; j < reportValues.size() ; j++)
                for(size_t i = 0 ; i < reportValues[j].size() ; i++)
                    std::cout << reportValues[j][i] << "\t" << std::flush ;
            std::cout << std::endl ;

        }
    }
    else
    {
        Vector x = getDisplacements() ;

        std::cout << "unknowns :" << x.size() << std::endl ;

        int npoints = get3DMesh()->begin()->getBoundingPoints().size() ;

        double volume = 0 ;

        double xavg = 0 ;

        for(auto k = get3DMesh()->begin() ; k != get3DMesh()->end() ; k++)
        {
            if(k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                double ar = k->area() ;
                volume += ar ;
                for(int l = 0 ; l < npoints ; l++)
                {
                    xavg += x[k->getBoundingPoint(l).getId()*3]*ar/npoints ;
                }
            }
        }

        xavg /= volume ;
        reportValues.push_back(Vector(41));
        std::pair<Vector, Vector> stempm = getFieldMinMax(REAL_STRESS_FIELD) ;
        std::pair<Vector, Vector> etempm = getFieldMinMax(STRAIN_FIELD) ;
        std::pair<Vector, Vector> vmm = getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
        Vector stemp = getAverageField(REAL_STRESS_FIELD) ;
        Vector etemp = getAverageField(STRAIN_FIELD) ;
        reportValues.back()[0] = x.max() ;
        reportValues.back()[1] = x.min() ;
        reportValues.back()[2] = xavg    ;
        reportValues.back()[3] = stempm.second[0] ;
        reportValues.back()[4] = stempm.first[0]  ;
        reportValues.back()[5] = stempm.second[1] ;
        reportValues.back()[6] = stempm.first[1]  ;
        reportValues.back()[7] = stempm.second[2] ;
        reportValues.back()[8] = stempm.first[2]   ;
        reportValues.back()[9] = stempm.second[3] ;
        reportValues.back()[10] = stempm.first[3]  ;
        reportValues.back()[11] = stempm.second[4] ;
        reportValues.back()[12] = stempm.first[4]  ;
        reportValues.back()[13] = stempm.second[5] ;
        reportValues.back()[14] = stempm.first[5]   ;

        reportValues.back()[15] =  etempm.second[0] ;
        reportValues.back()[16] =  etempm.first[0]  ;
        reportValues.back()[17] =  etempm.second[1] ;
        reportValues.back()[18] =  etempm.first[1]  ;
        reportValues.back()[19] =  etempm.second[2]  ;
        reportValues.back()[20] =  etempm.first[2]  ;
        reportValues.back()[21] =  etempm.second[3] ;
        reportValues.back()[22] =  etempm.first[3]  ;
        reportValues.back()[23] =  etempm.second[4] ;
        reportValues.back()[24] =  etempm.first[4]  ;
        reportValues.back()[25] =  etempm.second[5]  ;
        reportValues.back()[26] =  etempm.first[5]  ;

        reportValues.back()[27] =  vmm.second[0]    ;
        reportValues.back()[28] =  vmm.first[0]     ;
        reportValues.back()[29] =   stemp[0]        ;
        reportValues.back()[30] =   stemp[1]        ;
        reportValues.back()[31] =   stemp[2]        ;
        reportValues.back()[32] =   stemp[3]        ;
        reportValues.back()[33] =   stemp[4]        ;
        reportValues.back()[34] =   stemp[5]        ;
        reportValues.back()[35] =   etemp[0]        ;
        reportValues.back()[36] =   etemp[1]        ;
        reportValues.back()[37] =   etemp[2]      ;
        reportValues.back()[38] =   etemp[3]        ;
        reportValues.back()[39] =   etemp[4]        ;
        reportValues.back()[40] =   etemp[5]      ;
        std::cout << std::endl ;
        if(vertical)
        {
            if(printHeader)
                std::cout << "        max value\t" <<  std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][0] << "\t" << std::flush ;
            std::cout <<std::endl ;
            if(printHeader)
                std::cout << "        min value\t" <<  std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][1] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "        avg value\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][2] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][3] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][4] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][5] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][6] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][7] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][8] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][9] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][10] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][11] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][12] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      max sigma23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][13] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "      min sigma23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][14] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][15] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][16] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][17] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][18] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][19] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][20] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][21] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][22] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][23] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][24] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max epsilon23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][25] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min epsilon23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][26] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    max von Mises\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][27] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "    min von Mises\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][28] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][29] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][30] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][31] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][32] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][33] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "  average sigma23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][34] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon11\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][35] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon22\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][36] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon33\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][37] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon12\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][38] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon13\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][39] << "\t" << std::flush ;
            std::cout << std::endl ;
            if(printHeader)
                std::cout << "average epsilon23\t" << std::flush ;
            for(size_t i = 0 ; i < reportValues.size() ; i++)
                std::cout << reportValues[i][40] << "\t" << std::flush ;
            std::cout << std::endl ;
        }
        else
        {
            if(printHeader)
                std::cout << "max value        \tmin value        \tavg value        \tmax sigma11      \tmin sigma11      \tmax sigma22      \tmin sigma22      \tmax sigma33      \tmin sigma33      \tmax sigma12      \tmin sigma12      \tmax sigma13      \tmin sigma13      \tmax sigma23      \tmin sigma23      \tmax epsilon11    \tmin epsilon11    \tmax epsilon22    \tmin epsilon22    \tmax epsilon33    \tmin epsilon33    \tmax epsilon12    \tmin epsilon12    \tmax epsilon13    \tmin epsilon13    \tmax epsilon23    \tmin epsilon23    \tmax von Mises    \tmin von Mises    \taverage sigma11  \taverage sigma22  \taverage sigma33  \taverage sigma12  \taverage sigma13  \taverage sigma23  \taverage epsilon11\taverage epsilon22\taverage epsilon33\taverage epsilon12\taverage epsilon13\taverage epsilon23" << std::endl ;
            for(size_t j = 0 ; j < reportValues.size() ; j++)
            {
                for(size_t i = 0 ; i < reportValues[j].size() ; i++)
                    std::cout << reportValues[j][i] << "\t" << std::flush ;
                std::cout << std::endl ;
            }

        }
    }
}

void FeatureTree::printReport(const std::vector<FieldType> & fields, bool vertical)
{

    if(is2D())
    {

        if(vertical)
        {
            for(auto f : fields)
            {
                Vector vals = get2DMesh()->getField(f) ;
                for(size_t i = 0 ; i < vals.size() ;  i++)
                    std::cout << vals[i] << std::endl ;
            }

        }
        else
        {
            for(auto f : fields)
            {
                Vector vals = get2DMesh()->getField(f) ;
                for(size_t i = 0 ; i < vals.size() ;  i++)
                    std::cout << vals[i] <<  "  " << std::flush ;
            }
            std::cout << std::endl ;
        }
    }
    else
    {
        if(vertical)
        {
            for(auto f : fields)
            {
                Vector vals = get3DMesh()->getField(f) ;
                for(size_t i = 0 ; i < vals.size() ;  i++)
                    std::cout << vals[i] << std::endl ;
            }

        }
        else
        {
            for(auto f : fields)
            {
                Vector vals = get3DMesh()->getField(f) ;
                for(size_t i = 0 ; i < vals.size() ;  i++)
                    std::cout << vals[i] <<  "  " << std::flush ;
            }
            std::cout << std::endl ;
        }
    }
}

void FeatureTree::homothety ( double before, double now, double after )
{
    std::valarray<bool> nodes ( getDisplacements().size() /2 ) ;
    if ( nodes.size() == 0 )
    {
        return ;
    }
    nodes = false ;

    if ( dtree->begin()->timePlanes() != 2 )
    {
        return ;
    }
    double stepa = after/now ;
    for ( auto t = dtree->begin() ; t != dtree->end() ; t++ )
    {
        for ( size_t j = 0 ; j < t->getBoundingPoints().size() /2 ; j++ )
        {
            size_t id = t->getBoundingPoint ( j ).getId() ;
            if ( !nodes[id] )
            {
                nodes[id] = true ;
                t->getBoundingPoint ( j ).getX() = t->getBoundingPoint ( t->getBoundingPoints().size() /2+j ).getX() ;
                t->getBoundingPoint ( j ).getY() = t->getBoundingPoint ( t->getBoundingPoints().size() /2+j ).getY() ;
            }
        }
        for ( size_t j = t->getBoundingPoints().size() /2 ; j < t->getBoundingPoints().size() ; j++ )
        {
            size_t id = t->getBoundingPoint ( j ).getId() ;
            if ( !nodes[id] )
            {
                nodes[id] = true ;
                if ( t->getBoundingPoint ( j ).getX() != 0.5 && t->getBoundingPoint ( j ).getX() != -0.5 && t->getBoundingPoint ( j ).getY() != 0.5 && t->getBoundingPoint ( j ).getY() != -0.5 )
                {
                    double r = sqrt ( t->getBoundingPoint ( j ).getX() *t->getBoundingPoint ( j ).getX() + t->getBoundingPoint ( j ).getY() *t->getBoundingPoint ( j ).getY() ) ;
                    if ( r < now )
                    {
                        t->getBoundingPoint ( j ).getX() /= now ;
                        t->getBoundingPoint ( j ).getY() /= now ;
                        t->getBoundingPoint ( j ).getX() *= after ;
                        t->getBoundingPoint ( j ).getY() *= after ;
                    }
                    else
                    {
                        Point inter ( t->getBoundingPoint ( j ).getX(), t->getBoundingPoint ( j ).getY() ) ;
                        dynamic_cast<Rectangle *> ( tree[0] )->project ( &inter ) ;
                        double rl = sqrt ( inter.getX() *inter.getX() + inter.getY() *inter.getY() ) ;
                        t->getBoundingPoint ( j ).getX() *= 1.- ( 1.-stepa ) * ( rl-r ) / ( rl-now ) ;
                        t->getBoundingPoint ( j ).getY() *= 1.- ( 1.-stepa ) * ( rl-r ) / ( rl-now ) ;
                    }
                }
            }
        }
        t->clearElementaryMatrix() ;
    }
//	tri[0]->getBoundingPoint(0).print() ;
}

}
