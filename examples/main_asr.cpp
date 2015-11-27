
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/stiffness.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/enrichmentmanagers/gelmanager.h"
#include "../features/expansiveZone.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../utilities/writer/triangle_writer.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

FeatureTree *featureTree ;

double placed_area = 0 ;

double basesize = 0.07 ;

Rectangle baseGeometry( basesize, basesize, 0, 0 ) ;

std::vector<std::pair<ExpansiveZone *, Feature *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;
std::vector<Vector> asrStress ;


double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

int totit = 1 ;

double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;
GelManager * gelManager ;


void step(std::vector<Feature *> & inclusions, std::vector<Feature *> & blocks)
{
    int nsteps = 1600;

    int intermediateCount=0 ;
    featureTree->setMaxIterationsPerStep( 3200 ) ;

    std::vector<Feature *> inclusionsAndBlocks =inclusions ;
    if(!blocks.empty())
        inclusionsAndBlocks.insert(inclusionsAndBlocks.end(), blocks.begin(), blocks.end());
    for( int i = 0 ; i < nsteps ; i++ )
    {
        std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;
        bool go_on = featureTree->step() ;
//         if(i%5 != 0)
//             continue ;
        std::string filename( "triangles" ) ;
        if( !go_on )
        {
            filename = std::string( "intermediate-triangles" ) ;
            intermediateCount++ ;
        }
        else
            intermediateCount = 0 ;

        if( intermediateCount > 4)
        {
            intermediateCount = 0 ;
            go_on = true ;
        }

        if( !featureTree->solverConverged() )
            filename = std::string( "failed-triangles" ) ;

        filename.append( itoa( totit++, 10 ) ) ;
        std::cout << filename << std::endl ;

        TriangleWriter writer(filename.c_str(), featureTree) ;
        writer.getField(TWFT_PRINCIPAL_STRESS ) ;
        writer.getField(TWFT_PRINCIPAL_STRAIN ) ;
        writer.getField( TWFT_CRITERION ) ;
        writer.getField( TWFT_STIFFNESS ) ;
        writer.getField( TWFT_DAMAGE ) ;

        writer.write() ;

        if( go_on )
        {
            std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax(REAL_STRESS_FIELD) ;
            std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax(STRAIN_FIELD) ;
            Vector stemp = featureTree->getAverageField(REAL_STRESS_FIELD) ;
            Vector etemp = featureTree->getAverageField(STRAIN_FIELD) ;
            std::vector<double> macro_strain = featureTree->getMedianMacroscopicStrain(&baseGeometry) ;
            expansion_reaction.push_back( std::make_pair( gelManager->getReactedFraction(), macro_strain[0]) ) ;
            expansion_stress_xx.push_back( std::make_pair( etemp[0], stemp[0] ) ) ;
            expansion_stress_yy.push_back( std::make_pair( etemp[1], stemp[1] ) ) ;
            apparent_extension.push_back( std::make_pair(macro_strain[0], macro_strain[1]) ) ;
            cracked_volume.push_back( featureTree->averageDamageInFeatures(inclusions) ) ;
            damaged_volume.push_back( featureTree->averageDamageNotInFeatures(inclusionsAndBlocks) ) ;
            Vector buffer(3) ;
            featureTree->averageFieldInFeatures(REAL_STRESS_FIELD, blocks, buffer) ;
            asrStress.push_back(buffer) ;
        }

        std::cout << "reaction" << "   "<< "eps xx" << "\t" << "eps yy" << "\t" << "sig xx" << "\t"<< "sig yy" << "\t" << "   dx "  << "\t"<< "   dy "  << "\t"
                  << "cracks"  << "\t"<< "damage"  << "\t"<< "asr Stress x"  << "\t"<< "asr Stress y" << std::endl ;
        for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
            std::cout << expansion_reaction[i].first << "\t"<< expansion_stress_xx[i].first << "\t"<< expansion_stress_yy[i].first << "\t"<< expansion_stress_xx[i].second << "\t"<< expansion_stress_yy[i].second << "\t"<< apparent_extension[i].first  << "\t"<< apparent_extension[i].second  << "\t"
                      << cracked_volume[i]  << "\t"<< damaged_volume[i]  << "\t" << asrStress[i][0]  << "   "<< asrStress[i][1] << std::endl ;

    }
}

int main( int argc, char *argv[] )
{
    double restraintDepth = 0.01 ;

    int nzones   = atof( argv[1] )  ;
    double dmax  = atof( argv[2] )  ;
    double fact  = atof( argv[3] )*2.*restraintDepth/basesize  ;
    double fact0 = atof( argv[4] )*2.*restraintDepth/basesize  ;

    

    Sample sample( nullptr, basesize + restraintDepth, basesize + restraintDepth, 0, 0 ) ;

    FeatureTree F( &sample ) ;
    featureTree = &F ;

    double itzSize = 0.00005;
    int inclusionNumber = 1500 ; 2500 ; 1 ; 

    Rectangle placeGeometry( basesize, basesize, 0, 0 ) ;

    std::vector<Feature *> feats  = PSDGenerator::get2DConcrete(&F, new /*ElasticOnly*/AggregateBehaviour(),  inclusionNumber, dmax*0.5, itzSize, new PSDBolomeA(), nullptr, 100000, 0.8, &placeGeometry) ;

    double volume = 0 ;

    for( size_t i = 0 ; i < feats.size() ; i++ )
        volume += feats[i]->area() ;

    if( !feats.empty() )
        std::cout << "n = " << feats.size() 
                  << ", largest r = " <<  feats.front()->getRadius() - itzSize
                  << ", smallest r =" <<  feats.back()->getRadius()  - itzSize 
                  << ", ratio = "     << (feats.front()->getRadius() - itzSize) / (feats.back()->getRadius() - itzSize) <<std::endl ;

    sample.setBehaviour( new /*ElasticOnly*/PasteBehaviour()) ;

    std::vector<Feature *> blocks ;

    Sample *voidtop = new Sample( nullptr, restraintDepth * .5, 
                                           restraintDepth * .5, 
                                           sample.getCenter().getX() - basesize*.5 - restraintDepth * .25, 
                                           sample.getCenter().getY() + basesize*.5 + restraintDepth * .25) ;
//     voidtop->isVirtualFeature = true ;
    voidtop->setBehaviour( new VoidForm() );
    F.addFeature( &sample, voidtop );
//     F.setSamplingFactor(voidtop, 8.);

    Sample *voidbottom = new Sample( nullptr, restraintDepth * .5, 
                                              restraintDepth * .5, 
                                              sample.getCenter().getX() - basesize*.5 - restraintDepth * .25, 
                                              sample.getCenter().getY() - basesize*.5 - restraintDepth * .25 ) ;
//     voidbottom->isVirtualFeature = true ;
    voidbottom->setBehaviour( new VoidForm() );
    F.addFeature( &sample, voidbottom );
//     F.setSamplingFactor(voidbottom, 8.);

    Sample *voidleft = new Sample( nullptr, restraintDepth * .5, 
                                            restraintDepth * .5, 
                                            sample.getCenter().getX() + basesize*.5 + restraintDepth * .25, 
                                            sample.getCenter().getY() + basesize*.5 + restraintDepth * .25 ) ;
//     voidleft->isVirtualFeature = true ;
    voidleft->setBehaviour( new VoidForm() );
    F.addFeature( &sample, voidleft );
//     F.setSamplingFactor(voidleft, 8.);

    Sample *voidright = new Sample( nullptr, restraintDepth * .5, 
                                             restraintDepth * .5, 
                                             sample.getCenter().getX() + basesize*.5 + restraintDepth * .25, 
                                             sample.getCenter().getY() - basesize*.5 - restraintDepth * .25 ) ;
//     voidright->isVirtualFeature = true ;
    voidright->setBehaviour( new VoidForm() );
    F.addFeature( &sample, voidright );
//     F.setSamplingFactor(voidright, 8.);

    //width are  1100000000 3340000000  4360000000      done: 
    //length are 1220000000 2180000000  3400000000      next: 

    Sample *blocktop = new Sample( nullptr, basesize, restraintDepth * .5, sample.getCenter().getX(), sample.getCenter().getY() + basesize*.5 + restraintDepth * .25 ) ;
    double Emin = std::max(std::max(fact, fact0)*sqrt(POINT_TOLERANCE),1e5) ;
    fact = std::max(fact, Emin) ;
    fact0 = std::max(fact0, Emin) ;

    blocktop->setBehaviour( new OrthotropicStiffness(fact, fact0, (fact0+fact)*.5,  0., 0.))  ;
    blocks.push_back(blocktop);


    F.addFeature( &sample, blocktop );

    Sample *blockbottom = new Sample( nullptr, basesize, 
                                               restraintDepth * .5, 
                                               sample.getCenter().getX(), 
                                               sample.getCenter().getY() - basesize*.5 - restraintDepth * .25 ) ;
    blockbottom->setBehaviour(new VoidForm()) ;
    blocks.push_back(blockbottom);

    F.addFeature( &sample, blockbottom );

    Sample *blockleft = new Sample( nullptr, restraintDepth * .5, 
                                             basesize, 
                                             sample.getCenter().getX() - ( basesize )*.5 - restraintDepth * .25, 
                                             sample.getCenter().getY() ) ;

    blockleft->setBehaviour(new VoidForm()) ;
    blocks.push_back(blockleft);

    F.addFeature( &sample, blockleft );

    Sample *blockright = new Sample( nullptr, restraintDepth * .5, 
                                              basesize, 
                                              sample.getCenter().getX() + ( basesize )*.5 + restraintDepth * .25, 
                                              sample.getCenter().getY() ) ;
    blockright->setBehaviour( new OrthotropicStiffness(fact, fact0, (fact0+fact)*.5,  0., 0.)) ;
    blocks.push_back(blockright);

    F.addFeature( &sample, blockright );

    for( size_t i = 0 ; i < feats.size() ; i++ )
        placed_area += feats[i]->area() ;

    std::cout << ", filling = " << placed_area / baseGeometry.area() * 100. << "%" << std::endl ;

    if( !feats.empty() )
    {
        std::cout << "largest inclusion with r = " << ( *feats.begin() )->getRadius() << std::endl ;
        std::cout << "smallest inclusion with r = " << ( *feats.rbegin() )->getRadius() << std::endl ;
        std::cout << "placed area = " <<  placed_area << std::endl ;
    }

    gelManager = new GelManager(&F, nzones/baseGeometry.area(), feats, 0.5, 1e-5) ;
    F.addManager(gelManager) ;
    F.setSamplingNumber( 60 ) ;

    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)  ) ;
    F.addBoundaryCondition( new BoundingBoxAndRestrictionDefinedBoundaryCondition(FIX_ALONG_XI, RIGHT,   sample.getCenter().getX()-basesize*.5, 
                                                                                                         sample.getCenter().getX()+basesize*.5+restraintDepth*.5,
                                                                                                         sample.getCenter().getY()-basesize*.5, 
                                                                                                         sample.getCenter().getY()+basesize*.5)  ) ;    
    F.addBoundaryCondition( new BoundingBoxAndRestrictionDefinedBoundaryCondition(FIX_ALONG_ETA, TOP,    sample.getCenter().getX()-basesize*.5, 
                                                                                                         sample.getCenter().getX()+basesize*.5,
                                                                                                         sample.getCenter().getY()-basesize*.5, 
                                                                                                         sample.getCenter().getY()+basesize*.5+restraintDepth*.5)  ) ;
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM) ) ;
    

    F.setOrder( LINEAR ) ;

    step( feats, blocks) ;

    return 0 ;
}
