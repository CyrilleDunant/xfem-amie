// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/void_form.h"
#include "../features/sample.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"

#include <fstream>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <limits>
#include <time.h>

using namespace Amie ;




// double sampleLength = 3.9 ; //5.5 ;
// double sampleHeight = 1.2 ;
// double supportLever = 1.7 ;//2.5 ;
double supportMidPointToEndClearance = 0.25 ;
double platewidth = 0.15 ;
double plateHeight = 0.051 ;
double rebarDiametre = 0.025 ; //sqrt(506e-6);//sqrt( 4.*0.000506/M_PI ) ;
double rebarEndCover = 0.047 ;

std::vector< double > loads ;
std::vector< double > deltas ;
std::vector< double > displacements ;
std::vector< double > damages ;
Vector fracCrit ( 0 ) ;

Vector x ( 0 ) ;


MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, TOP, -platewidth, platewidth, -10, 10, 0. ) ;

// Rectangle bcbox ( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) ;
// GeometryDefinedBoundaryCondition selfload ( SET_VOLUMIC_STRESS_ETA, &bcbox , -9025.2 ) ;
// GeometryDefinedBoundaryCondition shrinkagey ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
// GeometryDefinedBoundaryCondition shrinkagex ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;

void step(FeatureTree * featureTree, double supportLever, double sampleHeight )
{
    size_t nsteps = 600*4 ; //16*10;
    double delta_d = 0.0175e-3 ;

    double target =  load->getData() ;            
    target -= delta_d*20. ;  
    load->setData ( target ) ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {
//         bool go_on = false ;
//         featureTree->stepToCheckPoint(10, 1e-6) ;
//         std::cout << load->getScale() << "   "<< featureTree->averageDamage<< std::endl ;
//         if (load->getScale() > .85 || v%30 == 0)
//         {
//             go_on = true ;
//             if(load->getScale() > .85)
//             {
//                 target -= delta_d ;   
//                 load->setData ( target ) ;
//             }
//         }
//         else
//         {
//             continue ;
//         }

        bool go_on = featureTree->stepToCheckPoint(40, 1e-8) ;
        
        x.resize ( featureTree->getDisplacements ( -1, false ).size() ) ;
        x = featureTree->getDisplacements ( -1, false ) ;

        int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;

        double delta = 0;
        double deltacount = 0;
        double volume = 0 ;

        for ( auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++ )
        {
            if ( k->getBehaviour()->type != VOID_BEHAVIOUR )
            {

                volume += k->area() ;
                for ( int p = 0 ; p < npoints ; p++ )
                {

                    if ( dist ( Point ( 0, sampleHeight*.5 + plateHeight ), k->getBoundingPoint ( p ) ) < .01 )
                    {
                        deltacount++ ;
                        delta += x[k->getBoundingPoint ( p ).getId() * 2] ;
                    }
                }

            }
        }


        delta /= deltacount ;

        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD ) ;

        std::cout << std::endl ;
        std::cout << "average sigma11 : " << stemp[0] << std::endl ;
        std::cout << "average sigma22 : " << stemp[1] << std::endl ;
        std::cout << "average sigma12 : " << stemp[2] << std::endl ;
        std::cout << std::endl ;
        std::fstream ldfile ( "ldn", std::ios::out )  ;

        
        if ( go_on )
        {
            displacements.push_back ( 1000.* ( load->getData()*load->getScale() ) );
            loads.push_back ( stemp[1]/1000. );
            deltas.push_back ( delta );
            damages.push_back ( featureTree->averageDamage );

//             writerc.reset ( featureTree ) ;
//             writerc.getField ( TWFT_PRINCIPAL_STRESS ) ;
//             writerc.getField ( TWFT_PRINCIPAL_STRAIN ) ;
//             writerc.getField ( TWFT_CRITERION ) ;
//             writerc.getField ( TWFT_STIFFNESS_X ) ;
//             writerc.getField ( TWFT_STIFFNESS_Y ) ;
//             writerc.getField ( TWFT_DAMAGE ) ;
//             writerc.append() ;
//             writerc.writeSvg ( 0. ) ;        
            for ( size_t j = 0 ; j < loads.size() ; j++ )
            {
                ldfile << displacements[j] << "   " << loads[j] << "   " << damages[j] << "   " << deltas[j] << "\n" ;
            }
        }
        else
        {
            for ( size_t j = 0 ; j < loads.size() ; j++ )
            {
                ldfile << displacements[j] << "   " << loads[j] << "   " << damages[j] << "   " << deltas[j] << "\n" ;

            }
            ldfile <<  1000.* ( load->getData()*load->getScale() ) << "   " << stemp[1]/1000. << "   " << featureTree->averageDamage << "   " << delta << "\n" ;
        }
        
        ldfile.close();

        writer.reset ( featureTree ) ;
        writer.getField ( TWFT_PRINCIPAL_STRESS ) ;
        writer.getField ( TWFT_PRINCIPAL_STRAIN ) ;
        writer.getField ( TWFT_CRITERION ) ;
        writer.getField ( TWFT_STIFFNESS_X ) ;
        writer.getField ( TWFT_STIFFNESS_Y ) ;
        writer.getField ( TWFT_DAMAGE ) ;
        writer.append() ;

    }
}


int main ( int argc, char *argv[] )
{

    double  samplingNumber = atof ( argv[1] ) ;
    double sampleLength = atof ( argv[2] ) ;
    double sampleHeight = atof ( argv[3] ) ;
    double supportLever = sampleLength*.5-.250 ;
    double halfSampleOffset = sampleLength*.25 ;
 
    std::cerr << sampleLength << "  " << supportLever << std::endl ;

    double compressionCrit = -34.2e6 ;

    double E_steel = 200e9 ;
    double nu_steel = 0.01 ;
    double nu = 0.3 ;
    double E_paste = 37e9 ;

 
    Matrix m0_steel = Tensor::cauchyGreen ( std::make_pair ( E_steel,nu_steel ), true, SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;

    Sample sample ( nullptr, sampleLength*.5, sampleHeight+2.*plateHeight, halfSampleOffset, 0 ) ; 
    sample.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit,PLANE_STRESS, UPPER_BOUND, SPACE_TWO_DIMENSIONAL,MIRROR_Y ) ) ;

    Sample topsupport ( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
    topsupport.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample toprightvoid ( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth ) *.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
    toprightvoid.setBehaviour ( new VoidForm() ) ;

    Sample baseright ( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
    baseright.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample bottomcentervoid ( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 ) *.5, -sampleHeight*.5 - plateHeight*.5 ) ;
    bottomcentervoid.setBehaviour ( new VoidForm() ) ;

    Sample rightbottomvoid ( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 ) *.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
    rightbottomvoid.setBehaviour ( new VoidForm() ) ;

    FeatureTree F ( &sample ) ;

    F.addFeature ( &sample,&baseright ) ;
    F.addFeature ( &baseright,&bottomcentervoid );
    F.addFeature ( &baseright,&rightbottomvoid ) ;
    F.addFeature ( &sample, &topsupport);
    F.addFeature ( &topsupport, &toprightvoid);

    Triangle fineZone ( Point ( 0.,sampleHeight*.5 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.15, -sampleHeight*.5 ) ) ;
    Triangle finerZone ( Point ( 0.,sampleHeight*.5 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.075, -sampleHeight*.5 ) ) ;
    F.addRefinementZone ( &fineZone );
     F.addRefinementZone ( &finerZone );



    F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
    F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;
    F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
    F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;

    F.setSamplingNumber ( samplingNumber ) ;

    F.setSamplingRestriction ( SAMPLE_NO_RESTRICTION );

    F.setMaxIterationsPerStep ( 1600 );


    F.addPoint ( new Point ( supportLever, -sampleHeight*.5-plateHeight ) ) ;
    F.addPoint ( new Point ( platewidth, sampleHeight*.5 ) ) ;
    F.addBoundaryCondition ( load ) ;
    load->setActive(true);
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT ) ) ;
    F.addBoundaryCondition ( new BoundingBoxNearestNodeDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM, Point ( supportLever, -sampleHeight*.5-plateHeight ) ) ) ;
    F.setOrder ( LINEAR ) ;

    step(&F, supportLever, sampleHeight) ;

    return 0 ;
}
