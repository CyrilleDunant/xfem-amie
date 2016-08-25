// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/mazars.h"
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"
#include "../polynomial/vm_function_extra.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/void_form.h"
#include "../features/sample.h"
#include "../utilities/itoa.h"
#include "../utilities/parser/command_line_parser.h"
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
double platewidth = 0.005 ;
double rebarDiametre = 0.025 ; //sqrt(506e-6);//sqrt( 4.*0.000506/M_PI ) ;
double rebarEndCover = 0.047 ;

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
double mindelta = 1e-6 ;


// Rectangle bcbox ( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) ;
// GeometryDefinedBoundaryCondition selfload ( SET_VOLUMIC_STRESS_ETA, &bcbox , -9025.2 ) ;
// GeometryDefinedBoundaryCondition shrinkagey ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
// GeometryDefinedBoundaryCondition shrinkagex ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;

void step(FeatureTree * featureTree, double supportLever, double sampleHeight, BoundaryCondition * load, double loadingSpeed )
{
    for ( size_t v = 0 ;  ; v++ )
    {

        featureTree->stepToCheckPoint(10, 1e-4) ;

        int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;

        double delta = 0;
        double deltacount = 0;
        double volume = 0 ;

        for ( auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++ )
        {
            if ( k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                volume += k->area() ;
                for ( int p = npoints/2 ; p < npoints ; p++ )
                {

                    if (  k->getBoundingPoint ( p ).getX()  < .01 && std::abs(k->getBoundingPoint ( p ).getY()-sampleHeight*.5) < 0.01 )
                    {
                        deltacount++ ;
                        Vector dummy(k->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
                        k->getState().getField(DISPLACEMENT_FIELD, k->getBoundingPoint ( p ), dummy, false); 
                        delta += dummy[1] ;
                    }
                }
            }
        }


        delta /= deltacount ;

        if( v%10 == 0)
        {
            Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD ) ;

            std::cout << std::endl ;
            std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
            std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
            std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
            std::cout << std::endl ;
            
            std::fstream ldfile ( "ldn_3", std::ios::out | std::ios::app)  ;
            ldfile << 1000.* ( VirtualMachine().eval(load->getData()*load->getScale(), 0,0,0,featureTree->getInitialTime()) ) << "   " 
                   << stemp[1]/1e6 << "   " 
                   << featureTree->averageDamage << "   " 
                   << delta*1000. << "   " 
                   << featureTree->get2DMesh()->begin()->getBoundingPoint(npoints-1).getT() << "   " 
                   << featureTree->get2DMesh()->begin()->getBoundingPoint(0).getT()<<"\n" ;
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
        
        load->setData(featureTree->getCurrentTime()*loadingSpeed);
    }
}


int main ( int argc, char *argv[] )
{

  
//     Function base("1 x - y -") ;
//     Function mone("-1") ;
//     base.setNumberOfDerivatives(2) ;
//     base.setDerivative(XI, mone);
//     base.setDerivative(ETA, mone);
//     double p = 0 ;
//     VirtualMachine vm ;
//     vm.print(mone);
//     for(size_t i = 1 ; i < 20000 ; i++)
//     {
//       p += vm.deval(base, XI, 0,0,0) ;
//     }
//     std::cout << p << std::endl ;
//     p = 0 ;
//     for(size_t i = 1 ; i < 20000 ; i++)
//     {
//       p += i ;
//     }
//     std::cout << p << std::endl ;
//     
//     exit(0) ;
    
    // Beton

    CommandLineParser parser("Make a tri-point bending test on an homogeneous concrete sample") ;
    parser.addArgument("length", .39*.5, "length of the sample (default 3.9)") ;
    parser.addArgument("height", .12*.5, "height of the sample (default 1.2)") ;
    parser.addArgument("speed", -0.00144, "loading speed (default -0.144 (1 mm / 60000 s))") ;
    parser.parseCommandLine(argc, argv) ;

    double sampleLength     = parser.getNumeralArgument( "length") ;
    double sampleHeight     = parser.getNumeralArgument( "height") ;
    double loadingSpeed     = parser.getNumeralArgument( "speed") ;
    double supportLever     = sampleLength*.5-.02 ;
    double halfSampleOffset = sampleLength*.25 ;

    double E_steel = 200e9 ;
    double nu_steel = 0.01 ;
    double cstrain = -2.e-3; 
    double cstress = -38.0e6; 
    planeType pt = PLANE_STRESS;
    double k_elas = 40.4e9;
    double nu_elas = 0.3 ;
    Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    std::vector<std::pair<Matrix, Matrix> > branches ;
    std::vector<double> K_chaine_cp = { 5.4e11,  3.9e11, 2.02e11, 5.1e10 } ;
    for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
    {
        double tau = 5.*std::pow(10., (double) i - 2 ); std::cerr << "TAU "<< tau << std::endl ;
        Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON )  ; 
        Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
        branches.push_back(std::make_pair(K_i, Am_i)) ;
    }    
    
    FractureCriterion * mcft = new NonLocalSpaceTimeMCFT(cstress,k_elas, 0.032,UPPER_BOUND) ;
    FractureCriterion * mazar = new NonLocalSpaceTimeMazars(4.52e-5, k_elas, nu_elas, 10, cstress , cstrain, 0.032, pt ) ;
    DamageModel * linear = new SpaceTimeIsotropicLinearDamage(1.0) ;
 
    Sample sample ( nullptr, sampleLength*.5, sampleHeight, halfSampleOffset, 0 ) ;   
    
    sample.setBehaviour (new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, mazar->getCopy(), linear->getCopy() )); 
//     sample.setBehaviour (new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, mcft->getCopy(), linear->getCopy())); 

    FeatureTree F ( &sample ) ;


    Triangle fineZone ( Point ( 0.,sampleHeight*.7 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.1, -sampleHeight*.5 ) ) ;
    Triangle finerZone ( Point ( 0.,sampleHeight*.7 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.05, -sampleHeight*.5 ) ) ;
    F.addRefinementZone ( &fineZone );
    F.addRefinementZone ( &finerZone );

    parser.setFeatureTree(&F) ;

    F.setMaxIterationsPerStep ( 1600 );


    F.addPoint ( new Point ( supportLever, -sampleHeight*.5 ) ) ;
    F.addPoint ( new Point ( platewidth, sampleHeight*.5 ) ) ;
    
    BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, TOP_AFTER, -platewidth*2., platewidth*2., -10, 10, 0 ) ;
    
    F.addBoundaryCondition ( load ) ;
    load->setActive(true);
    
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxNearestNodeDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM_AFTER, Point ( supportLever, -sampleHeight*.5 ) ) ) ;
    F.setOrder ( LINEAR_TIME_LINEAR ) ;
    F.setDeltaTime(0.002);
    F.setMinDeltaTime(1e-12);
    F.setSamplingNumber(4);

    step(&F, supportLever, sampleHeight, load, loadingSpeed) ;

    return 0 ;
}
