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

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
double mindelta = 1e-6 ;


// Rectangle bcbox ( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) ;
// GeometryDefinedBoundaryCondition selfload ( SET_VOLUMIC_STRESS_ETA, &bcbox , -9025.2 ) ;
// GeometryDefinedBoundaryCondition shrinkagey ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
// GeometryDefinedBoundaryCondition shrinkagex ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;

void step(FeatureTree * featureTree, double supportLever, double sampleHeight, BoundaryCondition * load )
{


    for ( size_t v = 0 ;  ; v++ )
    {

        featureTree->stepToCheckPoint(5, 1e-4) ;

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

                    if (  k->getBoundingPoint ( p ).getX()  < .01 && std::abs(k->getBoundingPoint ( p ).getY()-sampleHeight*.5 - plateHeight) < 0.01)
                    {
                        deltacount++ ;
                        Vector dummy(2) ;
                        k->getState().getField(DISPLACEMENT_FIELD, k->getBoundingPoint ( p ), dummy, false); 
                        delta +=dummy[1] ;
                    }
                }

            }
        }


        delta /= deltacount ;



        if( v%10 == 0)
        {
            Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD, -1, 1 ) ;

            std::cout << std::endl ;
            std::cout << "average sigma11 : " << stemp[0] << std::endl ;
            std::cout << "average sigma22 : " << stemp[1] << std::endl ;
            std::cout << "average sigma12 : " << stemp[2] << std::endl ;
            std::cout << std::endl ;
            
            std::fstream ldfile ( "ldn_2", std::ios::out | std::ios::app)  ;
            ldfile << 1000.* ( VirtualMachine().eval(load->getDataFunction()*load->getScale(), 0,0,0,featureTree->getInitialTime()) ) << "   " << stemp[1]/1000. << "   " << featureTree->averageDamage << "   " << delta*1000. << "   " << featureTree->get2DMesh()->begin()->getBoundingPoint(npoints-1).getT() << "   " << featureTree->get2DMesh()->begin()->getBoundingPoint(0).getT()<<"\n" ;
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
}


int main ( int argc, char *argv[] )
{

    // Beton

    double samplingNumber   = atof ( argv[1] ) ;
    double sampleLength     = atof ( argv[2] ) ;
    double sampleHeight     = atof ( argv[3] ) ;
    double supportLever     = sampleLength*.5-.250 ;
    double halfSampleOffset = sampleLength*.25 ;

    double E_steel = 200e9 ;
    double nu_steel = 0.01 ;
    double cstrain = -2.e-3; 
    double cstress = -38.0e6; 
    planeType pt = PLANE_STRESS;
    double k_elas = 40.4e9;
    double nu_elas = 0.3 ;
    Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
    std::vector<std::pair<Matrix, Matrix> > branches ;
    std::vector<double> K_chaine_cp = { 5.4e11,  3.9e11, 2.02e11, 5.1e10 } ;
    for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
    {
        double tau = 5.*std::pow(10., (double) i - 2 ); std::cerr << "TAU "<< tau << std::endl ;
        Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  true,  SPACE_TWO_DIMENSIONAL )  ; 
        Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
        branches.push_back(std::make_pair(K_i, Am_i)) ;
    }    
    
    FractureCriterion * mcft = new NonLocalSpaceTimeMCFT(cstress,k_elas, 0.064,UPPER_BOUND,MIRROR_Y) ;
    FractureCriterion * mazar = new NonLocalSpaceTimeMazars(4.52e-5, k_elas, nu_elas, 10, cstress , cstrain, 0.064, pt ) ;
    DamageModel * linear = new SpaceTimeIsotropicLinearDamage(1.0) ;
 
    Matrix m0_steel = Tensor::cauchyGreen ( std::make_pair ( E_steel,nu_steel ), true, SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;

    Sample sample ( nullptr, sampleLength*.5, sampleHeight+2.*plateHeight, halfSampleOffset, 0 ) ;   

    Sample topsupport ( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;

    Sample toprightvoid ( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth ) *.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
    toprightvoid.setBehaviour ( new VoidForm() ) ;

    Sample baseright ( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;


    Sample bottomcentervoid ( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 ) *.5, -sampleHeight*.5 - plateHeight*.5 ) ;
    bottomcentervoid.setBehaviour ( new VoidForm() ) ;

    Sample rightbottomvoid ( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 ) *.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
    rightbottomvoid.setBehaviour ( new VoidForm() ) ;
    
    
//     sample.setBehaviour (new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, mcft->getCopy(), linear->getCopy() )); 
    sample.setBehaviour (new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, mcft->getCopy(), linear->getCopy() )); 
    topsupport.setBehaviour ( new Viscoelasticity ( PURE_ELASTICITY, m0_steel ) ) ;
    baseright.setBehaviour ( new Viscoelasticity ( PURE_ELASTICITY, m0_steel ) ) ;
    
    //     baseright.setBehaviour ( new Stiffness ( m0_steel ) ) ;
//     topsupport.setBehaviour ( new Stiffness ( m0_steel ) ) ;
    //     sample.setBehaviour ( new ConcreteBehaviour ( k_elas, nu_elas, compressionCrit,PLANE_STRESS, UPPER_BOUND, SPACE_TWO_DIMENSIONAL,MIRROR_Y ) ) ;
    
    FeatureTree F ( &sample ) ;

    F.addFeature ( &sample,&baseright ) ;
    F.addFeature ( &baseright,&bottomcentervoid );
    F.addFeature ( &baseright,&rightbottomvoid ) ;
    F.addFeature ( &sample, &topsupport);
    F.addFeature ( &topsupport, &toprightvoid);

    Triangle fineZone ( Point ( 0.,sampleHeight*.7 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.1, -sampleHeight*.5 ) ) ;
    Triangle finerZone ( Point ( 0.,sampleHeight*.7 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.05, -sampleHeight*.5 ) ) ;
    F.addRefinementZone ( &fineZone );
    F.addRefinementZone ( &finerZone );



    F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
    F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;
    F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
    F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;

    F.setSamplingNumber ( samplingNumber ) ;

    F.setSamplingRestriction ( 0 );

    F.setMaxIterationsPerStep ( 1600 );


    F.addPoint ( new Point ( supportLever, -sampleHeight*.5-plateHeight ) ) ;
    F.addPoint ( new Point ( platewidth, sampleHeight*.5 ) ) ;
    
    Function loadfunc = Function("t");
    loadfunc *= -0.00016 ; /*-5e5 ;*/
    BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, TOP_AFTER, -platewidth*2., platewidth*2., -10, 10, loadfunc ) ;
    
    F.addBoundaryCondition ( load ) ;
    load->setActive(true);
    
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxNearestNodeDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM_AFTER, Point ( supportLever, -sampleHeight*.5-plateHeight ) ) ) ;
    F.setOrder ( LINEAR_TIME_LINEAR ) ;
    F.setDeltaTime(mindelta*1e3);
    F.setMinDeltaTime(1e-9);

    step(&F, supportLever, sampleHeight, load) ;

    return 0 ;
}
