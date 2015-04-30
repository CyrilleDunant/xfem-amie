// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../utilities/samplingcriterion.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/mazars.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/stiffness.h"
#include "../physics/materials/aggregate_behaviour.cpp"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/paste_behaviour.h"


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
#include <typeinfo>
#include <limits>
#include <time.h>

#define DEBUG



using namespace Amie ;


FeatureTree * featureTree ;


std::vector< double > loads ;
std::vector< double > displacements ;
std::vector< double > loadsx ;
std::vector< double > displacementsx ;
std::vector< double > times ;

Vector fracCrit ( 0 ) ;

Vector b ( 0 ) ;
Vector x ( 0 ) ;
Vector sigma ( 0 ) ;
Vector sigma11 ( 0 ) ;
Vector sigma22 ( 0 ) ;
Vector sigma12 ( 0 ) ;
Vector epsilon ( 0 ) ;
Vector epsilon11 ( 0 ) ;
Vector epsilon22 ( 0 ) ;
Vector epsilon12 ( 0 ) ;
Vector vonMises ( 0 ) ;
Vector angle ( 0 ) ;

// BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition ( SET_ALONG_XI, RIGHT,0 ) ;


double factor = 25 ;
MinimumAngle cri ( M_PI/6. ) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

void step ()
{

    double last_time = 0 ;
    while (featureTree->getCurrentTime() < 48 ) 
    {
        bool go_on = featureTree->step() ;
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD,-1.,1. ) ;
        Vector etemp = featureTree->getAverageField ( STRAIN_FIELD,-1.,1. ) ;

        std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2]*1e6 << std::endl ;

        std::cout << std::endl ;

        if(go_on && std::abs(last_time - featureTree->getCurrentTime()) >.1)
        {
            last_time = featureTree->getCurrentTime() ;
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[0] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[0] );
            times.push_back(featureTree->getCurrentTime());

            std::fstream ldfile  ;
            ldfile.open ( "ldn", std::ios::out ) ;
            for ( size_t j = 0 ; j < loads.size() ; j++ )
            {
                ldfile << times[j] << " " << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] <<  "\n" ;
            }
            ldfile.close();
        }

    }

}




int main ( int argc, char *argv[] )
{
    // Beton
     double k_elas = 40.4e9;
     double nu_elas = 0.2 ;
     Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, true, 
SPACE_TWO_DIMENSIONAL ) ;
     std::vector<std::pair<Matrix, Matrix> > branches ;
     double factor_k = 1.;
         std::vector<double> K_chaine_cp = {5.4e11/factor_k,  3.9e11/factor_k, 2.02e11/factor_k,5.1e10/factor_k} ;
     for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
     {
         double tau = 5.*std::pow(10., (double) i - 2 ); std::cerr << "TAU "<< tau << std::endl ;
         Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  true,  SPACE_TWO_DIMENSIONAL )  ; 
         Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
         branches.push_back(std::make_pair(K_i, Am_i)) ;
     }
    Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    
    double cstrain = -2.e-3; 
    double cstress = -38.0e6; 
    planeType pt = PLANE_STRESS;

    Sample samplef(0.3, 0.6,  0.15, 0.3) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;
    ViscoelasticityAndFracture * pasterupt = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1e-9, 1.) );
    
     //ELAS+DAMAGE
//     StiffnessAndFracture * spasterupt = new StiffnessAndFracture(E_cp_elas, new NonLocalMCFT(-40e6,40e9,1.), new FiberBasedIsotropicLinearDamage(0.001, 1.));
    ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas,new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1e-9,1.));
    //StiffnessAndFracture * spasterupt = new StiffnessAndFracture(E_cp_elas,new NonLocalMazars(1.0e-4, k_elas, nu_elas, 100, cstress , cstrain, 1., pt ), new IsotropicLinearDamage());
    //ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, new NonLocalMazars(1.0e-4, k_elas, nu_elas, 100, cstress , cstrain, 1., pt ), new IsotropicLinearDamage());
    samplef.setBehaviour ( pasterupt  ) ;

    Function loadfunc = Function("t 12 /")         *f_range("t", 0., 3) +
    Function("0.25 t 3 - 48 / -")*f_range("t", 3, 6) +
    Function("t 6 - 72 / 0.1866666 +")    *f_range("t", 6, 72) ;
    loadfunc *= 0.0002 ;
    
//     for(double t = 0 ; t < 72  ;t += .1)
//     {
//         std::cout << t << "  " << VirtualMachine().eval(loadfunc,0,0,0,t) << std::endl ;
//     }
//     exit(0)
    BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, loadfunc) ;
    F.addBoundaryCondition ( loadr );

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM_AFTER ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.setDeltaTime(.1);
    F.setOrder(LINEAR_TIME_LINEAR) ;

    F.setMaxIterationsPerStep ( 50 );

    step () ;

//    F.getAssembly()->print() ;


    return 0 ;
}
