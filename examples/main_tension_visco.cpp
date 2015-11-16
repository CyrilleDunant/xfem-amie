// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../physics/fracturecriteria/mcft.h"
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../physics/fracturecriteria/mazars.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
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
std::vector< double > damage;

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

void step (const Function & loadfunc)
{

    displacements.push_back (0 );
    displacementsx.push_back ( 0 );
    loads.push_back ( 0 );
    loadsx.push_back ( 0 );
    times.push_back(0);
    damage.push_back(0);
    double last_time = 0 ;
    int itcounter = 0 ;
    while ( true ) 
    {
        bool go_on = featureTree->stepToCheckPoint(1, 1e-4) ;
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD,-1.,1. ) ;
        Vector etemp = featureTree->getAverageField ( STRAIN_FIELD,-1.,1. ) ;
        Vector dtemp = featureTree->getAverageField (SCALAR_DAMAGE_FIELD,-1.,1.) ;
        itcounter++ ;

        if(go_on/* && std::abs(sqrt(etemp[1]*etemp[1]/16e-10+stemp[1]*stemp[1]/2.25e12) - sqrt(displacements.back()*displacements.back()/16e-10+loads.back()*loads.back()/2.25e12)) > .01*/)
        {
            itcounter = 0 ;
            last_time = featureTree->getCurrentTime() ;
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[0] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[0] );
            times.push_back(featureTree->getCurrentTime());
            damage.push_back(dtemp[0]);

            std::fstream ldfile  ;
            ldfile.open ( "ldn", std::ios::out | std::ios::app) ;

            ldfile << times.back()           << "   " 
                   << displacements.back()   << "   " 
                   << loads.back()           << "   " 
                   << displacementsx.back()  << "   " 
                   << loadsx.back()          << "   " 
                   << damage.back()          << "   "
                   << VirtualMachine().eval(loadfunc, 0,0,0,featureTree->getCurrentTime()) << "\n" ;

            ldfile.close();
            if(displacements.back()*1000. > 0.09)
                exit(0) ;
        }

    }

}

int main ( int argc, char *argv[] )
{
    // Beton
     double k_elas = 40.4e9;
     double nu_elas = 0.2 ;
     Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, true,  SPACE_TWO_DIMENSIONAL , PLANE_STRESS) ;
     std::vector<std::pair<Matrix, Matrix> > branches ;
     double factor_k = 1.;
         std::vector<double> K_chaine_cp = {5.4e11/factor_k,  3.9e11/factor_k, 2.02e11/factor_k,5.1e10/factor_k} ;
     for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
     {
         double tau = 5.*std::pow(10., (double) i - 2 ); std::cerr << "TAU "<< tau << std::endl ;
         Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  true,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS )  ; 
         Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas, true,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS ) ;
         branches.push_back(std::make_pair(K_i, Am_i)) ;
     }
    Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    
    double cstrain = -2.e-3; 
    double cstress = -38.0e6; 
    planeType pt = PLANE_STRESS;

    Sample samplef(0.3, 0.6,  0.15, 0.3) ;
    Pore pore(0.0010,  0.15, 0.3) ;

    FeatureTree F ( &samplef ) ;
//     F.addFeature(&samplef, &pore);
    featureTree = &F ;

    //ViscoelasticityAndFracture * pasterupt = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1e-9, 1.) );
    FractureCriterion * mcft = new NonLocalSpaceTimeMCFT(cstress,k_elas, .038) ;
    FractureCriterion * mazar = new NonLocalSpaceTimeMazars(4.52e-5, k_elas, nu_elas, 10, cstress , cstrain, .038, pt ) ;
    DamageModel * linear = new SpaceTimeIsotropicLinearDamage(1.0) ;
    
    ViscoelasticityAndFracture * pasterupt = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, mazar->getCopy(), linear->getCopy() ); 
    ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, mazar->getCopy(), linear->getCopy()); 
    ViscoelasticityAndFracture * pasteruptm = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, mcft->getCopy(), linear->getCopy() ); 
    ViscoelasticityAndFracture * spasteruptm = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, mcft->getCopy(), linear->getCopy()); 
    LogarithmicCreepWithExternalParameters * logcreep = new LogarithmicCreepWithExternalParameters("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 30e9, creep_poisson = 0.2, creep_characteristic_time = 2") ;

    samplef.setBehaviour ( logcreep  ) ;

    Function loadfunc = Function("t 12 /")         *f_range("t", 0., 3) +
    Function("0.25 t 3 - 48 / -")                 *f_range("t", 3, 6) +
    Function("t 6 - 72 / 0.1866666 +")            *f_range("t", 6, 72) ;
    loadfunc *= 0.0002 ;
    loadfunc = Function("t 12 /") * 0.0002 ;
    
    std::vector<LoadingCycle> cycles ;
    cycles.push_back(LoadingCycle(&F, ULTIMATE_STRAIN, 4.6e-5, ETA, 0., 1e-6)) ;
    cycles.push_back(LoadingCycle(&F, ULTIMATE_STRESS, 1e6, ETA, &(cycles[0]), 1e-6)) ;
    cycles.push_back(LoadingCycle(&F, ULTIMATE_STRESS, 1e6*.8, ETA, &(cycles[1]), -1e-6)) ;
    cycles.push_back(LoadingCycle(&F, ULTIMATE_STRESS, 1e5, ETA, &(cycles[2]),1e-6)) ;
    
    std::vector<LagrangeMultiplierType> bctypes = {SET_ALONG_ETA,SET_ALONG_ETA,SET_ALONG_ETA,SET_ALONG_ETA} ;
    std::vector<BoundingBoxPosition> positions = {TOP,TOP,TOP,TOP} ;
    
//     for(double t = 0 ; t < 72  ;t += .1)
//     {
//         std::cout << t << "  " << VirtualMachine().eval(loadfunc,0,0,0,t) << std::endl ;
//     }
//     exit(0)
//     BoundingBoxCycleDefinedBoundaryCondition * loadr = new BoundingBoxCycleDefinedBoundaryCondition(cycles, bctypes, positions) ;
    BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, loadfunc) ;
    loadr->setActive(true);
    F.addBoundaryCondition ( loadr );

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM_AFTER ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.setDeltaTime(.05);
    F.setMinDeltaTime(.00005);
    F.setOrder(LINEAR_TIME_LINEAR) ;

    F.setMaxIterationsPerStep ( 50 );


    step (loadfunc) ;

//    F.getAssembly()->print() ;


    return 0 ;
}
