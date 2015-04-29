// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../utilities/samplingcriterion.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/mcft.h"
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
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.05e-07 ;

double effectiveRadius = 0.13506098 ; //.5*.165*sqrt(M_PI*.25) ;//
double rebarDiametre = 0.0254*sqrt ( M_PI*.25 ) ;

double delta_displacement =  1e-5 ;
double displacement_tolerance = 0.05*delta_displacement ;
double softeningFactor = 1. ;

double percent = 0.01 ;
double displacement  = 0 ;
double prescribedDisplacement = 0;
double derror = 0 ;
double ierror = 0 ;
double preverror = 0 ;
bool firstRun = true ;

double sampleLength = 5.5 ;
double sampleHeight = 1.2 ;
double supportLever = 2.5 ;
double supportMidPointToEndClearance = 0.25 ;
double platewidth = 0.15 ;
double plateHeight = 0.051 ;
double rebarEndCover = 0.047 ;

std::vector< double > loads ;
std::vector< double > displacements ;
std::vector< double > loadsx ;
std::vector< double > displacementsx ;
std::vector< double > damages ;
std::vector< double > cmod ;
std::vector< double > idx ;

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
BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, 0.) ;

void step ( size_t nsteps )
{

    Function loadfunc = Function("t 12 /")         *f_range("t", 0., 12) +
                        Function("1 t 12 - 12 / -")*f_range("t", 12, 24) +
                        Function("t 24 - 24 /")    *f_range("t", 24, 72) ;
    loadfunc *= 0.0002 ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {

        bool go_on = featureTree->step() ;
        if(go_on)
            loadr->setData(VirtualMachine().eval(loadfunc, 0,0,0, featureTree->getCurrentTime()));

        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD,-1.,1. ) ;
        Vector etemp = featureTree->getAverageField ( STRAIN_FIELD,-1.,1. ) ;

        std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2]*1e6 << std::endl ;

        std::cout << std::endl ;

        if(go_on)
        {
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[0] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[0] );

            std::fstream ldfile  ;
            ldfile.open ( "ldn", std::ios::out ) ;
            for ( size_t j = 0 ; j < loads.size() ; j++ )
            {
                ldfile << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] <<  "\n" ;
            }
            ldfile.close();
        }

    }

}




int main ( int argc, char *argv[] )
{
    omp_set_num_threads(1) ;


    double E_rec = 1e10;
    double k_kv = 3.e10;
    double am_kv = 5.*k_kv;
    double k_mx = 2e30;
    double am_mx = 10.0*k_mx ;
    double nu_rec = 0.2 ;
    double delay = 0. ;
    int b = 0 ;
    int a = 0 ;
    Matrix E_kv = Tensor::cauchyGreen( k_kv, nu_rec, true,  SPACE_TWO_DIMENSIONAL ) ;
    Matrix C_kv = Tensor::cauchyGreen( am_kv, nu_rec, true,  SPACE_TWO_DIMENSIONAL ) ;
    Matrix C_tau = Tensor::cauchyGreen( am_mx, nu_rec, true,  SPACE_TWO_DIMENSIONAL ) ;
    Matrix C_eta = Tensor::cauchyGreen( k_mx, nu_rec, true,  SPACE_TWO_DIMENSIONAL ) ;
    
    double mradius = .1 ; 
    double nu = 0.2 ;
    double E_paste = 30e9 ;

	Sample samplef(0.3, 0.6,  0.15, 0.3) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;

    LogarithmicCreepWithExternalParameters * creep = new LogarithmicCreepWithExternalParameters("young_modulus = 1e9, poisson_ratio = 0.2, creep_modulus = 10e9, creep_poisson = 0.2, creep_characteristic_time = 1") ;
    Matrix EC = creep->C ;
    Matrix EE = creep->E ;
    Matrix ER = creep->R ;
    Viscoelasticity * paste = new Viscoelasticity(BURGER, ER, ER*1, EC, EE);
    ViscoelasticityAndFracture * pasterupt = new ViscoelasticityAndFracture(BURGER, ER, ER*1, EC, EE, new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1.)
    );
    StiffnessAndFracture * spasterupt = new StiffnessAndFracture(EE, new NonLocalMCFT(-40e6,40e9,1.), new FiberBasedIsotropicLinearDamage(0.001, 1.)
    );
//    Stiffness * paste = new Stiffness(C_kv);
    samplef.setBehaviour ( spasterupt  ) ;
 

    F.addBoundaryCondition ( loadr );

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM_AFTER ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.setDeltaTime(.1);

    F.setMaxIterationsPerStep ( 3400 );

    step ( 720 ) ;

//    F.getAssembly()->print() ;


    return 0 ;
}
