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
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"


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

double strain_peak = 0.6*1e-4*5.0/(4.5*0.25);

std::vector< double > loads ;
std::vector< double > displacements ;
std::vector< double > loadsx ;
std::vector< double > displacementsx ;
std::vector< double > damages ;
std::vector< double > cmod ;
std::vector< double > times ;

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, 0.) ;

void step ( size_t nsteps )
{
    
    Function loadfunc = Function("t 12 /")*f_range("t", 0., 3)+
                        Function("0.25 t 3 - 48 / -")*f_range("t", 3, 9) +
                        Function("0.25 1 8 / - t 9 - 48 / + ")*f_range("t", 9, 96) ;
    loadfunc *= strain_peak ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {
        bool go_on = featureTree->step() ;
        if(go_on)
        loadr->setData(VirtualMachine().eval(loadfunc, 0,0,0, featureTree->getCurrentTime()));
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD,-1.,1. ) ;
        Vector etemp = featureTree->getAverageField ( STRAIN_FIELD,-1.,1. ) ;
	Vector dtemp = featureTree->getAverageField (SCALAR_DAMAGE_FIELD) ;

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
	    times.push_back(featureTree->getCurrentTime());
	    damages.push_back(dtemp[0]);	    
            std::fstream ldfile  ;
            ldfile.open ( "ldn", std::ios::out ) ;
            for ( size_t j = 0 ; j < loads.size() ; j++ )
            {
                ldfile <<  times[j] <<" " << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] << " " << damages[j] <<  "\n" ;
            }
            ldfile.close();
	    if(times.back() > 95)
                exit(0) ;
        }

    }

}




int main ( int argc, char *argv[] )
{
    omp_set_num_threads(1) ;
    // Beton
    double k_elas = 40.4e9;
    double nu_elas = 0.2 ;
    Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
    std::vector<std::pair<Matrix, Matrix> > branches ;
    std::vector<std::pair<Matrix, Matrix> > branches_mx ;
    double factor_k = 1.;
        std::vector<double> K_chaine_cp = {5.4e11/factor_k,  3.9e11/factor_k, 2.02e11/factor_k,5.1e10/factor_k} ;
	for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
	{
		double tau = 5*std::pow(10., (double) i - 2 );
		Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas, true,  SPACE_TWO_DIMENSIONAL )  ;
		Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
		branches.push_back(std::make_pair(K_i, Am_i)) ;
	}
        std::vector<double> K_chaine_mx_cp = {1.861269605473859118073396e10, 2.8823932807429079935446e9, 3.4017749304038078634156e9,
	  5.3973666770458357270061e9, 1.01057690570688572352998e10} ;
	std::vector<double> tau_chaine_mx_cp = {0.04647627034132179101144053, 0.4553612767686235043851535, 4.242249572847878343173296, 32.07187552289942944803857} ;
	Matrix K_mx_0 = Tensor::cauchyGreen(K_chaine_mx_cp[0], nu_elas, true,  SPACE_TWO_DIMENSIONAL )  ;
	for(size_t i = 0 ; i < K_chaine_mx_cp.size() ; i++)
	{
		Matrix K_mx_i = Tensor::cauchyGreen(K_chaine_mx_cp[i+1], nu_elas, true,  SPACE_TWO_DIMENSIONAL )  ;
		Matrix Am_mx_i = Tensor::cauchyGreen( K_chaine_mx_cp[i+1]*tau_chaine_mx_cp[i], nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
		branches_mx.push_back(std::make_pair(K_mx_i, Am_mx_i)) ;
	}	
    Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    
    double cstrain = -2.e-3; 
    double cstress = -38.0e6; 
    planeType pt = PLANE_STRESS;
    
    Sample samplef(0.3, 0.6,  0.15, 0.3) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;
    DamageModel * linear = new SpaceTimeFiberBasedIsotropicLinearDamage(0.01, 1.0);

     //ELAS+DAMAGE
    //StiffnessAndFracture * spasterupt = new StiffnessAndFracture(E_cp_elas, new NonLocalMCFT(-40e6,40e9,1.), new FiberBasedIsotropicLinearDamage(0.001, 1.));
    //ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas,new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1e-6, 1.));
    ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75., cstress , cstrain, 1., pt ), linear);
    ViscoelasticityAndFracture * vpasterupt = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75.0, cstress , cstrain, 1., pt ), linear);
    ViscoelasticityAndFracture * vmxpasterupt = new ViscoelasticityAndFracture(GENERALIZED_MAXWELL, K_mx_0, branches_mx, new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75.0, cstress , cstrain, 1., pt ), linear);  
//     Viscoelasticity * vmxpasterupt = new Viscoelasticity(GENERALIZED_MAXWELL, K_mx_0, branches_mx);  
    //    Stiffness * paste = new Stiffness(C_kv);
    samplef.setBehaviour ( vmxpasterupt ) ;

 

    F.addBoundaryCondition ( loadr );

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM_AFTER ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.setDeltaTime(.0075); 
    F.setMinDeltaTime(.02);
//     F.setOrder(LINEAR_TIME_LINEAR) ;
    F.setMaxIterationsPerStep ( 50 );

    //step ( 240 ) ;
    step ( 12400 ) ;

//    F.getAssembly()->print() ;


    return 0 ; 
}
