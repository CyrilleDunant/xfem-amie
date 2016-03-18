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
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/materials/aggregate_behaviour.cpp"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"


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


void step ( size_t nsteps, std::string app = std::string() )
{

    Function loadfunc = Function("t 12 /")*f_range("t", 0., 3)+
                        Function("0.25 t 3 - 48 / -")*f_range("t", 3, 9) +
                        Function("0.25 1 8 / - t 9 - 48 / +")*f_range("t", 9.0000001, 96) ;
    loadfunc *= strain_peak ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {
        bool go_on = featureTree->stepToCheckPoint(1,1e-4) ;
        featureTree->setDeltaTime(0.01);
        double load = 0. ;
//         if(go_on)
// 	{
// //	    std::cout << "GO-ON!" << std::endl ;
// //            featureTree->setDeltaTime(0.2) ;
//             load = VirtualMachine().eval(loadfunc, 0,0,0, featureTree->getCurrentTime()) ;
//             loadr->setData(VirtualMachine().eval(loadfunc, 0,0,0, featureTree->getCurrentTime()));
// 	}

        Vector stemp = featureTree->getAverageField ( GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD ) ;
        Vector etemp = featureTree->getAverageField ( GENERALIZED_VISCOELASTIC_STRAIN_FIELD ) ;
        Vector dtemp = featureTree->getAverageField (SCALAR_DAMAGE_FIELD) ;

        std::cout << "current time :" << featureTree->getCurrentTime() << std::endl ;
        std::cout << "current delta time :" << featureTree->getDeltaTime() << std::endl ;
        std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
        std::cout << "average sigma44 : " << stemp[4]/1e6 << std::endl ;
        std::cout << "average sigma66 : " << stemp[7]/1e6 << std::endl ;
        std::cout << "average sigma88 : " << stemp[10]/1e6 << std::endl ;
        std::cout << "average sigma1010 : " << stemp[13]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2]*1e6 << std::endl ;
        std::cout << "average epsilon44 : " << etemp[4]*1e6 << std::endl ;
        std::cout << "average epsilon66 : " << etemp[7]*1e6 << std::endl ;
        std::cout << "average epsilon88 : " << etemp[10]*1e6 << std::endl ;
        std::cout << "average epsilon1010 : " << etemp[13]*1e6 << std::endl ;
        std::cout << std::endl ;

        if(go_on)
        {
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[4] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[4] );
            times.push_back(featureTree->getCurrentTime());
            damages.push_back(dtemp[0]);
//            if(damages.size() > 2 && damages[ damages.size()-1 ] == damages[ damages.size()-2])
//		featureTree->setDeltaTime( 0.2 ) ;
            if(true) //v%50 == 0)
            {
                std::fstream ldfile  ;
                ldfile.open ( "ldnviscodamage"+app, std::ios::app | std::ios::out ) ;
                ldfile <<  times.back() <<" " << displacements.back() << "   " << loads.back() << "   " <<  displacementsx.back() << "   " << loadsx.back() << " " << damages.back() << " " << load << std::endl ;
                ldfile.close();

//                 TriangleWriter trg("triangles", featureTree, 1.) ;
//                 trg.getField( SCALAR_DAMAGE_FIELD ) ;
//                 trg.getField( GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD ) ;
//                 trg.getField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD ) ;
//                 trg.write() ;

            }

            double time_n = featureTree -> getCurrentTime() ;
            if(time_n > 30)
            {
                featureTree->getAssembly()->print() ;

                exit(0) ;
            }
        }

    }

}

int main ( int argc, char *argv[] )
{
//    omp_set_num_threads(1) ;
    // Beton
    double k_elas = 40.4e9;
    double nu_elas = 0.2 ;
    Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    std::vector<std::pair<Matrix, Matrix> > branches ;
    std::vector<std::pair<Matrix, Matrix> > branches_mx ;
    double factor_k = 1.;
    std::vector<double> K_chaine_cp = {5.4e11/factor_k,  3.9e11/factor_k, 2.02e11/factor_k,5.1e10/factor_k} ;
    for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
    {
        double tau = 5*std::pow(10., (double) i - 2 );
        Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON )  ;
        Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
        branches.push_back(std::make_pair(K_i, Am_i)) ;
    }
    std::vector<double> K_chaine_mx_cp = {1.861269605473859118073396e10, 2.8823932807429079935446e9, 3.4017749304038078634156e9,
                                          5.3973666770458357270061e9, 1.01057690570688572352998e10
                                         } ;
    std::vector<double> tau_chaine_mx_cp = {0.04647627034132179101144053, 0.4553612767686235043851535, 4.242249572847878343173296, 32.07187552289942944803857} ;
    Matrix K_mx_0 = Tensor::cauchyGreen(K_chaine_mx_cp[0], nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS , YOUNG_POISSON)  ;
    for(size_t i = 0 ; i < K_chaine_mx_cp.size()-1 ; i++)
    {
        Matrix K_mx_i = Tensor::cauchyGreen(K_chaine_mx_cp[i+1], nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON )  ;
        Matrix Am_mx_i = Tensor::cauchyGreen( K_chaine_mx_cp[i+1]*tau_chaine_mx_cp[i], nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
        branches_mx.push_back(std::make_pair(K_mx_i, Am_mx_i)) ;
    }
    Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    Viscoelasticity * pastemx = new Viscoelasticity( GENERALIZED_MAXWELL, K_mx_0, branches_mx) ;
    LogarithmicCreepWithExternalParameters * logcreep = new LogarithmicCreepWithExternalParameters("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 30e9, creep_poisson = 0.2, creep_characteristic_time = 2") ;

    double cstrain = -2.e-3;
    double cstress = -38.0e6;
    planeType pt = PLANE_STRESS;

    Sample samplef(0.3, 0.6,  0.15, 0.3) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;
    DamageModel * linear = new SpaceTimeIsotropicLinearDamage(1.0);

    //ELAS+DAMAGE
    //StiffnessAndFracture * spasterupt = new StiffnessAndFracture(E_cp_elas, new NonLocalMCFT(-40e6,40e9,1.), new FiberBasedIsotropicLinearDamage(0.001, 1.));
    //ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas,new NonLocalSpaceTimeMCFT(-40e6,40e9,1.), new SpaceTimeFiberBasedIsotropicLinearDamage(0.001, 1e-6, 1.));
    ViscoelasticityAndFracture * spasterupt = new ViscoelasticityAndFracture(PURE_ELASTICITY, E_cp_elas, new NonLocalSpaceTimeMCFT(cstress,k_elas,1.)/*NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75., cstress , cstrain, 1., pt )*/, linear);
    ViscoelasticityAndFracture * vpasterupt = new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches, new NonLocalSpaceTimeMCFT(cstress,k_elas,1.)/*new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75.0, cstress , cstrain, 1., pt )*/, linear);
    ViscoelasticityAndFracture * vmxpasterupt = new ViscoelasticityAndFracture(GENERALIZED_MAXWELL, K_mx_0, branches_mx, /*new NonLocalSpaceTimeMCFT(cstress,k_elas,1.)*/new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75.0, cstress , cstrain, 1., pt ), linear);
    Point t1(0.0001, 1.2e6) ;
    Point t2(0.0005, 0.) ;
    std::vector<Point> ptension, pcompression ;
    ptension.push_back(t1) ;
    ptension.push_back(t2) ;
    LogarithmicCreepWithExternalParameters * logcreeprupt = new LogarithmicCreepWithExternalParameters("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 30e9, creep_poisson = 0.2, creep_characteristic_time = 2", new NonLocalSpaceTimeMazars(1.0e-4, k_elas, nu_elas, 75.0, cstress , cstrain, 1., pt ), linear) ;//new SpaceTimeIsotropicLinearDamage()) ;
//     Viscoelasticity * vmxpasterupt = new Viscoelasticity(GENERALIZED_MAXWELL, K_mx_0, branches_mx);  e
    //    Stiffness * paste = new Stiffness(C_kv);

    samplef.setBehaviour ( vpasterupt ) ;

    if(argc > 2)
    {
        bool damage = (argc > 3 && std::string(argv[3]) == "d") ;
        if(std::string(argv[2]) == "mx")
            samplef.setBehaviour( damage ? vmxpasterupt : pastemx ) ;
        else if(std::string(argv[2]) == "kv")
            samplef.setBehaviour( damage ? vpasterupt : paste ) ;
        else if(std::string(argv[2]) == "log")
            samplef.setBehaviour( damage ? logcreeprupt : logcreep ) ;
    }

    Function loadfunc = Function("t 12 /")*f_range("t", 0., 3)+
                        Function("0.25 t 3 - 48 / -")*f_range("t", 3, 9) +
                        Function("0.25 1 8 / - t 9 - 48 / +")*f_range("t", 9.0000001, 96) ;
    loadfunc *= strain_peak ;

    BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP_AFTER, loadfunc) ;
    F.addBoundaryCondition ( loadr );
    loadr->setActive(true);

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT_AFTER ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM_AFTER ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.setDeltaTime(0.01);

    F.setMinDeltaTime(1e-4);
//     F.setOrder(LINEAR_TIME_LINEAR) ;


    //step ( 240 ) ;
    step (124000, std::string( argc > 2 ? argv[2] : "") ) ;
//    step (12, std::string( argc > 2 ? argv[2] : "") ) ;

    /*    std::vector<Point *> nodes = F.getNodes() ;
        for(size_t i = 0 ; i < nodes.size() ; i++)
    	nodes[i]->print() ;*/



    return 0 ;
}
