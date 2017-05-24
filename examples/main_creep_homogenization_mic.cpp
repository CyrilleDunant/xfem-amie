// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/orthotropicstiffness.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/polygonSample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../features/microstructuregenerator.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../geometry/geometry_base.h"
#include "../features/microstructuregenerator.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/c3s_behaviour.h"
#include "../physics/materials/csh_behaviour.h"
#include "../physics/materials/ch_behaviour.h"
#include "../physics/homogenization/composite.h"


#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;


int main(int argc, char *argv[])
{
    //donnees sorties

  int nbr_test = 10 ;
  for (int k = 0; k < (nbr_test-7); k++)
  {
        std::stringstream ss;
	  ss << k;
        std::string strfile = ss.str();
  
	  double poreArea = 0 ;
  double incArea = 0 ;
  double C3SArea = 0 ;
  double CHArea = 0 ;
  double AftArea = 0 ;
  double MsArea = 0 ;
  double cpArea = 0 ;
  double tol_param = 0.1e9;
  ///////////////////////////////////////////////////
  // Modele et Materiaux      ///////////////////////
  ///////////////////////////////////////////////////
	//C-S-H
  double k_elas = 8.2654e9 + 4*0.5e9 + 2*0.5e9;
  double nu_elas = 0.24;  
  Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C3S
  double E_C3S = 135.0e9;
  double nu_C3S = 0.3;
  Matrix E_C3S_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C2S
  double E_C2S = 130.0e9;
  double nu_C2S = 0.3;
  Matrix E_C2S_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
   //C3A
  double E_C3A = 145.0e9;
  double nu_C3A = 0.3;
  Matrix E_C3A_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
   //C4AF
  double E_C4AF = 125.0e9;
  double nu_C4AF = 0.3;
  Matrix E_C4AF_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C$
  double E_CS = 30.0e9;
  double nu_CS = 0.3;
  Matrix E_CS_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //CH
  double E_CH = 38.0e9;
  double nu_CH = 0.305;
  Matrix E_CH_elas = Tensor::cauchyGreen( E_CH, nu_CH, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //Aft    
  double E_Aft = 22.4e9;
  double nu_Aft = 0.25;
  Matrix E_Aft_elas = Tensor::cauchyGreen( E_Aft, nu_Aft, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;    
  //Ms+Autre    
  double E_Ms = 42.3e9;
  double nu_Ms = 0.324;
  Matrix E_Ms_elas = Tensor::cauchyGreen( E_Ms, nu_Ms, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;    
 //water
    double nu = 0.2 ;
    double E = 1e5 ;
    Matrix m0 =  Tensor::cauchyGreen(E,nu, SPACE_THREE_DIMENSIONAL);
  
  double f_inc;
  

  // Elastique
  ///////////////////////////////////////////////////
  
  std::vector<std::pair<Matrix, Matrix> > branches ;
  double factor_k = 0.38 + 4*(0.1)/10;
  std::vector<double> K_chaine_cp = {2.e11*factor_k,1.2e11*factor_k, 0.97e11*factor_k,7.4e10*factor_k, 2.6e10*factor_k,  2.2e10*factor_k} ;
  for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
  {
    double tau = 0.02*std::pow(10., (double) i - 1 );
    Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  SPACE_TWO_DIMENSIONAL , PLANE_STRESS, YOUNG_POISSON)  ;
    Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    branches.push_back(std::make_pair(K_i, Am_i)) ;
  }
  
  //Valeur cible du béton  
  ///////////////////////////////////////////////////
  // Viscoelastic
  ///////////////////////////////////////////////////
    //std::cout << " I am a test "<< std::endl ;
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[6] = new Stiffness(m0)  ; //water
    behaviourMap[1] =new Stiffness( E_C3S_elas)  ;  // C3S
    behaviourMap[2] = new Stiffness(E_C2S_elas)  ;  // C2S
    behaviourMap[3] = new Stiffness(E_C3A_elas)   ;  // C3A
  /*  behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, E_C4AF_elas,0,0)  ;  // C4AF
    behaviourMap[5] = new Viscoelasticity(PURE_ELASTICITY, E_CS_elas,0,0)   ;  // C$  
    behaviourMap[8] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ; // inner C-S-H 
    behaviourMap[10] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;  // outer C-S-H
    behaviourMap[11] = new Viscoelasticity(PURE_ELASTICITY, E_CH_elas,0,0)  ;//CH
    behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Aft_elas,0,0)  ;//Aft
    behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other
    behaviourMap[13] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other*/
      
   /* behaviourMap[6] = new Viscoelasticity(PURE_ELASTICITY, m0, 0,0)  ; //water
    behaviourMap[1] =new Viscoelasticity(PURE_ELASTICITY, E_C3S_elas,0,0)  ;  // C3S
    behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, E_C2S_elas,0,0)  ;  // C2S
    behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, E_C3A_elas,0,0)   ;  // C3A
    behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, E_C4AF_elas,0,0)  ;  // C4AF
    behaviourMap[5] = new Viscoelasticity(PURE_ELASTICITY, E_CS_elas,0,0)   ;  // C$  
    behaviourMap[8] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ; // inner C-S-H 
    behaviourMap[10] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;  // outer C-S-H
    behaviourMap[11] = new Viscoelasticity(PURE_ELASTICITY, E_CH_elas,0,0)  ;//CH
    behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Aft_elas,0,0)  ;//Aft
    behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other
    behaviourMap[13] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other*/
     std::vector<double> timesd ;
     timesd.push_back(1.1);
    FeatureTree F(argv[1], behaviourMap) ;

    F.setOrder(LINEAR_TIME_LINEAR) ;

    
   // F.step();
  //std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, vaggregates, 10000, 2e-1, 0.0000000001, psd_b, CIRCLE, 1.0, M_PI, 1.0e5, 0.4) ;
  //LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;
  ///////////////////////////////////////////////////
  // CONDITIONS AUX LIMITES      ////////////////////
  // Prendre garde à bien mettre AFTER_ /////////////
  ///////////////////////////////////////////////////
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, -10.e6)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_AFTER)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_AFTER)) ;
  
  F.setMaxIterationsPerStep(256) ;
  
    
  //Chargement
  int i_kcp = 0 ;
  double kcp_inc = 1e9 ;
  double k_error = 1.;
  double k_error_0 = 1.;
  double k_homo ;
  double nu_homo ;
  
//   int cache_cp = F.get2DMesh()->generateCache((Feature *) &samplers) ;
  std::vector<bool> i_cp ;
  std::vector<bool> iC3S ;
  std::vector<bool> iCH ;
  std::vector<bool> iAft ;
  std::vector<bool> iMs ;
  double lame_C3S = E_C3S/(1+nu_C3S);
  double lame_CH = E_CH/(1+nu_CH);
  double lame_Aft = E_Aft/(1+nu_Aft);
  double lame_Ms = E_Ms/(1+nu_Ms);
  
 /* std::string test_data = "homogenization_elas"  ;
  std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::trunc);
  fichier_data << "f_C3S ; " << C3SArea  
               << "; f_CH ; " << CHArea  
	       << " ; f_Aft ; " << AftArea
	       << " ; f_Ms ; " << MsArea
      << "; Kv_param ; " << K_chaine_cp[0] << " ; " <<  K_chaine_cp[1] << " ; "<<   K_chaine_cp[2] << " ; " <<   K_chaine_cp[3] << " ; " << std::endl ;  
  fichier_data<< "Time ;"  << "S11 ; " << "S22 ; "  << " S12 ;  " << " epsi11 ;  " << " epsi22 ;  " << " epsi12 ;  " << " dep2 ;  "  << " nbre points calcul dep2;  "  
       << " dep1 ;  "  << " nbre points calcul dep1;  "     << " dep3 ;  "  << " nbre points calcul dep3;  "  << "S11 C3S ; " << "S22 C3S ; "  << " S12 C3S;  "
       << "S11 Aft ; " << "S22 Aft ; "  << " S12 Aft;  " << "S11 Ms ; " << "S22 Ms ; "  << " S12 Ms;  "<< "S11 CSH ; " << "S22 CSH ; "  << " S12 CSH;  " 
       << " epsi11 C3S;  " << " epsi22 C3S;  " << " epsi12 C3S;" << "epsi11 Aft ; " << "epsi22 Aft ; "  << " epsi12 Aft;  "
       << "epsi11 Ms ; " << "epsi22 Ms ; "  << " epsi12 Ms;  "<< " epsi11 CSH;  " << " epsi22 CSH;  " << " epsi12 CSH;" <<  std::endl ;
  
*/
  
  
  double time_step = 0.01 ;    
 
  while(F.getCurrentTime() < 0.0126)
  {
    time_step *= 1.25 ;
    F.setDeltaTime(time_step) ;
    F.step(); // we have a problem here
      std::cout << " I am a test "<< std::endl ;
    
  std::cout << " =============== " << F.getAverageField( REAL_STRESS_FIELD)[0]  << " =============== " << F.getAverageField( REAL_STRESS_FIELD)[1]  << " =============== " 
    << F.getAverageField( REAL_STRESS_FIELD)[2]    << " =============== " 
    << " Strain " << F.getAverageField( TOTAL_STRAIN_FIELD)[0]  << " =============== " << F.getAverageField( TOTAL_STRAIN_FIELD)[1]  
    << " =============== " << F.getAverageField( TOTAL_STRAIN_FIELD)[2]  << std::endl ;  
  /*  std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
   std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
    fichier_data << F.getCurrentTime()<< " ; " <<F.getAverageField( REAL_STRESS_FIELD)[0]<< " ; " << F.getAverageField( REAL_STRESS_FIELD)[1]<< " ; "<< F.getAverageField( REAL_STRESS_FIELD)[2] 
    << " ; " << F.getAverageField( STRAIN_FIELD)[0] << " ; " << F.getAverageField( STRAIN_FIELD)[1] <<  " ; " << F.getAverageField( STRAIN_FIELD)[2] 
    << std::endl ;
    fichier_data.close();*/
    
    
  std::string test = "creep_homo_mic" + strfile ;
  TriangleWriter writer(test.c_str(), &F, 1.) ;
  writer.getField(TOTAL_STRAIN_FIELD) ;
  writer.getField(TWFT_STIFFNESS) ;	 
  writer.getField(REAL_STRESS_FIELD) ;
  writer.write() ;
  }
  }
  
  return 0 ;
  
}


