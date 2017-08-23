// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../utilities/samplingcriterion.h"
#include "../../features/features.h"
#include "../../physics/physics_base.h"
#include "../../physics/stiffness.h"
#include "../../physics/dual_behaviour.h"
#include "../../physics/logarithmic_creep.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../../physics/fracturecriteria/mohrcoulomb.h"
#include "../../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
#include "../../physics/fracturecriteria/ruptureenergy.h"
#include "../../physics/weibull_distributed_stiffness.h"
#include "../../features/pore.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../physics/orthotropicstiffness.h"
#include "../../features/sample.h"
#include "../../features/sample3d.h"
#include "../../features/polygonSample.h"
#include "../../features/inclusion.h"
#include "../../features/expansiveZone.h"
#include "../../features/crack.h"
#include "../../features/enrichmentInclusion.h"
#include "../../features/microstructuregenerator.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
#include "../../utilities/granulo.h"
#include "../../utilities/placement.h"
#include "../../utilities/itoa.h"
#include "../../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../../physics/fracturecriteria/maxstrain.h"
#include "../../geometry/geometry_base.h"
#include "../../features/microstructuregenerator.h"
#include "../../polynomial/vm_function_extra.h"
#include "../../utilities/writer/voxel_writer.h"
#include "../../physics/homogenization/composite.h"


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

  int nbr_test = 8 ;
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
  Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas,  SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C3S
  double E_C3S = 135.0e9;
  double nu_C3S = 0.3;
  Matrix E_C3S_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C2S
  double E_C2S = 130.0e9;
  double nu_C2S = 0.3;
  Matrix E_C2S_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
   //C3A
  double E_C3A = 145.0e9;
  double nu_C3A = 0.3;
  Matrix E_C3A_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
   //C4AF
  double E_C4AF = 125.0e9;
  double nu_C4AF = 0.3;
  Matrix E_C4AF_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C$
  double E_CS = 30.0e9;
  double nu_CS = 0.3;
  Matrix E_CS_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //CH
  double E_CH = 38.0e9;
  double nu_CH = 0.305;
  Matrix E_CH_elas = Tensor::cauchyGreen( E_CH, nu_CH, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //Aft    
  double E_Aft = 22.4e9;
  double nu_Aft = 0.25;
  Matrix E_Aft_elas = Tensor::cauchyGreen( E_Aft, nu_Aft, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;    
  //Ms+Autre    
  double E_Ms = 42.3e9;
  double nu_Ms = 0.324;
  Matrix E_Ms_elas = Tensor::cauchyGreen( E_Ms, nu_Ms, SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;    
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
    Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  SPACE_THREE_DIMENSIONAL , PLANE_STRESS, YOUNG_POISSON)  ;
    Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas,  SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    branches.push_back(std::make_pair(K_i, Am_i)) ;
  }
  
  //Valeur cible du béton  
  ///////////////////////////////////////////////////
  // Viscoelastic
  ///////////////////////////////////////////////////
    //std::cout << " I am a test "<< std::endl ;
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY,m0, branches.size())  ; //water
    behaviourMap[6] = new Viscoelasticity(PURE_ELASTICITY,m0, branches.size())  ; //water
    behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, E_C3S_elas, branches.size())  ;  // C3S
    behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY,E_C2S_elas, branches.size())  ;  // C2S
    behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY,E_C3A_elas, branches.size())   ;  // C3A
    behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, E_C4AF_elas,branches.size())  ;  // C4AF
    behaviourMap[5] = new Viscoelasticity(PURE_ELASTICITY, E_CS_elas,branches.size())   ;  // C$  
    behaviourMap[7] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    behaviourMap[8] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ; // inner C-S-H 
    behaviourMap[9] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    behaviourMap[10] = new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;  // outer C-S-H*/
 
    FeatureTree F(argv[1], behaviourMap) ;

    F.setOrder(LINEAR_TIME_LINEAR) ;

    
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, -1.e6)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_AFTER)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_AFTER)) ;
  
  F.setMaxIterationsPerStep(100) ;
  
  
  
  double time_step = 1.0 ;    
 
  while(F.getCurrentTime() < 10)
  {
    time_step *= 1.25 ;
    F.setDeltaTime(time_step) ;
    F.step(); 
      
   std::cout << "temps " <<"  "<<  F.getCurrentTime()<< "  " << F.getAverageField( REAL_STRESS_FIELD)[0]  << "," << F.getAverageField( REAL_STRESS_FIELD)[1]  << "," 
    << F.getAverageField( REAL_STRESS_FIELD)[2]    << "  " 
    << " Strain "<< " "  << F.getAverageField(  TOTAL_STRAIN_FIELD)[0]  << "," << F.getAverageField(  TOTAL_STRAIN_FIELD)[1]  
    << "," << F.getAverageField(  TOTAL_STRAIN_FIELD)[2] << "  " << "Creep" << " "<< F.getAverageField( TOTAL_STRAIN_FIELD)[0]/F.getAverageField( REAL_STRESS_FIELD)[0] 
    << "," << F.getAverageField(  TOTAL_STRAIN_FIELD)[1]/F.getAverageField( REAL_STRESS_FIELD)[1]<< "," 
    << F.getAverageField(  TOTAL_STRAIN_FIELD)[2]/F.getAverageField( REAL_STRESS_FIELD)[2]<<std::endl ;  
  
  
  }
  }
  
  return 0 ;
  
}


