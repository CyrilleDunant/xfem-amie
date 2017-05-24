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
#include "../physics/ageing_logarithmic_creep.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/features.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"
#include "../utilities/parser/command_line_parser.h"
#include "../utilities/parser/config_parser.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
      CommandLineParser parser("2D viscoelastic simulation") ;
	parser.addString("--paste-behaviour", "", "path to a *.ini file containing the behaviour of the cement paste") ;
	parser.addString("--aggregate-behaviour", "", "path to a *.ini file containing the behaviour of the aggregates") ;
 	parser.addValue("--duration", 5, "duration of the creep simulation") ;
 	parser.addValue("--set-delta-time-increment", 0.3, "increases the time step at each step") ;
	parser.parseCommandLine( argc, argv ) ;
 
 	std::string iniPaste = parser.getString("--paste-behaviour") ;
	std::string iniAgg = parser.getString("--aggregate-behaviour") ;
	double maxTime = 10 ;
	double incr =0.1 ;
	
	std::fstream time_and_densities(argv[2])   ;
    
	std::vector<double> timesd ;
	std::vector<double> forces ;
	std::vector<double> densities ;                               
	std::fstream loads(argv[3])   ;  

	Sample box(nullptr, 0.08,0.08,0.,0.) ;
	  ///////////////////////////////////////////////////
      // Modele et Materiaux      ///////////////////////
      ///////////////////////////////////////////////////
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
	Matrix E_C2S_elas = Tensor::cauchyGreen( E_C2S, nu_C2S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
	//C3A
	double E_C3A = 145.0e9;
	double nu_C3A = 0.3;
	Matrix E_C3A_elas = Tensor::cauchyGreen( E_C3A, nu_C3A, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
	//C4AF
	double E_C4AF = 125.0e9;
	double nu_C4AF = 0.3;
	Matrix E_C4AF_elas = Tensor::cauchyGreen( E_C4AF, nu_C4AF, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
	//C$
	double E_CS = 30.0e9;
	double nu_CS = 0.3;
	Matrix E_CS_elas = Tensor::cauchyGreen( E_CS, nu_CS, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
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
	
	std::map<unsigned char,Form *> behaviourMap ;
	behaviourMap[6] = new Stiffness(m0)  ; //water
	behaviourMap[1] = new Stiffness( E_C3S_elas)  ;  // C3S
	behaviourMap[2] = new Stiffness(E_C2S_elas)  ;  // C2S
	behaviourMap[3] = new Stiffness(E_C3A_elas)   ;  // C3A
	behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, E_C4AF_elas,0,0)  ;  // C4AF
	behaviourMap[5] = new Viscoelasticity(PURE_ELASTICITY, E_CS_elas,0,0)   ;  // C$  
	behaviourMap[8] = new AgeingLogarithmicCreep(E_cp_elas, new RealTimeLogCreepAccumulator()) ; // inner C-S-H 
	behaviourMap[10] = new AgeingLogarithmicCreep(E_cp_elas, new RealTimeLogCreepAccumulator()) ;  // outer C-S-H
	behaviourMap[11] = new Viscoelasticity(PURE_ELASTICITY, E_CH_elas,0,0)  ;//CH
	behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Aft_elas,0,0)  ;//Aft
	behaviourMap[12] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other
	behaviourMap[13] = new Viscoelasticity(PURE_ELASTICITY, E_Ms_elas,0,0)  ;//Afm other
    
       Form * paste = new AgeingLogarithmicCreep(E_cp_elas, new RealTimeLogCreepAccumulator()) ;
    
	FeatureTree F(argv[1], behaviourMap, timesd) ;
	box.setBehaviour( paste ) ;
        
	/*F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setSamplingNumber(20) ;
	F.setDeltaTime(0.1) ;
	F.setMinDeltaTime(1e-9) ;
	F.setSamplingRestriction( 0.001 ) ;
// 	std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, aggregate, 6000, 0.008, 0.000001, new PSDBolomeA(), new InclusionGenerator(), 6000*10) ; 
// */
	parser.setFeatureTree( &F ) ;

	double d = F.getDeltaTime() ;

	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, BOTTOM_AFTER) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, TOP_AFTER ) ) ;

	while(F.getCurrentTime() < maxTime)
	{
		F.step() ;
		Vector strain = F.getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD ) ;
		Vector stress = F.getAverageField( REAL_STRESS_FIELD ) ;
		std::cout << F.getCurrentTime() << "\t" ;
		for(size_t j = 0 ; j < strain.size() ; j++)
			std::cout << "strain"<< strain[j] << "\t" ;
		for(size_t j = 0 ; j < stress.size() ; j++)
			std::cout << "stress"<<stress[j] << "\t" ;
		std::cout << std::endl ;

		d += incr ;
		F.setDeltaTime( d ) ;

	}
	return 0 ;
}

