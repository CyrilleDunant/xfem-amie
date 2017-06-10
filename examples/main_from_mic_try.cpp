// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/c3s_behaviour.h"
#include "../physics/homogenization/phase.h"
#include "../physics/materials/csh_behaviour.h"
#include "../physics/materials/ch_behaviour.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#define DEBUG


using namespace Amie ;


void step(FeatureTree * featureTree, const std::vector<double> & forces)
{

    int nsteps = 150;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;
    int count = 0 ;

    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;

        Vector stra = featureTree->getAverageField(TOTAL_STRAIN_FIELD) ;
	Vector str = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        //std::cout  << featureTree->getCurrentTime() << "  " << stra[0] << "  " << stra[1]<< "  " << stra[2]<<  std::endl ;
	std::cout  << featureTree->getCurrentTime() << " , " << stra[1] << " , " << str[1] << " , "   << stra[2] << " , " << str[2] << " , "  << stra[3] << " , " << str[3] << std::endl ;
	
      }
}

int main(int argc, char *argv[])
{
    double nu = 0.2 ;
    double E = 1e5 ;
    Matrix m0 =  Tensor::cauchyGreen(std::make_pair(E,nu),SPACE_THREE_DIMENSIONAL);
    
    std::fstream time_and_densities(argv[2])   ;
    
    std::vector<double> timesd ;
    std::vector<double> forces ;
    std::vector<double> densities ;                               
    std::fstream loads(argv[3])   ;                  
    
    while(!loads.eof())
    {
       double load;
        loads >> load ;
        forces.push_back(load) ;
    }
         
    while(!time_and_densities.eof())
    {
        double t ;
	// double correction = 1.0 ; 
       // double d040 ;
        double d032,d040,d048,d032_,d040_,d048_ ;
        double p032, p040, p048 ;
        double sc032, sc040, sc048 ;
	time_and_densities >> t >> d032 >> d040 >> d048 >> p032 >> p040 >> p048 >> d032_ >> d040_ >> d048_ >> sc032 >> sc040 >> sc048;
        //time_and_densities >> t >> d040 ;
	 
        timesd.push_back(t) ;
        densities.push_back(d040);
    }
  
    for(size_t i = 0 ; i < densities.size() ; i++)
        densities[i] *= 1.1 ;
  
    
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // C3S
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[5] = new CSHBehaviour(INNER_CSH, densities,forces,timesd) ; // inner C-S-H 
    behaviourMap[6] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ; // inner C-S-H 
    behaviourMap[7] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ;  // outer C-S-H
    behaviourMap[8] = new CHBehaviour() ; // CH 
    behaviourMap[13] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ;  // outer C-S-H 
   //behaviourMap[6] = new CSHBehaviour(INNER_CSH, densities,std::vector<double>(),timesd) ; // inner C-S-H 
    //behaviourMap[7] = new CSHBehaviour(OUTER_CSH, densities,std::vector<double>(),timesd) ;  // outer C-S-H
    
    
//     behaviourMap[0] = new Stiffness(m0) ; // C3S
//     behaviourMap[1] = new Stiffness(m0) ;  // C3S
//     behaviourMap[6] = new Stiffness(m0) ;// inner C-S-H 
//     behaviourMap[7] = new CSHBehaviour(OUTER_CSH, ageing) ;  // outer C-S-H
//     behaviourMap[8] = new Stiffness(m0) ; // CH 
      
//     behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // pores ?
//     behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // C3S
//     behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // CH ?
//     behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // outer C-S-H
//     behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // inner C-S-H
    FeatureTree F( argv[1], behaviourMap, timesd ) ;
    F.setDeltaTime(1.);
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step(&F, forces) ;

    return 0 ;
}
