// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/c3s_behaviour.h"
#include "../physics/materials/csh_behaviour.h"
#include "../physics/materials/ch_behaviour.h"
#include "../physics/homogenization/composite.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#define DEBUG


using namespace Amie ;



 
void step(FeatureTree * featureTree, const std::vector<double> & times)
{

    int nsteps = 150;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;

    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
//         if(i+2 < times.size() )
//             featureTree->setDeltaTime(times[i+2] - times[i+1]);
//         featureTree->printReport( ( i == 0 ) , false);
        Vector str = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        Vector stra = featureTree->getAverageField(TOTAL_STRAIN_FIELD) ;
        std::cout  << featureTree->getCurrentTime() << "  " << str[2]/stra[2] << std::endl ;
//         std::stringstream filename ;
//         filename << "sphere_stiffness_" << i ;
//         
//         VoxelWriter vw0(filename.str().c_str(), 50) ;
//         vw0.getField(featureTree, VWFT_STIFFNESS) ;
//         vw0.write();
    }
}



int main(int argc, char *argv[])
{
    double nu = 0.2 ;
    double E = 1e5 ;
    Matrix m0 =  Tensor::cauchyGreen(E,nu, SPACE_THREE_DIMENSIONAL);

    
    std::fstream time_and_densities(argv[2])   ;
    
    std::vector<double> timesd ;
    std::vector<double> timesh ;
    
    std::vector<double> densities ;                               
    
         
//     while(!time_and_densities.eof())
//     {
//         double correction = .77 ; // sc is 1
//         double t ;
//         double d032,d040,d048,d032_,d040_,d048_ ;
//         double p032, p040, p048 ;
//         double sc032, sc040, sc048 ;
//                            // 1      2       3       4       5       6      7        8       9        10        11       12       13
//         time_and_densities >> t >> d032 >> d040 >> d048 >> p032 >> p040 >> p048 >> d032_ >> d040_ >> d048_ >> sc032 >> sc040 >> sc048;
//         timesd.push_back(t) ;
//         if(atoi(argv[3]) == 0)
//         {
//             densities.push_back(d032_ * correction);
//         }
//         else if(atoi(argv[3]) == 1)
//         {
//             densities.push_back(d040_ * correction);
//         }
//         else
//         {
//             densities.push_back(d048_ * correction );
//         }
//     }
    while(!time_and_densities.eof())
    {
        double t ; 
        time_and_densities >> t ;
        timesh.push_back(t) ;
    }
    for(double t = 0 ; t < 350 ; t++)
    {
	    double phi = 0.8202*pow(t, -0.073) ;
	    densities.push_back((1+t/2400)*1.5/(1+exp(-100.*((1.-phi)-0.38))));
	    timesd.push_back(t) ;
	    
    }
    
    
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // C3S
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[6] = new CSHBehaviour(INNER_CSH, densities,std::vector<double>(),timesd) ; // inner C-S-H 
    behaviourMap[8] = new CSHBehaviour(OUTER_CSH, densities,std::vector<double>(),timesd) ;  // outer C-S-H
    behaviourMap[9] = new CHBehaviour() ; // CH 
    
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
    FeatureTree F( argv[1], behaviourMap, timesh, argv[1],0,new CSHBehaviour(OUTER_CSH, densities,std::vector<double>(),timesd) ) ;
    F.setDeltaTime(5);
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT, -1e-6)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step(&F, timesd) ;

    return 0 ;
}
