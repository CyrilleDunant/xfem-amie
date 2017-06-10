// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/c3s_behaviour.h"
#include "../physics/materials/c2s_behaviour.h"
#include "../physics/materials/c3a_behaviour.h"
#include "../physics/materials/c4af_behaviour.h"
#include "../physics/materials/gyp_behaviour.h"
#include "../physics/materials/ettr_behaviour.h"
#include "../physics/materials/filler_behaviour.h"
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
	// porosity calculation is volume ratio between gel water /outer_CSH 
	// Use the same relationship between porosity and stiffness of outer_CSH
        double d035;
        double p035;
        double sc035;
	time_and_densities >> t >> d035 >> p035 >> sc035 ;
	 
        timesd.push_back(t) ;
        densities.push_back(d035);
    }
  
    for(size_t i = 0 ; i < densities.size() ; i++)
        densities[i] *= 1.1 ;
  
    
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // water
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[2]= new C2SBehaviour() ; //C2S
    behaviourMap[3]= new C3ABehaviour(); //C3A
    behaviourMap[4]= new C4AFBehaviour(); //C4AF
    behaviourMap[5] = new GypBehaviour(); //gypsum
    behaviourMap[6] = new Stiffness(m0) ;  water
    behaviourMap[7] = new CSHBehaviour(INNER_CSH, densities,forces,timesd) ;
    behaviourMap[8] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ;
    behaviourMap[9] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ;
    behaviourMap[10] = new CHBehaviour() ; //CH
    behaviourMap[11] = new EttrBehaviour() ;
    behaviourMap[12] = new EttrBehaviour() ;
    behaviourMap[13] = new EttrBehaviour() ;
    behaviourMap[14] = new CSHBehaviour(OUTER_CSH, densities,forces,timesd) ;
    
    
    FeatureTree F( argv[1], behaviourMap, timesd ) ;
    F.setDeltaTime(1.);
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step(&F, forces) ;

    return 0 ;
}
