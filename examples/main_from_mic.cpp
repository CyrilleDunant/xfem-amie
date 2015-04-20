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
        Vector stra = featureTree->getAverageField(STRAIN_FIELD) ;
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
    Matrix m0 =  Tensor::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    
    std::fstream time_and_densities(argv[2])   ;
    
    std::vector<double> timesd ;
    std::vector<double> densities ;                               
    
         
    while(!time_and_densities.eof())
    {
        double correction = 0.8 ;
        double t ;
        double d032,d040,d048 ;
        double p032, p040, p048 ;
        double phi ;
        time_and_densities >> t >> d032 >> d040 >> d048 >>p032>> p040>> p048;
        timesd.push_back(t) ;
        if(atoi(argv[3]) == 0)
        {
            densities.push_back(d032 * correction);
            phi = p032 ;
        }
        else if(atoi(argv[3]) == 1)
        {
            densities.push_back(d040 * correction);
            phi = p040 ;
        }
        else
        {
            densities.push_back(d048 * correction );
            phi = p048 ;
        }
        
        if(t > 15)
        {
            Stiffness sp(Tensor::cauchyGreen(std::make_pair(0,.4997), true,SPACE_THREE_DIMENSIONAL)) ;
            Phase porosity(&sp, phi) ;
            Stiffness scsh(Tensor::cauchyGreen(std::make_pair(1,.25), true,SPACE_THREE_DIMENSIONAL)) ;
            Phase csh(&scsh, 1.-phi) ;
            BiphasicSelfConsistentComposite sc(porosity,csh) ;
            double fac = sc.getBehaviour()->getTensor(Point())[0][0]/scsh.param[0][0]  ;
            if(fac < densities.back())
                densities.back() = fac ;
        }
        
    }
    
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // C3S
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[6] = new CSHBehaviour(INNER_CSH, densities,timesd) ; // inner C-S-H 
    behaviourMap[7] = new CSHBehaviour(OUTER_CSH, densities,timesd) ;  // outer C-S-H
    behaviourMap[8] = new CHBehaviour() ; // CH 
    
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
    F.setDeltaTime(10);
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
