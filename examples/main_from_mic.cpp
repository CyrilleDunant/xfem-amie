// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/c3s_behaviour.h"
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



 
void step(FeatureTree * featureTree )
{

    int nsteps = 70;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;

    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
//         featureTree->printReport( ( i == 0 ) , false);
        Vector str = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        std::cout << -str[2]*1e-7 << std::endl ;
    }
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();

//     VoxelWriter vw("sphere_stress_25", 200) ;
//     vw.getField(featureTree, VWFT_STRESS) ;
//     vw.write();
//     VoxelWriter vw1("sphere_strain", 400) ;
//     vw1.getField(featureTree, VWFT_STRESS) ;
//     vw1.write();
// 	VoxelWriter vw0("sphere_stiffness_25", 200) ;
// 	vw0.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw0.write();
//     exit(0) ;
}



int main(int argc, char *argv[])
{

    double nu = 0.2 ;
    double E = 1e6 ;
    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // C3S
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[2] = new C3SBehaviour() ;  // C2S
    behaviourMap[3] = new C3SBehaviour() ;  // aluminate
    behaviourMap[4] = new C3SBehaviour() ;  // ferrite
    behaviourMap[5] = new C3SBehaviour() ;  // gypsum
    behaviourMap[6] = new CSHBehaviour(INNER_CSH) ; // inner C-S-H 
    behaviourMap[7] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
    behaviourMap[8] = new CHBehaviour() ; // CH 
    behaviourMap[9] = new CHBehaviour() ; // Ettringite 
    behaviourMap[10] = new CHBehaviour() ; // Ettringite 
    behaviourMap[11] = new CHBehaviour() ; // monosulfo 
    behaviourMap[12] = new CHBehaviour() ; // monosulfo 
    behaviourMap[13] = new CHBehaviour() ; // iron hydroxide 
    behaviourMap[14] = new CHBehaviour() ; // iron hydroxide 
    behaviourMap[15] = new CHBehaviour() ; // hydrogarnet
    behaviourMap[16] = new CHBehaviour() ; // hydrogarnet
    behaviourMap[17] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
    behaviourMap[18] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
   
   std::vector<double> times ;
   for(double i = 0 ; i < 24 ; i++)
       times.push_back(i) ;
   
//     behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // pores ?
//     behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // C3S
//     behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // CH ?
//     behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // outer C-S-H
//     behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // inner C-S-H
    FeatureTree F( argv[1], behaviourMap, times ) ;
    F.setDeltaTime(1);
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step(&F) ;

    return 0 ;
}
