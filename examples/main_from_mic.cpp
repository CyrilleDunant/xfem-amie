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


FeatureTree * featureTree ;

void step()
{

    int nsteps = 1;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;

    for(size_t i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
        featureTree->printReport();
    }
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();

    VoxelWriter vw("sphere_stress_25", 200) ;
    vw.getField(featureTree, VWFT_STRESS) ;
    vw.write();
//     VoxelWriter vw1("sphere_strain", 400) ;
//     vw1.getField(featureTree, VWFT_STRESS) ;
//     vw1.write();
	VoxelWriter vw0("sphere_stiffness_25", 200) ;
	vw0.getField(featureTree, VWFT_STIFFNESS) ;
	vw0.write();
    exit(0) ;
}



int main(int argc, char *argv[])
{

    double nu = 0.2 ;
    double E = 2e9 ;
    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new C3SBehaviour() ;  // C3S
    behaviourMap[1] = new Stiffness(m0) ; // pores ?
    behaviourMap[2] = new CSHBehaviour(INNER_CSH) ; // inner C-S-H 
    behaviourMap[3] = new CHBehaviour() ; // CH ?
    behaviourMap[4] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
    behaviourMap[5] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
   
   
   
//     behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // pores ?
//     behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // C3S
//     behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // CH ?
//     behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // outer C-S-H
//     behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // inner C-S-H
    FeatureTree F( argv[1], behaviourMap ) ;
    featureTree = &F ;
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ZETA, FRONT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step() ;

    return 0 ;
}
