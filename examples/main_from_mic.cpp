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

        Vector x = featureTree->getDisplacements() ;

        std::cout << "unknowns :" << x.size() << std::endl ;

        int npoints = featureTree->get3DMesh()->begin()->getBoundingPoints().size() ; ;

        double volume = 0 ;

        double xavg = 0 ;

        for(auto k = featureTree->get3DMesh()->begin() ; k != featureTree->get3DMesh()->end() ; k++)
        {
            if(k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                double ar = k->volume() ;
                volume += ar ;
                for(size_t l = 0 ; l < npoints ; l++)
                {
                    xavg += x[k->getBoundingPoint(l).getId()*3]*ar/npoints ;
                }
            }
        }

        xavg /= volume ;
        std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax(REAL_STRESS_FIELD) ;
        std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax(STRAIN_FIELD) ;
        std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
        Vector stemp = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        Vector etemp = featureTree->getAverageField(STRAIN_FIELD) ;

        std::cout << std::endl ;
        std::cout << "max value :" << x.max() << std::endl ;
        std::cout << "min value :" << x.min() << std::endl ;
        std::cout << "avg value :" << xavg << std::endl ;

        std::cout << "max sigma11 :" << stempm.second[0]  << std::endl ;
        std::cout << "min sigma11 :" << stempm.first[0]   << std::endl ;
        std::cout << "max sigma12 :" << stempm.second[3]  << std::endl ;
        std::cout << "min sigma12 :" << stempm.first[3]   << std::endl ;
        std::cout << "max sigma13 :" << stempm.second[4]  << std::endl ;
        std::cout << "min sigma13 :" << stempm.first[4]   << std::endl ;
        std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
        std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;
        std::cout << "max sigma23 :" << stempm.second[5]  << std::endl ;
        std::cout << "min sigma23 :" << stempm.first[5]   << std::endl ;
        std::cout << "max sigma33 :" << stempm.second[2]  << std::endl ;
        std::cout << "min sigma33 :" << stempm.first[2]   << std::endl ;

        std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
        std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
        std::cout << "max epsilon12 :" << etempm.second[3] << std::endl ;
        std::cout << "min epsilon12 :" << etempm.first[3]  << std::endl ;
        std::cout << "max epsilon13 :" << etempm.second[4] << std::endl ;
        std::cout << "min epsilon13 :" << etempm.first[4]  << std::endl ;
        std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
        std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;
        std::cout << "max epsilon23 :" << etempm.second[5] << std::endl ;
        std::cout << "min epsilon23 :" << etempm.first[5]  << std::endl ;
        std::cout << "max epsilon33 :" << etempm.second[2] << std::endl ;
        std::cout << "min epsilon33 :" << etempm.first[2]  << std::endl ;

        std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
        std::cout << "min von Mises :" << vmm.first[0] << std::endl ;

        std::cout << "average sigma11 : " << stemp[0] << std::endl ;
        std::cout << "average sigma22 : " << stemp[1] << std::endl ;
        std::cout << "average sigma33 : " << stemp[2] << std::endl ;
        std::cout << "average sigma12 : " << stemp[3] << std::endl ;
        std::cout << "average sigma13 : " << stemp[4] << std::endl ;
        std::cout << "average sigma23 : " << stemp[5] << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
        std::cout << "average epsilon33 : " << etemp[2] << std::endl ;
        std::cout << "average epsilon12 : " << etemp[3] << std::endl ;
        std::cout << "average epsilon13 : " << etemp[4] << std::endl ;
        std::cout << "average epsilon23 : " << etemp[5] << std::endl ;

    }
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();

//     VoxelWriter vw("sphere_stress_25", 200) ;
//     vw.getField(featureTree, VWFT_STRESS) ;
//     vw.write();
// //     VoxelWriter vw1("sphere_strain", 400) ;
// //     vw1.getField(featureTree, VWFT_STRESS) ;
// //     vw1.write();
// 	VoxelWriter vw0("sphere_stiffness_25", 200) ;
// 	vw0.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw0.write();
    exit(0) ;
}



int main(int argc, char *argv[])
{

    double nu = 0.2 ;
    double E = 2e9 ;
    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ; // pores ?
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[2] = new CHBehaviour() ; // CH ?
    behaviourMap[3] = new CSHBehaviour(OUTER_CSH) ;  // outer C-S-H
    behaviourMap[4] = new CSHBehaviour(INNER_CSH) ; // inner C-S-H
//     behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // pores ?
//     behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // C3S
//     behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // CH ?
//     behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // outer C-S-H
//     behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // inner C-S-H
    FeatureTree F( "../examples/data/Voxels_50.txt", behaviourMap ) ;
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
