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

    int nsteps = 150;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;

    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
//         featureTree->printReport( ( i == 0 ) , false);
        Vector str = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        Vector stra = featureTree->getAverageField(STRAIN_FIELD) ;
        std::cout << str[2]/stra[2] << std::endl ;
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
    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

//     Function ageing =1.-f_exp(Function("t -11.76 /"));
    Function ageing =Function("0.008") + f_positivity("t 12 -")*Function("t 12 - 1200 /");
    std::vector<double> densities = {0.0005522204, 0.003067648 , 0.0051838581, 0.0091409803, 0.0150820684, 0.0228759923, 0.0318640245, 0.0412745565, 0.0504521234, 0.0590429114,
                                     0.0669085528, 0.0738898316, 0.0802903061, 0.0860628559, 0.0910846804, 0.0958182743, 0.099462604 , 0.1038763866, 0.107060117 , 0.1108952613, 0.1136154733, 0.1168817815, 0.1184522229, 0.1225589833, 0.1238130176, 0.1270749807, 0.1275711473, 0.1311868434, 0.1312075722, 0.1332925131, 0.1333034615, 0.1360270816, 0.1360297467, 0.1370194198, 0.1384001087, 0.1408002509, 0.1439522145, 0.1470341199, 0.1490669162, 0.1509709139, 0.152818818 , 0.1545770311, 0.1563228192, 0.1582439685, 0.1599688396, 0.1616184296, 0.1632487674, 0.1648086438, 0.1664198001, 0.1680241413 } ;
                                     
     for(size_t i = 25 ; i < densities.size() ; i++)
         densities[i] += 0.2*(i-25)/120 ;
     std::vector<double> timesd ;
     for(double i = 0 ; i < densities.size() ; i++)
       timesd.push_back(i) ;
    
//     for(double t = 0 ; t < 70 ; t++)
//     {
//         std::cout << VirtualMachine().eval(ageing, 0,0,0,t) << std::endl ;
//     }
// 
//     exit(0) ;
    
    std::map<unsigned char,Form *> behaviourMap ;
    behaviourMap[0] = new Stiffness(m0) ;  // C3S
    behaviourMap[1] = new C3SBehaviour() ;  // C3S
    behaviourMap[6] = new CSHBehaviour(INNER_CSH) ; // inner C-S-H 
    behaviourMap[7] = new CSHBehaviour(OUTER_CSH, densities,timesd) ;  // outer C-S-H
    behaviourMap[8] = new CHBehaviour() ; // CH 
    
//     behaviourMap[0] = new Stiffness(m0) ; // C3S
//     behaviourMap[1] = new Stiffness(m0) ;  // C3S
//     behaviourMap[6] = new Stiffness(m0) ;// inner C-S-H 
//     behaviourMap[7] = new CSHBehaviour(OUTER_CSH, ageing) ;  // outer C-S-H
//     behaviourMap[8] = new Stiffness(m0) ; // CH 
   
   std::vector<double> times ;
   for(double i = 0 ; i < 160 ; i++)
       times.push_back(i) ;
   
//     behaviourMap[0] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // pores ?
//     behaviourMap[1] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // C3S
//     behaviourMap[2] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // CH ?
//     behaviourMap[3] = new Viscoelasticity(PURE_ELASTICITY, m0) ;  // outer C-S-H
//     behaviourMap[4] = new Viscoelasticity(PURE_ELASTICITY, m0) ; // inner C-S-H
    FeatureTree F( argv[1], behaviourMap, times ) ;
    F.setDeltaTime(1);
    F.setOrder(LINEAR) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT, -1e-6)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
    

    step(&F) ;

    return 0 ;
}
