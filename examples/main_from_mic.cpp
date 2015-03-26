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
    Matrix m0 =  Tensor::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

//     Function ageing =1.-f_exp(Function("t -11.76 /"));
    Function ageing =Function("0.008") + f_positivity("t 12 -")*Function("t 12 - 1200 /");
    std::vector<double> densities = {0.0030, 0.0037, 0.0044, 0.0051, 0.0058, 0.0065, 0.0072, 0.0079, 0.0086, 0.0093,
        0.009       
, //0.01 - hi
        0.01181818, 0.01363636, 0.01545455, 0.01727273, 0.01909091, 0.02090909, 0.02272727, 0.02454545, 0.02636364, 0.02818182,
        0.03
        , //0.05 - hi
        0.03363636, 0.03727273, 0.04090909, 0.04454545, 0.04818182, 0.05181818, 0.05545455, 0.05909091, 0.06272727, 0.06636364,
        0.07        
, //0.11 - hi
        0.07363636, 0.07727273, 0.08090909, 0.08454545, 0.08818182, 0.09181818, 0.09545455, 0.09909091, 0.10272727, 0.10636364,
        0.11
        ,  //0.12 - hi
        0.1136364, 0.1172727, 0.1209091, 0.1245455, 0.1281818, 0.1318182, 0.1354545, 0.1390909, 0.1427273, 0.1463636,
        0.15
,          //0.13 - low  ; 0.18255213 - hi
        0.1518182, 0.1536364, 0.1554545, 0.1572727, 0.1590909, 0.1609091, 0.1627273, 0.1645455, 0.1663636, 0.1681818,
        0.17
,          //0.14 - low ; 0.187814602 - hi
        0.1709091, 0.1718182, 0.1727273, 0.1736364, 0.1745455, 0.1754545, 0.1763636, 0.1772727, 0.1781818, 0.1790909,
        0.18        
,   //0.193393385 - hi
        0.19401703
,0.194550838
,0.195166033
,0.195663744
,0.196090483
,0.196489615
,0.196858146
,0.197209009
,0.197581075
,
        0.197934204
,
        0.198279303
,0.198618934
,0.198954894
,0.199314187
,0.199666172
,0.200006228
,0.200338256
,0.200668875
,0.200997652
,
        0.201325906
,
        0.20203677
,0.202745681
,0.203401935
,0.203999343
,0.204528733
,0.204981613
,0.205361634
,0.205695803
,0.205997172
,
        0.207364795
,
        0.210096611
,0.212544122
,0.214386231
,0.215744574
,0.2167839
,0.217604962
,0.218270125
,0.218820062
,0.219282427
,
        0.219676683
,
        0.22001692
,0.220313595
,0.220574629
,0.22080613
,0.221012889
,0.221198712
,0.221366662
,0.221519235
,0.221658479
,
        0.221786097
,
        0.221903513
,0.222011926
,0.222112356
,0.222205673
,0.222292628
,0.222373865
,0.222449949
,0.222521369
,0.222588556
,
        0.22265189
,
        0.222711705
,0.222768299
,0.222821936
,0.222872853
,0.222921262
,0.222967354
,0.2230113
,0.223053255
,0.22309336
,
        0.223131742
,
        0.223168518
,0.223203793
,0.223237665
,0.223270221 ,0.223301544
,0.223331707
,0.223360781
,0.223388828
,0.223415907
,
        0.223442072
,
        0.223467374
,0.22349186
,0.223515571
,0.22353855
,0.223560834 } ;
                                     
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
    F.setDeltaTime(10);
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
