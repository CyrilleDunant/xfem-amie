// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 505-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/stiffness.h"
#include "../features/sample3d.h"
#include "../features/polygonSample3d.h"
#include "../features/loftedPolygonSample3d.h"
#include "../features/inclusion3d.h"
#include "../polynomial/vm_function_extra.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/homogenization/composite.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>

using namespace Amie ;

GeometryAndFaceDefinedSurfaceBoundaryCondition * waterload ;

void step(FeatureTree * featureTree)
{
    int nsteps = 5;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;
    Point crest(0,-40,31) ;
    
    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
        Vector tmp(3) ;
        double size = 350 ;

        for(double z = -size*.5 ;  z <= size*.5 ; z += .5)
        {
            for(double y = -size*.5 ;  y <= size*.5 ; y += .5)
            {
                featureTree->get3DMesh()->getField(PRINCIPAL_REAL_STRESS_FIELD, Point(0,y,z),tmp) ;
                    std::cout <<  tmp[0]<< "  "<< std::flush ;
            }
            std::cout << std::endl ;
        }
/*        
        for(double z = -size*.5 ;  z <= size*.5 ; z += 1)
        {
            for(double y = -size*.5 ;  y <= size*.5 ; y += 1)
            {
                featureTree->get3DMesh()->getField(PRINCIPAL_REAL_STRESS_FIELD, Point(0,y,z),tmp) ;
                    std::cout <<  tmp[1]<< "  "<< std::flush ;
            }
            std::cout << std::endl ;
        }*/
        for(double z = -size*.5 ;  z <= size*.5 ; z += .5)
        {
            for(double y = -size*.5 ;  y <= size*.5 ; y += .5)
            {
                featureTree->get3DMesh()->getField(PRINCIPAL_REAL_STRESS_FIELD, Point(0,y,z),tmp) ;
                    std::cout <<  tmp[2]<< "  "<< std::flush ;
            }
            std::cout << std::endl ;
        }
        
//         VoxelWriter vw1("stiffness", 50) ;
//         vw1.getField(featureTree, VWFT_STIFFNESS) ;
//         vw1.write();
        exit(0) ;

        if(i == 0)
        {
            VoxelWriter vw("stress_full", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            waterload->setData(Function("-9810 z 0 - *")*f_negativity(Function("z 0 -")));
            
        }
        if(i == 1)
        {
            VoxelWriter vw("stress_intermediate_0", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            waterload->setData(Function("-9810 z -30 - *")*f_negativity(Function("z -30 -")));
        }
        if(i == 2)
        {
            VoxelWriter vw("stress_intermediate_1", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            waterload->setData(Function("-9810 z -70 - *")*f_negativity(Function("z -70 -")));
        }
        if(i == 3)
        {
            VoxelWriter vw("stress_empty", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            waterload->setData(Function("-9810 z 30 - *")*f_negativity(Function("z -30 -")));
        }
        if(i == 4 )
        {
            VoxelWriter vw("stress_full_2", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
        }
        Vector crestDisp = featureTree->getDisplacements(&crest) ;
        std::cout << crestDisp[0] << "  " << crestDisp[1] << "  " << crestDisp[2] << std::endl ;
        
//         featureTree->printReport();
    }

}

int main(int argc, char *argv[])
{
    
    
    
    double nu = 0.2 ;
    double E = 20 ;
    Sample3D samplers(nullptr, 350,350,350,0,0,0) ;

    FeatureTree F(&samplers) ;

    nu = 0.2 ;
    E = 37e9 ;
    Matrix m1 = Tensor::cauchyGreen(E, nu,SPACE_THREE_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON);
    
//     for(double phi = 1 ;  phi >=0  ; phi -=0.01)
//     {
//         Phase porosity(new Stiffness(Tensor::cauchyGreen(std::make_pair(0.03,0.00), true,SPACE_THREE_DIMENSIONAL)), 1.-phi) ;
//         Phase csh(new Stiffness(Tensor::cauchyGreen(std::make_pair(1,.25), true,SPACE_THREE_DIMENSIONAL)), phi) ;
//         BiphasicSelfConsistentComposite sc(csh,porosity) ;
//         double c00 = sc.getBehaviour()->getTensor(Point())[0][0] ;
//         double c55 = sc.getBehaviour()->getTensor(Point())[5][5] ;
//         double mu = c55 ;
//         double K = c00+mu*4./3. ;
//         E = 9.*K*(K-mu)/(3.*K-mu) ;
//         nu =mu/(3.*K-mu) ;
//         
//         std::cout << 1.-phi << "  "<< E << "  "<< nu << std::endl ;
//     }
//     exit(0) ;
    
//     Function loadfunc = Function("t 9 /")         *f_range("t", 0., 9) +
//                         Function("1 t 9 - 26 / -")*f_range("t", 9, 24) +
//                         Function("t 24 - 48 / 0.5 +")    *f_range("t", 24, 48) ;
//     for(double i = 0 ; i < 48 ; i+=0.1)
//     {
//         std::cout << i << "  " << VirtualMachine().eval(loadfunc, 0,0,0,i) << std::endl ;
//     }
//     exit(0) ;
    
    samplers.setBehaviour(new VoidForm()) ;

    std::valarray<Point *> damProfile(22) ;
    damProfile[0] =  new Point( 28, 148,    0) ;
    damProfile[1] =  new Point( 36, 148,    0) ;
    damProfile[2] =  new Point( 34, 134,    0) ;
    damProfile[3] =  new Point( 32, 114,    0) ;
    damProfile[4] =  new Point( 32,  91,    0) ;
    damProfile[5] =  new Point( 34,  75,    0) ;
    damProfile[6] =  new Point( 37,  61,    0) ;
    damProfile[7] =  new Point( 41,  44,    0) ;
    damProfile[8] =  new Point( 50,  27,    0) ;
    damProfile[9] =  new Point( 53,  27,    0) ;
    damProfile[10] = new Point( 53,  23,    0) ;
    damProfile[11] = new Point( 55,  18,    0) ;
    damProfile[12] = new Point( 61,  18,    0) ;
    damProfile[13] = new Point( 60,   9,    0) ;
    damProfile[14] = new Point( 13,   0,    0) ;
    damProfile[15] = new Point( 7 ,  13,     0) ;
    damProfile[16] = new Point( 10,  13,    0) ;
    damProfile[17] = new Point( 5 ,  37,     0) ;
    damProfile[18] = new Point( 5 ,  60,     0) ;
    damProfile[19] = new Point( 7 ,  76,     0) ;
    damProfile[20] = new Point( 11,  96,    0) ;
    damProfile[21] = new Point( 23, 134,    0) ;
    
    std::vector<Point> damArch ;
    for(double i = -1. ; i <= 1. ; i+=0.05)
    {
        damArch.push_back(Point(i*100., 40.*(1.-i*i) , 0)+Point(0, 0, -54));
    }
    
    std::valarray<Point *> groundProfile(8) ;  
    groundProfile[0] = new Point( -130,160 ) ;
    groundProfile[1] = new Point( -90, 160 ) ;
    groundProfile[2] = new Point( -90 , 55 ) ;
    groundProfile[3] = new Point( 90 , 55 ) ;
    groundProfile[4] = new Point( 90 ,160 ) ;
    groundProfile[5] = new Point( 130 ,160 ) ;
    groundProfile[6] = new Point( 130 ,0 ) ;
    groundProfile[7] = new Point( -130 ,0 ) ;
    
    
    PolygonalSample3D ground(&samplers,groundProfile ,Point(0, -160,0), Point(0, 100, -50)) ;
//     PolygonalSample3D dam(&samplers,damProfile ,Point(50, 0,0), Point(-100, 0, -54)) ;
    LoftedPolygonalSample3D dam(&ground, damProfile, damArch) ;
//     std::vector<Point> bb = dam.getBoundingBox() ;
//     
//     for(auto & p : bb)
//         p.print() ;
//     exit(0) ;
//     
    dam.setBehaviour(/*new Viscoelasticity( PURE_ELASTICITY, m1)*/new Stiffness(m1)) ;
    ground.setBehaviour(/*new Viscoelasticity( PURE_ELASTICITY, m1)*/new Stiffness(m1*.4)) ;
        
    F.addFeature(&samplers, &ground) ;
    F.addFeature(&ground, &dam) ;
//         F.addFeature(&samplers, &dam) ;
    F.setSamplingFactor(&ground, 2);

    F.setDeltaTime(1.);
//     F.addFeature(&dam, &gallery0) ;
//     F.addFeature(&dam, &gallery1) ;
    F.setSamplingNumber(atof(argv[1])) ;
    F.setPartition(1);

    //fixed at the bottom
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, RIGHT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, TOP)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;

//     //water load
    waterload = new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_NORMAL_STRESS, dam.getPrimitive(), Point(0,1,0) , Function("9810 30 z -  *")*f_negativity(Function("z 30 -")) ) ;
    F.addBoundaryCondition(waterload) ;
//     
    //selfweight
    F.addBoundaryCondition(new GlobalBoundaryCondition( SET_VOLUMIC_STRESS_ZETA, -23544 ) );
    
    step(&F) ;

    return 0 ;
}
