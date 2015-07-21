// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
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
    Point crest(-74,0,20) ;
    
    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
        Vector tmp(6) ;
        double size = 350 ;
        for(double x = -size*.5 ;  x <= size*.5 ; x += 1)
        {
            for(double y = -size*.5 ;  y <= size*.5 ; y += 1)
            {
                featureTree->get3DMesh()->getField(STRAIN_FIELD, Point(0,y,x),tmp) ;
                    std::cout <<  tmp[0]<< "  "<< std::flush ;
            }
            std::cout << std::endl ;
        }
        for(double x = -size*.5 ;  x <= size*.5 ; x += 1)
        {
            for(double y = -size*.5 ;  y <= size*.5 ; y += 1)
            {
                featureTree->get3DMesh()->getField(REAL_STRESS_FIELD, Point(0,y,x),tmp) ;
                    std::cout << tmp[0]<< "  "<< std::flush ;
            }
            std::cout << std::endl ;
        }
        
        
            VoxelWriter vw1("stiffness", 50) ;
            vw1.getField(featureTree, VWFT_STIFFNESS) ;
            vw1.write();
            exit(0) ;

        if(i == 0)
        {
            VoxelWriter vw("stress_full", 50) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            VoxelWriter vw1("stiffness", 50) ;
            vw1.getField(featureTree, VWFT_STIFFNESS) ;
            vw1.write();
            waterload->setData(Function("9810 x 41 + *")*f_positivity(Function("x 41 +")));
            exit(0) ;
            
        }
        if(i == 1)
        {

            waterload->setData(Function("9810 x 52 + *")*f_positivity(Function("x 52 +")));
        }
        if(i == 2)
        {
            VoxelWriter vw("stress_intermediate", 50) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            waterload->setData(Function("9810 x 63 + *")*f_positivity(Function("x 63 +")));
        }
        if(i == 3)
        {
            waterload->setData(Function("9810 x 74 + *")*f_positivity(Function("x 74 +")));
        }
        if(i == 4 )
        {
            VoxelWriter vw("stress_full", 50) ;
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
    Matrix m1 = Tensor::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);
    
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
    damProfile[0] =  new Point(-148,  28,  0) ;
    damProfile[1] =  new Point(-148,  36,  0) ;
    damProfile[2] =  new Point(-134,  34,  0) ;
    damProfile[3] =  new Point(-114,  32,  0) ;
    damProfile[4] =  new Point(- 91,  32,  0) ;
    damProfile[5] =  new Point(- 75,  34,  0) ;
    damProfile[6] =  new Point(- 61,  37,  0) ;
    damProfile[7] =  new Point(- 44,  41,  0) ;
    damProfile[8] =  new Point(- 27,  50,  0) ;
    damProfile[9] =  new Point(- 27,  53,  0) ;
    damProfile[10] = new Point(- 23,  53,  0) ;
    damProfile[11] = new Point(- 18,  55,  0) ;
    damProfile[12] = new Point(- 18,  61,  0) ;
    damProfile[13] = new Point(-  9,  60,  0) ;
    damProfile[14] = new Point(-  0,  13,  0) ;
    damProfile[15] = new Point(- 13,  7,   0) ;
    damProfile[16] = new Point(- 13,  10,  0) ;
    damProfile[17] = new Point(- 37,  5,   0) ;
    damProfile[18] = new Point(- 60,  5,   0) ;
    damProfile[19] = new Point(- 76,  7,   0) ;
    damProfile[20] = new Point(- 96,  11,  0) ;
    damProfile[21] = new Point(-134,  23,  0) ;
    
    std::vector<Point> damArch ;
    for(double i = -1. ; i <= 1. ; i+=0.05)
    {
        damArch.push_back(Point(i*100.,  20.*(1.-i*i),76 ));
    }
    
    std::valarray<Point *> groundProfile(4) ;  
    groundProfile[0] = new Point( -37,-22 ) ;
    groundProfile[1] = new Point( -37, 81 ) ;
    groundProfile[2] = new Point( 20 , 81 ) ;
    groundProfile[3] = new Point( 20 ,-22 ) ;
    
    std::valarray<Point *> galleryProfile0(5) ;
    galleryProfile0[0] = new Point(-1., -1.) ;
    galleryProfile0[1] = new Point(1., -1.) ;
    galleryProfile0[2] = new Point(2, 0) ;
    galleryProfile0[3] = new Point(1., 1.) ;
    galleryProfile0[4] = new Point(-1., 1.) ;
    
    std::valarray<Point *> galleryProfile1(5) ;
    galleryProfile1[0] = new Point(-1., -1.) ;
    galleryProfile1[1] = new Point(1., -1.) ;
    galleryProfile1[2] = new Point(2, 0) ;
    galleryProfile1[3] = new Point(1., 1.) ;
    galleryProfile1[4] = new Point(-1., 1.) ;
    
    std::vector<Point> galleryArch0 ;
    for(double i = -.9 ; i <= .9 ; i+=0.05)
    {
        galleryArch0.push_back(Point(i*100.,  20.*(1.-i*i),-50 ));
    }

    std::vector<Point> galleryArch1 ;
    for(double i = -.9 ; i <= .9 ; i+=0.05)
    {
        galleryArch1.push_back(Point(i*100.,  20.*(1.-i*i) ,0));
    }
    
    
    PolygonalSample3D ground(&samplers,groundProfile ,Point(250, 0,0), Point(-125, 0, -140)) ;
    PolygonalSample3D dam(&samplers,damProfile ,Point(200, 0,0), Point(-100, 0, -74)) ;
//     LoftedPolygonalSample3D dam(&ground, damProfile,damArch) ;
    
    dam.setBehaviour(/*new Viscoelasticity( PURE_ELASTICITY, m1)*/new Stiffness(m1)) ;
    ground.setBehaviour(/*new Viscoelasticity( PURE_ELASTICITY, m1)*/new Stiffness(m1*.4)) ;
    
    LoftedPolygonalSample3D gallery0(&dam, galleryProfile0, galleryArch0) ;
    gallery0.setBehaviour(new VoidForm()) ;
    
    LoftedPolygonalSample3D gallery1(&dam, galleryProfile1, galleryArch1) ;
    gallery1.setBehaviour(new VoidForm()) ;
    
    F.addFeature(&samplers, &ground) ;
    F.addFeature(&ground, &dam) ;
    F.setSamplingFactor(&ground, 1.5);
    F.setDeltaTime(1.);
//     F.addFeature(&dam, &gallery0) ;
//     F.addFeature(&dam, &gallery1) ;
    F.setSamplingNumber(atof(argv[1])) ;
    F.setPartition(1);

    //fixed at the bottom
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, RIGHT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, FRONT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, BACK)) ;

    //water load
    waterload = new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_NORMAL_STRESS, dam.getPrimitive(), Point(0,0,1) , Function("9810 z 30 +  *")*f_positivity(Function("z 30 +")) ) ;
    F.addBoundaryCondition(waterload) ;
    waterload->setData(Function("9810 z 74 + *")*f_positivity(Function("z 74 +")));
    
    //selfweight
    F.addBoundaryCondition(new GlobalBoundaryCondition( SET_VOLUMIC_STRESS_ZETA, 23544 ) );
    
    step(&F) ;

    return 0 ;
}
