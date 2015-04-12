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
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>

using namespace Amie ;

GeometryAndFaceDefinedSurfaceBoundaryCondition * waterload ;

void step(FeatureTree * featureTree, LoftedPolygonalSample3D * dam)
{
    Sample3D samplers(nullptr, 249.99,249.99,249.99,0,0,0) ;
    int nsteps = 5;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;
    Point crest(-74,0,20) ;
    
    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;

        if(i == 0)
        {
            VoxelWriter vw("stress_empty", 200) ;
            vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
            vw.write();
            VoxelWriter vw1("stiffness", 200) ;
            vw1.getField(featureTree, VWFT_STIFFNESS) ;
            vw1.write();
            waterload->setData(Function("9810 x 41 + *")*f_positivity(Function("x 41 +")));
            
        }
        if(i == 1)
        {

            waterload->setData(Function("9810 x 52 + *")*f_positivity(Function("x 52 +")));
        }
        if(i == 2)
        {
            VoxelWriter vw("stress_intermediate", 200) ;
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
    Sample3D samplers(nullptr, 250,250,250,0,0,0) ;

    FeatureTree F(&samplers) ;

    nu = 0.2 ;
    E = 37e9 ;
    Matrix m1 = Tensor::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);
    
//     for(double phi = 0 ;  phi <= 1 ; phi +=0.01)
//     {
//         Phase porosity(new Stiffness(Tensor::cauchyGreen(std::make_pair(0,.4997), true,SPACE_THREE_DIMENSIONAL)), 1.-phi) ;
//         Phase csh(new Stiffness(Tensor::cauchyGreen(std::make_pair(1,.25), true,SPACE_THREE_DIMENSIONAL)), phi) ;
//         BiphasicSelfConsistentComposite sc(porosity,csh) ;
//         std::cout << phi << "  "<< sc.getBehaviour()->getTensor(Point())[0][0]/((1+.25)*(1.-2*.25)/(1-.25)) << std::endl ;
//     }
//     exit(0) ;
    
    samplers.setBehaviour(new VoidForm()) ;

    std::valarray<Point *> damProfile(4) ;
    damProfile[0] = new Point(0, 8, 0) ;
    damProfile[1] = new Point(150, 8, 0) ;
    damProfile[2] = new Point(150, -8, 0) ;
    damProfile[3] = new Point(0, -17, 0) ;
    
    std::vector<Point> damArch ;
    for(double i = -1. ; i <= 1. ; i+=0.05)
    {
        damArch.push_back(Point(0,i*100.,  20.*(1.-i*i) ));
    }
    
    std::valarray<Point *> galleryProfile0(5) ;
    galleryProfile0[0] = new Point(-1., -1., 0) ;
    galleryProfile0[1] = new Point(1., -1., 0) ;
    galleryProfile0[2] = new Point(2, 0, 0) ;
    galleryProfile0[3] = new Point(1., 1., 0) ;
    galleryProfile0[4] = new Point(-1., 1., 0) ;
    
    std::valarray<Point *> galleryProfile1(5) ;
    galleryProfile1[0] = new Point(-1., -1., 0) ;
    galleryProfile1[1] = new Point(1., -1., 0) ;
    galleryProfile1[2] = new Point(2, 0, 0) ;
    galleryProfile1[3] = new Point(1., 1., 0) ;
    galleryProfile1[4] = new Point(-1., 1., 0) ;
    
    std::vector<Point> galleryArch0 ;
    for(double i = -.9 ; i <= .9 ; i+=0.05)
    {
        galleryArch0.push_back(Point(-50,i*100.,  20.*(1.-i*i) ));
    }

    std::vector<Point> galleryArch1 ;
    for(double i = -.9 ; i <= .9 ; i+=0.05)
    {
        galleryArch1.push_back(Point(0,i*100.,  20.*(1.-i*i) ));
    }
    
    LoftedPolygonalSample3D dam(&samplers, damProfile,damArch) ;
    
    dam.setBehaviour(/*new Viscoelasticity( PURE_ELASTICITY, m1)*/new Stiffness(m1)) ;
    
    LoftedPolygonalSample3D gallery0(&dam, galleryProfile0, galleryArch0) ;
    gallery0.setBehaviour(new VoidForm()) ;
    
    LoftedPolygonalSample3D gallery1(&dam, galleryProfile1, galleryArch1) ;
    gallery1.setBehaviour(new VoidForm()) ;
    
    F.addFeature(&samplers, &dam) ;
    F.setDeltaTime(1.);
    F.addFeature(&dam, &gallery0) ;
    F.addFeature(&dam, &gallery1) ;
    F.setSamplingNumber(atof(argv[1])) ;
    F.setPartition(1);

    //fixed at the bottom
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, RIGHT)) ;

    //fixed by the mountain
    std::pair<Point, Point> normals = dam.getEndNormals() ;
    F.addBoundaryCondition(new GeometryAndFaceDefinedSurfaceBoundaryCondition(FIX_ALONG_ALL, dam.getPrimitive(), normals.first)) ;
    F.addBoundaryCondition(new GeometryAndFaceDefinedSurfaceBoundaryCondition(FIX_ALONG_ALL, dam.getPrimitive(), normals.second)) ;

    //water load
    waterload = new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_NORMAL_STRESS, dam.getPrimitive(), Point(0,0,1) , Function("9810 x 30 +  *")*f_positivity(Function("x 30 +")) ) ;
    F.addBoundaryCondition(waterload) ;
    
    //selfweight
    F.addBoundaryCondition(new GlobalBoundaryCondition( SET_VOLUMIC_STRESS_XI, 23544 ) );
    
    step(&F, &dam) ;

    return 0 ;
}
