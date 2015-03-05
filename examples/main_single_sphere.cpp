// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/laplacian.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/polygonSample3d.h"
#include "../features/loftedPolygonSample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../utilities/writer/voxel_writer.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#define DEBUG


using namespace Amie ;




void step(FeatureTree * featureTree)
{

    int nsteps = 1;// number of steps between two clicks on the opengl thing
    featureTree->setMaxIterationsPerStep(50) ;

    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
//         featureTree->printReport();
    }
    VoxelWriter vw1("sphere_stiffness", 200) ;
    vw1.getField(featureTree, VWFT_STIFFNESS) ;
    vw1.write();

    VoxelWriter vw("sphere_stress", 200) ;
    vw.getField(featureTree, VWFT_PRINCIPAL_STRESS) ;
    vw.write();
// 	VoxelWriter vw0("sphere_strain", 50) ;
// 	vw0.getField(featureTree, VWFT_STRAIN) ;
// 	vw0.write();
    exit(0) ;
}

int main(int argc, char *argv[])
{

    double nu = 0.2 ;
    double E = 20 ;
    Sample3D samplers(nullptr, 250,250,250,0,0,0) ;

    FeatureTree F(&samplers) ;

    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    nu = 0.2 ;
    E = 40 ;
    Matrix m1 = Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);
    
    nu = 0.2 ;
    E = 80 ;
    Matrix m2 = Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    samplers.setBehaviour(new VoidForm()) ;

    std::valarray<Point *> damProfile(4) ;
    damProfile[0] = new Point(0, 40, 0) ;
    damProfile[1] = new Point(150, 20, 0) ;
    damProfile[2] = new Point(150, -20, 0) ;
    damProfile[3] = new Point(0, -20, 0) ;
    
    std::vector<Point> damArch ;
    for(double i = -1. ; i <= 1. ; i+=0.2)
    {
        damArch.push_back(Point(0,i*100.,  20.*(1.-i*i) ));
    }
    
     std::valarray<Point *> galleryProfile(5) ;
    galleryProfile[0] = new Point(-1.5, -1.5, 0) ;
    galleryProfile[1] = new Point(1.5, -1.5, 0) ;
    galleryProfile[2] = new Point(2.2, 0, 0) ;
    galleryProfile[3] = new Point(1.5, 1.5, 0) ;
    galleryProfile[4] = new Point(-1.5, 1.5, 0) ;
    
    std::vector<Point> galleryArch ;
    for(double i = -.9 ; i <= .9 ; i+=0.2)
    {
        galleryArch.push_back(Point(-50,i*100.,  20.*(1.-i*i) ));
    }

    LoftedPolygonalSample3D dam(&samplers, damProfile,damArch) ;
    dam.setBehaviour(new Stiffness(m1)) ;
    
    LoftedPolygonalSample3D gallery(&dam, galleryProfile, galleryArch) ;
    gallery.setBehaviour(new VoidForm()) ;
    
    F.addFeature(&samplers, &dam) ;
    F.addFeature(&dam, &gallery) ;
    F.setSamplingNumber(atof(argv[1])) ;


    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ALL, RIGHT)) ;

    std::pair<Point, Point> normals = dam.getEndNormals() ;
    F.addBoundaryCondition(new GeometryAndFaceDefinedSurfaceBoundaryCondition(FIX_ALONG_ALL, dam.getPrimitive(), normals.first)) ;
    F.addBoundaryCondition(new GeometryAndFaceDefinedSurfaceBoundaryCondition(FIX_ALONG_ALL, dam.getPrimitive(), normals.second)) ;
    F.addBoundaryCondition(new GeometryDefinedBoundaryCondition( SET_VOLUMIC_STRESS_XI, dam.getPrimitive() , -1 ) );
//     F.setProjectionOnBoundaries(false) ;
    F.setOrder(LINEAR) ;

    step(&F) ;

    return 0 ;
}
