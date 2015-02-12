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

    for(size_t i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
        featureTree->printReport();
    }
    VoxelWriter vw1("sphere_stiffness", 200) ;
    vw1.getField(featureTree, VWFT_STIFFNESS) ;
    vw1.write();

    VoxelWriter vw("sphere_stress", 200) ;
    vw.getField(featureTree, VWFT_STRESS) ;
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
    Sample3D samplers(nullptr, 300,300,300,0,0,0) ;

    FeatureTree F(&samplers) ;

    Matrix m0 =  Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    nu = 0.2 ;
    E = 40 ;
    Matrix m1 = Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);
    
    nu = 0.2 ;
    E = 80 ;
    Matrix m2 = Material::cauchyGreen(std::make_pair(E,nu), true,SPACE_THREE_DIMENSIONAL);

    samplers.setBehaviour(new Stiffness(m0)) ;

    std::valarray<Point *> pts(5) ;
    pts[0] = new Point(-20, 20, 0) ;
    pts[1] = new Point(20, 20, 0) ;
    pts[2] = new Point(5, 5, 0) ;
    pts[3] = new Point(20, -20, 0) ;
    pts[4] = new Point(-20, -20, 0) ;
    PolygonalSample3D inc(&samplers, pts, Point(200, 0, 200), Point(-100, 50, -100)) ;
    inc.setBehaviour(new Stiffness(m1)) ;
    
    std::valarray<Point *> pts0(5) ;
   
    pts0[1] = new Point(20, 20, 0) ;
    pts0[2] = new Point(20, -20, 0) ;
    pts0[3] = new Point(-20, -20, 0) ;
    pts0[4] = new Point(-5, 0, 0) ;
    pts0[0] = new Point(-20, 20, 0) ;
    PolygonalSample3D inc0(&samplers, pts0, Point(200, 0, -200), Point(-100, -50, 100)) ;
    inc0.setBehaviour(new Stiffness(m2)) ;


    F.addFeature(&samplers, &inc) ;
    F.addFeature(&samplers, &inc0) ;
    F.setSamplingNumber(atof(argv[1])) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, FRONT, -1.)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
//     F.setProjectionOnBoundaries(false) ;
    F.setOrder(QUADRATIC) ;

    step(&F) ;

    return 0 ;
}
