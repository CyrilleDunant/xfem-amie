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

    for(size_t i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
//         featureTree->printReport();
    }
    VoxelWriter vw1("sphere_stiffness", 150) ;
    vw1.getField(featureTree, VWFT_STIFFNESS) ;
    vw1.write();

    VoxelWriter vw("sphere_stress", 100) ;
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

    std::valarray<Point *> pts(4) ;
    pts[0] = new Point(0, -20, 0) ;
    pts[1] = new Point(150, -20, 0) ;
    pts[2] = new Point(150, 20, 0) ;
    pts[3] = new Point(0, 20, 0) ;
    std::vector<Point> ipts ;
    ipts.push_back(Point(-100,0 ,  0));
    ipts.push_back(Point(-70 ,8,  0));
    ipts.push_back(Point(-40 ,15,  0));
    ipts.push_back(Point(-20 ,18,  0));
    ipts.push_back(Point( 0  ,20,  0));
    ipts.push_back(Point( 20 ,18,  0));
    ipts.push_back(Point( 40 ,15,  0));
    ipts.push_back(Point( 70 ,8,  0));
    ipts.push_back(Point( 100,0 ,  0));
    LoftedPolygonalSample3D inc(&samplers, pts,ipts) ;
//     inc.isVirtualFeature = true ;
    inc.setBehaviour(new Stiffness(m1)) ;
    

    F.addFeature(&samplers, &inc) ;
    F.setSamplingNumber(atof(argv[1])) ;
    F.setSamplingFactor ( &inc, 1 ) ;
    

//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -1.)) ;
//     F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -1.)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM)) ;
    
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, LEFT)) ;
    
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, RIGHT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, RIGHT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, RIGHT)) ;
    Vector loadv(3, 0.) ; loadv[1] = 1 ;
//     F.addBoundaryCondition(new  GlobalForceBoundaryCondition(loadv)) ;
//     F.setProjectionOnBoundaries(false) ;
    F.setOrder(LINEAR) ;

    step(&F) ;

    return 0 ;
}
