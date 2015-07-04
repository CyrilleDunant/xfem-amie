// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../features/sample.h"
#include "../features/polygonSample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../utilities/writer/triangle_writer.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#define DEBUG


using namespace Amie ;

std::vector<Feature *> pockets ;
FeatureTree * featureTree ;

void step()
{

    int nsteps = 1 ;
    featureTree->setMaxIterationsPerStep(800) ;
    MultiTriangleWriter writerc( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
    for(int i = 0 ; i < nsteps ; )
    {
        bool converged = featureTree->step() ;

        if(converged)
        {
            i++ ;
            for(size_t j = 0 ; j < pockets.size() ; j++)
                dynamic_cast<ExpansiveZone *>(pockets[j])->setRadius(pockets[j]->getRadius()+0.000001) ;
        }
       featureTree->printReport();
        writerc.reset( featureTree ) ;
        writerc.getField( REAL_STRESS_FIELD ) ;
        writerc.getField( STRAIN_FIELD ) ;
//             writerc.getField( TWFT_DAMAGE ) ;
        writerc.getField( TWFT_STIFFNESS ) ;
        writerc.append() ;
//             writerc.writeSvg(0., true) ;
    }

    exit(0) ;
}

double partitionScore(const std::vector<int> & triplet)
{
    double avg = (triplet[0]+triplet[1]+triplet[2])/3. ;
    return (triplet[0]-avg)*(triplet[0]-avg)+(triplet[1]-avg)*(triplet[1]-avg)+(triplet[2]-avg)*(triplet[2]-avg) ;
}

// positions should be 
// 0 0 
// 0 0.004 + sym + quad
// 0 0.006 + sym + quad
// 0 0.008 + sym + quad

//tests should be 5 5 
//   0  0
//  -5 -5
//  -10 -10
//  -15 -15
//  -5 -10
//  -5 -15
//  -10 -15
int main(int argc, char *argv[])
{

    Sample samplers(nullptr, 200,200,0,0) ;

    FeatureTree F(&samplers) ;
    featureTree = &F ;
    int setup = -1 ; //atoi(argv[2]) ;
    int load = 1 ; //atoi(argv[3]) ;

    samplers.setBehaviour(new ElasticOnlyPasteBehaviour()) ;
    Vector a(0.,6) ;// a[0] = 1 ; a[1] = 1 ; a[2] = 1 ;
// 	ExpansiveZone3D inc(&samplers,100, 200, 200, 200, m1*4, a) ;
    ExpansiveZone inc(&samplers, 50, 0, 0,new ElasticOnlyAggregateBehaviour()) ;
//     std::valarray<Point *> pts(5) ;
    
//     pts[0] = new Point(-80, 50) ; 
//     pts[1] = new Point(0.00, 20) ; 
//     pts[2] = new Point(80, 80) ;
//     pts[3] = new Point(50, -80) ; 
//     pts[4] = new Point(-80, -80) ;
//     PolygonalSample inc(&samplers, pts) ;
    
// 	inc->setBehaviour(new StiffnessWithImposedDeformation(m1*4.,a)) ;
// 	inc.setBehaviour(new Stiffness(m1*4)) ;
    


    F.addFeature(&samplers, &inc) ;
// 	F.addFeature(&samplers, inc0) ;
    F.setSamplingNumber(atof(argv[1])) ;
    
    std::vector<Geometry *> inclusion ;
    inclusion.push_back(new Circle(0.007, 0, 0));
//     for (size_t i = 0 ; i < 1 ; i++)
    if(setup == 0)
    {
        ExpansiveZone * p = new ExpansiveZone(nullptr, 0.00000001, 0, 0.00, new GelBehaviour()) ;
        pockets.push_back(p);
        featureTree->addFeature(&inc, p);
    }
    if(setup == 1)
    {
        ExpansiveZone * p = new ExpansiveZone(nullptr, 0.00000001, 0, 0.004, new GelBehaviour()) ;
        pockets.push_back(p);
        featureTree->addFeature(&inc, p);
    }
    if(setup == 10)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.004, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.004, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
    }
    if(setup == 100)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.004, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.004, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
        ExpansiveZone * p2 = new ExpansiveZone(nullptr, 0.00000001, 0.004, 0.00, new GelBehaviour()) ;
        pockets.push_back(p2);
        featureTree->addFeature(&inc, p2);
        ExpansiveZone * p3 = new ExpansiveZone(nullptr, 0.00000001, 0.004, 0.00, new GelBehaviour()) ;
        pockets.push_back(p3);
        featureTree->addFeature(&inc, p3);
    }
    if(setup == 2)
    {
        ExpansiveZone * p = new ExpansiveZone(nullptr, 0.00000001, 0, 0.006, new GelBehaviour()) ;
        pockets.push_back(p);
        featureTree->addFeature(&inc, p);
    }
    if(setup == 20)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.006, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.006, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
    }
    if(setup == 200)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.006, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.006, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
        ExpansiveZone * p2 = new ExpansiveZone(nullptr, 0.00000001, 0.006, 0.00, new GelBehaviour()) ;
        pockets.push_back(p2);
        featureTree->addFeature(&inc, p2);
        ExpansiveZone * p3 = new ExpansiveZone(nullptr, 0.00000001, 0.006, 0.00, new GelBehaviour()) ;
        pockets.push_back(p3);
        featureTree->addFeature(&inc, p3);
    }
    if(setup == 3)
    {
        ExpansiveZone * p = new ExpansiveZone(nullptr, 0.00000001, 0, 0.008, new GelBehaviour()) ;
        pockets.push_back(p);
        featureTree->addFeature(&inc, p);
    }
    if(setup == 30)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.008, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.008, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
    }
    if(setup == 300)
    {
        ExpansiveZone * p0 = new ExpansiveZone(nullptr, 0.00000001, 0, 0.008, new GelBehaviour()) ;
        pockets.push_back(p0);
        featureTree->addFeature(&inc, p0);
        ExpansiveZone * p1 = new ExpansiveZone(nullptr, 0.00000001, 0, -0.008, new GelBehaviour()) ;
        pockets.push_back(p1);
        featureTree->addFeature(&inc, p1);
        ExpansiveZone * p2 = new ExpansiveZone(nullptr, 0.00000001, 0.008, 0.00, new GelBehaviour()) ;
        pockets.push_back(p2);
        featureTree->addFeature(&inc, p2);
        ExpansiveZone * p3 = new ExpansiveZone(nullptr, 0.00000001, 0.008, 0.00, new GelBehaviour()) ;
        pockets.push_back(p3);
        featureTree->addFeature(&inc, p3);
    }
    
//     placement2DInInclusions(samplers.getPrimitive(),inclusion, pockets ) ;
// 	F.setProjectionOnBoundaries(false) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_RIGHT_BACK)) ;
//
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP_LEFT_BACK)) ;
//
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_FRONT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_FRONT)) ;

    if(load == 1)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -5e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -5e6)) ;
    }
    if(load == 2)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -10e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -10e6)) ;
    }
    if(load == 3)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -15e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -15e6)) ;
    }
    if(load == 4)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -5e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -10e6)) ;
    }
    if(load == 5)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -5e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -15e6)) ;
    }
    if(load == 6)
    {
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -10e6)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -15e6)) ;
    }
//     F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(SET_NORMAL_STRESS, inc.getPrimitive(), -1)) ;

// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;

// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, TOP_LEFT_BACK)) ;

// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.setOrder(LINEAR) ;
//     F.setPartition(64);

    step() ;
// 	delete dt ;

    return 0 ;
}
