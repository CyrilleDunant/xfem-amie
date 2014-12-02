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
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
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

    int nsteps = 200 ;
    featureTree->setMaxIterationsPerStep(100) ;
    MultiTriangleWriter writerc( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
    for(size_t i = 0 ; i < nsteps ; )
    {
        double converged = featureTree->step() ;

        Vector x = featureTree->getDisplacements() ;

        if(converged)
        {
            i++ ;
            for(size_t j = 0 ; j < pockets.size() ; j++)
                dynamic_cast<ExpansiveZone *>(pockets[j])->setRadius(pockets[j]->getRadius()+0.000001) ;
        }
        std::cout << "unknowns :" << x.size() << std::endl ;

        int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;

        double volume = 0 ;

        double xavg = 0 ;

        for(auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++)
        {
            if(k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                double ar = k->area() ;
                volume += ar ;
                for(size_t l = 0 ; l < npoints ; l++)
                {
                    xavg += x[k->getBoundingPoint(l).getId()*2]*ar/npoints ;
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
        std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
        std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;
        std::cout << "max sigma12 :" << stempm.second[2]  << std::endl ;
        std::cout << "min sigma12 :" << stempm.first[2]   << std::endl ;

        std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
        std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
        std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
        std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;        
        std::cout << "max epsilon12 :" << etempm.second[2] << std::endl ;
        std::cout << "min epsilon12 :" << etempm.first[2]  << std::endl ;

        std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
        std::cout << "min von Mises :" << vmm.first[0] << std::endl ;

        std::cout << "average sigma11 : " << stemp[0] << std::endl ;
        std::cout << "average sigma22 : " << stemp[1] << std::endl ;
        std::cout << "average sigma12 : " << stemp[2] << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2] << std::endl ;

        if(converged)
        {
            writerc.reset( featureTree ) ;
            writerc.getField( REAL_STRESS_FIELD ) ;
            writerc.getField( STRAIN_FIELD ) ;
            writerc.getField( TWFT_DAMAGE ) ;
            writerc.append() ;
//             writerc.writeSvg(0., true) ;
        }
    }

    exit(0) ;
}

double partitionScore(const std::vector<int> & triplet)
{
    double avg = (triplet[0]+triplet[1]+triplet[2])/3. ;
    return (triplet[0]-avg)*(triplet[0]-avg)+(triplet[1]-avg)*(triplet[1]-avg)+(triplet[2]-avg)*(triplet[2]-avg) ;
}

int main(int argc, char *argv[])
{

    Sample samplers(nullptr, 0.02,0.02,0,0) ;

    FeatureTree F(&samplers) ;
    featureTree = &F ;
    

    samplers.setBehaviour(new PasteBehaviour()) ;
    Vector a(0.,6) ;// a[0] = 1 ; a[1] = 1 ; a[2] = 1 ;
// 	ExpansiveZone3D inc(&samplers,100, 200, 200, 200, m1*4, a) ;
    Inclusion inc( 0.008, 0, 0) ;

// 	inc->setBehaviour(new StiffnessWithImposedDeformation(m1*4.,a)) ;
// 	inc.setBehaviour(new Stiffness(m1*4)) ;
    inc.setBehaviour(new AggregateBehaviour()) ;


    F.addFeature(&samplers, &inc) ;
// 	F.addFeature(&samplers, inc0) ;
    F.setSamplingNumber(atof(argv[1])) ;
    
    std::vector<Geometry *> inclusion ;
    inclusion.push_back(new Circle(0.007, 0, 0));
    for (size_t i = 0 ; i < 1 ; i++)
    {
        ExpansiveZone * p = new ExpansiveZone(nullptr, 0.00000001, 0, 0.006, new GelBehaviour()) ;
        pockets.push_back(p);
        featureTree->addFeature(&inc, p);
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

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, RIGHT, -5e6)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP, -5e6)) ;
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
    F.setPartition(64);

    step() ;
// 	delete dt ;

    return 0 ;
}
