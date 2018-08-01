// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../features/sample.h"
#include "../features/crack.h"
#include "../features/incompatibleModes.h"

#include "../physics/materials/paste_behaviour.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_with_large_deformation.h"

#include "../utilities/writer/triangle_writer.h"

using namespace Amie ;

FeatureTree * featureTree ;


int main(int argc, char *argv[])
{
	int nsteps = 100 ;

	double width = 0.12;
	double height = 0.0025;
	Sample sample(width, height , 0,0) ;
    Point tb(width*.5-1e-5, 0) ;

	featureTree = new FeatureTree(&sample) ;
	
	Point * a = new Point(-0.025,0) ;
	Point * b = new Point(0.025,0) ;
	BranchedCrack * Crack1 = new BranchedCrack( a,  b);
	Crack1->setEnrichementRadius(0.01);

// 	featureTree->addFeature(&sample, Crack1);//MY

// 	sample.setBehaviour(new OrthotropicStiffness(E_paste, E_paste*.5,  E_paste*.5/(2.*1-nu*0.5),  nu, M_PI*.15)) ;
	sample.setBehaviour(new Stiffness/*WithLargeDeformation*/(10000, 0.4999)) ;// 0.499:0.000517251
    

//     BoundingBoxNearestNodeDefinedBoundaryCondition *force = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_STRESS_XI , RIGHT, tb, 0) ;
    
    BoundingBoxDefinedBoundaryCondition *force = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI , RIGHT, -0.00001) ;
    featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_RIGHT)) ;
    
    
    featureTree->addBoundaryCondition(force) ;
//     featureTree->addBoundaryCondition(new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_FORCE_XI , BOTTOM_RIGHT, tb, -1e-4)) ;
//     featureTree->addBoundaryCondition(new GlobalBoundaryCondition(SET_VOLUMIC_STRESS_ETA , -1e-6)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_LEFT)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
//     featureTree->addFeature(nullptr,new IncompatibleModes()) ;

	featureTree->setSamplingNumber(atof(argv[1])) ;
	featureTree->setOrder(LINEAR) ;
	featureTree->setMaxIterationsPerStep(10);

	
	MultiTriangleWriter writerm( "displacements_enrichment_head", "displacements_enrichment", nullptr ) ;
//     featureTree->largeStrains = true ;
    bool go_on = true ; 
    
    for(int v = 0 ; v < nsteps ; v++)
    {
       
        if(go_on)
           force->setData(-0.01*v/nsteps) ;// force->setData(-(v+0.05*((double)rand()/RAND_MAX-0.5))*0.00001) ; //
        
        do
        {
            go_on = featureTree->step() ;
            
        } while (!go_on) ;
        
        Vector dp = featureTree->getDisplacements(&tb) ;
        Vector str = featureTree->getAverageField(REAL_FINITE_STRESS_FIELD) ;
        Vector stri = featureTree->getAverageField(REAL_STRESS_FIELD) ;
        std::cout << featureTree->getIterationCount() << "  " << featureTree->getResidualError() << "  "<< featureTree->volumeVariation() << "  "
                  << dp[0]<<"  "<<  dp[1]<< "  " 
                  << str[0]<< "  "<< str[1]<< "  " 
                  << stri[0]<< "  "<< stri[1]<<"  "<< stri[2]<< "  "<< force->getData() << std::endl ;
        
                  

//         if(go_on)
//         {
            writerm.reset( featureTree ) ;
//             writerm.getField( TWFT_LARGE_DEFORMATION_TRANSFORM ) ;
            writerm.getField(REAL_FINITE_STRESS_FIELD) ;
            writerm.append() ;
//         }

        
//         if(!go_on)
//             break ;
    }
        
//     for(int i = 0 ; i < 10 ; i++)
//     {
//         force->setData(force->getData()*2.) ;
//         for(int v = 0 ; v < nsteps ; v++)
//         {
//             std::cout << featureTree->residualError << std::endl ;
//     // 		Crack1->print() ;
//             bool go_on = featureTree->step() ;
//             
//             if(v%40 == 0)
//             {
//                 writerm.reset( featureTree ) ;
//                 writerm.getField( TWFT_LARGE_DEFORMATION_TRANSFORM ) ;
//                 writerm.getField( TWFT_LARGE_DEFORMATION_ANGLE ) ;
//                 writerm.append() ;
//             }
// 
//             
//             if(!go_on)
//                 break ;
//         }            
// 
//     }
        

	return 0 ;
}
