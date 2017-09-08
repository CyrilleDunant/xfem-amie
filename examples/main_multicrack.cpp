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

#include "../utilities/writer/triangle_writer.h"

using namespace Amie ;

FeatureTree * featureTree ;


int main(int argc, char *argv[])
{
	int nsteps = 1 ;
	double nu = 0.2 ;
	double E_paste = 30e9 ;

	double width = 0.12;
	double height = 0.02;
	Sample sample(width, height , width*.5, height*.5) ;

	featureTree = new FeatureTree(&sample) ;
	
	Point * a = new Point(-0.025,0) ;
	Point * b = new Point(0.025,0) ;
	BranchedCrack * Crack1 = new BranchedCrack( a,  b);
	Crack1->setEnrichementRadius(0.01);

// 	featureTree->addFeature(&sample, Crack1);//MY

// 	sample.setBehaviour(new OrthotropicStiffness(E_paste, E_paste*.5,  E_paste*.5/(2.*1-nu*0.5),  nu, M_PI*.15)) ;
	sample.setBehaviour(new Stiffness(100, 0.499)) ;
    

    BoundingBoxDefinedBoundaryCondition *force = new BoundingBoxDefinedBoundaryCondition(SET_FORCE_ETA , RIGHT, -.0006) ;
    
//  	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , RIGHT)) ;
	featureTree->addBoundaryCondition(force) ;
//     featureTree->addBoundaryCondition(new GlobalBoundaryCondition(SET_VOLUMIC_STRESS_ETA , -1e6)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , LEFT)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
//     featureTree->addFeature(nullptr,new IncompatibleModes()) ;

	featureTree->setSamplingNumber(atof(argv[1])) ;
// 	featureTree->setOrder(QUADRATIC) ;
	featureTree->setMaxIterationsPerStep(1);

	
	MultiTriangleWriter writerm( "displacements_enrichment_head", "displacements_enrichment", nullptr ) ;
    featureTree->tgsolve = true ;

    for(int v = 0 ; v < nsteps ; v++)
    {
        std::cout << featureTree->residualError << std::endl ;
// 		Crack1->print() ;
        bool go_on = featureTree->step() ;

//         if(v%40 == 0 )
        {
            writerm.reset( featureTree ) ;
            writerm.getField( TWFT_LARGE_DEFORMATION_TRANSFORM ) ;
            writerm.getField( TWFT_LARGE_DEFORMATION_ANGLE ) ;
            writerm.getField( TWFT_LARGE_DEFORMATION_FORCES ) ;
            writerm.append() ;
        }
        
        if(!go_on)
            break ;
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
