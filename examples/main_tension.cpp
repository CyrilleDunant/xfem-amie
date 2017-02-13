// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/radialstiffnessgradient.h"
#include "../physics/physics_base.h"
#include "../physics/homogenization/composite.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/fractionmcft.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/fracturecriteria/druckerprager.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/damagemodels/plasticstrain.h"
#include "../physics/damagemodels/prandtlreussplasticity.h"
#include "../physics/stiffness.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/fraction_stiffness_and_fracture.h"
#include "../physics/void_form.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/layeredinclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/optimizer.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/orthotropicstiffness.h"


#include <fstream>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>

#define DEBUG



using namespace Amie ;




FeatureTree * featureTree ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.05e-07 ;

double effectiveRadius = 0.13506098 ; //.5*.165*sqrt(M_PI*.25) ;//
double rebarDiametre = 0.0254*sqrt ( M_PI*.25 ) ;

double delta_displacement =  1e-5 ;
double displacement_tolerance = 0.05*delta_displacement ;
double softeningFactor = 1. ;

double percent = 0.01 ;
double displacement  = 0 ;
double prescribedDisplacement = 0;
double derror = 0 ;
double ierror = 0 ;
double preverror = 0 ;
bool firstRun = true ;

double sampleLength = 5.5 ;
double sampleHeight = 1.2 ;
double supportLever = 2.5 ;
double supportMidPointToEndClearance = 0.25 ;
double platewidth = 0.15 ;
double plateHeight = 0.051 ;
double rebarEndCover = 0.047 ;

std::vector< double > loads ;
std::vector< double > displacements ;
std::vector< double > loadsx ;
std::vector< double > displacementsx ;
std::vector< double > damages ;
std::vector< double > cmod ;
std::vector< double > idx ;

Vector fracCrit ( 0 ) ;

Vector b ( 0 ) ;
Vector x ( 0 ) ;
Vector sigma ( 0 ) ;
Vector sigma11 ( 0 ) ;
Vector sigma22 ( 0 ) ;
Vector sigma12 ( 0 ) ;
Vector epsilon ( 0 ) ;
Vector epsilon11 ( 0 ) ;
Vector epsilon22 ( 0 ) ;
Vector epsilon12 ( 0 ) ;
Vector vonMises ( 0 ) ;
Vector angle ( 0 ) ;

// BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition ( SET_ALONG_XI, RIGHT,0 ) ;
BoundingBoxNearestNodeDefinedBoundaryCondition * loadr = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_XI, RIGHT,Point(.2, 0), 0) ;

double factor = 25 ;
MinimumAngle cri ( M_PI/6. ) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;
double pmax = 0 ;

MultiTriangleWriter writer ( "triangles_stiff_head", "triangles_stiff_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;
MultiTriangleWriter writerr ( "triangles_relaxed_head", "triangles_relaxed_layers", nullptr ) ;

void step ( size_t nsteps, Sample * samplef )
{

    size_t tries = 0 ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {

        bool go_on = featureTree->step() ;
	bool relaxed = tries%2 == 0 && tries > 6 ;
	if(relaxed)
	  featureTree->addBoundaryCondition( loadr );
	
        if ( go_on )
        {
	  tries++ ;
	  count++ ;

	  pmax = std::min(loadr->getData(), pmax) ;
	  loadr->setData(pmax-.25e-2) ;
	  if(tries%2 == 0 && tries > 6)
	    featureTree->removeBoundaryCondition( loadr );

        }
        else
            nsteps++ ;

        x.resize ( featureTree->getDisplacements().size() ) ;
        x = featureTree->getDisplacements() ;

        double volume = 0 ;
        double xavg = 0 ;
        double yavg = 0 ;
        std::fstream dfile ;
        if ( go_on )
        {
            dfile.open ( "dprofile", std::ios::out | std::ios::app ) ;
        }
        for ( auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++ )
        {
            if ( k->getBehaviour()->getDamageModel() && k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                if ( go_on )
                {
                    dfile <<  k->getCenter().getX() << "  "<< k->getBehaviour()->getDamageModel()->getState().max() << "  " ;
                }
                double ar = k->area() ;
                volume += ar ;
                for ( size_t l = 0 ; l < k->getBoundingPoints().size() ; l++ )
                {
                    xavg += x[k->getBoundingPoint ( l ).getId() *2]*ar/ k->getBoundingPoints().size() ;
                    yavg += x[k->getBoundingPoint ( l ).getId() *2+1]*ar/ k->getBoundingPoints().size() ;
                }
            }
        }
        if ( go_on )
        {
            dfile << "\n" ;
        }
        dfile.close();

        xavg /= volume ;
        yavg /= volume ;
//         std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax ( REAL_STRESS_FIELD ) ;
//         std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax ( STRAIN_FIELD ) ;
//         std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax ( VON_MISES_REAL_STRESS_FIELD ) ;
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD ) ;
        Vector etemp = featureTree->getAverageField ( TOTAL_STRAIN_FIELD ) ;

// 	Vector tmp(3) ;
// 	for(double x = 0 ;  x <= 0.3 ; x += .001)
//         {
//             for(double y = 0 ;  y <= 0.3 ; y += .001)
//             {
//                 featureTree->get2DMesh()->getField(PRINCIPAL_REAL_STRESS_FIELD, Point(x,y,0),tmp) ;
//                     std::cout <<  tmp[0]<< "  "<< std::flush ;
//             }
//             std::cout << std::endl ;
//         }
//         std::cout << std::endl ;
//         std::cout << "max value :" << x.max() << std::endl ;
//         std::cout << "min value :" << x.min() << std::endl ;
//         std::cout << "avg x value :" << xavg << std::endl ;
//         std::cout << "avg y value :" << xavg << std::endl ;
// 
//         std::cout << "max sigma11 :" << stempm.second[0]/1e6  << std::endl ;
//         std::cout << "min sigma11 :" << stempm.first[0]/1e6   << std::endl ;
//         std::cout << "max sigma12 :" << stempm.second[2]/1e6  << std::endl ;
//         std::cout << "min sigma12 :" << stempm.first[2]/1e6   << std::endl ;
//         std::cout << "max sigma22 :" << stempm.second[1]/1e6  << std::endl ;
//         std::cout << "min sigma22 :" << stempm.first[1]/1e6   << std::endl ;
// 
//         std::cout << "max epsilon11 :" << etempm.second[0]*1e6 << std::endl ;
//         std::cout << "min epsilon11 :" << etempm.first[0]*1e6  << std::endl ;
//         std::cout << "max epsilon12 :" << etempm.second[2]*1e6 << std::endl ;
//         std::cout << "min epsilon12 :" << etempm.first[2]*1e6  << std::endl ;
//         std::cout << "max epsilon22 :" << etempm.second[1]*1e6 << std::endl ;
//         std::cout << "min epsilon22 :" << etempm.first[1]*1e6  << std::endl ;
// 
//         std::cout << "max von Mises :" << vmm.second[0]/1e6 << std::endl ;
//         std::cout << "min von Mises :" << vmm.first[0]/1e6 << std::endl ;

        std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2]*1e6 << std::endl ;

        std::cout << std::endl ;


        if ( go_on )
        {
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[0] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[0] );
        }

        std::fstream ldfile  ;
        ldfile.open ( "ldn", std::ios::out ) ;
        for ( size_t j = 0 ; j < loads.size() ; j++ )
        {
            ldfile << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] <<  "\n" ;
        }
        ldfile.close();


	if (!relaxed)
	{
	  writer.reset ( featureTree ) ;
	  writer.getField ( TWFT_CRITERION ) ;
	  writer.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
	  writer.getField ( MECHANICAL_STRAIN_FIELD ) ;
	  writer.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
  // 	writer.getField ( IMPOSED_STRAIN_FIELD ) ;
  //             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
  //             writer.getField ( TWFT_DAMAGE ) ;
	  writer.append() ;
	}

        if ( go_on && !relaxed)
	{
	    writerc.reset ( featureTree ) ;
            writerc.getField ( TWFT_CRITERION ) ;
            writerc.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
            writerc.getField ( MECHANICAL_STRAIN_FIELD ) ;
	    writerc.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
// 	    writerc.getField ( IMPOSED_STRAIN_FIELD ) ;
//             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//             writer.getField ( TWFT_DAMAGE ) ;
            writerc.append() ;
	}
	if(relaxed)
	{
	    writerr.reset ( featureTree ) ;
            writerr.getField ( TWFT_CRITERION ) ;
            writerr.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
            writerr.getField ( MECHANICAL_STRAIN_FIELD ) ;
	    writerr.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
// 	    writerc.getField ( IMPOSED_STRAIN_FIELD ) ;
//             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//             writer.getField ( TWFT_DAMAGE ) ;
            writerr.append() ;
	}
//         if ( go_on )
//         {
//             writerc.reset ( featureTree ) ;
//             writerc.getField ( TWFT_CRITERION ) ;
//             writerc.getField ( TWFT_STIFFNESS_X ) ;
//             writerc.getField ( TWFT_STIFFNESS_Y ) ;
//             writerc.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//             writerc.getField ( TWFT_DAMAGE ) ;
//             writerc.append() ;
//         }
        //(1./epsilon11.getX())*( stressMoyenne.getX()-stressMoyenne.getY()*modulePoisson);
	    
// 	  samplef->setBehaviour(new Stiffness(30e9,0.1+0.1*v)) ;
    }

}




int main ( int argc, char *argv[] )
{

//   for(double soft = 0 ; soft < .6 ; soft += .05)
//   {
//     double E_agg=59e9 ;
//     double E_paste=12e9 ;
//     AggregateBehaviour * agg = new AggregateBehaviour(true, false,E_agg, 0.3) ;
//     PasteBehaviour * paste = new PasteBehaviour(true, false, E_paste*(1.-soft), 0.3*(1.-soft)) ;
//     
//     Phase matrix(paste, 0.3) ;
//     
//     Phase aggregate(agg,0.7 ) ;
//     
//     MoriTanakaMatrixInclusionComposite mt(matrix,aggregate) ;
//     std::cout << mt.getBehaviour()->getTensor(Point())[0][0] << std::endl ;
//     delete agg ;
//     delete paste ;
//   }
//   exit(0) ;
  

    double compressionCrit = -32.6e6 ;
    double mradius = .05 ; // .010 ;//

    // More or less a 5754 Al alloy
    double nu = 0.33 ;
    double E = 70e9 ;


	Sample samplef(0.05, 1.2,  0, 0) ;
//     Sample samplef ( 100, 100,  50, 50 ) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;


    double verticaloffset = 0. ;
    double minDist = std::min ( samplef.width(), samplef.height() ) ;
    TriangularInclusion t0 ( Point ( minDist*.5,-minDist*.5-mradius*.5 ),
                             Point ( minDist*.5,-minDist*.5+mradius*.5 ),
                             Point ( -minDist*.5,minDist*.5-mradius*.5 ) );
    TriangularInclusion t1 ( Point ( minDist*.5,-minDist*.5+mradius*.5 ),
                             Point ( -minDist*.5,minDist*.5+mradius*.5 ),
                             Point ( -minDist*.5,minDist*.5-mradius*.5 ) );
    t0.isVirtualFeature = true ;
    t1.isVirtualFeature = true ;

//     Sample r0 ( mradius, samplef.height(),samplef.getCenter().getX(), samplef.getCenter().getY() ) ;
//     r0.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*0.8,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
//     dynamic_cast<ConcreteBehaviour *> ( r0.getBehaviour() )->materialRadius = mradius ;
//     r0.isVirtualFeature = true ;
//     r0.setBehaviourSource ( &samplef );
//     F.addFeature ( &samplef, &r0 );

    PrandtlReussPlasticStrain * t0damagemodel = new PrandtlReussPlasticStrain() ;
    PrandtlReussPlasticStrain * t1damagemodel = new PrandtlReussPlasticStrain() ;
	t0.setBehaviour( new StiffnessAndFracture(E,nu, new NonLocalVonMises(30e6*.9, mradius),new PrandtlReussPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRESS));
	t1.setBehaviour( new StiffnessAndFracture(E,nu, new NonLocalVonMises(30e6*.9, mradius),new PrandtlReussPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRESS));
	samplef.setBehaviour(new StiffnessAndFracture(E,nu, new NonLocalVonMises(30e6, mradius),new PrandtlReussPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRESS));

//     t0.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*.96,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) );
//     t1.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*.96,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) );
    t0.setBehaviourSource ( &samplef );
    t1.setBehaviourSource ( &samplef );
//     F.addFeature ( &samplef, &t0 );
//     F.addFeature ( &samplef, &t1 );
    
    
    // 1:10 -> -0.063126
    // 1:5  -> -0.137859
    //             0         0.25       0.5       1         2      
    // 1:2  -> -0.341575  -0.558602 -0.630257 -0.689162 -0.784332 
    // 1:1  -> -0.412712
//     EllipsoidalInclusion inclusion(&samplef, samplef.getCenter(), Point(0.,samplef.width()*.1), Point(samplef.height()*.1/2, 0.)) ;
//     F.addFeature(&samplef, &inclusion);
// //     F.setSamplingFactor(&inclusion, 4);
//     inclusion.setBehaviour(new Stiffness(E_paste*4, nu)) ;
// //     samplef.setBehaviour(new Stiffness(E_paste,0)) ;
//     samplef.setBehaviour ( new PasteBehaviour(false, false, E_paste, nu,  3e6, 3.6e9,4e9, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, 1, 0.) );

//     dynamic_cast<ConcreteBehaviour *> ( samplef.getBehaviour() )->materialRadius = mradius ;

    F.addBoundaryCondition ( loadr );
// 	F.addBoundaryCondition(loadt);

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, BOTTOM ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, TOP ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,TOP ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;

    F.setOrder ( QUADRATIC ) ;
// F.addPoint(new Point(0, 0)) ;

    F.setMaxIterationsPerStep ( 5000 );

    step ( 300, &samplef ) ;


    return 0 ;
}
