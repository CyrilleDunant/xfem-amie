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
#include "../physics/damagemodels/prandtlgrauertplasticstrain.h"
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

BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition ( SET_ALONG_ETA, TOP,0 ) ;
// BoundingBoxNearestNodeDefinedBoundaryCondition * loadr = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_XI, RIGHT,Point(.02, 0), 0) ;

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

void step ( size_t nsteps, RectangularFeature * samplef )
{

    int tries = 0 ;
    int every = 200 ;
    bool relaxed = false ;
    bool go_on = true ;
    bool setbc = true ;
    


    std::fstream ldfile ;
    ldfile.open("loadStrain_m.txt", std::ios_base::in|std::ios_base::out|std::ios_base::trunc) ;
    featureTree->step() ;
    ldfile << 1 <<"  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 <<  "  " << 0 << "  " << 0 << "  "  << 0 << "  " << 1 <<  "  " << 0 <<  "  " << 0 << "  " << 0 << "  "  << 0 << "  " <<std::endl ;
    
    int centrecache = featureTree->get2DMesh()->generateCache(new Circle(0.0025, 0, 0)) ;
    Point p(0,0) ;
    DelaunayTriangle * centreElement = featureTree->get2DMesh()->getUniqueConflictingElement(&p) ;
    
    writer.reset ( featureTree ) ;
    writer.getField ( TWFT_CRITERION ) ;
    writer.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
    writer.getField ( TOTAL_STRAIN_FIELD ) ;
    writer.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
//     writer.getField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
// 	writer.getField ( IMPOSED_STRAIN_FIELD ) ;
//             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
    writer.getField ( TWFT_STIFFNESS ) ;

//             writer.getField ( TWFT_DAMAGE ) ;
    writer.append() ;

    writerc.reset ( featureTree ) ;
    writerc.getField ( TWFT_CRITERION ) ;
    writerc.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
    writerc.getField ( TOTAL_STRAIN_FIELD ) ;
    writerc.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
//     writerc.getField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
    // 	    writerc.getField ( IMPOSED_STRAIN_FIELD ) ;
    //             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//               writer.getField ( TWFT_DAMAGE ) ;
    writerc.getField ( TWFT_STIFFNESS ) ;
    writerc.append() ;

    for ( size_t v = 0 ; v < nsteps ; v++ )
    {

        if(false && go_on && tries%every == 0 && tries)
        {
            featureTree->removeBoundaryCondition( loadr );
            setbc = false ;
            relaxed = true ;
        }
        else if(go_on)
        {
//             if(v < 30)
                loadr->setData(loadr->getData()+1e-5/**.005*/) ; //e-4*.25 is mostly OK
//             else
//                 loadr->setData(loadr->getData()+2e-4/**.005*/) ; //e-4*.25 is mostly OK
            if(!setbc)
            {
                featureTree->addBoundaryCondition( loadr );
                setbc = true ;
            }
            relaxed = false ;
        }

        go_on = featureTree->step() ;

        if ( go_on )
            tries++ ;
        else
            nsteps++ ;

//         if(false &&go_on && (tries+1)%every == 0 && tries)
//             featureTree->thresholdScoreMet = 0.00001 ;
//         else
//             featureTree->thresholdScoreMet = 0.01 ;

        double volume = 0 ;
        double xavg = 0 ;
        double yavg = 0 ;

//         if ( go_on )
//         {
        Vector stemp(2) ;//= featureTree->getAverageField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
        Vector istemp(2); //= featureTree->getAverageField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
        Vector etemp(2) ;//= featureTree->getAverageField ( PRINCIPAL_MECHANICAL_STRAIN_FIELD ) ;
        Vector ietemp(2); //= featureTree->getAverageField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
        Vector astemp = featureTree->getAverageField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
        Vector aetemp = featureTree->getAverageField ( PRINCIPAL_TOTAL_STRAIN_FIELD ) ;
        stemp = featureTree->get2DMesh()->getSmoothedField ( PRINCIPAL_REAL_STRESS_FIELD, centrecache, centreElement,0. ) ;
        istemp = featureTree->get2DMesh()->getSmoothedField ( PRINCIPAL_IMPOSED_STRESS_FIELD, centrecache, centreElement,0.  ) ;
        etemp = featureTree->get2DMesh()->getSmoothedField ( PRINCIPAL_MECHANICAL_STRAIN_FIELD, centrecache, centreElement,0.  ) ;
        ietemp = featureTree->get2DMesh()->getSmoothedField ( PRINCIPAL_IMPOSED_STRAIN_FIELD, centrecache, centreElement,0.  ) ;

        displacements.push_back ( etemp[1] );
        displacementsx.push_back ( etemp[0] );
        loads.push_back ( stemp[1] );
        loadsx.push_back ( stemp[0] );
        std::cout << "\naverage sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << std::endl ;

        ldfile << go_on << "  "                                 //1
               << stemp[0]/1e6 << "  "                          //2
               << stemp[1]/1e6 << "  "                          //3
               << etemp[0]*1e6 << "  "                          //4
               << etemp[1]*1e6 <<  "  "                         //5
               << istemp[0]/1e6 << "  "                         //6
               << istemp[1]/1e6 << "  "                         //7
               << ietemp[0]*1e6 << "  "                         //8
               << ietemp[1]*1e6 << "  "                         //9
               << featureTree->getIterationCount() << "  "      //10
               << astemp[0]/1e6 << "  "                         //11
               << astemp[1]/1e6 << "  "                         //12
               << aetemp[0]*1e6 << "  "                         //13
               << aetemp[1]*1e6 << std::endl ;                  //14


//         }

               if(v%5 == 0)
               {
// 	if (!relaxed)
// 	{
        writer.reset ( featureTree ) ;
        writer.getField ( TWFT_CRITERION ) ;
        writer.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
        writer.getField ( TOTAL_STRAIN_FIELD ) ;
        writer.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
// 	  writer.getField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
        // 	writer.getField ( IMPOSED_STRAIN_FIELD ) ;
        //             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
        writer.getField ( TWFT_STIFFNESS ) ;
        writer.append() ;
// 	}

        if ( go_on && !relaxed)
        {
            writerc.reset ( featureTree ) ;
            writerc.getField ( TWFT_CRITERION ) ;
            writerc.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
            writerc.getField ( TOTAL_STRAIN_FIELD ) ;
            writerc.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
// 	    writerc.getField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
            writer.getField ( TWFT_STIFFNESS ) ;
// 	    writerc.getField ( IMPOSED_STRAIN_FIELD ) ;
//             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//             writer.getField ( TWFT_DAMAGE ) ;
            writerc.append() ;
        }
        if(go_on && relaxed)
        {
            writerr.reset ( featureTree ) ;
            writerr.getField ( TWFT_CRITERION ) ;
            writerr.getField ( PRINCIPAL_REAL_STRESS_FIELD ) ;
            writerr.getField ( TOTAL_STRAIN_FIELD ) ;
            writerr.getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD ) ;
// 	    writerr.getField ( PRINCIPAL_IMPOSED_STRESS_FIELD ) ;
            writer.getField ( TWFT_STIFFNESS ) ;
// 	    writerc.getField ( IMPOSED_STRAIN_FIELD ) ;
//             writer.getField ( PRINCIPAL_STRESS_ANGLE_FIELD ) ;
//             writer.getField ( TWFT_DAMAGE ) ;
            writerr.append() ;
        }
               }

    }
    ldfile.close() ;

}

int main ( int argc, char *argv[] )
{
//     double nu = 0.25 ;
//     double nu_mat = .01 ;
//     InclusionGeometryType t = INCLUSION_IS_SPHERE ;
//     Stiffness * agg = new Stiffness(3.2, nu, SPACE_THREE_DIMENSIONAL) ;
//     Stiffness * paste = new Stiffness(3.2e-8, nu_mat, SPACE_THREE_DIMENSIONAL) ;
//     
//     double a = 1. ; double b = 1./10  ; double c = 1./10 ;
//     
//     Phase matrix(paste,0, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
//     Phase aggregate(agg, 1.-0, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
//     
//     
//     BiphasicSelfConsistentComposite hint (matrix, aggregate) ;
//     hint.apply() ;
//     double dphi = .1 ;
//     Matrix dC = hint.C ;
//     Matrix dC1 = hint.C ;
//     for(double soft = .0 ; soft <= .6 ; soft += dphi)
//     {
// 
//         Phase matrix(paste,soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
//         Phase aggregate(agg, 1.-soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
// 
//         
//         Matrix C = BiphasicSelfConsistentComposite(matrix , aggregate, hint).getBehaviour()->getTensor(Point()) ;
//     //       Matrix C =   ReussMatrixInclusionComposite(aggregate,matrix).getBehaviour()->getTensor(Point())  ;
//         double nu_eff = C[0][1]/((C[0][0]+C[0][1])) ;
//         double E_eff = C[0][0]*(1.+nu_eff)*(1.-2.*nu_eff)/(1.-nu_eff) ;
//         std::cout << soft <<"  "<< E_eff*10 <<"  " <<nu_eff << std::endl ;
//         
//     //     C.print() ;
//     //     dC1 = dC ;
//     //     dC = C-hint.C ;
//         hint.C = C;
//     } 
//     dphi= 0.01 ;
//     for(double soft = .6 ; soft <= .85 ; soft += dphi)
//     {
// 
//         Phase matrix(paste,soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
//         Phase aggregate(agg, 1.-soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
// 
//         
//         Matrix C = BiphasicSelfConsistentComposite(matrix , aggregate, hint).getBehaviour()->getTensor(Point()) ;
//     //       Matrix C =   ReussMatrixInclusionComposite(aggregate,matrix).getBehaviour()->getTensor(Point())  ;
//         double nu_eff = C[0][1]/((C[0][0]+C[0][1])) ;
//         double E_eff = C[0][0]*(1.+nu_eff)*(1.-2.*nu_eff)/(1.-nu_eff) ;
//         std::cout << soft <<"  "<< E_eff*10 <<"  " <<nu_eff << std::endl ;
//         
//     //     C.print() ;
//     //     dC1 = dC ;
//     //     dC = C-hint.C ;
//         hint.C = C;
//     } 
//     dphi= 0.01 ;
//     for(double soft = .85 ; soft <= 1.01 ; soft += dphi)
//     {
// 
//         Phase matrix(paste,soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
//         Phase aggregate(agg, 1.-soft, SPACE_THREE_DIMENSIONAL, t, a, b, c) ;
// 
//         
//         Matrix C = BiphasicSelfConsistentComposite(matrix , aggregate, hint).getBehaviour()->getTensor(Point()) ;
//     //       Matrix C =   ReussMatrixInclusionComposite(aggregate,matrix).getBehaviour()->getTensor(Point())  ;
//         double nu_eff = C[0][1]/((C[0][0]+C[0][1])) ;
//         double E_eff = C[0][0]*(1.+nu_eff)*(1.-2.*nu_eff)/(1.-nu_eff) ;
//         std::cout << soft <<"  "<< E_eff*10 <<"  " <<nu_eff << std::endl ;
//         
//     //     C.print() ;
//     //     dC1 = dC ;
//     //     dC = C-hint.C ;
//         hint.C = C;
//     } 
//   
//   
//   delete agg ;
//   delete paste ;
//   exit(0) ;


    double compressionCrit = -32.6e6 ;

    //http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1517-70762010000200028
    double mradius = .0016 ; // .010 ;//

    // More or less a 5754 Al alloy
//     double nu = 0.2 ; //0.33 ;
    double E = 30e9 ; //70e9 ;


    RectangularFeature samplef(0.005, 0.021,  0, 0) ;
// 	Sample samplef(0.01, 0.01,  0, 0) ;
//     Sample samplef ( 100, 100,  50, 50 ) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;
//     F.setProjectionOnBoundaries(false) ;


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

    RectangularFeature r0 ( 0.0015, 0.009,+0.001+0.0015*.5, 0 ) ;
    r0.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &r0 );

    Inclusion i0(0.0015, samplef.width()*.5, 0.0045) ;
    i0.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &i0 );

    Inclusion i1(0.0015, samplef.width()*.5, -0.0045) ;
    i1.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &i1 );

    RectangularFeature r1 ( 0.0015, 0.009,-0.001-0.0015*.5, 0 ) ;
    r1.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &r1 );

    Inclusion i2(0.0015, -samplef.width()*.5, 0.0045) ;
    i2.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &i2 );

    Inclusion i3(0.0015, -samplef.width()*.5, -0.0045) ;
    i3.setBehaviour ( new  VoidForm() ) ;
    F.addFeature ( &samplef, &i3 );


    // More or less a 5754 Al alloy
    StiffnessAndFracture  * pr = new StiffnessAndFracture(69e9,0.33, new NonLocalDeviatoricVonMises(220e6, mradius),new PrandtlReussPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;
    StiffnessAndFracture  * pg = new StiffnessAndFracture(69e9,0.2, new NonLocalDeviatoricVonMises(220e6, mradius),new PrandtlGrauertPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;

    samplef.setBehaviour(pg);
//     samplef.setBehaviour(new ConcreteBehaviour());



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
    F.largeStrains = true ;
    loadr->setActive(true);
// 	F.addBoundaryCondition(loadt);

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM ) ) ;
//     F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA, TPO ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI,BOTTOM_LEFT ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, TOP_LEFT ) ) ;
//     F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,TOP_LEFT ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;


    F.setOrder ( LINEAR ) ;
    //Sample samplef(0.005, 0.021,  0, 0) ;    
    F.addRefinementZone(new Circle(0.003,0,0));
    F.addRefinementZone(new Circle(0.0015,0,0));
// F.addPoint(new Point(0, 0)) ;

    F.setMaxIterationsPerStep ( 1000 );
    F.thresholdScoreMet = 0.0001 ;


    step ( 60, &samplef ) ;


    return 0 ;
}
