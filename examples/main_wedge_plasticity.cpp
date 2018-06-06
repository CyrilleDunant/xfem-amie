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
#include "../physics/collisiondetectors/geometrybasedcontact.h"
#include "../physics/contactmodels/linearcontactforce.h"
#include "../physics/contactmodels/linearcontactdisplacement.h"
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


Rectangle rect(1, 1, 1.*.5+0.035*.5+0.00005, 1.*.5+0.00005) ;
BoundingBoxAndRestrictionDefinedBoundaryCondition * loadr = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_XI, RIGHT,-1, 1,0.00005, 1) ;
BoundingBoxAndRestrictionDefinedBoundaryCondition * loadf = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, RIGHT,-1, 1,0.00005, 1) ;
BoundingBoxDefinedBoundaryCondition * contact = new BoundingBoxDefinedBoundaryCondition(CONTACT_CONDITION, RIGHT, new GeometryBasedContact(&rect,0.024*1e-2), new LinearContactForce(&rect)) ;
// BoundingBoxDefinedBoundaryCondition * contact = new BoundingBoxDefinedBoundaryCondition(CONTACT_CONDITION, RIGHT, new GeometryBasedContact(&rect,0.024*1e-2), new LinearContactDisplacement(&rect)) ;
// BoundingBoxAndRestrictionDefinedBoundaryCondition * loadf = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, TOP,0.00005, 1,-1, 1) ;
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



void step ( size_t nsteps, Sample * samplef )
{

    int tries = 0 ;
    int every = 200 ;
    bool relaxed = false ;
    bool go_on = true ;
    bool setbc = true ;
    double tsteps = nsteps ;

    std::fstream ldfile ;
    ldfile.open("loadStrain_m.txt", std::ios_base::in|std::ios_base::out|std::ios_base::trunc) ;
    featureTree->step() ;
    ldfile << 1 <<"  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 <<  "  " << 0 << "  " << 0 << "  "  << 0 << "  " << 1 <<  std::endl ;
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

        /*if(false && go_on && tries%every == 0 && tries)
        {
            featureTree->removeBoundaryCondition( loadr );
            setbc = false ;
            relaxed = true ;
        }
        else */if(go_on)
        {
//             std::cout << "moving!" << std::endl ;
            rect.bottomLeft.print();
            transform(&rect, TRANSLATE, Point(-0.024/tsteps, 0.)) ;
            rect.bottomLeft.print();
//             rect.getCenter().print() ;
//             loadr->setData(loadr->getData()-/*1e9*1./nsteps*/0.024*1./tsteps) ;
//             if(!setbc)
//             {
//                 featureTree->addBoundaryCondition( loadr );
//                 setbc = true ;
//             }
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
        Point p(0.000,0.000) ;
        featureTree->get2DMesh()->getField ( PRINCIPAL_REAL_STRESS_FIELD, p, stemp ) ;
        featureTree->get2DMesh()->getField ( PRINCIPAL_IMPOSED_STRESS_FIELD, p, istemp ) ;
        featureTree->get2DMesh()->getField ( PRINCIPAL_MECHANICAL_STRAIN_FIELD,  p, etemp ) ;
        featureTree->get2DMesh()->getField ( PRINCIPAL_IMPOSED_STRAIN_FIELD,  p, ietemp ) ;

        displacements.push_back ( etemp[1] );
        displacementsx.push_back ( etemp[0] );
        loads.push_back ( stemp[1] );
        loadsx.push_back ( stemp[0] );
        std::cout << "\naverage sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "max score : " << featureTree->maxScore << std::endl ;
        std::cout << std::endl ;

        ldfile << go_on << "  " << stemp[0]/1e6 << "  " << stemp[1]/1e6 << "  " << etemp[0]*1e6 << "  " << etemp[1]*1e6 <<  "  "
               << istemp[0]/1e6 << "  " << istemp[1]/1e6 << "  "  << ietemp[0]*1e6 << "  " << ietemp[1]*1e6 << "  " << featureTree->getIterationCount() << std::endl ;


//         }

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
    ldfile.close() ;

}

int main ( int argc, char *argv[] )
{


    //http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1517-70762010000200028
    double mradius = .005 ; // .010 ;//

    // More or less a 5754 Al alloy
//     double nu = 0.2 ; //0.33 ;
    double E = 134e9 ; //70e9 ;
    double nu = .3 ;


    Sample samplef(0.035, 0.0279,  0, 0) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;


    StiffnessAndFracture  * pr = new StiffnessAndFracture(E, nu, new NonLocalDeviatoricVonMises(400e6, mradius),new PrandtlReussPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;
    StiffnessAndFracture  * pg = new StiffnessAndFracture(E, nu, new NonLocalDeviatoricVonMises(400e6, mradius),new PrandtlGrauertPlasticStrain(),SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;
    Stiffness  * sf = new Stiffness(E, nu) ;

    samplef.setBehaviour(sf);

    
//     F.addBoundaryCondition ( loadr );
//     F.addBoundaryCondition ( loadf );
    F.addBoundaryCondition ( contact );
    F.largeStrains = true ;
//     loadr->setActive(true);
// 	F.addBoundaryCondition(loadt);
    
    
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT ) ) ;
//     F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,TOP_LEFT ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;
    F.addPoint(new Point(0.035*.5, 0.00005)); 

    F.setOrder ( LINEAR ) ;
// F.addPoint(new Point(0, 0)) ;

    F.setMaxIterationsPerStep ( 10000 );
    F.thresholdScoreMet = 1e-2 ;


    step ( 100, &samplef ) ;


    return 0 ;
}
