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
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/fractionmcft.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "../physics/fracturecriteria/druckerprager.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/damagemodels/plasticstrain.h"
#include "../physics/stiffness.h"
#include "../physics/materials/aggregate_behaviour.cpp"
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
#include "../utilities/itoa.h"
#include "../utilities/random.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
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

BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition ( SET_ALONG_XI, RIGHT,0 ) ;
// BoundingBoxDefinedBoundaryCondition * loadr = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP,0) ;

double factor = 25 ;
MinimumAngle cri ( M_PI/6. ) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

void step ( size_t nsteps )
{

    size_t nit = 2 ;
    size_t ntries = 5;
    size_t dsteps = 60 ;
    size_t tries = 0 ;
    size_t dit = 0 ;
    int totit = 0 ;
    for ( size_t v = 0 ; v < nsteps ; v++ )
    {
        tries = 0 ;

        tries++ ;

        bool go_on = featureTree->step() ;
        double appliedForce = loadr->getData() *effectiveRadius*2.*rebarDiametre;
        if ( go_on )
        {
// 			loadr->setData(sin(double(count)/40.)*5e-5) ;
// 			if(count < 80)
// 				loadr->setData(loadr->getData()-1e-5) ;
// 			else
// 				loadr->setData(loadr->getData()+1e-7) ;
            count++ ;
            loadr->setData ( loadr->getData() +6e-6 ) ;
// 			loadr->setData(loadr->getData()-1e-5) ;

// 			loadt->setData(0) ;
        }

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
            if ( k->getBehaviour()->type != VOID_BEHAVIOUR )
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
        std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax ( REAL_STRESS_FIELD ) ;
        std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax ( STRAIN_FIELD ) ;
        std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax ( VON_MISES_REAL_STRESS_FIELD ) ;
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD ) ;
        Vector etemp = featureTree->getAverageField ( STRAIN_FIELD ) ;

        std::cout << std::endl ;
        std::cout << "max value :" << x.max() << std::endl ;
        std::cout << "min value :" << x.min() << std::endl ;
        std::cout << "avg x value :" << xavg << std::endl ;
        std::cout << "avg y value :" << xavg << std::endl ;

        std::cout << "max sigma11 :" << stempm.second[0]/1e6  << std::endl ;
        std::cout << "min sigma11 :" << stempm.first[0]/1e6   << std::endl ;
        std::cout << "max sigma12 :" << stempm.second[2]/1e6  << std::endl ;
        std::cout << "min sigma12 :" << stempm.first[2]/1e6   << std::endl ;
        std::cout << "max sigma22 :" << stempm.second[1]/1e6  << std::endl ;
        std::cout << "min sigma22 :" << stempm.first[1]/1e6   << std::endl ;

        std::cout << "max epsilon11 :" << etempm.second[0]*1e6 << std::endl ;
        std::cout << "min epsilon11 :" << etempm.first[0]*1e6  << std::endl ;
        std::cout << "max epsilon12 :" << etempm.second[2]*1e6 << std::endl ;
        std::cout << "min epsilon12 :" << etempm.first[2]*1e6  << std::endl ;
        std::cout << "max epsilon22 :" << etempm.second[1]*1e6 << std::endl ;
        std::cout << "min epsilon22 :" << etempm.first[1]*1e6  << std::endl ;

        std::cout << "max von Mises :" << vmm.second[0]/1e6 << std::endl ;
        std::cout << "min von Mises :" << vmm.first[0]/1e6 << std::endl ;

        std::cout << "average sigma11 : " << stemp[0]/1e6 << std::endl ;
        std::cout << "average sigma22 : " << stemp[1]/1e6 << std::endl ;
        std::cout << "average sigma12 : " << stemp[2]/1e6 << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0]*1e6 << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1]*1e6 << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2]*1e6 << std::endl ;

        std::cout << std::endl ;


        if ( go_on )
        {
            std::cout << appliedForce/1000. << std::endl ;
            displacements.push_back ( etemp[1] );
            displacementsx.push_back ( etemp[0] );
            loads.push_back ( stemp[1] );
            loadsx.push_back ( stemp[0] );
        }

        std::fstream ldfile  ;
        ldfile.open ( "ldn", std::ios::out ) ;
        for ( int j = 0 ; j < loads.size() ; j++ )
        {
            ldfile << displacements[j] << "   " << loads[j] << "   " <<  displacementsx[j] << "   " << loadsx[j] <<  "\n" ;
        }
        ldfile.close();


        if ( true )
        {
            writer.reset ( featureTree ) ;
            writer.getField ( TWFT_CRITERION ) ;
            writer.getField ( TWFT_STIFFNESS_X ) ;
            writer.getField ( TWFT_STIFFNESS_Y ) ;
            writer.getField ( PRINCIPAL_ANGLE_FIELD ) ;
            writer.getField ( TWFT_DAMAGE ) ;
            writer.append() ;
        }
        if ( go_on )
        {
            writerc.reset ( featureTree ) ;
            writerc.getField ( TWFT_CRITERION ) ;
            writerc.getField ( TWFT_STIFFNESS_X ) ;
            writerc.getField ( TWFT_STIFFNESS_Y ) ;
            writerc.getField ( PRINCIPAL_ANGLE_FIELD ) ;
            writerc.getField ( TWFT_DAMAGE ) ;
            writerc.append() ;
        }
        //(1./epsilon11.getX())*( stressMoyenne.getX()-stressMoyenne.getY()*modulePoisson);
    }

}




int main ( int argc, char *argv[] )
{


    double compressionCrit = -32.6e6 ;
    double mradius = .1 ; // .010 ;//

    double nu = 0.2 ;
    double E_paste = 30e9 ;


// 	Sample samplef(0.3, 0.6,  0.15, 0.3) ;
    Sample samplef ( 0.6, 0.3,  0.3, 0.15 ) ;

    FeatureTree F ( &samplef ) ;
    featureTree = &F ;


    double verticaloffset = 0 ; //1./14.*samplef.height()
    double minDist = std::min ( samplef.width(), samplef.height() ) ;
    TriangularInclusion t0 ( Point ( minDist*.5,-minDist*.5-mradius+verticaloffset ),
                             Point ( minDist*.5,-minDist*.5+mradius+verticaloffset ),
                             Point ( -minDist*.5,minDist*.5+mradius+verticaloffset ) );
    TriangularInclusion t1 ( Point ( minDist*.5,-minDist*.5-mradius+verticaloffset ),
                             Point ( -minDist*.5,minDist*.5+mradius+verticaloffset ),
                             Point ( -minDist*.5,minDist*.5-mradius+verticaloffset ) );
    t0.isVirtualFeature = true ;
    t1.isVirtualFeature = true ;

    Sample r0 ( mradius*.25, samplef.height(),samplef.getCenter().getX(), samplef.getCenter().getY() ) ;
    r0.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*0.9,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
    dynamic_cast<ConcreteBehaviour *> ( r0.getBehaviour() )->materialRadius = mradius ;
    r0.isVirtualFeature = true ;
    r0.setBehaviourSource ( &samplef );
//     F.addFeature ( &samplef, &r0 );

    PlasticStrain * t0damagemodel = new PlasticStrain() ;
    PlasticStrain * t1damagemodel = new PlasticStrain() ;

// 	t0.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-20e6*.95, -20e6*.95,E_paste,0.1 , mradius),t0damagemodel));
// 	t1.setBehaviour( new StiffnessAndFracture(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) , new DruckerPrager(-20e6*.95, -20e6*.95,E_paste,0.1 , mradius),t1damagemodel));

    t0.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*.96,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) );
    t1.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit*.96,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) );
    t0.setBehaviourSource ( &samplef );
    t1.setBehaviourSource ( &samplef );
//     F.addFeature ( &samplef, &t0 );
//     F.addFeature ( &samplef, &t1 );
    
    Pore pore(mradius*.25, samplef.getCenter()) ;
    F.addFeature(&samplef, &pore);
    F.setSamplingFactor(&pore, 4);
// 	samplef.setBehaviour(new Stiffness(Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL, PLANE_STRAIN))) ;
    samplef.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
    dynamic_cast<ConcreteBehaviour *> ( samplef.getBehaviour() )->materialRadius = mradius ;

    F.addBoundaryCondition ( loadr );
// 	F.addBoundaryCondition(loadt);

    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT ) ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_ETA,BOTTOM ) ) ;

    F.setSamplingNumber ( atof ( argv[1] ) ) ;

    F.setOrder ( LINEAR ) ;
// F.addPoint(new Point(0, 0)) ;

    F.setMaxIterationsPerStep ( 3400 );

    step ( 10 ) ;


    return 0 ;
}
