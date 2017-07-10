// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/void_form.h"
#include "../features/sample.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/concrete_behaviour.h"

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
#include <limits>
#include <time.h>

using namespace Amie ;




FeatureTree * featureTree ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.05e-07 ;

double delta_displacement =  1e-5 ;
double displacement_tolerance = 0.05 * delta_displacement ;
double softeningFactor = 1. ;

double percent = 0.01 ;
double displacement  = 0 ;
double prescribedDisplacement = 0;
double derror = 0 ;
double ierror = 0 ;
double preverror = 0 ;
bool firstRun = true ;

double sampleLength = 3.9 ; //5.5 ;
double sampleHeight = 1.2 ;
double supportLever = 1.7 ;//2.5 ;
double supportMidPointToEndClearance = 0.25 ;
double platewidth = 0.15 ;
double plateHeight = 0.051 ;
double rebarDiametre = 0.025 ; //sqrt(506e-6);//sqrt( 4.*0.000506/M_PI ) ;
double rebarEndCover = 0.047 ;
double phi = 3.*rebarDiametre  ;
double psi = 2.*0.0084261498  ;
bool haveStirrups = false ;

std::vector<std::pair<double, double> > load_displacement ;
std::vector< double > loads ;
std::vector< double > deltas ;
std::vector< double > displacements ;
std::vector< double > damages ;
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

MultiTriangleWriter writer ( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc ( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

Function loadFunction ( "0" ) ;
// BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_STRESS_ETA, TOP, -platewidth, platewidth, -10, 10, loadFunction ) ;
BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition ( SET_ALONG_ETA, TOP, -platewidth, platewidth, -10, 10, 0. ) ;

//  BoundingBoxNearestNodeDefinedBoundaryCondition * load = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_ETA, TOP, Point(0., sampleHeight*.5)) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP,0) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, 0) ;

Rectangle bcbox ( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) ;
GeometryDefinedBoundaryCondition selfload ( SET_VOLUMIC_STRESS_ETA, &bcbox , -9025.2 ) ;
GeometryDefinedBoundaryCondition shrinkagey ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
GeometryDefinedBoundaryCondition shrinkagex ( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
int count = 0 ;
double aggregateArea = 0;


void step()
{

    size_t nsteps = 600*4 ; //16*10;
    double delta_d = 5.*0.0175e-3 ;

    for ( size_t v = 0 ; v < nsteps ; v++ )
    {
        y_max = 0 ;
        x_max = 0 ;
        y_min = 0 ;
        x_min = 0 ;

        bool go_on = true ;

        go_on = featureTree->step() ;
        if ( go_on )
        {
            load->setData ( load->getData()-delta_d ) ;
        }

        x.resize ( featureTree->getDisplacements ( -1, false ).size() ) ;
        x = featureTree->getDisplacements ( -1, false ) ;

        Vector forces ( featureTree->getAssembly()->getForces() ) ;

        std::cerr << "unknowns :" << x.size() << std::endl ;


        int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;

        double delta = 0;
        double deltacount = 0;
        double volume = 0 ;
        double xavg = 0 ;
        double yavg = 0 ;

        for ( auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++ )
        {
            if ( k->getBehaviour()->type != VOID_BEHAVIOUR )
            {

                double ar = k->area() ;
                volume += ar ;
                for ( int p = 0 ; p < npoints ; p++ )
                {
                    xavg += x[k->getBoundingPoint ( p ).getId() *2]*ar/npoints ;
                    yavg += x[k->getBoundingPoint ( p ).getId() *2+1]*ar/npoints ;


                    if ( dist ( Point ( supportLever, -sampleHeight*.5 + 0.064 ), k->getBoundingPoint ( p ) ) < .1 )
                    {
                        deltacount++ ;
                        delta += x[k->getBoundingPoint ( p ).getId() * 2] ;
                    }
                }

            }
        }

        xavg /= volume ;
        yavg /= volume ;
        delta /= deltacount ;
        std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax ( REAL_STRESS_FIELD ) ;
        std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax ( TOTAL_STRAIN_FIELD ) ;
        std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax ( VON_MISES_REAL_STRESS_FIELD ) ;
        Vector stemp = featureTree->getAverageField ( REAL_STRESS_FIELD ) ;
        Vector etemp = featureTree->getAverageField ( TOTAL_STRAIN_FIELD ) ;

        std::cout << std::endl ;
        std::cout << "max value :" << x.max() << std::endl ;
        std::cout << "min value :" << x.min() << std::endl ;
        std::cout << "avg x value :" << xavg << std::endl ;
        std::cout << "avg y value :" << xavg << std::endl ;

        std::cout << "max sigma11 :" << stempm.second[0]  << std::endl ;
        std::cout << "min sigma11 :" << stempm.first[0]   << std::endl ;
        std::cout << "max sigma12 :" << stempm.second[2]  << std::endl ;
        std::cout << "min sigma12 :" << stempm.first[2]   << std::endl ;
        std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
        std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;

        std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
        std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
        std::cout << "max epsilon12 :" << etempm.second[2] << std::endl ;
        std::cout << "min epsilon12 :" << etempm.first[2]  << std::endl ;
        std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
        std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;

        std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
        std::cout << "min von Mises :" << vmm.first[0] << std::endl ;

        std::cout << "average sigma11 : " << stemp[0] << std::endl ;
        std::cout << "average sigma22 : " << stemp[1] << std::endl ;
        std::cout << "average sigma12 : " << stemp[2] << std::endl ;
        std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
        std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
        std::cout << "average epsilon12 : " << etemp[2] << std::endl ;

        std::cout << std::endl ;

        if ( go_on )
        {
            displacements.push_back ( 1000.* ( load->getData() +delta_d ) );
            loads.push_back ( stemp[1]/1000. );
            deltas.push_back ( delta/deltacount );
            damages.push_back ( featureTree->averageDamage );
        }


        if ( go_on )
        {
            std::cout << stemp[1]/1000. << "  " << displacements.back() << "  " << damages.back() << std::endl ;
        }


        std::fstream ldfile ( "ldn", std::ios::out )  ;
        for ( size_t j = 0 ; j < loads.size() ; j++ )
        {
            ldfile << displacements[j] << "   " << loads[j] << "   " << damages[j] << "   " << deltas[j] << "\n" ;

        }
        if ( !go_on )
        {
            ldfile <<  1000.* ( load->getData() ) << "   " << stemp[1]/1000. << "   " << featureTree->averageDamage << "   " << delta/deltacount << "\n" ;
        }
        ldfile.close();


        if ( true )
        {
            writer.reset ( featureTree ) ;
            writer.getField ( TWFT_PRINCIPAL_STRESS ) ;
            writer.getField ( TWFT_PRINCIPAL_STRAIN ) ;
            writer.getField ( TWFT_CRITERION ) ;
            writer.getField ( TWFT_STIFFNESS_X ) ;
            writer.getField ( TWFT_STIFFNESS_Y ) ;
            writer.getField ( TWFT_DAMAGE ) ;
            writer.append() ;
        }

        if ( go_on )
        {
            writerc.reset ( featureTree ) ;
            writerc.getField ( TWFT_PRINCIPAL_STRESS ) ;
            writerc.getField ( TWFT_PRINCIPAL_STRAIN ) ;
            writerc.getField ( TWFT_CRITERION ) ;
            writerc.getField ( TWFT_STIFFNESS_X ) ;
            writerc.getField ( TWFT_STIFFNESS_Y ) ;
            writerc.getField ( TWFT_DAMAGE ) ;
            writerc.append() ;
//             writerc.writeSvg ( 0. ) ;
        }
// 		if ( !go_on )
// 			break ;

    }
}


int main ( int argc, char *argv[] )
{

    double softeningFactor = M_PI*.24;

    double  samplingNumber = atof ( argv[1] ) ;
    haveStirrups = atof ( argv[2] ) ;
    sampleLength = atof ( argv[3] ) ;
    sampleHeight = atof ( argv[4] ) ;
    supportLever = sampleLength*.5-.250 ;

    std::cout << sampleLength << "  " << supportLever << std::endl ;

    double compressionCrit = -34.2e6 ;

    std::cout << "phi = "<< phi << ", psi = " << psi << std::endl ;
// 	double mradius = 0.1; //0.015 ;//0.055 ;//.11 ; // .015
// 	double nradius = mradius*2.5 ;

    double E_steel = 200e9 ;
    double nu_steel = 0.01 ;
    double nu = 0.3 ;
    double E_paste = 37e9 ;
    double E_steel_effective = M_PI*0.5*rebarDiametre*rebarDiametre*E_steel/(rebarDiametre*rebarDiametre)*.75/**.8*(phi/0.4)*/ ;

    double halfSampleOffset = sampleLength*.25 ;

    Matrix m0_paste = Tensor::cauchyGreen ( E_paste, nu, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON ) ;

// 	//redimensionned so that we get in shear the right moment of inertia
// 	Matrix m0_steel = Tensor::orthothropicCauchyGreen(E_steel, E_steel, E_steel*(1.-nu_steel)*.5*.13/(1.-nu_steel*nu_steel), nu_steel,PLANE_STRESS_FREE_G) ;
//
    Matrix m0_steel = Tensor::cauchyGreen ( E_steel,nu_steel, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON ) ;
    Matrix m0_steel_effective = Tensor::cauchyGreen ( E_steel_effective,nu_steel, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN, YOUNG_POISSON ) ;


    Sample sample ( nullptr, sampleLength*.5, sampleHeight+2.*plateHeight, halfSampleOffset, 0 ) ;
    Sample samplebulk ( nullptr, sampleLength*.5, sampleHeight+2.*plateHeight, halfSampleOffset, 0 ) ;
    Sample samplestirrupbulk ( nullptr, sampleLength*.5, sampleHeight+2.*plateHeight, halfSampleOffset, 0 ) ;

    Sample topsupport ( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
    topsupport.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample topsupportbulk ( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
    topsupportbulk.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample topsupportstirrupbulk ( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
    topsupportstirrupbulk.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample toprightvoid ( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth ) *.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
    toprightvoid.setBehaviour ( new VoidForm() ) ;

    Sample toprightvoidbulk ( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth ) *.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
    toprightvoidbulk.setBehaviour ( new VoidForm() ) ;
    
    Sample toprightvoidstirrupbulk ( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth ) *.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
    toprightvoidstirrupbulk.setBehaviour ( new VoidForm() ) ;

    Sample baseright ( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
    baseright.setBehaviour ( new Stiffness ( m0_steel ) ) ;
// 	baseright.isVirtualFeature = true ;

    Sample baserightbulk ( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
    baserightbulk.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample baserightstirrupbulk ( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
    baserightstirrupbulk.setBehaviour ( new Stiffness ( m0_steel ) ) ;

    Sample bottomcentervoid ( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 ) *.5, -sampleHeight*.5 - plateHeight*.5 ) ;
    bottomcentervoid.setBehaviour ( new VoidForm() ) ;
    bottomcentervoid.isVirtualFeature = true ;

    Sample rightbottomvoid ( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 ) *.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
    rightbottomvoid.setBehaviour ( new VoidForm() ) ;
    rightbottomvoid.isVirtualFeature = true ;

    Sample bottomcentervoidbulk ( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 ) *.5, -sampleHeight*.5 - plateHeight*.5 ) ;
    bottomcentervoidbulk.setBehaviour ( new VoidForm() ) ;


    Sample rightbottomvoidbulk ( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 ) *.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
    rightbottomvoidbulk.setBehaviour ( new VoidForm() ) ;


    double rebarcenter = ( sampleLength*.5 - rebarEndCover ) *.5 ;
    double rebarlength = ( sampleLength - rebarEndCover*2. ) *.5 ;
    Sample rebar0 ( &sample, rebarlength, rebarDiametre , rebarcenter,  -sampleHeight*.5 + 0.064 ) ;
    rebar0.setBehaviour ( new StiffnessAndFracture ( m0_steel_effective*softeningFactor, new VonMises ( 490e6 ) ) );
// 	rebar0.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

    Sample rebar1 ( &sample, rebarlength, rebarDiametre, rebarcenter,  -sampleHeight*.5 + 0.064 + 0.085 ) ;
    rebar1.setBehaviour ( new StiffnessAndFracture ( m0_steel_effective*softeningFactor, new VonMises ( 490e6 ) ) );
// 	rebar1.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

    Sample rebar2 ( &sample, rebarlength, rebarDiametre, rebarcenter,  sampleHeight*.5 - 0.064 ) ;
    rebar2.setBehaviour ( new StiffnessAndFracture ( m0_steel_effective*softeningFactor, new VonMises ( 490e6 ) ) );
// 	rebar2.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

    Sample rebar3 ( &sample, rebarlength, rebarDiametre, rebarcenter,  sampleHeight*.5 - 0.064 - 0.085 ) ;
    rebar3.setBehaviour ( new StiffnessAndFracture ( m0_steel_effective*softeningFactor, new VonMises ( 490e6 ) ) );
// 	rebar3.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );


    std::vector<Sample*> stirrups ;

    for ( size_t i = 0 ;  i < 7 ; i++ )
    {
        stirrups.push_back ( new Sample ( 0.0084261498, sampleHeight - 2.* ( 0.064 ), 0.175 + i*0.35, 0. ) );
        stirrups.back()->setBehaviour ( new StiffnessAndFracture ( m0_steel, new VonMises ( 490e6 ) ) );
        stirrups.back()->getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius ( 0.01 );
    }
// 	for ( size_t i = 0 ;  i < 7 ; i++ )
// 	{
// 		stirrups.push_back( new Sample( 0.0084261498, sampleHeight - 2.*( 0.064 ), -0.175 - i*0.35, 0. ) );
// 		stirrups.back()->setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
// 		stirrups.back()->getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( mradius );
// 	}

    FeatureTree F ( &samplebulk, .4-3.*rebarDiametre-haveStirrups*2.*0.0084261498) ;
// 	F.addFeature(&box, &samplebulk);
    featureTree = &F ;


// 	sample.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplebulk.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplestirrupbulk.setBehaviour( new Stiffness( m0_paste ) ) ;


    samplebulk.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
    dynamic_cast<ConcreteBehaviour *> ( samplebulk.getBehaviour() )->variability = 0.00 ;
    dynamic_cast<ConcreteBehaviour *> ( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar0.getCenter().getY(),rebarDiametre ) );
    dynamic_cast<ConcreteBehaviour *> ( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar1.getCenter().getY(),rebarDiametre ) );
// 	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().getY(),rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().getY(),rebarDiametre));
    samplebulk.getBehaviour()->setSource ( sample.getPrimitive() );

    sample.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
    sample.isVirtualFeature = true ;
    dynamic_cast<ConcreteBehaviour *> ( sample.getBehaviour() )->variability = 0.00 ;
    dynamic_cast<ConcreteBehaviour *> ( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar0.getCenter().getY(),rebarDiametre ) );
  	dynamic_cast<ConcreteBehaviour *> ( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar1.getCenter().getY(),rebarDiametre ) );
// 	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().getY(),rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().getY(),rebarDiametre));


    samplestirrupbulk.setBehaviour ( new ConcreteBehaviour ( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
    samplestirrupbulk.isVirtualFeature = true ;
    dynamic_cast<ConcreteBehaviour *> ( samplestirrupbulk.getBehaviour() )->variability = 0.00 ;
    dynamic_cast<ConcreteBehaviour *> ( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar0.getCenter().getY(),rebarDiametre ) );
    dynamic_cast<ConcreteBehaviour *> ( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back ( std::make_pair ( rebar1.getCenter().getY(),rebarDiametre ) );
// 	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().getY(),rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().getY(),rebarDiametre));
    samplestirrupbulk.getBehaviour()->setSource ( sample.getPrimitive() );

    int stirruplayer = 1 ;
    int rebarlayer = -1 ;

    F.addFeature ( nullptr, &sample, rebarlayer, phi ) ;
    F.addFeature ( &samplebulk,&baserightbulk );
    F.addFeature ( &sample,&baseright, rebarlayer, phi ) ;
    F.addFeature ( &baseright,&bottomcentervoid, rebarlayer, phi );
    F.addFeature ( &baseright,&rightbottomvoid, rebarlayer, phi ) ;
    F.addFeature ( &samplebulk, &topsupportbulk);
    F.addFeature ( &sample, &topsupport, rebarlayer, phi);
    F.addFeature ( &topsupportbulk, &toprightvoidbulk);
    F.addFeature ( &topsupport, &toprightvoid, rebarlayer, phi);

//     Triangle fineZone ( Point ( 0.,sampleHeight*.5 ), Point ( 0.,-sampleHeight*.5 ), Point ( sampleLength*.5, -sampleHeight*.5 ) ) ;
//     F.addRefinementZone ( &fineZone );

    F.addFeature ( &baserightbulk,&bottomcentervoidbulk );
    F.addFeature ( &baserightbulk,&rightbottomvoidbulk ) ;


    if ( haveStirrups )
    {
        F.addFeature ( nullptr, &samplestirrupbulk, stirruplayer, psi ) ;
        F.addFeature( &samplestirrupbulk, &topsupportstirrupbulk, stirruplayer, psi ) ;
        F.addFeature ( &sample, stirrups[0], stirruplayer, psi ) ;
        F.addFeature ( nullptr,&baserightstirrupbulk, stirruplayer, psi );
        F.setSamplingFactor ( stirrups[0], 3 ) ;

        int nstirrups = 7 ;

        if ( sampleLength < 5 )
        {
            nstirrups = 5 ;
        }

        for ( int i = 1 ;  i < nstirrups ; i++ )
        {
            F.addFeature ( stirrups[i-1], stirrups[i], stirruplayer, psi ) ;
            F.setSamplingFactor ( stirrups[i], 3 ) ;
        }

        F.addFeature ( stirrups.back(), &rebar0, rebarlayer, phi ) ;
        F.addFeature ( stirrups.back(), &rebar1, rebarlayer, phi ) ;
        F.addFeature ( stirrups.back(), &rebar2, rebarlayer, phi ) ;
        F.addFeature ( stirrups.back(), &rebar3, rebarlayer, phi ) ;
// 		F.addFeature( &sample, &vrebar0 ) ;
// 		F.addFeature( &sample, &vrebar1 ) ;
// 		F.addFeature( &sample, &vrebar2 ) ;
// 		F.addFeature( &sample, &vrebar3 ) ;
    }
    else
    {
        F.addFeature ( &samplebulk, &rebar0, rebarlayer, phi ) ;
        F.addFeature ( &samplebulk, &rebar1, rebarlayer, phi ) ;
        F.addFeature ( &samplebulk, &rebar2, rebarlayer, phi ) ;
        F.addFeature ( &samplebulk, &rebar3, rebarlayer, phi ) ;
// 		F.addFeature( &sample, &vrebar0 ) ;
// 		F.addFeature( &sample, &vrebar1 ) ;
// 		F.addFeature( &sample, &vrebar2 ) ;
// 		F.addFeature( &sample, &vrebar3 ) ;
    }


// 	F.setSamplingFactor( &samplebulk, 3 ) ;
//     F.setSamplingFactor ( &rebar0, 4 ) ;
//     F.setSamplingFactor ( &rebar1, 4 ) ;

//     F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
//     F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;
//     F.setSamplingFactor ( &bottomcentervoid, 1./3 ) ;
//     F.setSamplingFactor ( &rightbottomvoid, 1./3 ) ;

    F.setSamplingFactor ( &rebar2, 1./20 ) ;
    F.setSamplingFactor ( &rebar3, 1./20 ) ;
    F.setSamplingFactor ( &baseright, 2 ) ;
    F.setSamplingFactor ( &topsupport, 2 ) ;
    F.setSamplingNumber ( samplingNumber ) ;

    F.setSamplingRestriction ( 0 );


// 	F.addPoint( new Point( supportLever+platewidth*.02, -sampleHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.02, -sampleHeight*.5 ) ) ;



// 	F.addPoint( new Point(platewidth, sampleHeight*.5)) ;
    F.setMaxIterationsPerStep ( 900 );
    F.thresholdScoreMet = 0.001 ;


    F.addPoint ( new Point ( supportLever,                -sampleHeight*.5-plateHeight ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.5,  -sampleHeight*.5-plateHeight ) ) ;
// 	F.addPoint( new Point( supportLever+platewidth*.5,  -sampleHeight*.5-plateHeight ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.25, -sampleHeight*.5-plateHeight ) ) ;
// 	F.addPoint( new Point( supportLever+platewidth*.25, -sampleHeight*.5-plateHeight ) ) ;
//
// 	F.addPoint( new Point( supportLever,                -sampleHeight*.5-plateHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.5,  -sampleHeight*.5-plateHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever+platewidth*.5,  -sampleHeight*.5-plateHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.25, -sampleHeight*.5-plateHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever+platewidth*.25, -sampleHeight*.5-plateHeight*.5 ) ) ;
//
// 	F.addPoint( new Point( supportLever-platewidth*.25, -sampleHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever+platewidth*.25, -sampleHeight*.5 ) ) ;
    F.addPoint ( new Point ( platewidth, sampleHeight*.5 ) ) ;
    F.addBoundaryCondition ( load ) ;
    F.addBoundaryCondition ( new BoundingBoxDefinedBoundaryCondition ( FIX_ALONG_XI, LEFT ) ) ;
    F.addBoundaryCondition ( new BoundingBoxNearestNodeDefinedBoundaryCondition ( FIX_ALONG_ETA, BOTTOM, Point ( supportLever, -sampleHeight*.5-plateHeight ) ) ) ;
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM, Point( supportLever, -sampleHeight*.5-plateHeight)  )) ;
    F.setOrder ( LINEAR ) ;

    step() ;

    return 0 ;
}
