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
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/fracturecriteria/boundedvonmises.h"
#include "../physics/stiffness.h"
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
#include "../utilities/itoa.h"
#include "../utilities/random.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/homogenization/composite.h"

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
#include <GL/glut.h>
#include <time.h>



#define DEBUG

#define ID_QUIT 1
#define ID_ZOOM 5
#define ID_UNZOOM 6
#define ID_NEXT10 7
#define ID_NEXT100 3
#define ID_NEXT1000 4
#define ID_NEXT 2
#define ID_NEXT_TIME 0
#define ID_REFINE 8
#define ID_AMPLIFY 9
#define ID_DEAMPLIFY 10
#define ID_FRAC_CRIT 23

#define ID_DISP 11
#define ID_STRAIN_XX 12
#define ID_STRAIN_XY 13
#define ID_STRAIN_YY 14
#define ID_STRESS_XX 15
#define ID_STRESS_XY 16
#define ID_STRESS_YY 17
#define ID_STIFNESS 18
#define ID_ELEM 19
#define ID_VON_MISES 20
#define ID_ANGLE 22
#define ID_ENRICHMENT 21

#define DISPLAY_LIST_DISPLACEMENT 1
#define DISPLAY_LIST_ELEMENTS 2
#define DISPLAY_LIST_STRAIN_XX 3
#define DISPLAY_LIST_STRAIN_YY 4
#define DISPLAY_LIST_STRAIN_XY 5
#define DISPLAY_LIST_STRESS_XX 6
#define DISPLAY_LIST_STRESS_YY 7
#define DISPLAY_LIST_STRESS_XY 8
#define DISPLAY_LIST_CRACK 9
#define DISPLAY_LIST_STIFFNESS 10
#define DISPLAY_LIST_VON_MISES 11
#define DISPLAY_LIST_ANGLE 23
#define DISPLAY_LIST_ENRICHMENT 12
#define DISPLAY_LIST_STIFFNESS_DARK 24
#define DISPLAY_LIST_FRAC_CRIT 25


using namespace Mu ;




FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

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
double rebarDiametre = sqrt( 4*0.000506/M_PI ) ;
double rebarEndCover = 0.047 ;
double phi = 0. ;
double psi = 0. ;
bool haveStirrups = false ;

std::vector<DelaunayTriangle *> tris__ ;
double apriori_command = 0 ;
std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;
std::vector<std::pair<double, double> > load_displacement ;
std::vector< double > loads ;
std::vector< double > deltas ;
std::vector< double > displacements ;
std::vector< double > damages ;
Vector fracCrit( 0 ) ;

Vector b( 0 ) ;
Vector x( 0 ) ;
Vector sigma( 0 ) ;
Vector sigma11( 0 ) ;
Vector sigma22( 0 ) ;
Vector sigma12( 0 ) ;
Vector epsilon( 0 ) ;
Vector epsilon11( 0 ) ;
Vector epsilon22( 0 ) ;
Vector epsilon12( 0 ) ;
Vector vonMises( 0 ) ;
Vector angle( 0 ) ;

MultiTriangleWriter writer( "triangles_head", "triangles_layers", nullptr ) ;
MultiTriangleWriter writerc( "triangles_converged_head", "triangles_converged_layers", nullptr ) ;

Function loadFunction("0") ;
BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_STRESS_ETA, TOP, -platewidth, platewidth, -10, 10, loadFunction ) ;

//  BoundingBoxNearestNodeDefinedBoundaryCondition * load = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_ETA, TOP, Point(0., sampleHeight*.5+plateHeight)) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP,0) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, 0) ;
GeometryDefinedBoundaryCondition * selfload = new GeometryDefinedBoundaryCondition( SET_STRESS_ETA, new Rectangle( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) , -9025.2 ) ;
GeometryDefinedBoundaryCondition * shrinkagey = new GeometryDefinedBoundaryCondition( SET_STRESS_ETA, new Rectangle( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) , -2e-3 ) ;
GeometryDefinedBoundaryCondition * shrinkagex = new GeometryDefinedBoundaryCondition( SET_STRESS_XI, new Rectangle( sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5 ) , -2e-3 ) ;
size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 25 ;
MinimumAngle cri( M_PI / 6. ) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;

void computeDisplacement()
{
	x.resize( featureTree->getDisplacements().size() ) ;
	x = featureTree->getDisplacements() ;
	Circle C( .01, 0, 0.05 ) ;
	std::vector<DelaunayTriangle *> t = featureTree->getElements2D( &C ) ;
	std::vector<int> indices ;

	for ( size_t i = 0 ; i < t.size() ; i++ )
	{
		for ( size_t c = 0 ;  c < t[i]->getBoundingPoints().size() ; c++ )
		{
			if ( C.in( t[i]->getBoundingPoint( c ) ) )
				indices.push_back( t[i]->getBoundingPoint( c ).id ) ;
		}
	}

	std::sort( indices.begin(), indices.end() ) ;

	auto e = std::unique( indices.begin(), indices.end() ) ;
	displacement = 0 ;

	for ( auto i = indices.begin() ; i != e ; i++ )
	{
		displacement += x[( *i )*2+1] / ( e - indices.begin() ) ;
	}


// 	Circle C(.000001, -0.0025, -0.05) ;
// 	std::vector<DelaunayTriangle *> t = featureTree->get2DMesh()->getConflictingElements(&C) ;
// 	std::vector<int> indices ;
// 	for(size_t i = 0 ; i < t.size() ; i++)
// 	{
// 		for(size_t c = 0 ;  c < t[i]->getBoundingPoints().size() ; c++ )
// 		{
// 		if(C.in(t[i]->getBoundingPoint(c)))
// 			indices.push_back(t[i]->getBoundingPoint(c).id) ;
// 		}
// 	}
//
// 	std::sort(indices.begin(), indices.end()) ;
// 	std::vector<int>::iterator e = std::unique(indices.begin(), indices.end()) ;
// 	displacement = 0 ;
// 	for(std::vector<int>::iterator i = indices.begin() ; i != e ; i++)
// 	{
// 		displacement+=x[(*i)*2.]/(e-indices.begin()) ;
// 	}
//
// 	Circle C0(.000001, 0.0025, -0.05) ;
// 	std::vector<DelaunayTriangle *> t0 = featureTree->get2DMesh()->getConflictingElements(&C) ;
// 	std::vector<int> indices0 ;
// 	for(size_t i = 0 ; i < t0.size() ; i++)
// 	{
// 		for(size_t c = 0 ;  c < t0[i]->getBoundingPoints().size() ; c++ )
// 		{
// 		if(C0.in(t0[i]->getBoundingPoint(c)))
// 			indices0.push_back(t0[i]->getBoundingPoint(c).id) ;
// 		}
// 	}
//
// 	std::sort(indices0.begin(), indices0.end()) ;
// 	std::vector<int>::iterator e0 = std::unique(indices0.begin(), indices0.end()) ;
// 	double displacement0 = 0 ;
// 	for(std::vector<int>::iterator i = indices0.begin() ; i != e0 ; i++)
// 	{
// 		displacement0+=x[(*i)*2.]/(e0-indices0.begin()) ;
// 	}
//
// displacement = displacement-displacement0 ;
}

void step()
{
	
	size_t nsteps = 3000 ; //16*10;
	size_t nit = 2 ;
	size_t tries = 0 ;
	int totit = 0 ;

	for ( size_t v = 0 ; v < nsteps ; v++ )
	{
		y_max = 0 ;
		x_max = 0 ;
		y_min = 0 ;
		x_min = 0 ;
		double appliedForce = 2.*.4*load->getData()*platewidth / 1000. ;
		
		bool go_on = true ;

		go_on = featureTree->step() ;

		if ( go_on  && v > 0)
		{
			Function x("x") ;
			Function f = (x)/(platewidth) ;
			Function df = 3.*f*f-2.*f*f*f ;
			double l_real = 2e5*tries ;
			loadFunction = f_negativity(x-platewidth)*(-1e3*2-l_real) ; //5e4*52
			load->setData( loadFunction ) ;
			tries++ ;
		}
		else if (v == 0)
		{
			Function x("x") ;
			Function f = (x)/(platewidth) ;
			Function df = 3.*f*f-2.*f*f*f ;
			Function loadFunction = f_negativity(x-platewidth)*(-1e3*2) ;
			load->setData( loadFunction ) ;
			tries++ ;
		}

		triangles = featureTree->getActiveElements2D() ;
		x.resize( featureTree->getDisplacements(-1, false).size() ) ;
		x = featureTree->getDisplacements(-1, false) ;



		sigma.resize( triangles.size()*triangles[0]->getBoundingPoints().size()*3 ) ;
		epsilon.resize( triangles.size()*triangles[0]->getBoundingPoints().size()*3 ) ;
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrainInAllLayers(false) ;
		sigma.resize( sigma_epsilon.first.size() ) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize( sigma_epsilon.second.size() ) ;
		epsilon = sigma_epsilon.second ;
		sigma11.resize( sigma.size() / 3 ) ;
		sigma22.resize( sigma.size() / 3 ) ;
		sigma12.resize( sigma.size() / 3 ) ;
		epsilon11.resize( sigma.size() / 3 ) ;
		fracCrit.resize( sigma.size() / 3 ) ;
		epsilon22.resize( sigma.size() / 3 ) ;
		epsilon12.resize( sigma.size() / 3 ) ;
		vonMises.resize( sigma.size() / 3 ) ;
		angle.resize( sigma.size() / 3 ) ;

		Vector forces( featureTree->getAssembly()->getForces() ) ;

		std::cerr << "unknowns :" << x.size() << std::endl ;


		int npoints = triangles[0]->getBoundingPoints().size() ;

		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx_0 = 0 ;
		double e_xx_1 = 1 ;
		double ex_count_0 = 1 ;
		double ex_count_1 = 1 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;
		double forceCheck = 0 ;

		int tsize = 0 ;

		double deltacount = 0 ;

		double delta = 0 ;

		std::set<Point *> used ;

		for ( size_t k = 0 ; k < triangles.size() ; k++ )
		{
			bool in = false ;
// 			if(v == 8)
// 			{
// 				if(triangles[k]-> getCenter().x < 0.1)
// 				{
// 					Vector stre = triangles[k]->getState().getStress(triangles[k]-> getCenter()) ;
// 					std::cout << triangles[k]-> getCenter().y << "  "<< stre[0] << "  " << stre[1] << "  " << stre[2] << std::endl ;
// 				}
// 				if(k == triangles.size()-1)
// 					exit(0) ;
// 			}

			for ( size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++ )
			{
				if ( triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR )
				{
					if ( x[triangles[k]->getBoundingPoint( p ).id*2] > x_max )
						x_max = x[triangles[k]->getBoundingPoint( p ).id * 2] ;

					if ( x[triangles[k]->getBoundingPoint( p ).id*2] < x_min )
						x_min = x[triangles[k]->getBoundingPoint( p ).id * 2] ;

					if ( x[triangles[k]->getBoundingPoint( p ).id*2 + 1] > y_max )
						y_max = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] ;

					if ( x[triangles[k]->getBoundingPoint( p ).id*2 + 1] < y_min )
						y_min = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] ;

					if ( !triangles[k]->getBehaviour()->fractured() )
					{
						if ( used.find( &triangles[k]->getBoundingPoint( p ) ) == used.end() && triangles[k]->getBoundingPoint( p ).x <= .15 && triangles[k]->getBoundingPoint( p ).y > sampleHeight*.4999 )
						{
							used.insert( &triangles[k]->getBoundingPoint( p ) ) ;
							forceCheck += forces[triangles[k]->getBoundingPoint( p ).id * 2 + 1] ;
						}

						if ( true || triangles[k]->getBoundingPoint( p ).y > sampleHeight*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).x) < 0.0001 )
						{
							e_xx_0 = std::min(e_xx_0,x[triangles[k]->getBoundingPoint( p ).id * 2 + 1]) ;
						}

						if ( triangles[k]->getBoundingPoint( p ).y >= sampleHeight*.5 )
						{
							e_xx_1 = 0 ;
// 							ex_count_1++ ;
						}
					}

					if ( dist( Point( supportLever, -sampleHeight*.5 + 0.064 ), triangles[k]->getBoundingPoint( p ) ) < .1 )
					{
						deltacount++ ;
						delta += x[triangles[k]->getBoundingPoint( p ).id * 2] ;
					}
				}
			}

			area += triangles[k]->area() ;

			if ( triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				tsize++ ;

				if ( triangles[k]->getBehaviour()->param[0][0] > E_max )
					E_max = triangles[k]->getBehaviour()->param[0][0] ;

				if ( triangles[k]->getBehaviour()->param[0][0] < E_min )
					E_min = triangles[k]->getBehaviour()->param[0][0] ;
			}

			sigma11[k*npoints] = sigma[k*npoints*3];

			sigma22[k*npoints] = sigma[k*npoints*3+1];
			sigma12[k*npoints] = sigma[k*npoints*3+2];
			sigma11[k*npoints+1] = sigma[k*npoints*3+3];
			sigma22[k*npoints+1] = sigma[k*npoints*3+4];
			sigma12[k*npoints+1] = sigma[k*npoints*3+5];
			sigma11[k*npoints+2] = sigma[k*npoints*3+6];
			sigma22[k*npoints+2] = sigma[k*npoints*3+7];
			sigma12[k*npoints+2] = sigma[k*npoints*3+8];
			if ( npoints > 3 )
			{
				sigma11[k*npoints+3] = sigma[k*npoints*3+9];
				sigma22[k*npoints+3] = sigma[k*npoints*3+10];
				sigma12[k*npoints+3] = sigma[k*npoints*3+11];
				sigma11[k*npoints+4] = sigma[k*npoints*3+12];
				sigma22[k*npoints+4] = sigma[k*npoints*3+13];
				sigma12[k*npoints+4] = sigma[k*npoints*3+14];
				sigma11[k*npoints+5] = sigma[k*npoints*3+15];
				sigma22[k*npoints+5] = sigma[k*npoints*3+16];
				sigma12[k*npoints+5] = sigma[k*npoints*3+17];
			}

			epsilon11[k*npoints] = epsilon[k*npoints*3];

			epsilon22[k*npoints] = epsilon[k*npoints*3+1];
			epsilon12[k*npoints] = epsilon[k*npoints*3+2];
			epsilon11[k*npoints+1] = epsilon[k*npoints*3+3];
			epsilon22[k*npoints+1] = epsilon[k*npoints*3+4];
			epsilon12[k*npoints+1] = epsilon[k*npoints*3+5];
			epsilon11[k*npoints+2] = epsilon[k*npoints*3+6];
			epsilon22[k*npoints+2] = epsilon[k*npoints*3+7];
			epsilon12[k*npoints+2] = epsilon[k*npoints*3+8];

			if ( npoints > 3 )
			{
				epsilon11[k*npoints+3] = epsilon[k*npoints*3+9];
				epsilon22[k*npoints+3] = epsilon[k*npoints*3+10];
				epsilon12[k*npoints+3] = epsilon[k*npoints*3+11];
				epsilon11[k*npoints+4] = epsilon[k*npoints*3+12];
				epsilon22[k*npoints+4] = epsilon[k*npoints*3+13];
				epsilon12[k*npoints+4] = epsilon[k*npoints*3+14];
				epsilon11[k*npoints+5] = epsilon[k*npoints*3+15];
				epsilon22[k*npoints+5] = epsilon[k*npoints*3+16];
				epsilon12[k*npoints+5] = epsilon[k*npoints*3+17];
			}

// 			std::cout << "VM start" << std::endl ;
			for ( size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++ )
			{
				Vector vm0(0., 2) ;
				triangles[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, triangles[k]->getBoundingPoint(l), vm0, false) ;
				vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
// 				if(vonMises[k*triangles[k]->getBoundingPoints().size()+l] > 1000)
// 				{
// 					triangles[k]->print() ;
// 					std::cout << vonMises[k*triangles[k]->getBoundingPoints().size()+l] << std::endl ;
// 					std::cout << vm0[0] << "  " << vm0[1] << std::endl ;
// 					std::cout << x[triangles[k]->getBoundingPoint(0).id*2] << ", " << x[triangles[k]->getBoundingPoint(0).id*2+1] << "  "<< x[triangles[k]->getBoundingPoint(1).id*2] << ", " << x[triangles[k]->getBoundingPoint(1).id*2+1] << "  "<<  x[triangles[k]->getBoundingPoint(2).id*2] << ", " << x[triangles[k]->getBoundingPoint(2).id*2+1] << std::endl ;
// 					exit(0) ;
// 				}
				
				Vector agl(0., 1) ;
				triangles[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, triangles[k]->getBoundingPoint(l), agl, false) ;
				angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl[0] ;

				if ( triangles[k]->getBehaviour()->getFractureCriterion() )
				{
					fracCrit[k*triangles[k]->getBoundingPoints().size() + l] = triangles[k]->getBehaviour()->getFractureCriterion()->grade( triangles[k]->getState() ) ;
				}
			}
// 			std::cout << "VM end" << std::endl ;
			double ar = triangles[k]->area() ;

// 			if(!haveStirrups)
// 			{
// 				if(k < triangles.size()/2)
// 						ar *= 1.-phi ;
// 				else
// 					ar *= phi ;
// 			}
// 			else
// 			{
// 				if(k < triangles.size()/3)
// 						ar *= 1.-phi-psi ;
// 				else if( k < 2*triangles.size()/3)
// 					ar *= phi ;
// 				else
// 					ar *= psi ;
// 			}
				
			Vector avgsig(3) ;
			Vector avgeps(3) ;
			triangles[k]->getState().getAverageField(REAL_STRESS_FIELD, avgsig);
			triangles[k]->getState().getAverageField(STRAIN_FIELD, avgeps);

			avg_e_xx += avgeps[0] * ar;
			avg_e_yy += avgeps[1] * ar;
			avg_e_xy += avgeps[2] * ar;
			avg_s_xx += avgsig[0] * ar;
			avg_s_yy += avgsig[1] * ar;
			avg_s_xy += avgsig[2] * ar;

// 			if ( triangles[k]->getEnrichmentFunctions().size() == 0 )
// 			{
// 				for ( int l = 0 ; l < npoints ;l++ )
// 				{
// 					avg_e_xx_nogel += ( epsilon11[k*npoints+l] / npoints ) * ar;
// 					avg_e_yy_nogel += ( epsilon22[k*npoints+l] / npoints ) * ar;
// 					avg_e_xy_nogel += ( epsilon12[k*npoints+l] / npoints ) * ar;
// 					avg_s_xx_nogel += ( sigma11[k*npoints+l] / npoints ) * ar;
// 					avg_s_yy_nogel += ( sigma22[k*npoints+l] / npoints ) * ar;
// 					avg_s_xy_nogel += ( sigma12[k*npoints+l] / npoints ) * ar;
// 
// 				}
// 
// 				nogel_area += ar ;
// 			}
		}

		if ( go_on )
		{
			displacements.push_back( 1000.*e_xx_0 );
			loads.push_back( avg_s_yy/1000. );
			deltas.push_back( delta/deltacount );
			damages.push_back( featureTree->averageDamage );
		}

		if ( v % 5 == 0 )
		{

			std::cout << std::endl ;
			std::cout << "load :" << avg_s_yy/1000. << std::endl ;
			std::cout << "load check :" << .4*forceCheck / 1000 << std::endl ;
			std::cout << "delta :" << delta*1000./deltacount << std::endl ;
			std::cout << "displacement :" << 1000.*e_xx_0  << std::endl ;
			std::cout << "max value :" << x_max << std::endl ;
			std::cout << "min value :" << x_min << std::endl ;
			std::cout << "max sigma11 :" << sigma11.max() / 1000000. << std::endl ;
			std::cout << "min sigma11 :" << sigma11.min() / 1000000. << std::endl ;
			std::cout << "max sigma12 :" << sigma12.max() / 1000000. << std::endl ;
			std::cout << "min sigma12 :" << sigma12.min() / 1000000. << std::endl ;
			std::cout << "max sigma22 :" << sigma22.max() / 1000000. << std::endl ;
			std::cout << "min sigma22 :" << sigma22.min() / 1000000. << std::endl ;

			std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
			std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
			std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
			std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
			std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
			std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;

			std::cout << "max von Mises :" << vonMises.max() / 1000000. << std::endl ;
			std::cout << "min von Mises :" << vonMises.min() / 1000000. << std::endl ;

			std::cout << "average sigma11 : " << ( avg_s_xx / area ) / 1000000. << std::endl ;
			std::cout << "average sigma22 : " << ( avg_s_yy / area ) / 1000000. << std::endl ;
			std::cout << "average sigma12 : " << ( avg_s_xy / area ) / 1000000. << std::endl ;
			std::cout << "average epsilon11 : " << avg_e_xx / area << std::endl ;
			std::cout << "average epsilon22 : " << avg_e_yy / area << std::endl ;
			std::cout << "average epsilon12 : " << avg_e_xy / area << std::endl ;

		}
		else
		{
			std::cout << " ( " << avg_s_yy/1000. << "  ,  " <<  1000.*e_xx_0  << " ) "<< std::endl ;
		}


		if ( go_on )
			std::cout << appliedForce << "  " << displacements.back() << "  " << damages.back() << std::endl ;

		std::fstream ldfile( "ldn", std::ios::out )  ;

		for ( int j = 0 ; j < loads.size() ; j++ )
		{
			ldfile << displacements[j] << "   " << loads[j] << "   " << damages[j] << "   " << deltas[j] << "\n" ;
			
		}
		ldfile <<  1000.*e_xx_0 << "   " << avg_s_yy/1000. << "   " << featureTree->averageDamage << "   " << delta/deltacount << "\n" ;
		ldfile.close();
		if ( true )
		{
			writer.reset( featureTree ) ;
			writer.getField( TWFT_PRINCIPAL_STRESS ) ;
			writer.getField( TWFT_PRINCIPAL_STRAIN ) ;
			writer.getField( TWFT_CRITERION ) ;
			writer.getField( TWFT_PRINCIPAL_ANGLE ) ;
			writer.getField( TWFT_STIFFNESS_X ) ;
			writer.getField( TWFT_STIFFNESS_Y ) ;
			writer.getField( TWFT_CRACKS ) ;
			writer.getField( TWFT_DAMAGE ) ;
			writer.append() ;
		}
		if(go_on)
		{
			writerc.reset( featureTree ) ;
			writerc.getField( TWFT_PRINCIPAL_STRESS ) ;
			writerc.getField( TWFT_PRINCIPAL_STRAIN ) ;
			writerc.getField( TWFT_CRITERION ) ;
			writerc.getField( TWFT_PRINCIPAL_ANGLE ) ;
			writerc.getField( TWFT_STIFFNESS_X ) ;
			writerc.getField( TWFT_STIFFNESS_Y ) ;
			writerc.getField( TWFT_CRACKS ) ;
			writerc.getField( TWFT_DAMAGE ) ;
			writerc.append() ;
		}
// 		if ( !go_on )
// 			break ;

		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);

	}
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZones( int n, std::vector<Inclusion * > & incs , FeatureTree & F )
{
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;

	for ( size_t i = 0 ; i < incs.size() ; i++ )
	{
		aggregateArea += incs[i]->area() ;

		for ( int j = 0 ; j < n ; j++ )
		{
			double radius = 0.000001 ;
			double rangle = ( 2.*rand() / RAND_MAX - 1. ) * M_PI ;
			double rradius = ( double )rand() / RAND_MAX * ( incs[i]->getRadius() - 2.*radius * 64 ) ;
			Point center = incs[i]->getCenter() + Point( rradius * cos( rangle ), rradius * sin( rangle ) ) ;

			bool alone  = true ;

			for ( size_t k = 0 ; k < ret.size() ; k++ )
			{
				if ( squareDist( center, ret[k].first->Circle::getCenter() ) < 2.*( 2.*radius )*( 2.*radius )*64*64 )
				{
					alone = false ;
					break ;
				}
			}

			if ( alone )
			{
				ExpansiveZone * z = new ExpansiveZone( incs[i], radius, center.x, center.y, gel ) ;
				ret.push_back( std::make_pair( z, incs[i] ) ) ;
				F.addFeature( incs[i], z ) ;
			}
		}
	}

	std::cout << "initial Reacted Area = " << M_PI*0.000001*0.000001*ret.size() << " in " << ret.size() << " zones" << std::endl ;

	return ret ;
}

void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
	int i;
	double f, p, q, t;

	if ( s == 0 )
	{
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}

	h /= 60.;                        // sector 0 to 5

	i = ( int )floor( h );
	f = h - i;                      // factorial part of h
	p = v * ( 1. - s );
	q = v * ( 1. - s * f );
	t = v * ( 1. - s * ( 1. - f ) );

	switch ( i )
	{

		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;

		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;

		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;

		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;

		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;

		default:                // case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}

void init( void )
{
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_SMOOTH ); // Enables Smooth Shading
	glEnable( GL_LINE_SMOOTH ) ;
	glEnable( GL_POLYGON_SMOOTH );
	glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST ) ;

// 	glPointSize(std::max(0.4*(double)width()/(double)columns, 1.));
	glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );                                 // Black Background

	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE );
}

void Menu( int selection )
{
	switch ( selection )
	{

		case ID_NEXT:
		{
			step() ;
			dlist = false ;
			break ;
		}

		case ID_NEXT_TIME:
		{
			timepos += 0.0001 ;
			break ;
		}

		case ID_DISP :
		{
			current_list = DISPLAY_LIST_DISPLACEMENT ;
			break ;
		}

		case ID_STIFNESS :
		{
			current_list = DISPLAY_LIST_STIFFNESS ;
			break ;
		}

		case ID_STRAIN_XX :
		{
			current_list = DISPLAY_LIST_STRAIN_XX ;
			break ;
		}

		case ID_STRAIN_YY :
		{
			current_list = DISPLAY_LIST_STRAIN_YY ;
			break ;
		}

		case ID_STRAIN_XY :
		{
			current_list = DISPLAY_LIST_STRAIN_XY ;
			break ;
		}

		case ID_STRESS_XX :
		{
			current_list = DISPLAY_LIST_STRESS_XX ;
			break ;
		}

		case ID_STRESS_YY :
		{
			current_list = DISPLAY_LIST_STRESS_YY ;
			break ;
		}

		case ID_STRESS_XY :
		{
			current_list = DISPLAY_LIST_STRESS_XY ;
			break ;
		}

		case ID_ELEM :
		{
			current_list = DISPLAY_LIST_ELEMENTS ;
			break ;
		}

		case ID_VON_MISES:
		{
			current_list = DISPLAY_LIST_VON_MISES ;
			break ;
		}

		case ID_ANGLE:
		{
			current_list = DISPLAY_LIST_ANGLE ;
			break ;
		}

		case ID_ENRICHMENT:
		{
			current_list = DISPLAY_LIST_ENRICHMENT ;
			break ;
		}

		case ID_FRAC_CRIT:
		{
			current_list = DISPLAY_LIST_FRAC_CRIT ;
			break ;
		}

		case ID_QUIT :
			exit( 0 ) ;

		case ID_ZOOM :
		{
			factor *= 1.5 ;
			break ;
		}

		case ID_UNZOOM :
		{
			factor /= 1.5 ;
			break ;
		}

		case ID_AMPLIFY :
		{
			x *= 10 ;
// 			sigma11 *= 1.5 ;
// 			sigma22 *= 1.5 ;
// 			sigma12 *= 1.5 ;

			for ( size_t k = 0 ; k < triangles.size() ; k++ )
			{
				/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
				bool in = false ;

				for ( size_t m = 0 ; m < tris__.size() ; m++ )
				{
					if ( triangles[k] == tris__[m] )
					{
						in = true ;
						break ;
					}
				}

				cracked.push_back( in ) ;



				if ( !in && !triangles[k]->getBehaviour()->fractured() )
				{

					for ( size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++ )
					{
						if ( x[triangles[k]->getBoundingPoint( p ).id*2] > x_max )
							x_max = x[triangles[k]->getBoundingPoint( p ).id * 2];

						if ( x[triangles[k]->getBoundingPoint( p ).id*2] < x_min )
							x_min = x[triangles[k]->getBoundingPoint( p ).id * 2];

						if ( x[triangles[k]->getBoundingPoint( p ).id*2 + 1] > y_max )
							y_max = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1];

						if ( x[triangles[k]->getBoundingPoint( p ).id*2 + 1] < y_min )
							y_min = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1];
					}
				}
			}

			dlist = false ;

			break ;
		}

		case ID_DEAMPLIFY :
		{
			x /= 1.5 ;
// 			sigma11 /= 1.5 ;
// 			sigma22 /= 1.5 ;
// 			sigma12 /= 1.5 ;
			dlist = false ;
			break ;
		}
	}
}

void reshape( int w, int h )
{
	if ( h == 0 )                                      // Prevent A Divide By Zero By
		h = 1;                                         // Making Height Equal One

	glViewport( 0, 0, ( int )( w*factor ), ( int )( h*factor ) );

	gluPerspective(( double )h / ( double )w, 1., 1.f, 45.0f );
}

void Display( void )
{

	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ) ;
	glFlush();
	glutSwapBuffers();
	glMatrixMode( GL_PROJECTION ) ;
	glLoadIdentity() ;
	glOrtho( -4.5 / factor, 4.5 / factor, -4.5 / factor, 4.5 / factor, -4.5, 4.5 );
// 	glEnable( GL_POLYGON_OFFSET_FILL );
// 	glPolygonOffset( 0.5, 0.5 );

	//std::cout << x.max() << std::endl ;
	//std::cout << x.min() << std::endl ;

// 	double x_max = std::abs(x).min() ;
// 	double y_max = std::abs(x).min() ;
//
// 	for(size_t k = 0 ; k < x.size()/2 ; k++)
// 	{
// 		if(x[k*2]*x[k*2]+x[k*2+1]*x[k*2+1] > x_max*x_max+y_max*y_max )
// 		{
// 			x_max = x[k*2] ;
// 			y_max = x[k*2+1] ;
// 		}
// 	}

	if ( !dlist )
	{


		glNewList( DISPLAY_LIST_DISPLACEMENT,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{
			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->getBoundingPoint( 0 ).id * 2];
				double vy = x[triangles[j]->getBoundingPoint( 0 ).id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - sqrt((( vx - x_min )*( vx - x_min ) + ( vy - y_min )*( vy - y_min ) ) / (( x_max - x_min )*( x_max - x_min ) + ( y_max - y_min )*( y_max - y_min ) ) )*300., 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;

				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt((( vx - x_min )*( vx - x_min ) + ( vy - y_min )*( vy - y_min ) ) / (( x_max - x_min )*( x_max - x_min ) + ( y_max - y_min )*( y_max - y_min ) ) )*300., 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;

					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double crit_max = fracCrit.max() ;
		double crit_min = fracCrit.min() ;
		glNewList( DISPLAY_LIST_FRAC_CRIT,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{
			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR  && triangles[j]->getBehaviour()->getFractureCriterion() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->getBoundingPoint( 0 ).id * 2];
				double vy = x[triangles[j]->getBoundingPoint( 0 ).id * 2 + 1];
				double s = 1. ;

				if ( fracCrit[j*triangles[j]->getBoundingPoints().size()] < 1e-12 )
					s = .2 ;

				glBegin( GL_TRIANGLE_FAN );

				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( fracCrit[j*triangles[j]->getBoundingPoints().size()] - crit_min ) / ( crit_max - crit_min ), s, 1. ) ;

				glColor3f( c1, c2, c3 ) ;

				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 + 0 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( fracCrit[j*triangles[j]->getBoundingPoints().size() + k] - crit_min ) / ( crit_max - crit_min ), s, 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;


		double sigma11_min = sigma11.min() ;
		double sigma11_max = sigma11.max() ;
		glNewList( DISPLAY_LIST_STRAIN_XX,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{
			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma11[j*triangles[j]->getBoundingPoints().size()] - sigma11_min ) / ( sigma11_max - sigma11_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma11[j*triangles[j]->getBoundingPoints().size()] - sigma11_min ) / ( sigma11_max - sigma11_min ), 1., 1. ) ;

				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma11[j*triangles[j]->getBoundingPoints().size() + k] - sigma11_min ) / ( sigma11_max - sigma11_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double vonMises_max = vonMises.max() ;
		double vonMises_min = vonMises.min() ;

		glNewList( DISPLAY_LIST_VON_MISES,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{
			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( vonMises[j*triangles[j]->getBoundingPoints().size()] - vonMises_min ) / ( vonMises_max - vonMises_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;

				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( vonMises[j*triangles[j]->getBoundingPoints().size() + k] - vonMises_min ) / ( vonMises_max - vonMises_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;


		double angle_max = angle.max() ;
		double angle_min = angle.min() ;
		glNewList( DISPLAY_LIST_ANGLE,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( angle[j*triangles[j]->getBoundingPoints().size()] - angle_min ) / ( angle_max - angle_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;

				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( angle[j*triangles[j]->getBoundingPoints().size() + k] - angle_min ) / ( angle_max - angle_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double sigma22_min = sigma22.min() ;
		double sigma22_max = sigma22.max() ;

		glNewList( DISPLAY_LIST_STRAIN_YY,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma22[j*triangles[j]->getBoundingPoints().size()] - sigma22_min ) / ( sigma22_max - sigma22_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma22[j*triangles[j]->getBoundingPoints().size() + k] - sigma22_min ) / ( sigma22_max - sigma22_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double sigma12_min = sigma12.min() ;
		double sigma12_max = sigma12.max() ;
		glNewList( DISPLAY_LIST_STRAIN_XY,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma12[j*triangles[j]->getBoundingPoints().size()] - sigma12_min ) / ( sigma12_max - sigma12_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( sigma12[j*triangles[j]->getBoundingPoints().size() + k] - sigma12_min ) / ( sigma12_max - sigma12_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		glNewList( DISPLAY_LIST_STIFFNESS,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );

				Point a = triangles[j]->inLocalCoordinates( triangles[j]->getBoundingPoint( 0 ) ) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( triangles[j]->getBehaviour()->getTensor( a )[0][0] - E_min ) / ( E_max - E_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];
					a = triangles[j]->inLocalCoordinates( triangles[j]->getBoundingPoint( k ) ) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( triangles[j]->getBehaviour()->getTensor( a )[0][0] - E_min ) / ( E_max - E_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		glNewList( DISPLAY_LIST_STIFFNESS_DARK,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );

				Point a = triangles[j]->inLocalCoordinates( triangles[j]->getBoundingPoint( 0 ) ) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( triangles[j]->getBehaviour()->getTensor( a )[0][0] - E_min ) / ( E_max - E_min ), 1., .2 ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];
					a = triangles[j]->inLocalCoordinates( triangles[j]->getBoundingPoint( k ) ) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( triangles[j]->getBehaviour()->getTensor( a )[0][0] - E_min ) / ( E_max - E_min ), 1., .2 ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double epsilon11_min = epsilon11.min() ;
		double epsilon11_max = epsilon11.max() ;
		glNewList( DISPLAY_LIST_STRESS_XX,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon11[j*triangles[j]->getBoundingPoints().size()] - epsilon11_min ) / ( epsilon11_max - epsilon11_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon11[j*triangles[j]->getBoundingPoints().size() + k] - epsilon11_min ) / ( epsilon11_max - epsilon11_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;


		double epsilon22_min = epsilon22.min() ;
		double epsilon22_max =  epsilon22.max() ;

		glNewList( DISPLAY_LIST_STRESS_YY,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon22[j*triangles[j]->getBoundingPoints().size()] - epsilon22_min ) / ( epsilon22_max - epsilon22_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon22[j*triangles[j]->getBoundingPoints().size() + k] - epsilon22_min ) / ( epsilon22_max - epsilon22_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		double epsilon12_min = epsilon12.min() ;
		double epsilon12_max =  epsilon12.max() ;

		glNewList( DISPLAY_LIST_STRESS_XY,  GL_COMPILE ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !triangles[j]->getBehaviour()->fractured() )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				double vx = x[triangles[j]->first->id * 2];
				double vy = x[triangles[j]->first->id * 2 + 1];

				glBegin( GL_TRIANGLE_FAN );
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon12[j*triangles[j]->getBoundingPoints().size()] - epsilon12_min ) / ( epsilon12_max - epsilon12_min ), 1., 1. ) ;
				glColor3f( c1, c2, c3 ) ;
				glVertex2f( double( triangles[j]->getBoundingPoint( 0 ).x + vx ) , double( triangles[j]->getBoundingPoint( 0 ).y + vy ) );

				for ( size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					vx = x[triangles[j]->getBoundingPoint( k ).id * 2];
					vy = x[triangles[j]->getBoundingPoint( k ).id * 2 + 1];

					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*( epsilon12[j*triangles[j]->getBoundingPoints().size() + k] - epsilon12_min ) / ( epsilon12_max - epsilon12_min ), 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd() ;
			}
		}

		glEndList() ;

		glNewList( DISPLAY_LIST_ENRICHMENT,  GL_COMPILE ) ;
		glBegin( GL_TRIANGLES );

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{

			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				double c1 ;
				double c2 ;
				double c3 ;

				int enrichment = triangles[j]->getEnrichmentFunctions().size() ;
				//HSVtoRGB( &c1, &c2, &c3, 180. + 180.*(sigma12[j]-sigma12.min())/(sigma12.max()-sigma12.min()), 1., 1. )


				if ( enrichment )
				{
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*enrichment / 20., 1., 1. ) ;
					glColor3f( c1, c2, c3 ) ;
					double vx = x[triangles[j]->first->id * 2];
					double vy = x[triangles[j]->first->id * 2 + 1];

					glVertex2f( double( triangles[j]->first->x + vx ) ,
					            double( triangles[j]->first->y + vy ) );

					vx = x[triangles[j]->second->id * 2];
					vy = x[triangles[j]->second->id * 2 + 1];

					glVertex2f( double( triangles[j]->second->x + vx ) ,
					            double( triangles[j]->second->y + vy ) );


					vx = x[triangles[j]->third->id * 2];
					vy = x[triangles[j]->third->id * 2 + 1];


					glVertex2f( double( triangles[j]->third->x + vx ) ,
					            double( triangles[j]->third->y + vy ) );
				}

// 				else
// 				{
// 					double vx = x[triangles[j]->first->id*2];
// 					double vy = x[triangles[j]->first->id*2+1];
// 					Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(0)) ;
//
// 					HSVtoRGB( &c1, &c2, &c3, 0,0, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min)) ;
// 					glColor3f(c1, c2, c3) ;
//
// 					glVertex2f( double(triangles[j]->first->x + vx) ,
// 					            double(triangles[j]->first->y + vy) );
//
// 					vx = x[triangles[j]->second->id*2];
// 					vy = x[triangles[j]->second->id*2+1];
//
// 					glVertex2f( double(triangles[j]->second->x + vx) ,
// 					            double(triangles[j]->second->y + vy) );
//
//
// 					vx = x[triangles[j]->third->id*2];
// 					vy = x[triangles[j]->third->id*2+1];
//
//
// 					glVertex2f( double(triangles[j]->third->x + vx) ,
// 					            double(triangles[j]->third->y + vy) );
//
//
// 				}
			}
		}

		glEnd();

		glEndList() ;


		glNewList( DISPLAY_LIST_ELEMENTS,  GL_COMPILE ) ;
		glColor3f( 1, 1, 1 ) ;

		for ( unsigned int j = 0 ; j < triangles.size() ; j++ )
		{
			if ( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				if ( triangles[j]->getBehaviour()->fractured() )
					glColor3f( 1, 0, 0 ) ;
				else
					glColor3f( 1, 1, 1 ) ;

				glBegin( GL_LINE_LOOP );

				for ( size_t k = 0 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					double vx = 0 ;//x[triangles[j]->getBoundingPoint(k).id*2];
					double vy = 0 ;//x[triangles[j]->getBoundingPoint(k).id*2+1];

					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd();
			}
			else
			{
				glColor3f( 0, 0, 1 ) ;
				glBegin( GL_LINE_LOOP );

				for ( size_t k = 0 ; k < triangles[j]->getBoundingPoints().size() ; k++ )
				{
					double vx = 0 ;//x[triangles[j]->getBoundingPoint(k).id*2];
					double vy = 0 ;//x[triangles[j]->getBoundingPoint(k).id*2+1];

					glVertex2f( double( triangles[j]->getBoundingPoint( k ).x + vx ) ,  double( triangles[j]->getBoundingPoint( k ).y + vy ) );

				}

				glEnd();
			}

			glColor3f( 1, 1, 1 ) ;
		}

		glEndList() ;


		glNewList( DISPLAY_LIST_CRACK,  GL_COMPILE ) ;
		glLineWidth( 4 ) ;
// 		for(size_t k  = 0 ; k < crack.size() ; k++)
// 		{
// 			glColor3f(1, 0, 0) ;
// // 			for(unsigned int j=0 ; j< tris__.size() ; j++ )
// // 			{
// // 				glBegin(GL_LINE_LOOP);
// // 				double vx = x[tris__[j]->first->id*2];
// // 				double vy = x[tris__[j]->first->id*2+1];
// //
// // 				glVertex2f( double(tris__[j]->first->x/*+ vx*/) ,
// // 				            double(tris__[j]->first->y/*+ vy*/) );
// //
// // 				vx = x[tris__[j]->second->id*2];
// // 				vy = x[tris__[j]->second->id*2+1];
// //
// // 				glVertex2f( double(tris__[j]->second->x/*+ vx*/) ,
// // 				            double(tris__[j]->second->y/*+ vy*/) );
// //
// // 				vx = x[tris__[j]->third->id*2];
// // 				vy = x[tris__[j]->third->id*2+1];
// //
// // 				glVertex2f( double(tris__[j]->third->x/*+ vx*/) ,
// // 				            double(tris__[j]->third->y/*+ vy*/) );
// // 				glEnd();
// // 			}
// //
// // 			glColor3f(0, 1, 1) ;
// 			glBegin(GL_LINES) ;
// 			for(size_t j=0 ; j< crack[k]->getBoundingPoints().size()-1 ; j++ )
// 			{
// 				glVertex2f( double(crack[k]->getBoundingPoint(j).x) ,
// 				            double(crack[k]->getBoundingPoint(j).y) );
// 				glVertex2f( double(crack[k]->getBoundingPoint(j+1).x) ,
// 				            double(crack[k]->getBoundingPoint(j+1).y) );
// 			}
// 			glEnd();
// 		}

// 		for(unsigned int j=0 ; j< triangles.size() ; j++ )
// 		{
// 			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
// 			{
//
//
// 				Vector t = triangles[j]->getState()->getPrincipalStresses(triangles[j]->getCenter()) ;
// 				glBegin(GL_LINE_LOOP);
//
// 				glColor3f(1, 1, 1) ;
// 				glVertex2f( triangles[j]->getCenter().x ,  triangles[j]->getCenter().y  );
// 				glColor3f(1, 1, 1) ;
// 				glVertex2f( triangles[j]->getCenter().x +5.*t[0],  triangles[j]->getCenter().y +5.*t[1] );
//
// 				glEnd();
// 			}
//
// 			glColor3f(1, 1, 1) ;
// 		}
		glLineWidth( 1 ) ;
		glEndList() ;

		dlist = true ;
		glCallList( current_list ) ;
	}
	else
	{
		//glCallList(DISPLAY_LIST_DISPLACEMENT) ;
		//glCallList(DISPLAY_LIST_STRAIN) ;
		double c1, c2, c3 = 0;
		HSVtoRGB( &c1, &c2, &c3, 180. + 0, 1., 1. ) ;
// 		glBegin(GL_LINE) ;
// 		glVertex2f(3.5 ,
// 		           3. );
// 		glVertex2f(3.5 ,
// 		           -3. );
// 		glEnd() ;

		if ( current_list != DISPLAY_LIST_ENRICHMENT )
			glCallList( current_list ) ;

		if ( current_list == DISPLAY_LIST_ENRICHMENT )
		{
			glCallList( DISPLAY_LIST_STIFFNESS_DARK ) ;
			glCallList( current_list ) ;

		}

		glCallList( DISPLAY_LIST_CRACK ) ;

// 		if(current_list == DISPLAY_LIST_ELEMENTS)
// 			glCallList(DISPLAY_LIST_CRACK) ;

		glColor3f( 1, 1, 1 ) ;


	}

	glColor3f( 1, 0, 0 ) ;

	glFlush();
	glutSwapBuffers();
}

int main( int argc, char *argv[] )
{
	
	double softeningFactor = 1 ; .85 ;
	
	sampleLength = atof( argv[3] ) ;
	sampleHeight = atof( argv[4] ) ;
	supportLever = sampleLength*.5-.250 ;

	std::cout << sampleLength << "  " << supportLever << std::endl ;

	double compressionCrit = -37.0e6*softeningFactor ;
	phi =  3.*rebarDiametre/.4  ;
	
	psi = 2.*0.0084261498/.4  ;
	std::cout << "phi = "<< phi << ", psi = " << psi << std::endl ; 
// 	double mradius = 0.1; //0.015 ;//0.055 ;//.11 ; // .015
// 	double nradius = mradius*2.5 ;

	double E_steel = 200e9 ; 
	double nu_steel = 0.2 ;
	double nu = 0.3 ;
	double E_paste = 37e9*softeningFactor ;

	double halfSampleOffset = sampleLength*.25 ;
	
	Matrix m0_paste = Material::cauchyGreen(std::make_pair(E_paste,nu), true,SPACE_TWO_DIMENSIONAL,PLANE_STRAIN) ;
	
// 	//redimensionned so that we get in shear the right moment of inertia
// 	Matrix m0_steel = Material::orthothropicCauchyGreen(E_steel, E_steel, E_steel*(1.-nu_steel)*.5*.13/(1.-nu_steel*nu_steel), nu_steel,PLANE_STRESS_FREE_G) ;
// 		
	Matrix m0_steel = Material::cauchyGreen(std::make_pair(E_steel,nu_steel), true,SPACE_TWO_DIMENSIONAL,PLANE_STRAIN) ;

	Sample box( nullptr, sampleLength*.5+ plateHeight*2, sampleHeight + plateHeight*2, halfSampleOffset, 0 ) ;
	box.setBehaviour( new VoidForm() ) ;

	Sample sample( nullptr, sampleLength*.5, sampleHeight, halfSampleOffset, 0 ) ;
	Sample samplebulk( nullptr, sampleLength*.5, sampleHeight, halfSampleOffset, 0 ) ;
	Sample samplestirrupbulk( nullptr, sampleLength*.5, sampleHeight, halfSampleOffset, 0 ) ;
	
	Sample topsupport( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
	topsupport.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;

	Sample topsupportbulk( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
	topsupportbulk.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;
	Sample topsupportstirrupbulk( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
	topsupportstirrupbulk.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;
	Sample toprightvoid( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth )*.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
	toprightvoid.setBehaviour( new VoidForm() ) ;
	
	Sample baseright( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baseright.setBehaviour( new Stiffness( m0_steel )) ;
	Sample baserightbulk( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baserightbulk.setBehaviour( new Stiffness( m0_steel )) ;
	Sample baserightstirrupbulk( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baserightstirrupbulk.setBehaviour(new Stiffness( m0_steel )) ;

	Sample bottomcentervoid( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 )*.5, -sampleHeight*.5 - plateHeight*.5 ) ;
	bottomcentervoid.setBehaviour( new VoidForm() ) ;
	Sample rightbottomvoid( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 )*.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
	rightbottomvoid.setBehaviour( new VoidForm() ) ;

	
	Sample rebar0(&sample, sampleLength*.5 - rebarEndCover, rebarDiametre, (sampleLength*.5 - rebarEndCover)*.5,  -sampleHeight*.5 + 0.064 ) ;
	rebar0.setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
	rebar0.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

	Sample rebar1(&sample, sampleLength*.5 - rebarEndCover, rebarDiametre, (sampleLength*.5 - rebarEndCover)*.5,  -sampleHeight*.5 + 0.064 + 0.085 ) ;
	rebar1.setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
	rebar1.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );
	
	Sample rebar2(&sample, sampleLength*.5 - rebarEndCover, rebarDiametre, (sampleLength*.5 - rebarEndCover)*.5,  sampleHeight*.5 - 0.064 ) ;
	rebar2.setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
	rebar2.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

	Sample rebar3(&sample, sampleLength*.5 - rebarEndCover, rebarDiametre, (sampleLength*.5 - rebarEndCover)*.5,  sampleHeight*.5 - 0.064 - 0.085 ) ;
	rebar3.setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
	rebar3.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );
	

	std::vector<Sample*> stirrups ;

	for ( size_t i = 0 ;  i < 7 ; i++ )
	{
		stirrups.push_back( new Sample( 0.0084261498, sampleHeight - 2.*( 0.064 ), 0.175 + i*0.35, 0. ) );
		stirrups.back()->setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
		stirrups.back()->getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );
	}
// 	for ( size_t i = 0 ;  i < 7 ; i++ )
// 	{
// 		stirrups.push_back( new Sample( 0.0084261498, sampleHeight - 2.*( 0.064 ), -0.175 - i*0.35, 0. ) );
// 		stirrups.back()->setBehaviour( new StiffnessAndFracture( m0_steel, new VonMises( 490e6 ) ) );
// 		stirrups.back()->getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( mradius );
// 	}

	FeatureTree F( &box ) ;
	F.addFeature(&box, &samplebulk);
	featureTree = &F ;

	
// 	sample.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplebulk.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplestirrupbulk.setBehaviour( new Stiffness( m0_paste ) ) ;

	samplebulk.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064+0.085,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064-0.085,rebarDiametre));
	
	sample.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
	sample.isVirtualFeature = true ;
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064+0.085,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064-0.085,rebarDiametre));
	
	
	samplestirrupbulk.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL ) ) ;
	samplestirrupbulk.isVirtualFeature = true ;
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(-0.6+0.064+0.085,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(0.6-0.064-0.085,rebarDiametre));


	F.addBoundaryCondition( load ) ;
//	F.addBoundaryCondition(selfload) ;
//	F.addBoundaryCondition(shrinkagex) ;
//	F.addBoundaryCondition(shrinkagey) ;
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM, Point( -supportLever, -sampleHeight*.5 - plateHeight) ) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT ) ) ;
	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM, Point( supportLever, -sampleHeight*.5-plateHeight) ) ) ;

	int stirruplayer = 1 ;
	int rebarlayer = 0 ;

	F.addFeature( nullptr, &topsupport, rebarlayer, phi ) ;
	F.addFeature( nullptr, &topsupportbulk ) ;
	F.addFeature( nullptr, &toprightvoid ) ;
	F.addFeature( nullptr, &sample, rebarlayer, phi ) ;
	F.addFeature( nullptr,&baseright, rebarlayer, phi ) ;
	F.addFeature( nullptr,&baserightbulk);

	F.addFeature( nullptr,&bottomcentervoid);
	F.addFeature( nullptr,&rightbottomvoid) ;


	if ( atoi( argv[2] ) )
	{
		haveStirrups = true ;
		F.addFeature( nullptr, &samplestirrupbulk, stirruplayer, psi ) ;
		F.addFeature( nullptr, &topsupportstirrupbulk, stirruplayer, psi ) ;
		F.addFeature( &sample, stirrups[0], stirruplayer, psi ) ;
		F.addFeature( nullptr,&baserightstirrupbulk, stirruplayer, psi);
		F.setSamplingFactor( stirrups[0], 3 ) ;

		int nstirrups = 7 ;

		if ( sampleLength < 5 )
			nstirrups = 5 ;

		for ( size_t i = 1 ;  i < nstirrups ; i++ )
		{
			F.addFeature( stirrups[i-1], stirrups[i], stirruplayer, psi ) ;
			F.setSamplingFactor( stirrups[i], 3 ) ;
		}

		F.addFeature( stirrups.back(), &rebar0, rebarlayer, phi ) ;
		F.addFeature( stirrups.back(), &rebar1, rebarlayer, phi ) ;
		F.addFeature( stirrups.back(), &rebar2, rebarlayer, phi ) ;
		F.addFeature( stirrups.back(), &rebar3, rebarlayer, phi ) ;
// 		F.addFeature( &sample, &vrebar0 ) ;
// 		F.addFeature( &sample, &vrebar1 ) ;
// 		F.addFeature( &sample, &vrebar2 ) ;
// 		F.addFeature( &sample, &vrebar3 ) ;
	}
	else
	{
		F.addFeature( &samplebulk, &rebar0, rebarlayer, phi ) ;
		F.addFeature( &samplebulk, &rebar1, rebarlayer, phi ) ;
		F.addFeature( &samplebulk, &rebar2, rebarlayer, phi ) ;
		F.addFeature( &samplebulk, &rebar3, rebarlayer, phi ) ;
// 		F.addFeature( &sample, &vrebar0 ) ;
// 		F.addFeature( &sample, &vrebar1 ) ;
// 		F.addFeature( &sample, &vrebar2 ) ;
// 		F.addFeature( &sample, &vrebar3 ) ;
	}


	F.setSamplingFactor( &samplebulk, 3 ) ;
	F.setSamplingFactor( &rebar0, 3 ) ;
	F.setSamplingFactor( &rebar1, 3 ) ;
	
	F.setSamplingFactor( &rebar2, 3 ) ;
	F.setSamplingFactor( &rebar3, 3 ) ;
	F.setSamplingNumber( atoi( argv[1] ) ) ;
	F.setOrder( LINEAR ) ;

	
// 	F.addPoint( new Point( supportLever+platewidth*.02, -sampleHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.02, -sampleHeight*.5 ) ) ;
	F.addPoint( new Point( supportLever, -sampleHeight*.5 - plateHeight ) ) ;
// 	F.addPoint( new Point(platewidth, sampleHeight*.5)) ;
	F.setMaxIterationsPerStep( 3200 );
	
	triangles = F.getElements2D() ;
	step() ;

	glutInit( &argc,argv ) ;
	glutInitDisplayMode( GLUT_RGBA ) ;
	glutInitWindowSize( 600, 600 ) ;
	glutReshapeFunc( reshape ) ;
	glutCreateWindow( "coucou !" ) ;

	int submenu = glutCreateMenu( Menu ) ;

	glutAddMenuEntry( " Displacements ", ID_DISP );
	glutAddMenuEntry( " Stress (s) xx ", ID_STRAIN_XX );
	glutAddMenuEntry( " Stress (s) yy ", ID_STRAIN_YY );
	glutAddMenuEntry( " Stress (s) xy ", ID_STRAIN_XY );
	glutAddMenuEntry( " Strain (e) xx ", ID_STRESS_XX );
	glutAddMenuEntry( " Strain (e) yy ", ID_STRESS_YY );
	glutAddMenuEntry( " Strain (e) xy ", ID_STRESS_XY );
	glutAddMenuEntry( " Elements      ", ID_ELEM );
	glutAddMenuEntry( " Stiffness     ", ID_STIFNESS );
	glutAddMenuEntry( " Von Mises     ", ID_VON_MISES );
	glutAddMenuEntry( " Princ. angle  ", ID_ANGLE );
	glutAddMenuEntry( " Frac. crit    ", ID_FRAC_CRIT );
	glutAddMenuEntry( " Enrichment    ", ID_ENRICHMENT );

	glutCreateMenu( Menu ) ;

	glutAddMenuEntry( " Step          ", ID_NEXT );
	glutAddMenuEntry( " Step time     ", ID_NEXT_TIME );
	glutAddMenuEntry( " Zoom in       ", ID_ZOOM );
	glutAddMenuEntry( " Zoom out      ", ID_UNZOOM );
	glutAddMenuEntry( " Amplify       ", ID_AMPLIFY );
	glutAddMenuEntry( " Deamplify     ", ID_DEAMPLIFY );
	glutAddSubMenu( " Display       ", submenu );
	glutAddMenuEntry( " Quit          ", ID_QUIT ) ;


	glutAttachMenu( GLUT_RIGHT_BUTTON ) ;

	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_SMOOTH );

	glutDisplayFunc( Display ) ;
	glutMainLoop() ;

// 	delete dt ;

	return 0 ;
}
