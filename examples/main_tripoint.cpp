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

using namespace Mu ;




FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;

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
double phi = 3.*rebarDiametre/.4  ;
double psi = 2.*0.0084261498/.4  ;
bool haveStirrups = false ;

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
// BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_STRESS_ETA, TOP, -platewidth, platewidth, -10, 10, loadFunction ) ;
BoundingBoxAndRestrictionDefinedBoundaryCondition * load = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_ETA, TOP, -platewidth, platewidth, -10, 10, 0. ) ;

//  BoundingBoxNearestNodeDefinedBoundaryCondition * load = new BoundingBoxNearestNodeDefinedBoundaryCondition(SET_ALONG_ETA, TOP, Point(0., sampleHeight*.5+plateHeight)) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP,0) ;
// BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, 0) ;

Rectangle bcbox(sampleLength*.5001, sampleHeight*1.001, sampleLength*.25, sampleHeight*.5) ;
GeometryDefinedBoundaryCondition selfload( SET_VOLUMIC_STRESS_ETA, &bcbox , -9025.2 ) ;
GeometryDefinedBoundaryCondition shrinkagey( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
GeometryDefinedBoundaryCondition shrinkagex( SET_VOLUMIC_STRESS_ETA, &bcbox , -2e-3 ) ;
int count = 0 ;
double aggregateArea = 0;


void step()
{
	
	size_t nsteps = 600*4 ; //16*10;
	size_t nit = 2 ;
	size_t tries = 0 ;
	int totit = 0 ;
	double delta_d = 0.0175e-3 ;

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
			load->setData( load->getData()-delta_d ) ;
		}
// 		if ( go_on  && v > 0)
// 		{
// 			Function x("x") ;
// 			Function f = (x)/(platewidth) ;
// 			Function df = 3.*f*f-2.*f*f*f ;
// 			double l_real = 2e5*tries ;
// 			loadFunction = f_negativity(x-platewidth)*(-1e3*2-l_real) ; //5e4*52
// 			load->setData( loadFunction ) ;
// 			tries++ ;
// 		}
// 		else if (v == 0)
// 		{
// 			Function x("x") ;
// 			Function f = (x)/(platewidth) ;
// 			Function df = 3.*f*f-2.*f*f*f ;
// 			Function loadFunction = f_negativity(x-platewidth)*(-1e3*2) ;
// 			load->setData( loadFunction ) ;
// 			tries++ ;
// 		}

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

						if (triangles[k]->getBoundingPoint( p ).y > sampleHeight*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).x) < 0.0001 )
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

			if(!haveStirrups)
			{
				if(k < triangles.size()/2)
					ar *= 1.-phi ;
				else
					ar *= phi ;
			}
			else
			{
				if(k < triangles.size()/3)
						ar *= 1.-phi-psi ;
				else if( k < 2*triangles.size()/3)
					ar *= phi ;
				else
					ar *= psi ;
			}
				
			Vector avgsig(3) ;
			Vector avgeps(3) ;
			triangles[k]->getState().getAverageField(REAL_STRESS_FIELD, avgsig,0,0);
			triangles[k]->getState().getAverageField(STRAIN_FIELD, avgeps,0,0);
			
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
			displacements.push_back( 1000.*(load->getData()+delta_d));
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
			std::cout << "displacement :" << 1000.*e_xx_0  << " ; " << 1000.*(load->getData()) << std::endl ;
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
			std::cout << " ( " << avg_s_yy/1000. << "  ,  " <<  1000.*(load->getData())  << " ) "<< std::endl ;
		}

		if ( go_on )
			std::cout << avg_s_yy/1000. << "  " << displacements.back() << "  " << damages.back() << std::endl ;

		
		std::fstream ldfile( "ldn", std::ios::out )  ;
		for ( int j = 0 ; j < loads.size() ; j++ )
		{
			ldfile << displacements[j] << "   " << loads[j] << "   " << damages[j] << "   " << deltas[j] << "\n" ;
			
		}
		if(!go_on)
		  ldfile <<  1000.*(load->getData()) << "   " << avg_s_yy/1000. << "   " << featureTree->averageDamage << "   " << delta/deltacount << "\n" ;
		ldfile.close();
		
		
		if ( true )
		{
			writer.reset( featureTree ) ;
			writer.getField( TWFT_PRINCIPAL_STRESS ) ;
			writer.getField( TWFT_PRINCIPAL_STRAIN ) ;
			writer.getField( TWFT_CRITERION ) ;
			writer.getField( TWFT_STIFFNESS_X ) ;
			writer.getField( TWFT_STIFFNESS_Y ) ;
			writer.getField( TWFT_DAMAGE ) ;
			writer.append() ;
		}
		
		if ( go_on )
		{
			writerc.reset( featureTree ) ;
			writerc.getField( TWFT_DAMAGE ) ;
			writerc.append() ;
			writerc.writeSvg(0.) ;
		}
// 		if ( !go_on )
// 			break ;

		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);

	}
}


int main( int argc, char *argv[] )
{

	double softeningFactor = 1 ; .85 ;
	
	sampleLength = atof( argv[3] ) ;
	sampleHeight = atof( argv[4] ) ;
	supportLever = sampleLength*.5-.250 ;

	std::cout << sampleLength << "  " << supportLever << std::endl ;

	double compressionCrit = -34.2e6*softeningFactor ;

	std::cout << "phi = "<< phi << ", psi = " << psi << std::endl ; 
// 	double mradius = 0.1; //0.015 ;//0.055 ;//.11 ; // .015
// 	double nradius = mradius*2.5 ;

	double E_steel = 200e9 * M_PI *.25; 
	double nu_steel = 0.3 ;
	double nu = 0.3 ;
	double E_paste = 32.4e9*softeningFactor ;

	double halfSampleOffset = sampleLength*.25 ;
	
	Matrix m0_paste = Material::cauchyGreen(std::make_pair(E_paste,nu), true, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;
	
// 	//redimensionned so that we get in shear the right moment of inertia
// 	Matrix m0_steel = Material::orthothropicCauchyGreen(E_steel, E_steel, E_steel*(1.-nu_steel)*.5*.13/(1.-nu_steel*nu_steel), nu_steel,PLANE_STRESS_FREE_G) ;
// 		
	Matrix m0_steel = Material::cauchyGreen(std::make_pair(E_steel,nu_steel), true, SPACE_TWO_DIMENSIONAL, PLANE_STRAIN) ;
	

	Sample sample( nullptr, sampleLength*.5, sampleHeight+plateHeight, halfSampleOffset, -plateHeight*.5 ) ;
	Sample samplebulk( nullptr, sampleLength*.5, sampleHeight+plateHeight, halfSampleOffset, -plateHeight*.5 ) ;
	Sample samplestirrupbulk( nullptr, sampleLength*.5, sampleHeight+plateHeight, halfSampleOffset, -plateHeight*.5 ) ;
	
// 	Sample topsupport( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
// 	topsupport.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;
// 	topsupport.isVirtualFeature = true ;
// 
// 	Sample topsupportbulk( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
// 	topsupportbulk.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;
// 	
// 	Sample topsupportstirrupbulk( nullptr, platewidth, plateHeight, platewidth*.5, sampleHeight*.5 + plateHeight*.5 ) ;
// 	topsupportstirrupbulk.setBehaviour( new VoidForm()/*Stiffness( m0_steel )*/ ) ;
// 	
// 	Sample toprightvoid( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth )*.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
// 	toprightvoid.setBehaviour( new VoidForm() ) ;
// 	toprightvoid.isVirtualFeature = true ;
// 	
// 	Sample toprightvoidbulk( nullptr, sampleLength*.5 - platewidth, plateHeight, ( sampleLength*.5 - platewidth )*.5 + platewidth, sampleHeight*.5 + plateHeight*.5 ) ;
// 	toprightvoidbulk.setBehaviour( new VoidForm() ) ;
	
	Sample baseright( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baseright.setBehaviour(  new ConcreteBehaviour( E_steel, nu_steel, 1000.*compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL )/*new Stiffness( m0_steel )*/) ;
// 	baseright.isVirtualFeature = true ;
	
	Sample baserightbulk( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baserightbulk.setBehaviour(  new ConcreteBehaviour( E_steel, nu_steel, 1000.*compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL )/*new Stiffness( m0_steel )*/) ;
	
	Sample baserightstirrupbulk( platewidth, plateHeight, supportLever, -sampleHeight*.5 - plateHeight*.5 ) ;
	baserightstirrupbulk.setBehaviour( new ConcreteBehaviour( E_steel, nu_steel, 1000.*compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL )/*new Stiffness( m0_steel )*/) ;

	Sample bottomcentervoid( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 )*.5, -sampleHeight*.5 - plateHeight*.5 ) ;
	bottomcentervoid.setBehaviour( new VoidForm() ) ;
	bottomcentervoid.isVirtualFeature = true ;
	
	Sample rightbottomvoid( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 )*.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
	rightbottomvoid.setBehaviour( new VoidForm() ) ;
	rightbottomvoid.isVirtualFeature = true ;
	
	Sample bottomcentervoidbulk( supportLever - platewidth*.5, plateHeight, ( supportLever - platewidth*.5 )*.5, -sampleHeight*.5 - plateHeight*.5 ) ;
	bottomcentervoidbulk.setBehaviour( new VoidForm() ) ;
	
	
	Sample rightbottomvoidbulk( supportMidPointToEndClearance - platewidth*.5, plateHeight, sampleLength*.5 - ( supportMidPointToEndClearance - platewidth*.5 )*.5,  -sampleHeight*.5 - plateHeight*.5 ) ;
	rightbottomvoidbulk.setBehaviour( new VoidForm() ) ;
	
	
	double rebarcenter = (sampleLength*.5 - rebarEndCover)*.5 ;
	double rebarlength = (sampleLength - rebarEndCover*2.)*.5 ;
	Sample rebar0(&sample, rebarlength, rebarDiametre, rebarcenter,  -sampleHeight*.5 + 0.064 ) ;
	rebar0.setBehaviour( new Stiffness/*AndFracture*/( m0_steel/*, new VonMises( 490e6 )*/ ) );
// 	rebar0.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

	Sample rebar1(&sample, rebarlength, rebarDiametre, rebarcenter,  -sampleHeight*.5 + 0.064 + 0.085 ) ;
	rebar1.setBehaviour( new Stiffness/*AndFracture*/( m0_steel/*, new VonMises( 490e6 )*/ ) );
// 	rebar1.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );
	
	Sample rebar2(&sample, rebarlength, rebarDiametre, rebarcenter,  sampleHeight*.5 - 0.064 ) ;
	rebar2.setBehaviour( new Stiffness/*AndFracture*/( m0_steel/*, new VonMises( 490e6 )*/ ) );
// 	rebar2.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );

	Sample rebar3(&sample, rebarlength, rebarDiametre, rebarcenter,  sampleHeight*.5 - 0.064 - 0.085 ) ;
	rebar3.setBehaviour( new Stiffness/*AndFracture*/( m0_steel/*, new VonMises( 490e6 )*/ ) );
// 	rebar3.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius( 0.01 );
	

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

	FeatureTree F( &samplebulk ) ;
// 	F.addFeature(&box, &samplebulk);
	featureTree = &F ;

	
// 	sample.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplebulk.setBehaviour( new Stiffness( m0_paste ) ) ;
// 	samplestirrupbulk.setBehaviour( new Stiffness( m0_paste ) ) ;

	
	samplebulk.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL,MIRROR_Y ) ) ;
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar0.getCenter().y,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar1.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplebulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().y,rebarDiametre));
	samplebulk.getBehaviour()->setSource(sample.getPrimitive());
	
	sample.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL,MIRROR_Y ) ) ;
	sample.isVirtualFeature = true ;
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar0.getCenter().y,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar1.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( sample.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().y,rebarDiametre));
	
	
	samplestirrupbulk.setBehaviour( new ConcreteBehaviour( E_paste, nu, compressionCrit,PLANE_STRAIN, UPPER_BOUND, SPACE_TWO_DIMENSIONAL,MIRROR_Y ) ) ;
	samplestirrupbulk.isVirtualFeature = true ;
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->variability = 0.00 ;
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar0.getCenter().y,rebarDiametre));
	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar1.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar2.getCenter().y,rebarDiametre));
// 	dynamic_cast<ConcreteBehaviour *>( samplestirrupbulk.getBehaviour() )->rebarLocationsAndDiameters.push_back(std::make_pair(rebar3.getCenter().y,rebarDiametre));
	samplestirrupbulk.getBehaviour()->setSource(sample.getPrimitive());

	int stirruplayer = 1 ;
	int rebarlayer = 0 ;

	F.addFeature( nullptr, &sample, rebarlayer, phi ) ;
	F.addFeature( &samplebulk,&baserightbulk);
	F.addFeature( &sample,&baseright, rebarlayer, phi ) ;
	F.addFeature( &baseright,&bottomcentervoid, rebarlayer, phi);
	F.addFeature( &baseright,&rightbottomvoid, rebarlayer, phi) ;
	
	Triangle fineZone(Point(0.,sampleHeight*.5), Point(0.,-sampleHeight*.5), Point(sampleLength*.5, -sampleHeight*.5)) ;
	F.addRefinementZone(&fineZone);
	
// 	F.addFeature( nullptr, &topsupportbulk ) ;
// 	F.addFeature( nullptr, &toprightvoid ) ;
	

	F.addFeature( &baserightbulk,&bottomcentervoidbulk);
	F.addFeature( &baserightbulk,&rightbottomvoidbulk) ;


	if ( atoi( argv[2] ) )
	{
		haveStirrups = true ;
		F.addFeature( nullptr, &samplestirrupbulk, stirruplayer, psi ) ;
// 		F.addFeature( nullptr, &topsupportstirrupbulk, stirruplayer, psi ) ;
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


// 	F.setSamplingFactor( &samplebulk, 3 ) ;
	F.setSamplingFactor( &rebar0, 4 ) ;
	F.setSamplingFactor( &rebar1, 4 ) ;
	
	F.setSamplingFactor( &bottomcentervoid, 1./3 ) ;
	F.setSamplingFactor( &rightbottomvoid, 1./3 ) ;
	F.setSamplingFactor( &bottomcentervoid, 1./3 ) ;
	F.setSamplingFactor( &rightbottomvoid, 1./3 ) ;
	
	F.setSamplingFactor( &rebar2, 4 ) ;
	F.setSamplingFactor( &rebar3, 4 ) ;
	F.setSamplingNumber( atoi( argv[1] ) ) ;
	F.setOrder( LINEAR ) ;
	F.setSamplingRestriction(SAMPLE_NO_RESTRICTION);

	
// 	F.addPoint( new Point( supportLever+platewidth*.02, -sampleHeight*.5 ) ) ;
// 	F.addPoint( new Point( supportLever-platewidth*.02, -sampleHeight*.5 ) ) ;


	
// 	F.addPoint( new Point(platewidth, sampleHeight*.5)) ;
	F.setMaxIterationsPerStep( 3200 );
	
	
	F.addPoint( new Point( supportLever,                -sampleHeight*.5-plateHeight ) ) ;
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
	F.addPoint( new Point( platewidth, sampleHeight*.5 ) ) ;
	F.addBoundaryCondition( load ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, LEFT ) ) ;
	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM, Point( supportLever, -sampleHeight*.5-plateHeight)  )) ;
	triangles = F.getElements2D() ;
	step() ;
// 	delete dt ;

	return 0 ;
}
