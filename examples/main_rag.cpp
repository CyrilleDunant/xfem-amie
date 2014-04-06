
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
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../physics/orthotropicstiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/itoa.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../physics/void_form.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../utilities/writer/triangle_writer.h"


#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>



using namespace Mu ;

FeatureTree *featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
int nzones = 0.10 ;
double placed_area = 0 ;

double stress = 15e6 ;

double restraintDepth = 0.01 ;

Sample sample( nullptr, 0.07 + restraintDepth, 0.07 + restraintDepth, 0, 0 ) ;
Rectangle baseGeometry( 0.07, 0.07, 0, 0 ) ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;

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

double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

int totit = 1 ;

double factor = 200 ;
bool nothingToAdd = false ;
int count = 0 ;
double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;


void step()
{
	int nsteps = 320;
	int nstepstot = 320;
	featureTree->setMaxIterationsPerStep( 400 ) ;
// 	fastForward(4, 10) ;

	for( size_t i = 0 ; i < nsteps ; i++ )
	{
		std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;
		bool go_on = featureTree->step() ;

		if( featureTree->solverConverged() )
		{
			cracked_volume.push_back( featureTree->crackedVolume ) ;
			damaged_volume.push_back( featureTree->damagedVolume ) ;
		}

		//
		//
		triangles = featureTree->getElements2D() ;

		x.resize( featureTree->getDisplacements().size() ) ;
		x = featureTree->getDisplacements() ;
		sigma.resize( triangles.size()*triangles[0]->getBoundingPoints().size() * 3 ) ;
		epsilon.resize( triangles.size()*triangles[0]->getBoundingPoints().size() * 3 ) ;

		// 	sigma = F.strainFromDisplacements() ;
		// 	epsilon = F.stressFromDisplacements() ;
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
		sigma.resize( sigma_epsilon.first.size() ) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize( sigma_epsilon.second.size() ) ;
		epsilon = sigma_epsilon.second ;

		sigma11.resize( sigma.size() / 3 ) ;
		sigma22.resize( sigma.size() / 3 ) ;
		sigma12.resize( sigma.size() / 3 ) ;
		epsilon11.resize( sigma.size() / 3 ) ;
		epsilon22.resize( sigma.size() / 3 ) ;
		epsilon12.resize( sigma.size() / 3 ) ;
		vonMises.resize( sigma.size() / 3 ) ;
		angle.resize( sigma.size() / 3 ) ;

		std::cout << "unknowns :" << x.size() << std::endl ;

		cracked.clear() ;

		int npoints = triangles[0]->getBoundingPoints().size() ;

		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx_max = 0 ;
		double e_xx_min = 0 ;
		double e_yy_max = 0 ;
		double e_yy_min = 0 ;
		double ex_count = 0 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;
		
		double e_xx_max_count = 0 ;
		double e_xx_min_count = 0 ;
		double e_yy_max_count = 0 ;
		double e_yy_min_count = 0 ;

		for( size_t k = 0 ; k < triangles.size() ; k++ )
		{
			/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
			bool in = false ;

			for( size_t m = 0 ; m < tris__.size() ; m++ )
			{
				if( triangles[k] == tris__[m] )
				{
					in = true ;
					break ;
				}
			}

			cracked.push_back( in ) ;



			if( !in /*&& !triangles[k]->getBehaviour()->fractured()*/ && baseGeometry.in( triangles[k]->getCenter() ) )
			{

				for( size_t p = 0 ; p < triangles[k]->getBoundingPoints().size() ; p++ )
				{
					if( x[triangles[k]->getBoundingPoint( p ).id * 2] > x_max )
						x_max = x[triangles[k]->getBoundingPoint( p ).id * 2];

					if( x[triangles[k]->getBoundingPoint( p ).id * 2] < x_min )
						x_min = x[triangles[k]->getBoundingPoint( p ).id * 2];

					if( x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] > y_max )
						y_max = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1];

					if( x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] < y_min )
						y_min = x[triangles[k]->getBoundingPoint( p ).id * 2 + 1];

					if( triangles[k]->getBoundingPoint( p ).x > baseGeometry.width()*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).y) < .01 )
					{
//						if( e_xx_max < x[triangles[k]->getBoundingPoint( p ).id * 2] )
						e_xx_max += x[triangles[k]->getBoundingPoint( p ).id * 2] ;
 						e_xx_max_count++ ;
					}

					if( triangles[k]->getBoundingPoint( p ).x < -baseGeometry.width()*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).y) < .01 )
					{
//						if( e_xx_min > x[triangles[k]->getBoundingPoint( p ).id * 2] )
						e_xx_min += x[triangles[k]->getBoundingPoint( p ).id * 2] ;
 						e_xx_min_count++ ;

// 						ex_count++ ;
					}

					if( triangles[k]->getBoundingPoint( p ).y > baseGeometry.height()*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).x) < .01  )
					{
//						if( e_yy_max < x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] )
						e_yy_max += x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] ;
 						e_yy_max_count++ ;
// 						ex_count++ ;
					}

					if( triangles[k]->getBoundingPoint( p ).y < -baseGeometry.height()*.4999 && std::abs(triangles[k]->getBoundingPoint( p ).x) < .01  )
					{
// 						if( e_yy_min > x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] )
						e_yy_min += x[triangles[k]->getBoundingPoint( p ).id * 2 + 1] ;
 						e_yy_min_count++ ;

// 						ex_count++ ;
					}
				}

				area += triangles[k]->area() ;

				if( triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR )
				{
					if( !triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] > E_max )
						E_max = triangles[k]->getBehaviour()->param[0][0] ;

					if( !triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] < E_min )
						E_min = triangles[k]->getBehaviour()->param[0][0] ;
				}

				sigma11[k * npoints] = sigma[k * npoints * 3];
				sigma22[k * npoints] = sigma[k * npoints * 3 + 1];
				sigma12[k * npoints] = sigma[k * npoints * 3 + 2];
				sigma11[k * npoints + 1] = sigma[k * npoints * 3 + 3];
				sigma22[k * npoints + 1] = sigma[k * npoints * 3 + 4];
				sigma12[k * npoints + 1] = sigma[k * npoints * 3 + 5];
				sigma11[k * npoints + 2] = sigma[k * npoints * 3 + 6];
				sigma22[k * npoints + 2] = sigma[k * npoints * 3 + 7];
				sigma12[k * npoints + 2] = sigma[k * npoints * 3 + 8];

				if( npoints > 3 )
				{
					sigma11[k * npoints + 3] = sigma[k * npoints * 3 + 9];
					sigma22[k * npoints + 3] = sigma[k * npoints * 3 + 10];
					sigma12[k * npoints + 3] = sigma[k * npoints * 3 + 11];
					sigma11[k * npoints + 4] = sigma[k * npoints * 3 + 12];
					sigma22[k * npoints + 4] = sigma[k * npoints * 3 + 13];
					sigma12[k * npoints + 4] = sigma[k * npoints * 3 + 14];
					sigma11[k * npoints + 5] = sigma[k * npoints * 3 + 15];
					sigma22[k * npoints + 5] = sigma[k * npoints * 3 + 16];
					sigma12[k * npoints + 5] = sigma[k * npoints * 3 + 17];
				}

				epsilon11[k * npoints] = epsilon[k * npoints * 3];
				epsilon22[k * npoints] = epsilon[k * npoints * 3 + 1];
				epsilon12[k * npoints] = epsilon[k * npoints * 3 + 2];
				epsilon11[k * npoints + 1] = epsilon[k * npoints * 3 + 3];
				epsilon22[k * npoints + 1] = epsilon[k * npoints * 3 + 4];
				epsilon12[k * npoints + 1] = epsilon[k * npoints * 3 + 5];
				epsilon11[k * npoints + 2] = epsilon[k * npoints * 3 + 6];
				epsilon22[k * npoints + 2] = epsilon[k * npoints * 3 + 7];
				epsilon12[k * npoints + 2] = epsilon[k * npoints * 3 + 8];

				if( npoints > 3 )
				{
					epsilon11[k * npoints + 3] = epsilon[k * npoints * 3 + 9];
					epsilon22[k * npoints + 3] = epsilon[k * npoints * 3 + 10];
					epsilon12[k * npoints + 3] = epsilon[k * npoints * 3 + 11];
					epsilon11[k * npoints + 4] = epsilon[k * npoints * 3 + 12];
					epsilon22[k * npoints + 4] = epsilon[k * npoints * 3 + 13];
					epsilon12[k * npoints + 4] = epsilon[k * npoints * 3 + 14];
					epsilon11[k * npoints + 5] = epsilon[k * npoints * 3 + 15];
					epsilon22[k * npoints + 5] = epsilon[k * npoints * 3 + 16];
					epsilon12[k * npoints + 5] = epsilon[k * npoints * 3 + 17];
				}

				for( size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++ )
				{
					Vector vm0(0., 3) ;
					triangles[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, triangles[k]->getBoundingPoint(l), vm0, false) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					Vector agl(0., 1) ;
					triangles[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, triangles[k]->getBoundingPoint(l), agl, false) ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl[0] ;
				}

				double ar = triangles[k]->area() ;

				for( size_t l = 0 ; l < npoints ; l++ )
				{
					avg_e_xx += ( epsilon11[k * npoints + l] / npoints ) * ar;
					avg_e_yy += ( epsilon22[k * npoints + l] / npoints ) * ar;
					avg_e_xy += ( epsilon12[k * npoints + l] / npoints ) * ar;
					avg_s_xx += ( sigma11[k * npoints + l] / npoints ) * ar;
					avg_s_yy += ( sigma22[k * npoints + l] / npoints ) * ar;
					avg_s_xy += ( sigma12[k * npoints + l] / npoints ) * ar;
				}

				if( triangles[k]->getEnrichmentFunctions().size() == 0 )
				{
					for( size_t l = 0 ; l < npoints ; l++ )
					{
						avg_e_xx_nogel += ( epsilon11[k * npoints + l] / npoints ) * ar;
						avg_e_yy_nogel += ( epsilon22[k * npoints + l] / npoints ) * ar;
						avg_e_xy_nogel += ( epsilon12[k * npoints + l] / npoints ) * ar;
						avg_s_xx_nogel += ( sigma11[k * npoints + l] / npoints ) * ar;
						avg_s_yy_nogel += ( sigma22[k * npoints + l] / npoints ) * ar;
						avg_s_xy_nogel += ( sigma12[k * npoints + l] / npoints ) * ar;

					}

					nogel_area += ar ;
				}

			}
			else
			{
				sigma11[k * npoints] = 0 ;
				sigma22[k * npoints] = 0 ;
				sigma12[k * npoints] = 0 ;
				sigma11[k * npoints + 1] = 0 ;
				sigma22[k * npoints + 1] = 0 ;
				sigma12[k * npoints + 1] = 0 ;
				sigma11[k * npoints + 2] = 0 ;
				sigma22[k * npoints + 2] = 0 ;
				sigma12[k * npoints + 2] = 0 ;

				if( npoints > 3 )
				{
					sigma11[k * npoints + 3] = 0 ;
					sigma22[k * npoints + 3] = 0 ;
					sigma12[k * npoints + 3] = 0 ;
					sigma11[k * npoints + 4] = 0 ;
					sigma22[k * npoints + 4] = 0 ;
					sigma12[k * npoints + 4] = 0 ;
					sigma11[k * npoints + 5] = 0 ;
					sigma22[k * npoints + 5] = 0 ;
					sigma12[k * npoints + 5] = 0 ;
				}

				epsilon11[k * npoints] = 0 ;
				epsilon22[k * npoints] = 0 ;
				epsilon12[k * npoints] = 0 ;
				epsilon11[k * npoints + 1] = 0 ;
				epsilon22[k * npoints + 1] = 0 ;
				epsilon12[k * npoints + 1] = 0 ;
				epsilon11[k * npoints + 2] = 0 ;
				epsilon22[k * npoints + 2] = 0 ;
				epsilon12[k * npoints + 2] = 0 ;

				if( npoints > 3 )
				{
					epsilon11[k * npoints + 3] = 0 ;
					epsilon22[k * npoints + 3] = 0 ;
					epsilon12[k * npoints + 3] = 0 ;
					epsilon11[k * npoints + 4] = 0 ;
					epsilon22[k * npoints + 4] = 0 ;
					epsilon12[k * npoints + 4] = 0 ;
					epsilon11[k * npoints + 5] = 0 ;
					epsilon22[k * npoints + 5] = 0 ;
					epsilon12[k * npoints + 5] = 0 ;
				}

				for( size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++ )
				{
					vonMises[k * triangles[k]->getBoundingPoints().size() + l]  = 0 ;
					angle[k * triangles[k]->getBoundingPoints().size() + l]  = 0 ;
				}
			}
		}

		int tsize = 0 ;

		for( size_t j = 0 ; j < triangles.size() ; j++ )
		{
			if( triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR )
				tsize++ ;
		}

		std::string filename( "triangles" ) ;

		if( !go_on )
			filename = std::string( "intermediate-triangles" ) ;

		if( !featureTree->solverConverged() )
			filename = std::string( "failed-triangles" ) ;

		filename.append( itoa( totit++, 10 ) ) ;
		std::cout << filename << std::endl ;

		TriangleWriter writer(filename.c_str(), featureTree) ;
// 		writer.reset( featureTree ) ;
// 		writer.getField(TWFT_STRAIN_AND_STRESS ) ;
		writer.getField(TWFT_PRINCIPAL_STRESS ) ;
		writer.getField(TWFT_PRINCIPAL_STRAIN ) ;
		writer.getField( TWFT_VON_MISES ) ;
		writer.getField( TWFT_STIFFNESS ) ;
		writer.getField( TWFT_DAMAGE ) ;

		writer.write() ;
		e_xx_max /= e_xx_max_count ;
		e_xx_min /= e_xx_min_count ;
		e_yy_max /= e_yy_max_count ;
		e_yy_min /= e_yy_min_count ;

		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "max sigma11 :" << sigma11.max() << std::endl ;
		std::cout << "min sigma11 :" << sigma11.min() << std::endl ;
		std::cout << "max sigma12 :" << sigma12.max() << std::endl ;
		std::cout << "min sigma12 :" << sigma12.min() << std::endl ;
		std::cout << "max sigma22 :" << sigma22.max() << std::endl ;
		std::cout << "min sigma22 :" << sigma22.min() << std::endl ;

		std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
		std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
		std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
		std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
		std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
		std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;

		std::cout << "max von Mises :" << vonMises.max() << std::endl ;
		std::cout << "min von Mises :" << vonMises.min() << std::endl ;

		std::cout << "average sigma11 : " << avg_s_xx / area << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy / area << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy / area << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx / area << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy / area << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy / area << std::endl ;

		std::cout << "average sigma11 (no gel): " << avg_s_xx_nogel / nogel_area << std::endl ;
		std::cout << "average sigma22 (no gel): " << avg_s_yy_nogel / nogel_area << std::endl ;
		std::cout << "average sigma12 (no gel): " << avg_s_xy_nogel / nogel_area << std::endl ;
		std::cout << "average epsilon11 (no gel): " << avg_e_xx_nogel / nogel_area << std::endl ;
		std::cout << "average epsilon22 (no gel): " << avg_e_yy_nogel / nogel_area << std::endl ;
		std::cout << "average epsilon12 (no gel): " << avg_e_xy_nogel / nogel_area << std::endl ;

		std::cout << "apparent extension (x) " << (e_xx_max - e_xx_min)/sample.width() << std::endl ;
		std::cout << "apparent extension (y) " << (e_yy_max - e_yy_min)/sample.width() << std::endl ;

		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
		if( go_on )
		{
			featureTree->forceEnrichmentChange();
			double delta_r = sqrt( aggregateArea * 0.03 / ( ( double )zones.size() * M_PI ) ) / ( double )nstepstot ;
			double reactedArea = 0 ;

			Inclusion *current = nullptr ;

			if( !zones.empty() )
				current = zones[0].second ;

			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;

			for( size_t z = 0 ; z < zones.size() ; z++ )
			{

				zones[z].first->setRadius( zones[z].first->getRadius() + delta_r ) ;

				// 		zones[z].first->reset() ;
				if( zones[z].second == current )
				{
					current_area += zones[z].first->area() ;
					current_number++ ;
				}
				else
				{
					if( current_area / zones[z - 1].second->area() > 0.03 )
					{
						stopped_reaction++ ;

						for( size_t m = 0 ; m < current_number ; m++ )
						{
							reactedArea -= zones[z - 1 - m].first->area() ;
							zones[z - 1 - m].first->setRadius( zones[z].first->getRadius() - delta_r ) ;
							reactedArea += zones[z - 1 - m].first->area() ;
						}
					}

					current_area = zones[z].first->area() ;
					current_number = 1 ;
					current = zones[z].second ;
				}

				reactedArea += zones[z].first->area() ;
			}

			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in " << stopped_reaction << " aggs." << std::endl ;


			if( go_on )
			{
				expansion_reaction.push_back( std::make_pair( reactedArea / placed_area, avg_e_xx / area ) ) ;
				expansion_stress_xx.push_back( std::make_pair( ( avg_e_xx_nogel ) / ( nogel_area ), ( avg_s_xx_nogel ) / ( nogel_area ) ) ) ;
				expansion_stress_yy.push_back( std::make_pair( ( avg_e_yy_nogel ) / ( nogel_area ), ( avg_s_yy_nogel ) / ( nogel_area ) ) ) ;
				apparent_extension.push_back( std::make_pair((e_xx_max - e_xx_min)/sample.width(), (e_yy_max - e_yy_min)/sample.width() ) ) ;
			}

		}
			std::cout << "reaction" << "   "
// 			          << expansion_reaction[i].second << "   "
			          << "eps xx" << "\t"
			          << "eps yy" << "\t" 
								<< "sig xx" << "\t"
			          << "sig yy" << "\t"
			          << "   dx "  << "\t"
			          << "   dy "  << "\t"
			          << "cracks"  << "\t"
			          << "damage"  << "\t"
			          << std::endl ;
		for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
			std::cout << expansion_reaction[i].first << "\t"
// 			          << expansion_reaction[i].second << "\t"
			          << expansion_stress_xx[i].first << "\t"
			          << expansion_stress_yy[i].first << "\t" 
								<< expansion_stress_xx[i].second << "\t"
			          << expansion_stress_yy[i].second << "\t"
			          << apparent_extension[i].first  << "\t"
			          << apparent_extension[i].second  << "\t"
			          << cracked_volume[i]  << "\t"
			          << damaged_volume[i]  << "\t"
			          << std::endl ;

	}

	for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
		std::cout << expansion_reaction[i].first << "   "
		          << expansion_reaction[i].second << "   "
		          << expansion_stress_xx[i].first << "   "
		          << expansion_stress_xx[i].second << "   "
		          << expansion_stress_yy[i].first << "   "
		          << expansion_stress_yy[i].second << "   "
		          << apparent_extension[i].first  << "   "
		          << apparent_extension[i].second  << "   "
		          << cracked_volume[i]  << "   "
		          << damaged_volume[i]  << "   "
		          << std::endl ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZonesHomogeneously( int n, std::vector<Inclusion * > & incs , FeatureTree &F )
{
	double radiusFraction = 10 ;
	double radius = 0.00001 ;
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;

	std::vector<ExpansiveZone *> zonesToPlace ;

	srand(1);
	for( size_t i = 0 ; i < n ; i++ )
	{
		Point pos( ( ( double )rand() / RAND_MAX - .5 ) * ( sample.width() - radius * radiusFraction ), ( ( double )rand() / RAND_MAX - .5 ) * ( sample.height() - radius * radiusFraction ) ) ;
		bool alone  = true ;

		for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
		{
			if( squareDist( pos, zonesToPlace[j]->Circle::getCenter() ) < ( radius * radiusFraction + radius * radiusFraction ) * ( radius * radiusFraction + radius * radiusFraction ) )
			{
				alone = false ;
				break ;
			}
		}

		if( alone )
			zonesToPlace.push_back( new ExpansiveZone( nullptr, radius/15., pos.x, pos.y, gel ) ) ;
		else
			i-- ;
	}

	std::map<Inclusion *, int> zonesPerIncs ;

	for( size_t i = 0 ; i < zonesToPlace.size() ; i++ )
	{
		bool placed = false ;

		for( int j = 0 ; j < incs.size() ; j++ )
		{
			if( dist( zonesToPlace[i]->getCenter(), incs[j]->getCenter() ) < incs[j]->getRadius() - radius * radiusFraction /*&& incs[j]->getRadius() <= 0.008 && incs[j]->getRadius() > 0.004*/ && baseGeometry.in( zonesToPlace[i]->getCenter() ) )
			{
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature( incs[j], zonesToPlace[i] ) ;
				ret.push_back( std::make_pair( zonesToPlace[i], incs[j] ) ) ;
				placed = true ;
				break ;
			}
		}

		if( !placed )
			delete zonesToPlace[i] ;
	}

	int count = 0 ;

	for( auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i )
	{
		aggregateArea += i->first->area() ;
		count += i->second ;
// 		std::cout << aggregateArea << "  " << count << std::endl ;
	}

	std::cout << "initial Reacted Area = " << M_PI *radius *radius *ret.size() << " in " << ret.size() << " zones" << std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;
}

int main( int argc, char *argv[] )
{

	nzones = atof( argv[1] ) ;
	double dmax = atof( argv[2] ) ;

	Matrix m0_agg( 3, 3 ) ;
	m0_agg[0][0] = E_agg / ( 1 - nu * nu ) ;
	m0_agg[0][1] = E_agg / ( 1 - nu * nu ) * nu ;
	m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg / ( 1 - nu * nu ) * nu ;
	m0_agg[1][1] = E_agg / ( 1 - nu * nu ) ;
	m0_agg[1][2] = 0 ;
	m0_agg[2][0] = 0 ;
	m0_agg[2][1] = 0 ;
	m0_agg[2][2] = E_agg / ( 1 - nu * nu ) * ( 1. - nu ) / 2. ;

	Matrix m0_paste( 3, 3 ) ;
	m0_paste[0][0] = E_paste / ( 1 - nu * nu ) ;
	m0_paste[0][1] = E_paste / ( 1 - nu * nu ) * nu ;
	m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste / ( 1 - nu * nu ) * nu ;
	m0_paste[1][1] = E_paste / ( 1 - nu * nu ) ;
	m0_paste[1][2] = 0 ;
	m0_paste[2][0] = 0 ;
	m0_paste[2][1] = 0 ;
	m0_paste[2][2] = E_paste / ( 1 - nu * nu ) * ( 1. - nu ) / 2. ;

	Matrix m0_support( 3, 3 ) ;
	m0_support[0][0] = 1. / ( 1 - 0.2 * 0.2 ) ;
	m0_support[0][1] = 1. / ( 1 - 0.2 * 0.2 ) * 0.2 ;
	m0_support[0][2] = 0 ;
	m0_support[1][0] = 1. / ( 1 - 0.2 * 0.2 ) * 0.2 ;
	m0_support[1][1] = 1. / ( 1 - 0.2 * 0.2 ) ;
	m0_support[1][2] = 0 ;
	m0_support[2][0] = 0 ;
	m0_support[2][1] = 0 ;
	m0_support[2][2] = 1. / ( 1 - 0.2 * 0.2 ) * ( 1. - 0.2 ) / 2. ;

// 	Sample reinforcement0(nullptr, 8,.15,0,.5) ;
// 	reinforcement0.setBehaviour(new Stiffness(m0*5)) ;
//
// 	Sample reinforcement1(nullptr, 8,.15,0,-.5) ;
// 	reinforcement1.setBehaviour(new Stiffness(m0*5)) ;

	FeatureTree F( &sample ) ;
	featureTree = &F ;


	double itzSize = 0.00002;
//	int inclusionNumber = 1 ;
 	int inclusionNumber = 4096 ;
// 	std::vector<Inclusion *> inclusions = GranuloBolome(4.79263e-07, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);
//
// // 	if(inclusionNumber)
// // 		itzSize = inclusions[inclusions.size()/4]->getRadius() ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		delete inclusions[i] ;

	double masseInitiale = 1.06366e-05 * .9;
	double densite = 1.;

//	std::vector<Inclusion *> inclusions = GranuloBolome( masseInitiale, densite, BOLOME_A )( dmax, .0001, inclusionNumber, itzSize );
// 	std::vector<Inclusion *> inclusions = GranuloBolome(0.0000416, 1, BOLOME_D)(.0025, .1, inclusionNumber, itzSize);
// 	std::vector<Inclusion *> inclusions = ParticleSizeDistribution::get2DConcrete(dmax*0.5) ;//, masseInitiale, BOLOME_A, PSDEndCriteria(-1, 0.001, inclusionNumber)) ;

	
	
	std::vector<Feature *> feats  = ParticleSizeDistribution::get2DConcrete(&F, nullptr,  800, dmax*0.5, itzSize, BOLOME_A, CIRCLE, 1., M_PI, 100000, 0.8, &baseGeometry) ;
	std::vector<Inclusion *> inclusions ;
	
	for( size_t i = 0; i < feats.size() ; i++ )
		inclusions.push_back( static_cast<Inclusion *>( feats[i] ) ) ;

	Rectangle placeGeometry( 0.07, 0.07, 0, 0 ) ;
	int nAgg = 1 ;
	feats = placement( &placeGeometry, feats, &nAgg, 1, 6400 );
	double volume = 0 ;

	for( size_t i = 0 ; i < feats.size() ; i++ )
		volume += feats[i]->area() ;

	if( !feats.empty() )
		std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() - itzSize
		          << ", smallest r =" << feats.back()->getRadius() - itzSize << std::endl ;

	sample.setBehaviour( new PasteBehaviour()) ;
	

// 		sample.setBehaviour(new Stiffness(m0_paste)) ;
	Vector setExpansion(0., 3) ; setExpansion[0] = -0.0015 ; setExpansion[1] = -0.0015 ; setExpansion[2] = 0 ;
// 	sample.setBehaviour(new StiffnessWithImposedDeformation(m0_paste, setExpansion)) ;
	if( restraintDepth > 0 )
	{
		Sample *voidtop = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().x - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().y + ( sample.height() - restraintDepth )*.5 + 0.0025 ) ;
		voidtop->isVirtualFeature = true ;
		voidtop->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidtop );
		Sample *voidbottom = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().x - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().y - ( sample.height() - restraintDepth )*.5 - 0.0025 ) ;
		voidbottom->isVirtualFeature = true ;
		voidbottom->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidbottom );
		Sample *voidleft = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().x + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().y + ( sample.height() - restraintDepth )*.5 + 0.0025 ) ;
		voidleft->isVirtualFeature = true ;
		voidleft->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidleft );
		Sample *voidright = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().x + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().y - ( sample.height() - restraintDepth )*.5 - 0.0025 ) ;
		voidright->isVirtualFeature = true ;
		voidright->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidright );

		//width are 6544984695	10226538586	14726215564
		//length are 5113269293	26179938780	40906154344


		double fact = atof( argv[3] ) ;
		double fact0 = atof(argv[4]) ;


		Sample *blocktop = new Sample( nullptr, sample.width() - restraintDepth, restraintDepth * .5, sample.getCenter().x, sample.getCenter().y + ( sample.height() - restraintDepth )*.5 + restraintDepth * .25 ) ;
		blocktop->setBehaviour(/* new VoidForm()*/new OrthotropicStiffness(fact0*1e-4, fact0, fact0*1e-4*fact0/(fact0+fact0),  0.1, 0.) ) ;
		
		F.addFeature( nullptr, blocktop );
		F.setSamplingFactor(blocktop, 0.5);

		Sample *blockbottom = new Sample( nullptr, sample.width() - restraintDepth, restraintDepth * .5, sample.getCenter().x, sample.getCenter().y - ( sample.height() - restraintDepth )*.5 - restraintDepth * .25 ) ;
		blockbottom->setBehaviour( /*new VoidForm()*/new OrthotropicStiffness(fact0*1e-4, fact0, fact0*1e-4*fact0/(fact0+fact0),  0.1, 0.) ) ;
		F.addFeature( nullptr, blockbottom );
		F.setSamplingFactor(blockbottom, 0.5);

		Sample *blockleft = new Sample( nullptr, restraintDepth * .5, sample.height() - restraintDepth, sample.getCenter().x - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().y ) ;
		blockleft->setBehaviour( new OrthotropicStiffness(fact, fact*1e-4, fact*1e-4*fact/(fact+fact),  0.1, 0.)) ;
		F.addFeature( nullptr, blockleft );
		F.setSamplingFactor(blockleft, 0.5);

		Sample *blockright = new Sample( nullptr, restraintDepth * .5, sample.height() - restraintDepth, sample.getCenter().x + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().y ) ;
		blockright->setBehaviour(new  OrthotropicStiffness(fact, fact*1e-4, fact*1e-4*fact/(fact+fact),  0.1, 0.) ) ;
		F.addFeature( nullptr, blockright );
		F.setSamplingFactor(blockright, 0.5);
	}

	std::vector<Inclusion *> placedinclusions ;

	for( size_t i = 0 ; i < feats.size() ; i++ )
	{
		Point a( inclusions[i]->getCenter() + Point( 0, inclusions[i]->getRadius() ) ) ;
		Point b( inclusions[i]->getCenter() + Point( 0, -inclusions[i]->getRadius() ) ) ;
		Point c( inclusions[i]->getCenter() + Point( inclusions[i]->getRadius(), 0 ) ) ;
		Point d( inclusions[i]->getCenter() + Point( -inclusions[i]->getRadius(), 0 ) ) ;

		if( !( !baseGeometry.in( a ) && !baseGeometry.in( b ) && !baseGeometry.in( c ) && !baseGeometry.in( d ) ) )
		{
			inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
			AggregateBehaviour *stiff = new AggregateBehaviour() ;
// 			StiffnessWithImposedDeformation * stiff = new StiffnessWithImposedDeformation(m0_agg, setExpansion) ;
// 			Stiffness * stiff = new Stiffness(m0_agg) ;
			// 		stiff->variability = .5 ;
			inclusions[i]->setBehaviour( stiff ) ;
			F.addFeature( &sample, inclusions[i] ) ;
			F.setSamplingFactor(inclusions[i], 4.);
			placed_area += inclusions[i]->area() ;
			placedinclusions.push_back( inclusions[i] );
		}
	}

	inclusions.erase( inclusions.begin() + feats.size(), inclusions.end() ) ;
	std::cout << ", filling = " << placed_area / baseGeometry.area() * 100. << "%" << std::endl ;

	if( !inclusions.empty() )
	{
		std::cout << "largest inclusion with r = " << ( *inclusions.begin() )->getRadius() << std::endl ;
		std::cout << "smallest inclusion with r = " << ( *inclusions.rbegin() )->getRadius() << std::endl ;
		std::cout << "placed area = " <<  placed_area << std::endl ;
	}

	Circle cercle( .5, 0, 0 ) ;

	zones = generateExpansiveZonesHomogeneously(nzones, placedinclusions, F ) ;
	F.setSamplingNumber( 128 ) ;

	if( restraintDepth > 0 )
	{
// 		F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(FIX_ALONG_XI, new Rectangle(0.035+restraintDepth*.5, 0.07+restraintDepth*1.1, -(0.035+restraintDepth*.5)*.5, 0))) ;
// 		F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(FIX_ALONG_ETA, new Rectangle(0.07+restraintDepth*1.1,0.035+restraintDepth*.5, 0,-(0.035+restraintDepth*.5)*.5))) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , LEFT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , RIGHT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , BOTTOM ) ) ;
 		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , TOP ) ) ;
	}
	else
	{
// 		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA , TOP, -15e6)) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , BOTTOM_LEFT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , BOTTOM ) ) ;
	}
	F.setSamplingFactor(&sample, 2.);
	F.setOrder( LINEAR ) ;
//
	step() ;



	return 0 ;
}
