
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../features/sample3d.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone3d.h"
#include "../utilities/itoa.h"
#include "../utilities/granulo.h"
#include "../utilities/writer/voxel_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <random>


using namespace Amie ;

FeatureTree *featureTree ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;
double scale = 100 ;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double placed_area = 0 ;

double stress = 15e6 ;

double restraintDepth = 0 ; //0.01 ;

Sample3D sample( nullptr, 0.15*scale,0.15*scale,0.15*scale, 0.075*scale, 0.075*scale, 0.075*scale ) ;

bool firstRun = true ;

std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > expansion_stress_zz ;
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

int totit = 1 ;

double factor = 200 ;
MinimumAngle cri( M_PI / 6. ) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

void step()
{
	int nsteps = 10 ;
	int nstepstot = 10 ;
	featureTree->setMaxIterationsPerStep( 10000 ) ;
	
	for( int s = 0 ; s < nsteps ; s++ )
	{
		std::cout << "\r iteration " << s << "/" << nsteps << std::flush ;
		bool go_on = featureTree->step() ;
		
		if( featureTree->solverConverged() )
		{
			cracked_volume.push_back( featureTree->crackedVolume ) ;
			damaged_volume.push_back( featureTree->damagedVolume ) ;
		}
		

		Vector strain11(4*featureTree->get3DMesh()->begin().size()) ;
		Vector strain22(4*featureTree->get3DMesh()->begin().size()) ;
		Vector strain33(4*featureTree->get3DMesh()->begin().size()) ;
		Vector stress11(4*featureTree->get3DMesh()->begin().size()) ;
		Vector stress22(4*featureTree->get3DMesh()->begin().size()) ;
		Vector stress33(4*featureTree->get3DMesh()->begin().size()) ;
		{
			Vector stress(24*featureTree->get3DMesh()->begin().size()) ;
			Vector strain(24*featureTree->get3DMesh()->begin().size()) ;
			{
				std::pair<Vector,Vector> stress_strain ;
				stress_strain.first.resize(24*featureTree->get3DMesh()->begin().size()) ;
				stress_strain.second.resize(24*featureTree->get3DMesh()->begin().size()) ;
				stress_strain = featureTree->getStressAndStrain() ;		
				stress = stress_strain.first ;
				strain = stress_strain.second ;
			}
			
			std::cout << "get stress and strain..." << std::endl ;		
			for(auto i = featureTree->get3DMesh()->begin() ; i != featureTree->get3DMesh()->end() ; i++)
			{
				for(size_t j = 0 ; j < 4 ; j++)
				{
					stress11[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+0] ;
					stress22[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+1] ;
					stress33[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+2] ;
					strain11[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+0] ;
					strain22[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+1] ;
					strain33[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+2] ;
				}
			}
		}
		
		std::cout << "averaging stress and strain..." << std::endl ;		
		Vector average_stress(3) ;
		Vector average_strain(3) ;
		double total_volume = 0. ;
		for(auto i = featureTree->get3DMesh()->begin() ; i != featureTree->get3DMesh()->end() ; i++)
		{
			double volume = i->volume() ;
			total_volume += volume ;
			for(size_t j = 0 ; j < 4 ; j++)
			{
				average_stress[0] += stress11[4*i.getPosition()+j]*volume/4 ;
				average_stress[1] += stress22[4*i.getPosition()+j]*volume/4 ;
				average_stress[2] += stress33[4*i.getPosition()+j]*volume/4 ;
				average_strain[0] += strain11[4*i.getPosition()+j]*volume/4 ;
				average_strain[1] += strain22[4*i.getPosition()+j]*volume/4 ;
				average_strain[2] += strain33[4*i.getPosition()+j]*volume/4 ;
			}
		}
		
		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "max stress11 :" << stress11.max() << std::endl ;
		std::cout << "min stress11 :" << stress11.min() << std::endl ;
		std::cout << "max stress22 :" << stress22.max() << std::endl ;
		std::cout << "min stress22 :" << stress22.min() << std::endl ;
		std::cout << "max stress33 :" << stress33.max() << std::endl ;
		std::cout << "min stress33 :" << stress33.min() << std::endl ;
		
		std::cout << "max strain11 :" << strain11.max() << std::endl ;
		std::cout << "min strain11 :" << strain11.min() << std::endl ;
		std::cout << "max strain22 :" << strain22.max() << std::endl ;
		std::cout << "min strain22 :" << strain22.min() << std::endl ;
		std::cout << "max strain33 :" << strain33.max() << std::endl ;
		std::cout << "min strain33 :" << strain33.min() << std::endl ;
		
		std::cout << "average stress11 : " << average_stress[0]/total_volume << std::endl ;
		std::cout << "average stress22 : " << average_stress[1]/total_volume << std::endl ;
		std::cout << "average stress33 : " << average_stress[2]/total_volume << std::endl ;
		std::cout << "average strain11 : " << average_strain[0]/total_volume << std::endl ;
		std::cout << "average strain22 : " << average_strain[1]/total_volume << std::endl ;
		std::cout << "average strain33 : " << average_strain[2]/total_volume << std::endl ;

		VoxelWriter writer_stiffness("rag3d_out_stiffness_"+itoa(s), 100) ;
		writer_stiffness.getField(featureTree, VWFT_STIFFNESS) ;
		writer_stiffness.write();
		
		VoxelWriter writer_strain("rag3d_out_strain+itoa(i)_"+itoa(s), 100) ;
		writer_strain.getField(featureTree, VWFT_STRAIN) ;
		writer_strain.write();

		VoxelWriter writer_stress("rag3d_out_stress_"+itoa(s), 100) ;
		writer_stress.getField(featureTree, VWFT_STRESS) ;
		writer_stress.write();
		
		VoxelWriter writer_damage("rag3d_out_damage_"+itoa(s), 100) ;
		writer_damage.getField(featureTree, VWFT_DAMAGE) ;
		writer_damage.write();
				
		if( go_on )
		{
			featureTree->forceEnrichmentChange();
			double delta_r = std::pow( aggregateArea * 0.03 / ( ( double )zones.size() * 1.333333333 * M_PI ), 0.333333 ) / ( double )nstepstot ;
			double reactedArea = 0 ;
			
			Inclusion3D * current = nullptr ;
			
			if( !zones.empty() )
				current = zones[0].second ;
			
			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;
			
			for( size_t z = 0 ; z < zones.size() ; z++ )
			{
				
				zones[z].first->setRadius( zones[z].first->getRadius() + delta_r ) ;
				
				if( zones[z].second == current )
				{
					current_area += zones[z].first->volume() ;
					current_number++ ;
				}
				else
				{
					if( current_area / zones[z - 1].second->volume() > 0.03 )
					{
						stopped_reaction++ ;
						
						for( int m = 0 ; m < current_number ; m++ )
						{
							reactedArea -= zones[z - 1 - m].first->volume() ;
							zones[z - 1 - m].first->setRadius( zones[z].first->getRadius() - delta_r ) ;
							reactedArea += zones[z - 1 - m].first->volume() ;
						}
					}
					
					current_area = zones[z].first->volume() ;
					current_number = 1 ;
					current = zones[z].second ;
				}
				
				reactedArea += zones[z].first->volume() ;
			}
			
			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in " << stopped_reaction << " aggs." << std::endl ;
			
			
			if( go_on )
			{
				expansion_reaction.push_back( std::make_pair( reactedArea / placed_area, (average_strain[0]+average_strain[1]+average_strain[2]) / (3.*total_volume) ) ) ;
				expansion_stress_xx.push_back( std::make_pair( average_strain[0] / total_volume , average_stress[0] / total_volume ) ) ;
				expansion_stress_yy.push_back( std::make_pair( average_strain[1] / total_volume , average_stress[1] / total_volume ) ) ;
				expansion_stress_zz.push_back( std::make_pair( average_strain[2] / total_volume , average_stress[2] / total_volume ) ) ;
//				apparent_extension.push_back( std::make_pair( e_xx_max - e_xx_min, e_yy_max - e_yy_min ) ) ;
			}
			
		}
		std::cout << "reaction" << "   "
		// 			          << expansion_reaction[i].second << "   "
		<< "eps xx" << "\t"
		<< "eps yy" << "\t" 
		<< "eps zz" << "\t" 
		<< "sig xx" << "\t"
		<< "sig yy" << "\t"
		<< "sig zz" << "\t"
		<< "cracks"  << "\t"
		<< "damage"  << "\t"
		<< std::endl ;
		for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
			std::cout << expansion_reaction[i].first << "\t"
			<< expansion_stress_xx[i].first << "\t"
			<< expansion_stress_yy[i].first << "\t" 
			<< expansion_stress_zz[i].first << "\t" 
			<< expansion_stress_xx[i].second << "\t"
			<< expansion_stress_yy[i].second << "\t"
			<< expansion_stress_zz[i].second << "\t" 
			<< cracked_volume[i]  << "\t"
			<< damaged_volume[i]  << "\t"
			<< std::endl ;
			
	}
	for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
		std::cout << expansion_reaction[i].first << "\t"
		<< expansion_stress_xx[i].first << "\t"
		<< expansion_stress_yy[i].first << "\t" 
		<< expansion_stress_zz[i].first << "\t" 
		<< expansion_stress_xx[i].second << "\t"
		<< expansion_stress_yy[i].second << "\t"
		<< expansion_stress_zz[i].second << "\t" 
		<< cracked_volume[i]  << "\t"
		<< damaged_volume[i]  << "\t"
		<< std::endl ;
	
/*	
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
		<< std::endl ;*/
}

std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > generateExpansiveZonesHomogeneously( int n, std::vector<Inclusion3D *> incs, FeatureTree & F)
{
	double radiusFraction = 10 ;
	double radius = 0.00001*scale ;
	std::vector<std::pair<ExpansiveZone3D *, Inclusion3D *> > ret ;
	aggregateArea = 0 ;
	std::vector<ExpansiveZone3D *> zonesToPlace ;
	std::default_random_engine generator ;
	std::uniform_real_distribution< double > xpos( sample.getXSize()*0.01, sample.getXSize()*0.99 ) ;
	std::uniform_real_distribution< double > ypos( sample.getYSize()*0.01, sample.getYSize()*0.99 ) ;
	std::uniform_real_distribution< double > zpos( sample.getZSize()*0.01, sample.getZSize()*0.99 ) ;
	GelBehaviour * gel = new GelBehaviour(false, 22e9, 0.3, 0.5, SPACE_THREE_DIMENSIONAL) ;
	
	srand(1);
	for( int i = 0 ; i < n ; i++ )
	{
		Point pos( xpos(generator), ypos(generator), zpos(generator) ) ;
		pos.print() ;
		bool alone  = true ;
		
		for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
		{
			Sphere test(radius*radiusFraction, zonesToPlace[j]->getCenter() ) ;
			if(test.in(pos))
			{
				alone = false ;
				break ;
			}
		}
		
		if( alone )
			zonesToPlace.push_back( new ExpansiveZone3D( nullptr, radius, pos.getX(), pos.getY(), pos.getZ(), gel ) ) ;
/*		else
			i-- ;*/
	}
	
	std::map<Inclusion3D *, int> zonesPerIncs ;
	
	for( size_t i = 0 ; i < zonesToPlace.size() ; i++ )
	{
		bool placed = false ;
		
		for( size_t j = 0 ; j < incs.size() ; j++ )
		{
			if( squareDist3D( zonesToPlace[i]->getCenter(), incs[j]->getCenter() ) < std::pow(incs[j]->getRadius() - radius * radiusFraction, 2. ) )
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
		aggregateArea += i->first->volume() ;
		count += i->second ;
		std::cout << aggregateArea << "  " << count << std::endl ;
	}
	
	std::cout << "initial Reacted Volume = " << M_PI *radius *radius * radius*1.33333333333333333*ret.size() << " in " << ret.size() << " zones" << std::endl ;
	std::cout << "Reactive aggregate Volume = " << aggregateArea << std::endl ;
	return ret ;
}

int main( int argc, char *argv[] )
{
		int n = atof(argv[1]) ;
		int nsample = atof(argv[2]) ;
		
		FeatureTree F( &sample ) ;
		featureTree = &F ;

		std::vector<Inclusion3D *> inclusions ;
		if(n == 1)
			inclusions.push_back(new Inclusion3D(0.0623*scale,0.075*scale,0.075*scale,0.075*scale)) ;
		else
		{
			std::string file = "sphere_2024.getT()xt" ;
			std::vector<std::string> columns ;
			columns.push_back("center_x") ;
			columns.push_back("center_y") ;
			columns.push_back("center_z") ;
			columns.push_back("radius") ;
			GranuloFromFile spheres(file, columns) ;
			inclusions = spheres.getInclusion3D(n, scale) ;			
		}
		
		double volume = 0 ;		
		for( size_t i = 0 ; i < inclusions.size() ; i++ )
			volume += inclusions[i]->volume() ;
		
		sample.setBehaviour( new /*ElasticOnly*/PasteBehaviour(false, false, 12e9,0.3,2.9e6,SPACE_THREE_DIMENSIONAL) ) ;
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
//			inclusions[i]->setBehaviour(new GelBehaviour(59e9,0.3,0.1,SPACE_THREE_DIMENSIONAL)) ;
			inclusions[i]->setBehaviour(new /*ElasticOnly*/AggregateBehaviour(false, false, 59e9,0.3,5.7e6,SPACE_THREE_DIMENSIONAL)) ;
			F.addFeature(&sample, inclusions[i]) ;
			placed_area += inclusions[i]->volume() ;
		}
				
		std::cout << "placing zones" << std::endl ;
		int numberOfZones = 10*10*10 ;
		if(n==1)
			numberOfZones = 20 ;
		zones = generateExpansiveZonesHomogeneously( numberOfZones , inclusions, F ) ;
		
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , LEFT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , BOTTOM ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ZETA , BACK ) ) ;
		F.setSamplingNumber(nsample) ;
		F.setOrder( LINEAR ) ;
		step() ;
		
		return 0 ;
}
