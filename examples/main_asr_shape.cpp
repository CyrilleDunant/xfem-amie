
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/writer/triangle_writer.h"

#include <fstream>
#include <cmath>
#include <random>
#include <typeinfo>
#include <limits>
#include <time.h> 
#include <sys/time.h> 
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

double timepos = 0.00 ;
double percent = 0.10 ;
double placed_area = 0 ;

double stress = 15e6 ;

RectangularFeature sample(nullptr, 0.07, 0.07, 0.0, 0.0) ;
RectangularFeature placing(nullptr, 0.07, 0.07, 0.0, 0.0) ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma12(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon12(0) ; 
Vector vonMises(0) ; 
Vector angle(0) ; 

double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

int totit = 5 ;

double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;

double shape ;
double orientation ;
double spread ;

struct Zone
{
	GeometryType type ;
	
	ExpansiveZone * zone ;
	Inclusion * inc ;
	EllipsoidalInclusion * ell ;
	TriangularInclusion * tri ;
	RectangularInclusion * rec ;

	Zone() { type = CIRCLE;  zone = nullptr ; inc = nullptr ; ell = nullptr; tri = nullptr; rec = nullptr; }
	Zone(ExpansiveZone * z, Inclusion * i) { type = CIRCLE; zone = z ; inc = i ; ell = nullptr; tri = nullptr; rec = nullptr; }
	Zone(ExpansiveZone * z, EllipsoidalInclusion * i) { type = ELLIPSE; zone = z ; inc = nullptr ; ell = i; tri = nullptr; rec = nullptr; }
	Zone(ExpansiveZone * z, TriangularInclusion * i) { type = TRIANGLE; zone = z ; inc = nullptr ; ell = nullptr; tri = i; rec = nullptr; }
	Zone(ExpansiveZone * z, RectangularInclusion * i) { type = RECTANGLE; zone = z ; inc = nullptr ; ell = nullptr; tri = nullptr; rec = i; }
	Zone(const Zone & z) { type = z.type; zone = z.zone ; inc = z.inc ; ell = z.ell ; tri = z.tri ; rec = z.rec ; }
	
	Feature * feature() 
	{
		switch(type)
		{
			case CIRCLE:
				return inc ;
			case ELLIPSE:
				return ell ;
			case TRIANGLE:
				return tri ;
			case RECTANGLE:
				return rec ;
            default:
               return nullptr ; 
		}
		return nullptr ;
	}
	
	bool is(Feature * f) { return f == feature() ; }	
	
	double areaOfInclusion()
	{
		switch(type)
		{
			case CIRCLE:
				return inc->area() ;
			case ELLIPSE:
				return ell->area() ;
			case TRIANGLE:
				return tri->area() ;
			case RECTANGLE:
				return rec->area() ;
            default:
                std::cout << "unhandled geometry type" << std::endl ;
                exit(0) ;
                return 0 ;
		}
		return 0. ;
	}
	
} ;

std::vector<Zone> zones ;

void step(GeometryType ref, int samplingNumber)
{

	int nsteps = 1;
	int nstepstot = 20;
	int maxtries = 400 ;
	int tries = 0 ;
	featureTree->setMaxIterationsPerStep(40000) ;
	
	for(int i = 0 ; i < nsteps ; i++)
	{

		bool go_on = featureTree->step() ;
		
		if(featureTree->solverConverged())
		{
			cracked_volume.push_back(featureTree->crackedVolume) ;
			damaged_volume.push_back(featureTree->damagedVolume) ;
		}	
	
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
		sigma.resize(featureTree->get2DMesh()->begin().size()*featureTree->get2DMesh()->begin()->getBoundingPoints().size()*3) ;
		epsilon.resize(featureTree->get2DMesh()->begin().size()*featureTree->get2DMesh()->begin()->getBoundingPoints().size()*3) ;
		
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		
		sigma11.resize(sigma.size()/3) ;
		sigma22.resize(sigma.size()/3) ;
		sigma12.resize(sigma.size()/3) ;
		epsilon11.resize(sigma.size()/3) ;
		epsilon22.resize(sigma.size()/3) ;
		epsilon12.resize(sigma.size()/3) ;
		vonMises.resize(sigma.size()/3) ;
		angle.resize(sigma.size()/3) ;
		
		std::cout << "unknowns :" << x.size() << std::endl ;
		
		cracked.clear() ;
		
		int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;
		
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

		for(auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++)
		{
			bool in = false ;
			for(size_t m = 0 ; m < tris__.size() ; m++)
			{
				if(k == tris__[m])
				{
					in = true ;
					break ;
				}
			}
			cracked.push_back(in) ;
			
			
			
			if(!in /*&& !k->getBehaviour()->fractured()*/)
			{
				
				for(size_t p = 0 ;p < k->getBoundingPoints().size() ; p++)
				{
					if(x[k->getBoundingPoint(p).getId()*2] > x_max)
						x_max = x[k->getBoundingPoint(p).getId()*2];
					if(x[k->getBoundingPoint(p).getId()*2] < x_min)
						x_min = x[k->getBoundingPoint(p).getId()*2];
					if(x[k->getBoundingPoint(p).getId()*2+1] > y_max)
						y_max = x[k->getBoundingPoint(p).getId()*2+1];
					if(x[k->getBoundingPoint(p).getId()*2+1] < y_min)
						y_min = x[k->getBoundingPoint(p).getId()*2+1];
					if( k->getBoundingPoint( p ).getX() > sample.width()*.4999 && std::abs(k->getBoundingPoint( p ).getY()) < .01 )
					{
//						if( e_xx_max < x[k->getBoundingPoint( p ).getId() * 2] )
						e_xx_max += x[k->getBoundingPoint( p ).getId() * 2] ;
 						e_xx_max_count++ ;
					}

					if( k->getBoundingPoint( p ).getX() < -sample.width()*.4999 && std::abs(k->getBoundingPoint( p ).getY()) < .01 )
					{
//						if( e_xx_min > x[k->getBoundingPoint( p ).getId() * 2] )
						e_xx_min += x[k->getBoundingPoint( p ).getId() * 2] ;
 						e_xx_min_count++ ;

// 						ex_count++ ;
					}

					if( k->getBoundingPoint( p ).getY() > sample.height()*.4999 && std::abs(k->getBoundingPoint( p ).getX()) < .01  )
					{
//						if( e_yy_max < x[k->getBoundingPoint( p ).getId() * 2 + 1] )
						e_yy_max += x[k->getBoundingPoint( p ).getId() * 2 + 1] ;
 						e_yy_max_count++ ;
// 						ex_count++ ;
					}

					if( k->getBoundingPoint( p ).getY() < -sample.height()*.4999 && std::abs(k->getBoundingPoint( p ).getX()) < .01  )
					{
// 						if( e_yy_min > x[k->getBoundingPoint( p ).getId() * 2 + 1] )
						e_yy_min += x[k->getBoundingPoint( p ).getId() * 2 + 1] ;
 						e_yy_min_count++ ;

// 						ex_count++ ;
					}
				}
				area += k->area() ;
				if(k->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(!k->getBehaviour()->param.isNull() && k->getBehaviour()->param[0][0] > E_max)
						E_max = k->getBehaviour()->param[0][0] ;
					if(!k->getBehaviour()->param.isNull() && k->getBehaviour()->param[0][0] < E_min)
						E_min = k->getBehaviour()->param[0][0] ;
				}
					
				sigma11[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3];
				sigma22[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3+1];
				sigma12[k.getPosition()*npoints] = sigma[k.getPosition()*npoints*3+2];
				sigma11[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+3];
				sigma22[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+4];
				sigma12[k.getPosition()*npoints+1] = sigma[k.getPosition()*npoints*3+5];
				sigma11[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+6];
				sigma22[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+7];
				sigma12[k.getPosition()*npoints+2] = sigma[k.getPosition()*npoints*3+8];
				
				if(npoints >3)
				{
					sigma11[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+9];
					sigma22[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+10];
					sigma12[k.getPosition()*npoints+3] = sigma[k.getPosition()*npoints*3+11];
					sigma11[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+12];
					sigma22[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+13];
					sigma12[k.getPosition()*npoints+4] = sigma[k.getPosition()*npoints*3+14];
					sigma11[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+15];
					sigma22[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+16];
					sigma12[k.getPosition()*npoints+5] = sigma[k.getPosition()*npoints*3+17];
				}
				
				epsilon11[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3];
				epsilon22[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3+1];
				epsilon12[k.getPosition()*npoints] = epsilon[k.getPosition()*npoints*3+2];
				epsilon11[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+3];
				epsilon22[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+4];
				epsilon12[k.getPosition()*npoints+1] = epsilon[k.getPosition()*npoints*3+5];
				epsilon11[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+6];
				epsilon22[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+7];
				epsilon12[k.getPosition()*npoints+2] = epsilon[k.getPosition()*npoints*3+8];
				
				if(npoints > 3)
				{
					epsilon11[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+9];
					epsilon22[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+10];
					epsilon12[k.getPosition()*npoints+3] = epsilon[k.getPosition()*npoints*3+11];
					epsilon11[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+12];
					epsilon22[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+13];
					epsilon12[k.getPosition()*npoints+4] = epsilon[k.getPosition()*npoints*3+14];
					epsilon11[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+15];
					epsilon22[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+16];
					epsilon12[k.getPosition()*npoints+5] = epsilon[k.getPosition()*npoints*3+17];
				}  
				
				for(size_t l = 0 ; l < k->getBoundingPoints().size() ; l++)
				{
					Vector vm0(0., 3) ;
					k->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, k->getBoundingPoint(l), vm0, false) ;
					vonMises[k.getPosition()*k->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					Vector agl(0., 1) ;
					k->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, k->getBoundingPoint(l), agl, false) ;
					angle[k.getPosition()*k->getBoundingPoints().size()+l]  = agl[0] ;
				}
				
				double ar = k->area() ;
				for(int l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k.getPosition()*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k.getPosition()*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k.getPosition()*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k.getPosition()*npoints+l]/npoints)*ar;
				}
				
				if(k->getEnrichmentFunctions().size() == 0)
				{
					for(int l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k.getPosition()*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k.getPosition()*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k.getPosition()*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k.getPosition()*npoints+l]/npoints)*ar;
						
					}
					nogel_area+= ar ;
				}

			}
			else
			{
				sigma11[k.getPosition()*npoints] = 0 ;
				sigma22[k.getPosition()*npoints] = 0 ;
				sigma12[k.getPosition()*npoints] = 0 ;
				sigma11[k.getPosition()*npoints+1] = 0 ;
				sigma22[k.getPosition()*npoints+1] = 0 ;
				sigma12[k.getPosition()*npoints+1] = 0 ;
				sigma11[k.getPosition()*npoints+2] = 0 ;
				sigma22[k.getPosition()*npoints+2] = 0 ;
				sigma12[k.getPosition()*npoints+2] = 0 ;
				
				if(npoints >3)
				{
					sigma11[k.getPosition()*npoints+3] = 0 ;
					sigma22[k.getPosition()*npoints+3] = 0 ;
					sigma12[k.getPosition()*npoints+3] = 0 ;
					sigma11[k.getPosition()*npoints+4] = 0 ;
					sigma22[k.getPosition()*npoints+4] = 0 ;
					sigma12[k.getPosition()*npoints+4] = 0 ;
					sigma11[k.getPosition()*npoints+5] = 0 ;
					sigma22[k.getPosition()*npoints+5] = 0 ;
					sigma12[k.getPosition()*npoints+5] =0 ;
				}
				
				epsilon11[k.getPosition()*npoints] = 0 ;
				epsilon22[k.getPosition()*npoints] = 0 ;
				epsilon12[k.getPosition()*npoints] = 0 ;
				epsilon11[k.getPosition()*npoints+1] = 0 ;
				epsilon22[k.getPosition()*npoints+1] = 0 ;
				epsilon12[k.getPosition()*npoints+1] = 0 ;
				epsilon11[k.getPosition()*npoints+2] = 0 ;
				epsilon22[k.getPosition()*npoints+2] = 0 ;
				epsilon12[k.getPosition()*npoints+2] = 0 ;
				
				if(npoints > 3)
				{
					epsilon11[k.getPosition()*npoints+3] = 0 ;
					epsilon22[k.getPosition()*npoints+3] = 0 ;
					epsilon12[k.getPosition()*npoints+3] =0 ;
					epsilon11[k.getPosition()*npoints+4] = 0 ;
					epsilon22[k.getPosition()*npoints+4] = 0 ;
					epsilon12[k.getPosition()*npoints+4] =0 ;
					epsilon11[k.getPosition()*npoints+5] = 0 ;
					epsilon22[k.getPosition()*npoints+5] =0 ;
					epsilon12[k.getPosition()*npoints+5] = 0 ;
				}  
				
				for(size_t l = 0 ; l < k->getBoundingPoints().size() ; l++)
				{
					vonMises[k.getPosition()*k->getBoundingPoints().size()+l]  = 0 ;
					angle[k.getPosition()*k->getBoundingPoints().size()+l]  = 0 ;
				}
			}
		}
		std::string filename("triangles_") ;
		switch(ref)
		{
			case CIRCLE:
				filename.append("circle_") ;
				break ;
			case ELLIPSE:
				filename.append("ellipse_") ;
				break ;
			case TRIANGLE:
				filename.append("triangle_") ;
				break ;
			case RECTANGLE:
				filename.append("square_") ;
				break ;
            default:
                std::cout << "unhadled inclusino geometry type" << std::endl ;
                exit(0) ;
		}
		filename.append(itoa(samplingNumber, 10)) ;
		filename.append("_") ;
		filename.append(itoa(i, 10)) ;
		std::cout << filename << std::endl ;

		TriangleWriter writer(filename, featureTree) ;
//		writer.getField(TWFT_STRAIN_AND_STRESS) ;
//		writer.getField(TWFT_VON_MISES) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.write() ;
//		exit(0) ;

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
		
		std::cout << "average sigma11 : " << avg_s_xx/area << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy/area << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy/area << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx/area << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;
		
		std::cout << "apparent extension (y) " << (e_yy_max - e_yy_min)/sample.height() << std::endl ;
		std::cout << "apparent extension (x) " << (e_xx_max - e_xx_min)/sample.width() << std::endl ;

		std::cout << tries << std::endl ;

        if (tries < maxtries && featureTree->solverConverged() && go_on)
		{
			double delta_r = sqrt(aggregateArea*.05/((double)zones.size()*M_PI))/(double)nstepstot ;
            std::cout << "delta_r => " << delta_r << std::endl ;
			if(!featureTree->solverConverged())
				delta_r *= .01 ;

			double reactedArea = 0 ;
			
			Feature * current = nullptr ;
			if(!zones.empty())
				current = zones[0].feature() ;
			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;
			for(size_t z = 0 ; z < zones.size() ; z++)
			{
				zones[z].zone->setRadius(zones[z].zone->getRadius()+delta_r) ;	
				if(zones[z].is(current))
				{
					current_area += zones[z].zone->area() ;
					current_number++ ;
				}
				else
				{
					if(current_area/zones[z-1].areaOfInclusion() > 0.05)
					{
						stopped_reaction++ ;
						for(int m = 0 ; m < current_number ; m++)
						{
							reactedArea -= zones[z-1-m].zone->area() ;
							zones[z-1-m].zone->setRadius(zones[z].zone->getRadius()-delta_r) ;
							reactedArea += zones[z-1-m].zone->area() ;
						}
					}
					current_area = zones[z].zone->area() ;
					current_number = 1 ;
					current = zones[z].feature() ;
				}
				reactedArea += zones[z].zone->area() ;
			}
			
			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in "<< stopped_reaction << " aggs."<< std::endl ;

		
			if(featureTree->solverConverged())
			{
				expansion_reaction.push_back(std::make_pair(reactedArea/placed_area, avg_e_xx/area)) ;
				expansion_stress_xx.push_back(std::make_pair((avg_e_xx)/(area), (avg_s_xx)/(area))) ;
				expansion_stress_yy.push_back(std::make_pair((avg_e_yy)/(area), (avg_s_yy)/(area))) ;
				apparent_extension.push_back( std::make_pair((e_xx_max - e_xx_min)/sample.width(), (e_yy_max - e_yy_min)/sample.width() ) ) ;
			}
			
//			if (tries >= maxtries)
//				break ;
		}
		for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
		{
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
	}

}

std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	double radius = 0.0000005 ;
	std::default_random_engine generator ;
	double w = sample.width()*0.5-radius*1000. ;
	double h = sample.width()*0.5-radius*1000. ;
	std::uniform_real_distribution< double > xpos(-w,w) ;
	std::uniform_real_distribution< double > ypos(-h,h) ;
  	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(int i = 0 ; i < n ; i++)
	{
		Point pos( xpos(generator), ypos(generator) ) ;
		pos += sample.getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*1000.)*(radius*1000.))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(size_t j = 0 ; j < incs.size() ; j++)
		{
			Circle circle(incs[j]->getRadius() - radius*1000., incs[j]->getCenter()) ;
			if(circle.in(zonesToPlace[i]->getCenter()))
			{
				std::cout << dist(zonesToPlace[i]->getCenter(), incs[j]->getCenter() ) << std::endl ;
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(Zone(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
		
		if((int)ret.size() == max)
		  break ;
	}
//	exit(0) ;
	
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}


std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<EllipsoidalInclusion * > & incs , FeatureTree & F)
{
	std::default_random_engine generator ;
	double radius = 0.0000005 ;
	double w = sample.width()*0.5-radius*1000. ;
	double h = sample.width()*0.5-radius*1000. ;
	std::uniform_real_distribution< double > xpos(-w,w) ;
	std::uniform_real_distribution< double > ypos(-h,h) ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(int i = 0 ; i < n ; i++)
	{
		Point pos( xpos(generator), ypos(generator) ) ;
		pos += sample.getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*1000.+radius*1000.)*(radius*1000.+radius*1000.))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<EllipsoidalInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(size_t j = 0 ; j < incs.size() ; j++)
		{
			Ellipse ellipse(incs[j]->getCenter(), incs[j]->getMajorAxis()*0.9, incs[j]->getMinorAxis()*0.9) ;
			if(ellipse.in(zonesToPlace[i]->getCenter()))
			{
				if(!incs[j]->in(zonesToPlace[i]->getCenter())) {
					std::cout << i << ";" << j << "||" ;
				}
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(Zone(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
		if((int)ret.size() == max)
		  break ;
	}
	
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<TriangularInclusion * > & incs , FeatureTree & F)
{
	std::default_random_engine generator ;
	double radius = 0.0000005 ;
	double w = sample.width()*0.5-radius*1000. ;
	double h = sample.width()*0.5-radius*1000. ;
	std::uniform_real_distribution< double > xpos(-w,w) ;
	std::uniform_real_distribution< double > ypos(-h,h) ;
	std::vector<Zone> ret ;
	return ret ;
	aggregateArea = 0 ;
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(int i = 0 ; i < n ; i++)
	{
		Point pos( xpos(generator), ypos(generator) ) ;
		pos += sample.getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*1000.+radius*1000.)*(radius*1000.+radius*1000.))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), gel)) ;
	}
	std::cout << zonesToPlace.size() << " zones generated" << std::endl ;
	std::map<TriangularInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(size_t j = 0 ; j < incs.size() ; j++)
		{
			Triangle triangle(incs[j]->getBoundingPoint(0), incs[j]->getBoundingPoint(1), incs[j]->getBoundingPoint(2)) ;
			Point c = zonesToPlace[i]->getCenter() ;
			triangle.project(&c) ;
			if(triangle.in(zonesToPlace[i]->getCenter()) && dist(zonesToPlace[i]->getCenter(), c) > radius*1000.)
			{
				if(!incs[j]->in(zonesToPlace[i]->getCenter())) {
					std::cout << i << ";" << j << "||" ;
				}
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(Zone(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
		if((int)ret.size() == max)
		  break ;
	}
	
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<RectangularInclusion * > & incs , FeatureTree & F)
{
	std::default_random_engine generator ;
	double radius = 0.0000005 ;
	double w = sample.width()*0.5-radius*1000. ;
	double h = sample.width()*0.5-radius*1000. ;
	std::uniform_real_distribution< double > xpos(-w,w) ;
	std::uniform_real_distribution< double > ypos(-h,h) ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(int i = 0 ; i < n ; i++)
	{
		Point pos(xpos(generator),ypos(generator)) ;
		pos += sample.getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*60.+radius*60.)*(radius*60.+radius*60.))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<RectangularInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(size_t j = 0 ; j < incs.size() ; j++)
		{
			OrientedRectangle rectangle(incs[j]->getBoundingPoint(0), incs[j]->getBoundingPoint(1), incs[j]->getBoundingPoint(2), incs[j]->getBoundingPoint(3)) ;
			Point c = zonesToPlace[i]->getCenter() ;
			rectangle.project(&c) ;
			if(rectangle.in(zonesToPlace[i]->getCenter()) && dist(zonesToPlace[i]->getCenter(), c) > radius*60.)
			{
				if(!incs[j]->in(zonesToPlace[i]->getCenter())) {
					std::cout << i << ";" << j << "||" ;
				}
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(Zone(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
		if((int)ret.size() == max)
		  break ;
	}
	
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

TriangularInclusion * rotate(TriangularInclusion * tri, double alpha)
{
	size_t n = tri->getBoundingPoints().size() ;
	Point A = tri->getBoundingPoint(0) ;
	Point B = tri->getBoundingPoint(n/3) ;
	Point C = tri->getBoundingPoint(n*2/3) ;
	
	Point center = tri->getCenter() ;
	
	A -= center ;
	B -= center ;
	C -= center ;
	
	double a = A.angle()+alpha ;
	double b = B.angle()+alpha ;
	double c = C.angle()+alpha ;
	
	double ax = A.norm() ;
	double bx = B.norm() ;
	double cx = C.norm() ;
	
	A = Point(std::cos(a), std::sin(a))*ax ;
	B = Point(std::cos(b), std::sin(b))*bx ;
	C = Point(std::cos(c), std::sin(c))*cx ;

	A += center ;
	B += center ;
	C += center ;
	
	return new TriangularInclusion(tri->getFather(), A, B, C) ;
	
}

RectangularInclusion * rotate(RectangularInclusion * rec, double alpha)
{
	Point A = rec->getBoundingPoint(0) ;
	Point B = rec->getBoundingPoint(1) ;
	Point C = rec->getBoundingPoint(2) ;
	Point D = rec->getBoundingPoint(3) ;
	
	Point center = rec->getCenter() ;
	
	A -= center ;
	B -= center ;
	C -= center ;
	D -= center ;
	
	double a = A.angle()+alpha ;
	double b = B.angle()+alpha ;
	double c = C.angle()+alpha ;
	double d = D.angle()+alpha ;
	
	double ax = A.norm() ;
	double bx = B.norm() ;
	double cx = C.norm() ;
	double dx = D.norm() ;
	
	A = Point(std::cos(a), std::sin(a))*ax ;
	B = Point(std::cos(b), std::sin(b))*bx ;
	C = Point(std::cos(c), std::sin(c))*cx ;
	D = Point(std::cos(d), std::sin(d))*dx ;
	
	A += center ;
	B += center ;
	C += center ;
	D += center ;
	
	return new RectangularInclusion(rec->getFather(), A, B, C, D) ;
	
}

EllipsoidalInclusion * rotate(EllipsoidalInclusion * ell, double alpha)
{
	Point A = ell->getMajorAxis() ;
	Point B = ell->getMinorAxis() ;
	
	double a = A.angle()+alpha ;
	double b = B.angle()+alpha ;
	
	double ax = A.norm() ;
	double bx = B.norm() ;
	
	A = Point(std::cos(a), std::sin(a))*ax ;
	B = Point(std::cos(b), std::sin(b))*bx ;
	
	return new EllipsoidalInclusion(ell->getFather(), ell->getCenter(), A, B) ;
	
}

bool tintersects(std::vector<TriangularInclusion *> triinc, int index, RectangularFeature * box)
{
	if(triinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(int i = 0 ; i < index ; i++)
	{
		if(triinc[index]->getPrimitive()->intersects(triinc[i]->getPrimitive()))
		{
			return true ;
		}
	}

	for(size_t i = index+1 ; i < triinc.size() ; i++)
	{
		if(triinc[index]->getPrimitive()->intersects(triinc[i]->getPrimitive()))
		{
			return true ;
		}
	}
	return false ;
}

bool tintersects(std::vector<RectangularInclusion *> recinc, int index, RectangularFeature * box)
{
	if(recinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(int i = 0 ; i < index ; i++)
	{
		if(recinc[index]->getPrimitive()->intersects(recinc[i]->getPrimitive()))
		{
			return true ;
		}
	}
	
	for(size_t i = index+1 ; i < recinc.size() ; i++)
	{
		if(recinc[index]->getPrimitive()->intersects(recinc[i]->getPrimitive()))
		{
			return true ;
		}
	}
	return false ;
}

bool tintersects(std::vector<EllipsoidalInclusion *> ellinc, int index, RectangularFeature * box)
{
	if(ellinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(int i = 0 ; i < index ; i++)
	{
		if(ellinc[index]->getPrimitive()->intersects(ellinc[i]->getPrimitive()))
		{
			return true ;
		}
	}
	
	for(size_t i = index+1 ; i < ellinc.size() ; i++)
	{
		if(ellinc[index]->getPrimitive()->intersects(ellinc[i]->getPrimitive()))
		{
			return true ;
		}
	}
	return false ;
}

bool rotateUntilNoIntersection(std::vector<TriangularInclusion *> & triinc, int index, RectangularFeature * box)
{
	std::default_random_engine generator ;
	std::uniform_real_distribution<double> alpha(0.001*M_PI, 0.01*M_PI) ;
	if(tintersects(triinc, index, box))
	{
		int i = 0 ;
		while(i < 150)
		{
			triinc[index] = rotate(triinc[index],alpha(generator)) ;
			if(!tintersects(triinc, index, box))
			{
				return true ;
			}
			i++ ;
		}
		return false ;
	}
	else
		return true ;
}

bool rotateUntilNoIntersection(std::vector<RectangularInclusion *> & recinc, int index, RectangularFeature * box)
{
	std::default_random_engine generator ;
	std::uniform_real_distribution<double> alpha(0.15*M_PI, 0.95*M_PI) ;
	if(tintersects(recinc, index, box))
	{
		int i = 0 ;
		while(i < 500)
		{
			recinc[index] = rotate(recinc[index],alpha(generator)) ;
			if(!tintersects(recinc, index, box))
			{
				return true ;
			}
			i++ ;
		}
		return false ;
	}
	else
		return true ;
}

bool rotateUntilNoIntersection(std::vector<EllipsoidalInclusion *> & ellinc, int index, RectangularFeature * box)
{
	std::default_random_engine generator ;
	std::uniform_real_distribution<double> alpha(0.15*M_PI, 0.95*M_PI) ;
	if(tintersects(ellinc, index, box))
	{
		int i = 0 ;
		while(i < 50)
		{
			ellinc[index] = rotate(ellinc[index],alpha(generator)) ;
			if(!tintersects(ellinc, index, box))
			{
				return true ;
			}
			i++ ;
		}
		return false ;
	}
	else
		return true ;
}

double getAverageDamage( FeatureTree & F)
{
	double area = 0. ;
	double damage = 0. ;
	for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
	{
		double ar = i->area() ;
		double d = 0. ;
		if(i->getBehaviour() && i->getBehaviour()->getDamageModel())
		{
			d = i->getBehaviour()->getDamageModel()->getState().max() ;
		}
		damage += d*ar ;
		area += ar ;
	}
	return damage/area ;
}

int main(int argc, char *argv[])
{
	timeval time0, time1 ;
	gettimeofday(&time0, nullptr);
  	
	FeatureTree F(&sample) ;
	featureTree = &F ;

 	int inclusionNumber = 500 ;
	sample.setBehaviour( new PasteBehaviour(true) ) ;

	AggregateBehaviour * agg = new AggregateBehaviour(true) ;
	std::vector<Feature *> inclusions = PSDGenerator::get2DConcrete(&F, agg, inclusionNumber, 0.008, 0.00001,new PSDBolomeA(), new PolygonalInclusionGenerator( 3 ) ,inclusionNumber*100) ;

	std::vector<TriangularInclusion *> ellinc ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		ellinc.push_back(dynamic_cast<TriangularInclusion *>(inclusions[i])) ;

	std::cout << inclusions.size() << "\t" << ellinc.size() << std::endl ;

	zones = generateExpansiveZonesHomogeneously(1000, 100, ellinc, F) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, -0.00005*1e-3) ;
	F.addBoundaryCondition(stress) ;

	int nSampling = atof(argv[1]) ;
	
	F.setSamplingNumber(nSampling) ;
	F.setOrder(LINEAR) ;
//	F.setMaxIterationsPerStep(4000) ;
	

	std::string hop = "rag_triangle_out_" ;
	hop.append(std::string(argv[1])) ;
	std::fstream out ;
	out.open(hop.c_str(), std::ios::out) ;

	for(size_t i = 0 ; i < 30 ; i++)
	{
		F.step() ;
		Vector strain = F.getAverageField(TOTAL_STRAIN_FIELD) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD) ;
		double damage = getAverageDamage(F) ;
		std::cout << strain[0] << "\t" << strain[1] << std::endl ;

		std::vector<double> u_right ;
		std::vector<double> u_top ;
		Vector disp = F.getDisplacements() ;
		std::valarray<bool> done(disp.size()/2) ;
		done = false ;
		for(auto j = featureTree->get2DMesh()->begin() ; j != featureTree->get2DMesh()->end() ; j++)
		{
			for(size_t k = 0 ; k < j->getBoundingPoints().size() ; k++)
			{
				Point p = j->getBoundingPoint(k) ;
				if(!done[p.getId()])
				{
					if(p.getX() > 0.035*0.999)
						u_right.push_back(disp[p.getId()*2]) ;
					if(p.getY() > 0.035*0.999)
						u_top.push_back(disp[p.getId()*2+1]) ;
					done[p.getId()] = true ;
				}
			}
		}
		
		std::stable_sort( u_right.begin(), u_right.end() ) ;
		std::stable_sort( u_top.begin(), u_top.end() ) ;
		
		std::string filename = "rag_triangle_" ;
		filename.append(std::string(argv[1])) ;
		filename.append("_") ;
		filename.append(itoa(i)) ;
		TriangleWriter writer(filename, &F) ;
		writer.getField(REAL_STRESS_FIELD ) ;
		writer.getField(TOTAL_STRAIN_FIELD ) ;
		writer.getField(TWFT_VON_MISES) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.write() ;

		double delta_r = sqrt(aggregateArea*.005/((double)zones.size()*M_PI))/(double)30 ;
	        std::cout << "delta_r => " << delta_r << std::endl ;
		if(!F.solverConverged())
				delta_r *= .01 ;

		double reactedArea = 0 ;
			
		Feature * current = nullptr ;
		if(!zones.empty())
			current = zones[0].feature() ;
		double current_area = 0 ;
		int current_number = 0 ;
		int stopped_reaction = 0 ;
		for(size_t z = 0 ; z < zones.size() ; z++)
		{
			zones[z].zone->setRadius(zones[z].zone->getRadius()+delta_r) ;	
			if(zones[z].is(current))
			{
				current_area += zones[z].zone->area() ;
				current_number++ ;
			}
			else
			{
				if(current_area/zones[z-1].areaOfInclusion() > 0.05)
				{
					stopped_reaction++ ;
					for(int m = 0 ; m < current_number ; m++)
					{
						reactedArea -= zones[z-1-m].zone->area() ;
						zones[z-1-m].zone->setRadius(zones[z].zone->getRadius()-delta_r) ;
						reactedArea += zones[z-1-m].zone->area() ;
					}
				}
				current_area = zones[z].zone->area() ;
				current_number = 1 ;
				current = zones[z].feature() ;
			}
			reactedArea += zones[z].zone->area() ;
		}
			
		out << reactedArea/aggregateArea << "\t" << u_right[u_right.size()/2] << "\t" << u_top[u_top.size()/2] << "\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" << damage << std::endl ;
		std::cout << "reacted Area : " << reactedArea << ", reaction stopped in "<< stopped_reaction << " aggs."<< std::endl ;

		
//			if (tries >= maxtries)
//				break ;
	


//		step(reference, nSampling+i) ;
//		stress->setData( -0.005*(i+2)*1e-3 ) ;
	}
	
	double area = 0 ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
	    area += inclusions[i]->area() ;
	
	gettimeofday(&time1, nullptr);
	double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
	std::cout << delta*1e-6 << std::endl ;
	
	return 0 ;
}
