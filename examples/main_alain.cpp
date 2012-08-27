
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
#include "../physics/kelvinvoight.h"
#include "../physics/maxwell.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/microstructuregenerator.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/itoa.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/random.h"
#include "../utilities/writer/triangle_writer.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#include <sys/time.h>
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

double timepos = 0.00 ;
double percent = 0.10 ;
double placed_area = 0 ;

double stress = 15e6 ;

Sample sample(nullptr, 0.01, 0.01, 0.0, 0.0) ;

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

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

double shape ;
double orientation ;
double spread ;

MultiTriangleWriter writer( "triangles_head", "triangles_layers", nullptr ) ;

GelBehaviour * gel = new GelBehaviour() ;

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
		}
		return 0. ;
	}
	
} ;

std::vector<Zone> zones ;

void step(GeometryType ref, int samplingNumber, int stressXI, int stressETA)
{

	int nsteps = 10; // step 16 = 0.023
	int nstepstot = 10;
	int maxtries = 5 ;
	int tries = 0 ;
	featureTree->setMaxIterationsPerStep(1600) ;
	featureTree->setDeltaTime(24*60*3) ;
	
	for(size_t i = 0 ; i < nsteps ; i++)
	{
		if(i == 2)
			featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -5*1e6)) ;
  
		bool go_on = featureTree->step() ;
		
		if(featureTree->solverConverged())
		{
			cracked_volume.push_back(featureTree->crackedVolume) ;
			damaged_volume.push_back(featureTree->damagedVolume) ;
		}	
		triangles = featureTree->getElements2D() ;
		x.resize(featureTree->getDisplacements(-1, false).size()) ;
		x = featureTree->getDisplacements(-1, false) ;
		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain(-1, false) ;
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
		
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
			bool in = false ;
			for(size_t m = 0 ; m < tris__.size() ; m++)
			{
				if(triangles[k] == tris__[m])
				{
					in = true ;
					break ;
				}
			}
			cracked.push_back(in) ;
			
			
			
			if(!in && !triangles[k]->getBehaviour()->fractured())
			{
				
				for(size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++)
				{
					if(x[triangles[k]->getBoundingPoint(p).id*2] > x_max)
						x_max = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2] < x_min)
						x_min = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] > y_max)
						y_max = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] < y_min)
						y_min = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(triangles[k]->getBoundingPoint(p).x > sample.width()*.4999)
					{
						if(e_xx_max < x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_max=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).x < -sample.width()*.4999)
					{
						if(e_xx_min > x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_min=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).y > sample.height()*.4999)
					{
						if(e_yy_max < x[triangles[k]->getBoundingPoint(p).id*2+1])
							e_yy_max=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).y < -sample.height()*.4999)
					{
						if(e_yy_min > x[triangles[k]->getBoundingPoint(p).id*2+1])
							e_yy_min=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
// 						ex_count++ ;
					}
				}
				area += triangles[k]->area() ;
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] > E_max)
						E_max = triangles[k]->getBehaviour()->param[0][0] ;
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] < E_min)
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
				
//				std::cout << sigma22[k*npoints] << "\t" << sigma22[k*npoints+1] << "\t" << sigma22[k*npoints+2] << std::endl ;
				
				if(npoints >3)
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
				
				if(npoints > 3)
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
				
				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					Vector vm0(0., 3) ;
					triangles[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, triangles[k]->getBoundingPoint(l), vm0, false) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					Vector agl(0., 1) ;
					triangles[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, triangles[k]->getBoundingPoint(l), agl, false) ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl[0] ;
				}
				
				double ar = triangles[k]->area() ;
				for(size_t l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
				}
				
				if(triangles[k]->getEnrichmentFunctions().size() == 0)
				{
					for(size_t l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k*npoints+l]/npoints)*ar;
						
					}
					nogel_area+= ar ;
				}

			}
			else
			{
				sigma11[k*npoints] = 0 ;
				sigma22[k*npoints] = 0 ;
				sigma12[k*npoints] = 0 ;
				sigma11[k*npoints+1] = 0 ;
				sigma22[k*npoints+1] = 0 ;
				sigma12[k*npoints+1] = 0 ;
				sigma11[k*npoints+2] = 0 ;
				sigma22[k*npoints+2] = 0 ;
				sigma12[k*npoints+2] = 0 ;
				
				if(npoints >3)
				{
					sigma11[k*npoints+3] = 0 ;
					sigma22[k*npoints+3] = 0 ;
					sigma12[k*npoints+3] = 0 ;
					sigma11[k*npoints+4] = 0 ;
					sigma22[k*npoints+4] = 0 ;
					sigma12[k*npoints+4] = 0 ;
					sigma11[k*npoints+5] = 0 ;
					sigma22[k*npoints+5] = 0 ;
					sigma12[k*npoints+5] =0 ;
				}
				
				epsilon11[k*npoints] = 0 ;
				epsilon22[k*npoints] = 0 ;
				epsilon12[k*npoints] = 0 ;
				epsilon11[k*npoints+1] = 0 ;
				epsilon22[k*npoints+1] = 0 ;
				epsilon12[k*npoints+1] = 0 ;
				epsilon11[k*npoints+2] = 0 ;
				epsilon22[k*npoints+2] = 0 ;
				epsilon12[k*npoints+2] = 0 ;
				
				if(npoints > 3)
				{
					epsilon11[k*npoints+3] = 0 ;
					epsilon22[k*npoints+3] = 0 ;
					epsilon12[k*npoints+3] =0 ;
					epsilon11[k*npoints+4] = 0 ;
					epsilon22[k*npoints+4] = 0 ;
					epsilon12[k*npoints+4] =0 ;
					epsilon11[k*npoints+5] = 0 ;
					epsilon22[k*npoints+5] =0 ;
					epsilon12[k*npoints+5] = 0 ;
				}  
				
				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
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
		}
		filename.append(itoa(samplingNumber, 10)) ;
		filename.append("_") ;
		filename.append(itoa(stressXI, 10)) ;
		filename.append("_") ;
		filename.append(itoa(stressETA, 10)) ;
		filename.append("_") ;
		filename.append(itoa(i, 10)) ;
		std::cout << filename << std::endl ;

		TriangleWriter wrt(filename, featureTree) ;
		wrt.getField(TWFT_STRAIN_AND_STRESS) ;
		wrt.getField(TWFT_STIFFNESS) ;
		wrt.getField(TWFT_DAMAGE) ;
		wrt.write() ;
// 		if(go_on)
// 		{
// 			writer.reset( featureTree ) ;
// 			writer.getField( TWFT_PRINCIPAL_STRESS ) ;
// 			writer.getField( TWFT_PRINCIPAL_STRAIN ) ;
// 			writer.getField( TWFT_CRITERION ) ;
// 			writer.getField( TWFT_STIFFNESS ) ;
// 			writer.getField( TWFT_DAMAGE ) ;
// 			writer.append() ;
// 		}

		
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
		
		std::cout << "average sigma11 : " << avg_s_xx_nogel/nogel_area << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy_nogel/nogel_area << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy_nogel/nogel_area << std::endl ;
		
		std::cout << "apparent extension X " << e_xx_max-e_xx_min << std::endl ;
		std::cout << "apparent extension Y " << e_yy_max-e_yy_min << std::endl ;

		std::cout << tries << std::endl ;

		if (go_on && featureTree->solverConverged())
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
						for(size_t m = 0 ; m < current_number ; m++)
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
				apparent_extension.push_back(std::make_pair(e_xx_max-e_xx_min, e_yy_max-e_yy_min)) ;
			}
			tries = 0 ;
			
//			if (tries >= maxtries)
//				break ;
		}
		else
		{
			std::cout << "maximum tries reached..." << std::endl ;
		}
		for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
		{
			std::cout << /*expansion_reaction[i].first << "   " 
			<< expansion_reaction[i].second << "   " 
			<< expansion_stress_xx[i].first << "   " 
			<< expansion_stress_xx[i].second << "   " 
			<<*/ expansion_stress_yy[i].first << "\t" 
			<< expansion_stress_yy[i].second << "   " 
			/*<< apparent_extension[i].first  << "   "
			<< apparent_extension[i].second  << "   "
			<< cracked_volume[i]  << "   "
			<< damaged_volume[i]  << "   " */
			<< std::endl ;
		}
	}

}

std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
  
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample.width()*0.5-radius*60 ;
		double h = sample.width()*0.5-radius*60 ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
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
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			Circle circle(incs[j]->getRadius() - radius*60, incs[j]->getCenter()) ;
			if(circle.in(zonesToPlace[i]->getCenter()))
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
		
		if(ret.size() == max)
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

std::vector<Zone> generateExpansiveZonesAtCenter(std::vector<EllipsoidalInclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
 	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		Point pos = incs[i]->getCenter() ;
		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
		F.addFeature(incs[i],zonesToPlace[i]) ;
		ret.push_back(Zone(zonesToPlace[i],incs[i])) ;
		aggregateArea+= incs[i]->area() ;
	}

	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesAtCenter(std::vector<TriangularInclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		Point pos = incs[i]->getCenter() ;
		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x+incs[i]->getRadius()*0.5, pos.y, gel)) ;
		F.addFeature(incs[i],zonesToPlace[i*2]) ;
		F.addFeature(incs[i],zonesToPlace[i*2+1]) ;
		ret.push_back(Zone(zonesToPlace[i*2],incs[i])) ;
		ret.push_back(Zone(zonesToPlace[i*2+1],incs[i])) ;
		aggregateArea+= incs[i]->area() ;
	}

	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesAtCenter(std::vector<RectangularInclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		Point pos = incs[i]->getCenter() ;
		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
		F.addFeature(incs[i],zonesToPlace[i]) ;
		ret.push_back(Zone(zonesToPlace[i],incs[i])) ;
		aggregateArea+= incs[i]->area() ;
	}

	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesAtCenter(std::vector<Inclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		Point pos = incs[i]->getCenter() ;
		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
//		zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x+incs[i]->getRadius()*0.5, pos.y, gel)) ;
		F.addFeature(incs[i],zonesToPlace[i]) ;
//		F.addFeature(incs[i],zonesToPlace[i*2+1]) ;
		ret.push_back(Zone(zonesToPlace[i],incs[i])) ;
//		ret.push_back(Zone(zonesToPlace[i*2+1],incs[i])) ;
		aggregateArea+= incs[i]->area() ;
	}

	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<Zone> generateExpansiveZonesHomogeneously(int n, int max, std::vector<EllipsoidalInclusion * > & incs , FeatureTree & F)
{
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample.width()*0.5-radius*60 ;
		double h = sample.width()*0.5-radius*60 ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
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
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<EllipsoidalInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			Ellipse ellipse(incs[j]->getCenter(), incs[j]->getMajorAxis()*0.95, incs[j]->getMinorAxis()*0.95) ;
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
		if(ret.size() == max)
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
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample.width()*0.5-radius*60 ;
		double h = sample.width()*0.5-radius*60 ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
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
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<TriangularInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			Triangle triangle(incs[j]->getBoundingPoint(0), incs[j]->getBoundingPoint(1), incs[j]->getBoundingPoint(2)) ;
			Point c = zonesToPlace[i]->getCenter() ;
			triangle.project(&c) ;
			if(triangle.in(zonesToPlace[i]->getCenter()) && dist(zonesToPlace[i]->getCenter(), c) > radius*60.)
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
		if(ret.size() == max)
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
	RandomNumber gen ;
	std::vector<Zone> ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;

	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample.width()*0.5-radius*60 ;
		double h = sample.width()*0.5-radius*60 ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
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
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, gel)) ;
	}
	std::cout << zonesToPlace.size() << std::endl ;
	std::map<RectangularInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
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
		if(ret.size() == max)
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


bool tintersects(std::vector<TriangularInclusion *> triinc, int index, Sample * box)
{
	if(triinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(size_t i = 0 ; i < index ; i++)
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

bool tintersects(std::vector<RectangularInclusion *> recinc, int index, Sample * box)
{
	if(recinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(size_t i = 0 ; i < index ; i++)
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

bool tintersects(std::vector<EllipsoidalInclusion *> ellinc, int index, Sample * box)
{
	if(ellinc[index]->getPrimitive()->intersects(box->getPrimitive()))
	{
		return true ;
	}
	
	for(size_t i = 0 ; i < index ; i++)
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


bool rotateUntilNoIntersection(std::vector<TriangularInclusion *> & triinc, int index, Sample * box)
{
	if(tintersects(triinc, index, box))
	{
		TriangularInclusion * next ;
		int i = 0 ;
		while(i < 150)
		{
			double alpha = UniformDistribution(0.001*M_PI,0.01*M_PI).draw() ;
			triinc[index] = rotate(triinc[index],alpha) ;
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

bool rotateUntilNoIntersection(std::vector<RectangularInclusion *> & recinc, int index, Sample * box)
{
	if(tintersects(recinc, index, box))
	{
		RectangularInclusion * next ;
		int i = 0 ;
		while(i < 500)
		{
			double alpha = UniformDistribution(0.15*M_PI,0.95*M_PI).draw() ;
			recinc[index] = rotate(recinc[index],alpha) ;
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

bool rotateUntilNoIntersection(std::vector<EllipsoidalInclusion *> & ellinc, int index, Sample * box)
{
	if(tintersects(ellinc, index, box))
	{
		EllipsoidalInclusion * next ;
		int i = 0 ;
		while(i < 50)
		{
			double alpha = UniformDistribution(0.15*M_PI,0.95*M_PI).draw() ;
			ellinc[index] = rotate(ellinc[index],alpha) ;
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

struct RowMatrix
{
  Vector * values ;
  const size_t col ;
  
  RowMatrix(size_t n) : col(n) { values = new Vector(0., n*n) ; }
  
  double & operator() (size_t i, size_t j) 
  {
    return (*values)[i*col + j] ;
  }

  
  const double & operator() (size_t i, size_t j)  const
  {
    return (*values)[i*col + j] ;
  }
  
  Vector operator*(const Vector & v)
  {
    Vector ret(0.,col) ;
//    #pragma omp parallel for
    for(size_t r = 0 ; r < col ; r++)
    {
      size_t iterator = r*col ;
      for(size_t c = 0 ; c < col ; c++)
      {
	ret[r] += (*values)[iterator++]*v[c] ;
      }
    }
    return ret ;
  }

  void clear()
  {
   values->resize(0) ;
  }
  
  void print()
  {
    for(size_t r = 0 ; r < col ; r++)
    {
      for(size_t c = 0 ; c < col ; c++)
      {
	std::cout << (*this)(r,c) << "\t" ;
      }
    std::cout << std::endl ;
    }
  }
  
} ;

struct ColumnMatrix
{
  Vector * values ;
  const size_t dim ;
  
  ColumnMatrix(size_t n) : dim(n) { values = new Vector(0., dim*dim) ;}
  
  double & operator() (size_t i, size_t j) 
  {
    return (*values)[i + j*dim] ;
  }

  
  const double & operator() (size_t i, size_t j)  const
  {
    return (*values)[i + j*dim] ;
  }
  
  Vector operator*(const Vector & v)
  {
    Vector ret(0.,dim) ;
    for(size_t c = 0 ; c < dim ; ++c)
    {
      double vc = v[c] ;
      double * retc_iter = &ret[0] ;
      size_t iterator = c*dim ;
      for(size_t r = 0 ; r < dim ; ++r)
      {
	*(retc_iter++) += (*values)[iterator++]*vc ;
      }
    }
    return ret ;
  }
  
  void clear()
  {
   values->resize(0) ;
  }

  void multiply(const Vector & v, Vector & ret)
  {
    for(size_t c = 0 ; c < dim ; ++c)
    {
      double vc = v[c] ;
      double * retc_iter = &ret[0] ;
      size_t iterator = c*dim ;
      for(size_t r = 0 ; r < dim ; ++r)
      {
	*(retc_iter++) += (*values)[iterator++]*vc ;
      }
    }
  }
  
  void print()
  {
    for(size_t r = 0 ; r < dim ; r++)
    {
      for(size_t c = 0 ; c < dim ; c++)
      {
	std::cout << (*this)(r,c) << "\t" ;
      }
    std::cout << std::endl ;
    }
  }
  
  
} ;

struct BlockColumnMatrix
{
  ColumnMatrix AA ;
  ColumnMatrix AB ;
  ColumnMatrix BA ;
  ColumnMatrix BB ;
  size_t dim ;
  
  BlockColumnMatrix(size_t dim) : dim(dim), AA(dim/2), BB(dim/2), AB(dim/2), BA(dim/2) { }
  
  double & operator() (size_t i, size_t j) 
  {
    if(i < dim/2 && j < dim/2)
      return AA(i,j) ;
    if(i < dim/2 && j >= dim/2)
      return AB(i,j-dim/2) ;
    if(i >= dim/2 && j < dim/2)
      return BA(i-dim/2,j) ;
    if(i >= dim/2 && j >= dim/2)
      return BB(i-dim/2,j-dim/2) ;
  }

  
  const double & operator() (size_t i, size_t j)  const
  {
    if(i < dim/2 && j < dim/2)
      return AA(i,j) ;
    if(i < dim/2 && j >= dim/2)
      return AB(i,j-dim/2) ;
    if(i >= dim/2 && j < dim/2)
      return BA(i-dim/2,j) ;
    if(i >= dim/2 && j >= dim/2)
      return BB(i-dim/2,j-dim/2) ;
  }
  
  Vector operator*(const Vector & v)
  {
    Vector rAA(0., dim/2) ;
    Vector rAB(0., dim/2) ;
    Vector rBA(0., dim/2) ;
    Vector rBB(0., dim/2) ;

    Vector vAA(&v[0], dim/2) ;
    Vector vAB(&v[dim/2], dim/2) ;
    Vector vBA(&v[0], dim/2) ;
    Vector vBB(&v[dim/2], dim/2) ;
    
#pragma omp sections 
    {
#pragma omp section 
{
    AA.multiply(vAA, rAA) ;
}
#pragma omp section  
{
    AB.multiply(vAB, rAB) ;
}
#pragma omp section  
{
    BA.multiply(vBA, rBA) ;
}
#pragma omp section  
{
    BB.multiply(vBB, rBB) ;
}
    }
    
    Vector r(0., dim) ;
//     for(size_t i = 0 ; i < dim/2 ; i++)
//     {
//       r[i] = rAA[i] + rAB[i] ;
//       r[i+dim/2] = rBA[i] + rBB[i] ;
//     }
    return r ;
  }
  
  void clear()
  {
    AA.clear() ;
    AB.clear() ;
    BA.clear() ;
    BB.clear() ;
  }
} ;


const size_t dim = 40000 ;

int main(int argc, char *argv[])
{
  
  
  
/*	GeometryType reference = CIRCLE ;
	
	if(std::string(argv[2]) == std::string("--ellipse"))
		reference = ELLIPSE ;
	if(std::string(argv[2]) == std::string("--triangle"))
		reference = TRIANGLE ;
	if(std::string(argv[2]) == std::string("--rectangle"))
		reference = RECTANGLE ;*/
	
	
	FeatureTree F(&sample) ;
	featureTree = &F ;
	
	Matrix C(3,3) ;
	Matrix C0(3,3) ;
	Matrix E(3,3) ;
	
	double young = 1e9 ;
	double nu = 0.3 ;
	double eta = 1 ;
	double day = 1. ;
	
	C[0][0] = 1. ; C[1][0] = nu ; C[2][0] = 0. ;
	C[0][1] = nu ; C[1][1] = 1. ; C[2][1] = 0. ;
	C[0][2] = 0. ; C[1][2] = 0. ; C[2][2] = 1.-nu ;
	
	C *= young/(1.-nu*nu) ;
	E = C*eta*day ;
	C0 = C*3 ;
	
	sample.setBehaviour( new StandardLinearSolid(C0, C,E) ) ;

/*	std::vector<Feature *> feats ;
	std::vector<Inclusion *> inclusions ;
//	inclusions.push_back(new Inclusion(0.002,0.,0.)) ;
	std::vector<EllipsoidalInclusion *> ellinc ;
	std::vector<TriangularInclusion *> triinc ;
	std::vector<RectangularInclusion *> recinc ;

	if(reference != CIRCLE)
	{
		InclusionConverter cnv(reference, new ConstantDistribution(1.), new ConstantDistribution(1.), new ConstantDistribution(0.)) ;
		if(reference == ELLIPSE)
		{
			cnv.setAspectRatio(0.5) ;
		}
		feats = cnv.convert(inclusions) ;
	}
	
	switch(reference)
	{
		case ELLIPSE:
			ellinc.push_back(dynamic_cast<EllipsoidalInclusion *>(feats[0])) ;
			break ;
		case TRIANGLE:
			triinc.push_back(dynamic_cast<TriangularInclusion *>(feats[0])) ;
			std::cout << inclusions[0]->area() << "\t" << triinc[0]->area() << std::endl ;
			break ;
		case RECTANGLE:
			recinc.push_back(dynamic_cast<RectangularInclusion *>(feats[0])) ;
			break ;
	}
	
		
	sample.setBehaviour(new ViscoElasticOnlyPasteBehaviour()) ;
// 	ElasticOnlyAggregateBehaviour * agg = new ElasticOnlyAggregateBehaviour() ;
// 	placed_area = 0. ;
// 	for(size_t i = 0 ; i < inclusions.size() ; i++)
// 	{
// 		switch(reference)
// 		{
// 			case CIRCLE:
// 				inclusions[i]->setBehaviour(agg) ;
// 				F.addFeature(&sample, inclusions[i]) ;
// 				placed_area += inclusions[i]->area() ;
// 				break ;
// 			case ELLIPSE:
// 				ellinc[i]->setBehaviour(agg) ;
// 				F.addFeature(&sample, ellinc[i]) ;
// 				placed_area += ellinc[i]->area() ;
// 				break ;
// 			case TRIANGLE:
// 				triinc[i]->setBehaviour(agg) ;
// 				F.addFeature(&sample, triinc[i]) ;
// 				placed_area += triinc[i]->area() ;
// 				break ;
// 			case RECTANGLE:
// 				recinc[i]->setBehaviour(agg) ;
// 				F.addFeature(&sample, recinc[i]) ;
// 				placed_area += recinc[i]->area() ;
// 				break ;
// 		}
// 	}
    
	switch(reference)
	{
		case CIRCLE:
//			zones = generateExpansiveZonesHomogeneously(100, 3, inclusions, F) ;
//			zones = generateExpansiveZonesAtCenter(inclusions, F) ;
			break ;
		case ELLIPSE:
// 			zones = generateExpansiveZonesHomogeneously(100, 3, ellinc, F) ;
//			zones = generateExpansiveZonesAtCenter(ellinc, F) ;
			break ;
		case TRIANGLE:
//			zones = generateExpansiveZonesHomogeneously(100, 3, triinc, F) ;
// 			zones = generateExpansiveZonesAtCenter(triinc, F) ;
			break ;
		case RECTANGLE:
//			zones = generateExpansiveZonesHomogeneously(100, 3, recinc, F) ;
// 			zones = generateExpansiveZonesAtCenter(recinc, F) ;
			break ;
	}*/


	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -1e6)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BEFORE)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BEFORE)) ;

	int nSampling = 100 ;
	
	F.setSamplingNumber(nSampling) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;

//	step(reference, nSampling, 0,0) ;

	F.step() ;
	
	Vector x = F.getDisplacements() ;
	std::cout << x.max() << std::endl ;
	
	return 0 ;
}

