// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

#include <fstream>

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

Sample sample(NULL, 0.07, 0.07, 0, 0) ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<double> apparent_extension ;
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


std::string itoa(int value, int base) {

	enum { kMaxDigits = 35 };
	std::string buf;
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.

	// check that the base if valid
	if (base < 2 || base > 16) return buf;

	int quotient = value;
	
	
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += '-';
	
	std::reverse( buf.begin(), buf.end() );
	
	return buf;
	
}

void fastForward (int steps, int nstepstot)
{
	for(size_t i = 0 ; i < steps ; i++)
	{
		double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/(double)nstepstot ;

		double reactedArea = 0 ;
		
		Inclusion * current = NULL ;
		if(!zones.empty())
			current = zones[0].second ;
		double current_area = 0 ;
		int current_number = 0 ;
		int stopped_reaction = 0 ;
		for(size_t z = 0 ; z < zones.size() ; z++)
		{
			
			zones[z].first->setRadius(zones[z].first->getGeometry()->getRadius()+delta_r) ;	
			// 		zones[z].first->reset() ;
			if(zones[z].second == current)
			{
				current_area += zones[z].first->area() ;
				current_number++ ;
			}
			else
			{
				if(current_area/zones[z-1].second->area() > 0.03)
				{
					stopped_reaction++ ;
					for(size_t m = 0 ; m < current_number ; m++)
					{
						reactedArea -= zones[z-1-m].first->area() ;
						zones[z-1-m].first->setRadius(zones[z].first->getGeometry()->getRadius()-delta_r) ;
						reactedArea += zones[z-1-m].first->area() ;
					}
				}
				current_area = zones[z].first->area() ;
				current_number = 1 ;
				current = zones[z].second ;
			}
			reactedArea += zones[z].first->area() ;
		}
		
		std::cout << "reacted Area : " << reactedArea << ", reaction stopped in "<< stopped_reaction << " aggs."<< std::endl ;
	}
}

void setBC()
{
	triangles = featureTree->getTriangles() ;
	
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			if (triangles[k]->getBoundingPoint(c).y < -.499*sample.height() && triangles[k]->getBoundingPoint(c).x < -.499*sample.width())
			{
				featureTree->getAssembly()->setPoint( 0,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			else if (triangles[k]->getBoundingPoint(c).x < -.499*sample.width() && triangles[k]->getBoundingPoint(c).y > .499*sample.height())
			{
				featureTree->getAssembly()->setPointAlong( XI,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			else if(triangles[k]->getBoundingPoint(c).x > .499*sample.width() && triangles[k]->getBoundingPoint(c).y < -.499*sample.height())
			{
				featureTree->getAssembly()->setPointAlong( ETA,0 ,triangles[k]->getBoundingPoint(c).id) ;

//				for(size_t l = c+1 ;  l < triangles[k]->getBoundingPoints().size() ; l++)
//				{
//					if(triangles[k]->getBoundingPoint(l).x > .499*sample.width())
//					{
//						double d = std::abs(triangles[k]->getBoundingPoint(l).y-triangles[k]->getBoundingPoint(c).y) ;
// 						featureTree->getAssembly()->setPointAlong(ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
// 						featureTree->getAssembly()->setPointAlong(ETA, 0, triangles[k]->getBoundingPoint(l).id) ;
//						featureTree->getAssembly()->setForceOn( XI, -stress*d*.25 ,triangles[k]->getBoundingPoint(c).id) ;
//						featureTree->getAssembly()->setForceOn( XI, -stress*d*.25 ,triangles[k]->getBoundingPoint(l).id) ;
//						break ;
//					}
//				}
//				break ;
			}
		}

	}

}


void step()
{

	int nsteps = 20;
	int nstepstot = 40;
	int maxtries = 200 ;
	int tries = 0 ;
	
// 	fastForward(4, 10) ;
	
	for(size_t i = 0 ; i < nsteps ; i++)
	{
		std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;
		tries = !(nsteps < maxtries) ;
		bool go_on = true ;
		while(go_on && tries < maxtries)
		{
			setBC() ;
			featureTree->step(timepos) ;
			go_on = featureTree->solverConverged() &&  (featureTree->meshChanged() || featureTree->enrichmentChanged());
			std::cout << "." << std::flush ;
// 			timepos-= 0.0001 ;
			
			tries++ ;
		}
		std::cout << " " << tries << " tries." << std::endl ;
		
		if(featureTree->solverConverged())
		{
			cracked_volume.push_back(featureTree->crackedVolume) ;
			damaged_volume.push_back(featureTree->damagedVolume) ;
		}
	// 		
	// 		
		timepos+= 0.0001 ;
	
	
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		
	// 	sigma = F.strainFromDisplacements() ;
	// 	epsilon = F.stressFromDisplacements() ;
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
	/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
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
					Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l)) ;
					agl = (vm0[0] <= 0 && vm0[1] <= 0)*180. ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl ;
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
		std::string filename("triangles") ;
		filename.append(itoa(totit++, 10)) ;
		std::cout << filename << std::endl ;
		std::fstream outfile  ;
		outfile.open(filename.c_str(), std::ios::out) ;
		
		outfile << "TRIANGLES" << std::endl ;
		outfile << triangles.size() << std::endl ;
		outfile << 3 << std::endl ;
		outfile << 8 << std::endl ;
		
		for(size_t j = 0 ; j < triangles.size() ;j++)
		{
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile << triangles[j]->getBoundingPoint(l).x << " " << triangles[j]->getBoundingPoint(l).y << " ";
			}

			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  epsilon11[j*3+l] << " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  epsilon22[j*3+l] << " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<   epsilon12[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma11[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma22[j*3+l]<< " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  sigma12[j*3+l] << " ";
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile << vonMises[j*3+l]<< " " ;
			}
			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
			{
				outfile <<  triangles[j]->getBehaviour()->getTensor(Point(.3, .3))[0][0] << " ";
			}
			outfile << "\n" ;
		}
		
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
		
		std::cout << "average sigma11 (no gel): " << avg_s_xx_nogel/nogel_area << std::endl ;
		std::cout << "average sigma22 (no gel): " << avg_s_yy_nogel/nogel_area << std::endl ;
		std::cout << "average sigma12 (no gel): " << avg_s_xy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon11 (no gel): " << avg_e_xx_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon22 (no gel): " << avg_e_yy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon12 (no gel): " << avg_e_xy_nogel/nogel_area << std::endl ;
		
		std::cout << "apparent extension " << e_xx_max-e_xx_min << std::endl ;
		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
		if (tries < maxtries)
		{
			double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/(double)nstepstot ;
			if(!featureTree->solverConverged())
				delta_r *= .01 ;
			double reactedArea = 0 ;
			
			Inclusion * current = NULL ;
			if(!zones.empty())
				current = zones[0].second ;
			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;
			for(size_t z = 0 ; z < zones.size() ; z++)
			{
				
				zones[z].first->setRadius(zones[z].first->getGeometry()->getRadius()+delta_r) ;	
		// 		zones[z].first->reset() ;
				if(zones[z].second == current)
				{
					current_area += zones[z].first->area() ;
					current_number++ ;
				}
				else
				{
					if(current_area/zones[z-1].second->area() > 0.03)
					{
						stopped_reaction++ ;
						for(size_t m = 0 ; m < current_number ; m++)
						{
							reactedArea -= zones[z-1-m].first->area() ;
							zones[z-1-m].first->setRadius(zones[z].first->getGeometry()->getRadius()-delta_r) ;
							reactedArea += zones[z-1-m].first->area() ;
						}
					}
					current_area = zones[z].first->area() ;
					current_number = 1 ;
					current = zones[z].second ;
				}
				reactedArea += zones[z].first->area() ;
			}
			
			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in "<< stopped_reaction << " aggs."<< std::endl ;

		
			if(featureTree->solverConverged())
			{
				expansion_reaction.push_back(std::make_pair(reactedArea/placed_area, avg_e_xx/area)) ;
				expansion_stress_xx.push_back(std::make_pair((avg_e_xx_nogel)/(nogel_area), (avg_s_xx_nogel)/(nogel_area))) ;
				expansion_stress_yy.push_back(std::make_pair((avg_e_yy_nogel)/(nogel_area), (avg_s_yy_nogel)/(nogel_area))) ;
				apparent_extension.push_back(e_xx_max-e_xx_min) ;
			}
			
			if (tries >= maxtries)
				break ;
		}
	for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
		std::cout << expansion_reaction[i].first << "   " 
		<< expansion_reaction[i].second << "   " 
		<< expansion_stress_xx[i].first << "   " 
		<< expansion_stress_xx[i].second << "   " 
		<< expansion_stress_yy[i].first << "   " 
		<< expansion_stress_yy[i].second << "   " 
		<< apparent_extension[i]  << "   " 
		<< cracked_volume[i]  << "   " 
		<< damaged_volume[i]  << "   " 
		<< std::endl ;

	}

	for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
	std::cout << expansion_reaction[i].first << "   " 
	<< expansion_reaction[i].second << "   " 
	<< expansion_stress_xx[i].first << "   " 
	<< expansion_stress_xx[i].second << "   " 
	<< expansion_stress_yy[i].first << "   " 
	<< expansion_stress_yy[i].second << "   " 
	<< apparent_extension[i]  << "   " 
	<< cracked_volume[i]  << "   " 
	<< damaged_volume[i]  << "   " 
	<< std::endl ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZonesHomogeneously(int n, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	double nu_incompressible = .499997 ;
	
	double E = percent*E_csh ;
	double nu = nu_csh ; //nu_incompressible ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	Vector a(double(0), 3) ;
	a[0] = 0.5 ;
	a[1] = 0.5 ;
	a[2] = 0.00 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		Point pos(((double)rand()/RAND_MAX-.5)*(sample.width()-radius*60),((double)rand()/RAND_MAX-.5)*(sample.height()-radius*60)) ;
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
			zonesToPlace.push_back(new ExpansiveZone(incs[i], radius, pos.x, pos.y, m0, a)) ;
		else
			i-- ;
	}
	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			if(dist(zonesToPlace[i]->getCenter(), incs[j]->getCenter()) < incs[j]->getRadius()-radius*60)
			{
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
	}
	
	int count = 0 ;
	for(std::map<Inclusion *, int>::iterator i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
		std::cout << aggregateArea << "  " << count << std::endl ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZones(int n, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	double nu_incompressible = .499997 ;
	
	double E = percent*E_csh ;
	double nu = nu_incompressible ;

	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;
	double radius = 0.000005 ;
	for(size_t i = 0 ; i < incs.size() ; i++)
	{
		if( incs[i]->getRadius() < .001)
		{
			aggregateArea += incs[i]->area() ;
			for(int j = 0 ; j < n ; j++)
			{
				
				Point pos((2.*rand()/RAND_MAX-1.),(2.*rand()/RAND_MAX-1.)) ;
				pos /= pos.norm() ;
				pos *= (2.*rand()/RAND_MAX-1.)*(incs[i]->getRadius() - 0.00003) ;
				Point center = incs[i]->getCenter()+pos ; 
				
				bool alone  = true ;
				
				for(size_t k = 0 ; k < ret.size() ; k++ )
				{
					if (squareDist(center, ret[k].first->Circle::getCenter()) < (radius*60.+radius*60.)*(radius*60.+radius*60.))
					{
						alone = false ;
						break ;
					}
				}
				if (alone)
				{
					Vector a(double(0), 3) ;
					a[0] = 0.5 ;
					a[1] = 0.5 ;
					a[2] = 0.00 ;
					
					ExpansiveZone * z = new ExpansiveZone(incs[i], radius, center.x, center.y, m0, a) ;
					ret.push_back(std::make_pair(z, incs[i])) ;
				}
			}
		}
	}

	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i].first->setRadius(radius) ;
		F.addFeature(ret[i].second, ret[i].first) ;
	}
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > generateInclusionsAndPores(size_t n, double fraction, Matrix * tensor, Feature * father, FeatureTree * F)
{
// 	srandom(time(NULL)) ;
	size_t nombre_de_pores = static_cast<size_t>(round(n*fraction)) ;
	size_t nombre_d_inclusions = static_cast<size_t>(round(n*(1. - fraction))) ;
	
	std::pair<std::vector<Inclusion * >, std::vector<Pore * > > ret ;
	ret.first = std::vector<Inclusion * >() ;
	ret.second = std::vector<Pore * >() ;
	double v = 0 ;
	std::vector<Circle *> cercles ;
	for(size_t j =0 ; j < n ; j++)
	{
		
		double radius = .0005 + .0025*rand()/RAND_MAX ;
		
		Point center = Point(
		                      (2.*rand()/RAND_MAX-1.)*(.08-2.*radius-0.00001),
		                      (2.*rand()/RAND_MAX-1.)*(.02-2.*radius-0.00001)
		                    ); 
		bool alone  = true ;
		
		for(size_t k = 0 ; k < cercles.size() ; k++ )
		{
			if (squareDist(center, cercles[k]->getCenter()) < (radius+cercles[k]->getRadius()+0.00001)*(radius+cercles[k]->getRadius()+0.00001))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			cercles.push_back(new Circle(radius, center)) ;
			v+= M_PI*radius*radius ;
		}
		else
			j-- ;
		
	}
	for(size_t j =0 ; j < nombre_d_inclusions ; j++)
	{
		Vector imp(double(0),3) ;
		imp[0] = 0.01 ;
		imp[1] = 0.01 ;
		Inclusion * temp = new Inclusion(cercles[j]->getRadius(), cercles[j]->getCenter()) ;
		ret.first.push_back(temp) ;
// 		(*ret.first.rbegin())->setBehaviour(new StiffnessAndFracture(*tensor, new MohrCoulomb(1000000, -10000000))) ;
		(*ret.first.rbegin())->setBehaviour(new WeibullDistributedStiffness(*tensor, 1000000)) ;
		F->addFeature(father, temp) ;
	}
	
	for(size_t j =0 ; j < nombre_de_pores ; j++)
	{
		Pore * temp = new Pore(cercles[j+nombre_d_inclusions]->getRadius(), cercles[j+nombre_d_inclusions]->getCenter()) ;
		ret.second.push_back(temp) ;
		F->addFeature(father, temp) ;
	}
	
	for(size_t k = 0 ; k < cercles.size() ; k++ )
	{
		delete cercles[k] ;
	}
	
	std::cout << "initial aggregate volume was : " << v << std::endl ;
	aggregateArea = v ;
	return ret ;
}

void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
	int i;
	double f, p, q, t;
	if( s == 0 ) {
                // achromatic (grey)
		*r = *g = *b = v;
		return;
	}
	h /= 60.;                        // sector 0 to 5
	i = (int)floor( h );
	f = h - i;                      // factorial part of h
	p = v * ( 1. - s );
	q = v * ( 1. - s * f );
	t = v * ( 1. - s * ( 1. - f ) );
	switch( i ) {
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

void init(void) 
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);   // Enables Smooth Shading
	glEnable(GL_LINE_SMOOTH) ;
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) ;
	
// 	glPointSize(std::max(0.4*(double)width()/(double)columns, 1.));
	glClearColor(0.0f,0.0f,0.0f,0.0f);                                      // Black Background

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE); 
}

void Menu(int selection)
{
	switch (selection)
	{
	case ID_NEXT:
		{
			step() ;
			dlist = false ;
			break ;
		}
	case ID_NEXT_TIME:
		{
			timepos +=0.0001 ;
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

	case ID_QUIT : exit(0) ;
		
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
			
			for(size_t k = 0 ; k < triangles.size() ; k++)
			{
/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
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

void reshape(int w, int h)
{
	if (h==0)                                          // Prevent A Divide By Zero By
		h=1;                                           // Making Height Equal One
	
	glViewport(0, 0, (int)(w*factor), (int)(h*factor));
	gluPerspective((double)h/(double)w,1.,1.f,45.0f);
}

void Display(void)
{
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
	glFlush();
	glutSwapBuffers();
	glMatrixMode(GL_PROJECTION) ;
	glLoadIdentity() ;
	glOrtho(-4.5/factor, 4.5/factor, -4.5/factor, 4.5/factor, -4.5, 4.5);
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
	
	if(!dlist)
	{
		
		
		glNewList( DISPLAY_LIST_DISPLACEMENT,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j] && !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(0).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(0).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min)))*300., 1., 1. ) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double sigma11_min = sigma11.min() ;
		double sigma11_max = sigma11.max() ;
		glNewList(  DISPLAY_LIST_STRAIN_XX,  GL_COMPILE ) ;
		
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()+k]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double vonMises_max = vonMises.max() ;
		double vonMises_min = vonMises.min() ;
		
		glNewList(  DISPLAY_LIST_VON_MISES,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()+k]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		
		double angle_max = angle.max() ;
		double angle_min = angle.min() ;
		glNewList(  DISPLAY_LIST_ANGLE,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()]-angle_min)/(angle_max-angle_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()+k]-angle_min)/(angle_max-angle_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double sigma22_min = sigma22.min() ;
		double sigma22_max = sigma22.max() ;
		
		glNewList(  DISPLAY_LIST_STRAIN_YY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()+k]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		double sigma12_min = sigma12.min() ;
		double sigma12_max = sigma12.max() ;
		glNewList(  DISPLAY_LIST_STRAIN_XY,  GL_COMPILE ) ;
		
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()+k]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_STIFFNESS,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(0)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(k)) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_STIFFNESS_DARK,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(0)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., .2) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(k)) ;
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., .2) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		double epsilon11_min = epsilon11.min() ;
		double epsilon11_max = epsilon11.max() ;
		glNewList(  DISPLAY_LIST_STRESS_XX,  GL_COMPILE ) ;
		
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()+k]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		
		double epsilon22_min = epsilon22.min() ;
		double epsilon22_max =  epsilon22.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_YY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()+k]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		
		glEndList() ;
		
		double epsilon12_min = epsilon12.min() ;
		double epsilon12_max =  epsilon12.max() ;
		
		glNewList(  DISPLAY_LIST_STRESS_XY,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->first->id*2]; 
				double vy = x[triangles[j]->first->id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(0).x + vx) , double(triangles[j]->getBoundingPoint(0).y + vy) );
				
				for(size_t k = 1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()+k]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x + vx) ,  double(triangles[j]->getBoundingPoint(k).y + vy) );
					
				}
				glEnd() ;
			}
		}
		glEndList() ;
		
		glNewList(  DISPLAY_LIST_ENRICHMENT,  GL_COMPILE ) ;
		glBegin(GL_TRIANGLES);
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				int enrichment = triangles[j]->getEnrichmentFunctions().size() ;
				//HSVtoRGB( &c1, &c2, &c3, 180. + 180.*(sigma12[j]-sigma12.min())/(sigma12.max()-sigma12.min()), 1., 1. ) 
				
				
				if(enrichment)
				{
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*enrichment/20., 1., 1.) ;
					glColor3f(c1, c2, c3) ;
					double vx = x[triangles[j]->first->id*2]; 
					double vy = x[triangles[j]->first->id*2+1]; 
					
					glVertex2f( double(triangles[j]->first->x + vx) ,
					            double(triangles[j]->first->y + vy) );
					
					vx = x[triangles[j]->second->id*2];
					vy = x[triangles[j]->second->id*2+1]; 
					
					glVertex2f( double(triangles[j]->second->x + vx) ,
					            double(triangles[j]->second->y + vy) );
					
					
					vx = x[triangles[j]->third->id*2]; 
					vy = x[triangles[j]->third->id*2+1]; 
					
					
					glVertex2f( double(triangles[j]->third->x + vx) ,
					            double(triangles[j]->third->y + vy) );
				}
			}
		}
		glEnd();
		glEndList() ;
		
		
		glNewList(  DISPLAY_LIST_ELEMENTS,  GL_COMPILE ) ;
		glColor3f(1, 1, 1) ;
		for(unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				if(triangles[j]->getBehaviour()->fractured())
					glColor3f(1, 0, 0) ;
				else
					glColor3f(1, 1, 1) ;
				
				glBegin(GL_LINE_LOOP);
				for(size_t k = 0 ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					double vx = x[triangles[j]->getBoundingPoint(k).id*2]; 
					double vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					glVertex2f( double(triangles[j]->getBoundingPoint(k).x+vx) ,  double(triangles[j]->getBoundingPoint(k).y+vy) );
					
				}
				glEnd();
			}
			
			glColor3f(1, 1, 1) ;
		}
		glEndList() ;
		
		
		glNewList(  DISPLAY_LIST_CRACK,  GL_COMPILE ) ;
		glLineWidth(4) ;
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
		glLineWidth(1) ;
		glEndList() ;
		
		dlist = true ;
		glCallList(current_list) ;
	}
	else
	{
		//glCallList(DISPLAY_LIST_DISPLACEMENT) ;
		//glCallList(DISPLAY_LIST_STRAIN) ;
		double c1, c2, c3 = 0;
		HSVtoRGB( &c1, &c2, &c3, 180. + 0, 1., 1.) ;
// 		glBegin(GL_LINE) ;
// 		glVertex2f(3.5 ,
// 		           3. );
// 		glVertex2f(3.5 ,
// 		           -3. );
// 		glEnd() ;
		
		if(current_list != DISPLAY_LIST_ENRICHMENT)
			glCallList(current_list) ;
		if(current_list == DISPLAY_LIST_ENRICHMENT)
		{
			glCallList(DISPLAY_LIST_STIFFNESS_DARK) ;
			glCallList(current_list) ;
			
		}
		
		glCallList(DISPLAY_LIST_CRACK) ;
// 		if(current_list == DISPLAY_LIST_ELEMENTS)
// 			glCallList(DISPLAY_LIST_CRACK) ;
		
		glColor3f(1, 1, 1) ;
		
		
	}
	glColor3f(1, 0, 0) ;
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char *argv[])
{

	percent = atof(argv[1]) ;

	Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ; 
	
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ; 
	
// 	Sample reinforcement0(NULL, 8,.15,0,.5) ;
// 	reinforcement0.setBehaviour(new Stiffness(m0*5)) ;
// 	
// 	Sample reinforcement1(NULL, 8,.15,0,-.5) ;
// 	reinforcement1.setBehaviour(new Stiffness(m0*5)) ;
	
	FeatureTree F(&sample) ;
	featureTree = &F ;


	double itzSize = 0.00005;
	int inclusionNumber = 8000 ;
// 	int inclusionNumber = 4096 ;
// 	std::vector<Inclusion *> inclusions = GranuloBolome(4.79263e-07, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);
// 
// // 	if(inclusionNumber)
// // 		itzSize = inclusions[inclusions.size()/4]->getRadius() ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		delete inclusions[i] ;

	double masseInitiale = .00000743;
	double densite = 1.;
	std::vector<Inclusion *> inclusions = GranuloBolome(masseInitiale, densite, BOLOME_A)(.008, .0001, inclusionNumber, itzSize);
// 	std::vector<Inclusion *> inclusions = GranuloBolome(0.0000416, 1, BOLOME_D)(.0025, .1, inclusionNumber, itzSize);

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;


	inclusions.clear() ;
	for(size_t i = 0; i < feats.size() ; i++)
		inclusions.push_back(static_cast<Inclusion *>(feats[i])) ;

	
	int nAgg = 1 ;
	feats=placement(sample.getPrimitive(), feats, &nAgg, 6400);
	double volume = 0 ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		volume += feats[i]->area() ;
	if(!feats.empty())
		std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
		<< ", smallest r =" << feats.back()->getRadius() 
		<< ", filling = " << volume/sample.area()*100.<< "%"<< std::endl ; 

	sample.setBehaviour(new WeibullDistributedStiffness(m0_paste, 13500000)) ;

	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
		WeibullDistributedStiffness * stiff = new WeibullDistributedStiffness(m0_agg,57000000) ;
// 		stiff->variability = .5 ;
		inclusions[i]->setBehaviour(stiff) ;
		F.addFeature(&sample,inclusions[i]) ;
		placed_area += inclusions[i]->area() ;
	}
	inclusions.erase(inclusions.begin()+feats.size(), inclusions.end()) ;
	
	if(!inclusions.empty())
	{
		std::cout << "largest inclusion with r = " << (*inclusions.begin())->getRadius() << std::endl ;
		std::cout << "smallest inclusion with r = " << (*inclusions.rbegin())->getRadius() << std::endl ;
		std::cout << "placed area = " <<  placed_area << std::endl ;
	}
	Circle cercle(.5, 0,0) ;

	zones = generateExpansiveZones(20, inclusions, F) ;

	F.sample(1200) ;

	F.setOrder(LINEAR) ;

	F.generateElements() ;
// 	
	step() ;
	
	glutInit(&argc, argv) ;	
	glutInitDisplayMode(GLUT_RGBA) ;
	glutInitWindowSize(600, 600) ;
	glutReshapeFunc(reshape) ;
	glutCreateWindow("coucou !") ;
	
	int submenu = glutCreateMenu(Menu) ;
	
	glutAddMenuEntry(" Displacements ", ID_DISP);
	glutAddMenuEntry(" Strain (s) xx ", ID_STRAIN_XX);
	glutAddMenuEntry(" Strain (s) yy ", ID_STRAIN_YY);
	glutAddMenuEntry(" Strain (s) xy ", ID_STRAIN_XY);
	glutAddMenuEntry(" Stress (e) xx ", ID_STRESS_XX);
	glutAddMenuEntry(" Stress (e) yy ", ID_STRESS_YY);
	glutAddMenuEntry(" Stress (e) xy ", ID_STRESS_XY);
	glutAddMenuEntry(" Elements      ", ID_ELEM);
	glutAddMenuEntry(" Stiffness     ", ID_STIFNESS);
	glutAddMenuEntry(" Von Mises     ", ID_VON_MISES);
	glutAddMenuEntry(" Princ. angle  ", ID_ANGLE);
	glutAddMenuEntry(" Enrichment    ", ID_ENRICHMENT);
	
	glutCreateMenu(Menu) ;

 	glutAddMenuEntry(" Step          ", ID_NEXT);
	glutAddMenuEntry(" Step time     ", ID_NEXT_TIME);
	glutAddMenuEntry(" Zoom in       ", ID_ZOOM);
	glutAddMenuEntry(" Zoom out      ", ID_UNZOOM);
	glutAddMenuEntry(" Amplify       ", ID_AMPLIFY);
	glutAddMenuEntry(" Deamplify     ", ID_DEAMPLIFY);
	glutAddSubMenu(  " Display       ", submenu);
	glutAddMenuEntry(" Quit          ", ID_QUIT) ;
	
	
	glutAttachMenu(GLUT_RIGHT_BUTTON) ;
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	
	glutDisplayFunc(Display) ;
	glutMainLoop() ;
	
// 	delete dt ;

	for(size_t i = 0 ; i < inclusions.size() ; i++)
		delete inclusions[i] ;
	
	return 0 ;
}
