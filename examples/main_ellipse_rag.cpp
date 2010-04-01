// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
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
#include "../physics/kelvinvoight.h"
#include "../physics/fracturecriteria/vonmises.h"
#include "../physics/spatially_distributed_stiffness.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../features/expansiveZone.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

#include <iostream>  // I/O 
#include <fstream>   // file I/O
#include <iomanip>   // format manipulation

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
#define ID_FRAC_CRIT 23

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
using namespace std;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;
std::vector<BranchedCrack *> crack ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
double percent = 0.70 ;
double placed_area = 0 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, EllipsoidalInclusion *> > zones ;

double width = 0.07;
double height = 0.07;
Sample sample(NULL, width, height, 0.035, 0.035) ;
	
std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;
std::vector<std::pair<double, double> > young_modulus ;


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
Vector fracCrit(0) ;

Vector g_count(0) ;

double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;
double E_stiff = E_agg*10. ;//stiffer
double E_soft = E_agg/10.; //stiffest

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

int totit = 5 ;

std::vector<double> energy ;


std::string itoa(int value, int base) {

	enum { kMaxDigits = 35 };
	std::string buf;
	buf.reserve( kMaxDigits ); 
	if (base < 2 || base > 16) return buf;
	int quotient = value;
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	if ( value < 0 && base == 10) buf += '-';	
	std::reverse( buf.begin(), buf.end() );
	return buf;
	
}



void setBC()
{
	triangles = featureTree->getTriangles() ;
	
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			if (triangles[k]->getBoundingPoint(c).y < 0.001*sample.height() && triangles[k]->getBoundingPoint(c).x < 0.001*sample.width())
			{
				featureTree->getAssembly()->setPoint( 0,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			else if (triangles[k]->getBoundingPoint(c).x < 0.001*sample.width())
			{
				featureTree->getAssembly()->setPointAlong( XI,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			else if(triangles[k]->getBoundingPoint(c).y < 0.001*sample.height())
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
 	/*
 	bool cracks_did_not_touch = true;
 	size_t max_growth_steps = 1;
 	size_t countit = 0;	
// 	
 	while ( (cracks_did_not_touch) && (countit < max_growth_steps) )
 	{
 		countit++;
 		std::cout << "\r iteration " << countit << "/" << max_growth_steps << std::flush ;
 		setBC() ;
//       
 		int limit = 0 ;
 		while(!featureTree->step(timepos) && limit < 50)//as long as we can update the features
 		{
 			std::cout << "." << std::flush ;
 // 			timepos-= 0.0001 ;
 			setBC() ;
 			limit++ ;
 	  // check if the two cracks did not touch
 			cracks_did_not_touch = true;
 			for(size_t j = 0 ; j < crack.size() ; j++)
 			{
 				for(size_t k = j+1 ; k < crack.size() ; k++)
 				{
 		  
 					if (static_cast<SegmentedLine *>(crack[j])->intersects(static_cast<SegmentedLine *>(crack[k])))
 					{
 						cracks_did_not_touch = false;
 						break;
 					}
 		  //Circle headj(radius, crack[j]->getHead());
 		  
 		  //	      head0.intersects(crack[1]);
 				}
			}
 		}
 	
//   // Prints the crack geo to a file for each crack
// 		if (cracks_did_not_touch == false) // if cracks touched
// 		{
// 			std::cout << "** Cracks touched exporting file **" << endl;	
// 				// Print the state of the cracks to a file  
// 				std::string filename = "crackGeo.txt";
// 				fstream filestr;
// 				filestr.open (filename.c_str(), fstream::in | fstream::out | fstream::app);
// 				filestr << "Crack vertices" << std::endl;
// 				filestr << "x" << " " << "y" << std::endl ; 
// 		}
   
 
// 		timepos+= 0.0001 ;
 		double da = 0 ;
 
		triangles = featureTree->getTriangles() ;
 		x.resize(featureTree->getDisplacements().size()) ;
 		x = featureTree->getDisplacements() ;
 
 		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
 		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
 	
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
 		fracCrit.resize(sigma.size()/3) ;
 		g_count.resize(sigma.size()/3) ;
 		std::cout << "unknowns :" << x.size() << std::endl ;
 		
 		if(crack.size() > 0)
 			tris__ = crack[0]->getElements(featureTree->get2DMesh()) ;
 		
 		for(size_t k = 1 ; k < crack.size() ; k++)
 		{
 			std::vector<DelaunayTriangle *> temp = crack[k]->getElements(featureTree->get2DMesh()) ;
 			if(tris__.empty())
 				tris__ = temp ;
 			else if(!temp.empty())
 				tris__.insert(tris__.end(), temp.begin(), temp.end() ) ;
 		}
 		cracked.clear() ;
 		
 		int npoints = triangles[0]->getBoundingPoints().size() ;
 	
 		double area = 0 ;
 		double avg_e_xx = 0;
 		double avg_e_yy = 0;
 		double avg_e_xy = 0;
 		double avg_s_xx = 0;
 		double avg_s_yy = 0;
 		double avg_s_xy = 0;
 		double e_xx = 0 ;
 		double ex_count = 0 ;
 		double enr = 0 ;
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
 		
 		
 		
 			if(!in && !triangles[k]->getBehaviour()->fractured() && triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
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
 					if(triangles[k]->getBoundingPoint(p).x > 0.0799)
 					{
 						e_xx+=x[triangles[k]->getBoundingPoint(p).id*2] ;
 						ex_count++ ;
 					}
 				}
 				area += triangles[k]->area() ;
 				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
 				{
 					if(triangles[k]->getBehaviour()->param[0][0] > E_max)
 						E_max = triangles[k]->getBehaviour()->param[0][0] ;
 					if(triangles[k]->getBehaviour()->param[0][0] < E_min)
 						E_min = triangles[k]->getBehaviour()->param[0][0] ;
 				}
 				
 				g_count[k*npoints]++ ;
 				g_count[k*npoints+1]++ ;
 				g_count[k*npoints+2]++ ;
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
 					g_count[k*npoints+3]++ ;
 					g_count[k*npoints+4]++ ;
 					g_count[k*npoints+5]++ ;
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
 				
 				for(size_t l = 0 ; l < npoints ; l++)
 				{
 					Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
 					vonMises[k*npoints+l]  += sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
 	
 					double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l)) ;
 					angle[k*npoints+l]  += agl ;
 					if(triangles[k]->getBehaviour()->getFractureCriterion())
 					{
 						fracCrit[k*npoints+l] = triangles[k]->getBehaviour()->getFractureCriterion()->grade(triangles[k]->getState()) ;
 						enr += triangles[k]->getState().elasticEnergy() ;
 					}
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
 	
 			}
 			else
 			{
 				sigma11[k*npoints] += 0 ;
 				sigma22[k*npoints] += 0 ;
 				sigma12[k*npoints] += 0 ;
 				sigma11[k*npoints+1] += 0 ;
 				sigma22[k*npoints+1] += 0 ;
 				sigma12[k*npoints+1] += 0 ;
 				sigma11[k*npoints+2] += 0 ;
 				sigma22[k*npoints+2] += 0 ;
 				sigma12[k*npoints+2] += 0 ;
 				
 				if(npoints >3)
 				{
 					sigma11[k*npoints+3] += 0 ;
 					sigma22[k*npoints+3] += 0 ;
 					sigma12[k*npoints+3] += 0 ;
 					sigma11[k*npoints+4] += 0 ;
 					sigma22[k*npoints+4] += 0 ;
 					sigma12[k*npoints+4] += 0 ;
 					sigma11[k*npoints+5] += 0 ;
 					sigma22[k*npoints+5] += 0 ;
 					sigma12[k*npoints+5] += 0 ;
 				}
 				
 				epsilon11[k*npoints] += 0 ;
				epsilon22[k*npoints] += 0 ;
 				epsilon12[k*npoints] += 0 ;
 				epsilon11[k*npoints+1] += 0 ;
 				epsilon22[k*npoints+1] += 0 ;
 				epsilon12[k*npoints+1] += 0 ;
 				epsilon11[k*npoints+2] += 0 ;
 				epsilon22[k*npoints+2] += 0 ;
 				epsilon12[k*npoints+2] += 0 ;
 				
 				if(npoints > 3)
 				{
 					epsilon11[k*npoints+3] += 0 ;
 					epsilon22[k*npoints+3] += 0 ;
 					epsilon12[k*npoints+3] += 0 ;
 					epsilon11[k*npoints+4] += 0 ;
 					epsilon22[k*npoints+4] += 0 ;
 					epsilon12[k*npoints+4] += 0 ;
 					epsilon11[k*npoints+5] += 0 ;
 					epsilon22[k*npoints+5] += 0 ;
 					epsilon12[k*npoints+5] += 0 ;
 				}  
 				
 				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
 				{
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  += 0 ;
 					angle[k*triangles[k]->getBoundingPoints().size()+l]  += 0 ;
 				}
 			}
 		}
// 	
// 		std::string filename("triangles") ;
// 		filename.append(itoa(totit++, 10)) ;
// 		std::cout << filename << std::endl ;
// 		std::fstream outfile  ;
// 		outfile.open(filename.c_str(), std::ios::out) ;
// 		
// 		outfile << "TRIANGLES" << std::endl ;
// 		outfile << triangles.size() << std::endl ;
// 		outfile << 3 << std::endl ;
// 		outfile << 8 << std::endl ;
// 		
// 		for(size_t j = 0 ; j < triangles.size() ;j++)
// 		{
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile << triangles[j]->getBoundingPoint(l).x << " " << triangles[j]->getBoundingPoint(l).y << " ";
// 			}
// 
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  epsilon11[j*3+l] << " ";
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  epsilon22[j*3+l] << " " ;
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<   epsilon12[j*3+l]<< " " ;
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  sigma11[j*3+l]<< " " ;
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  sigma22[j*3+l]<< " ";
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  sigma12[j*3+l] << " ";
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
//				outfile << vonMises[j*3+l]<< " " ;
// 			}
// 			for(size_t l = 0 ; l < triangles[j]->getBoundingPoints().size() ; l++)
// 			{
// 				outfile <<  triangles[j]->getBehaviour()->getTensor(Point(.3, .3))[0][0] << " ";
// 			}
// 			outfile << "\n" ;
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
 		
 		std::cout << "average sigma11 : " << avg_s_xx/area << std::endl ;
 		std::cout << "average sigma22 : " << avg_s_yy/area << std::endl ;
 		std::cout << "average sigma12 : " << avg_s_xy/area << std::endl ;
 		std::cout << "average epsilon11 : " << avg_e_xx/area << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
 		std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;
 			
 		std::cout << "energy index :" << enr << std::endl ;
 		energy.push_back(enr) ;
 
 		if(limit < 2)
 			break ;
 	}
 	for(size_t i = 0 ; i < energy.size() ; i++)
 		std::cout << energy[i] << std::endl ;*/
 }

void stepOLD()
{
/*
	int nsteps = 5;
	int nstepstot = 10;
	int maxtries = 20 ;
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
		double avg_Emodulus_xx = 0 ;
		double avg_Emodulus_yy = 0 ;
		double avg_Emodulus_xy = 0 ;
		
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
	//		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;
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
					if(triangles[k]->getBoundingPoint(p).x > sample.width()*.9999)
					{
						if(e_xx_max < x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_max=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).x < sample.width()*.0001)
					{
						if(e_xx_min > x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_min=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).y > sample.height()*.9999)
					{
						if(e_yy_max < x[triangles[k]->getBoundingPoint(p).id*2+1])
							e_yy_max=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).y < sample.height()*.0001)
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
					avg_Emodulus_xx += ar*(sigma11[k*npoints+l]/npoints)/(epsilon11[k*npoints+l]/npoints) ;
					avg_Emodulus_yy += ar*(sigma22[k*npoints+l]/npoints)/(epsilon22[k*npoints+l]/npoints) ;
					avg_Emodulus_xy += ar*(sigma12[k*npoints+l]/npoints)/(epsilon12[k*npoints+l]/npoints) ;
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
		
		std::cout << "apparent extension on the X axis: " << e_xx_max-e_xx_min << std::endl ;
		std::cout << "apparent extension on the Y axis: " << e_yy_max-e_yy_min << std::endl ;
		std::cout << "anisotropy factor: " << (e_yy_max-e_yy_min)/e_xx_max-e_xx_min << std::endl ;
		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
		if (tries < maxtries)
		{
			double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/(double)nstepstot ;
			if(!featureTree->solverConverged())
				delta_r *= .01 ;
			double reactedArea = 0 ;
			
			EllipsoidalInclusion * current = NULL ;
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
				apparent_extension.push_back(std::make_pair((e_xx_max-e_xx_min),(e_yy_max-e_yy_min))) ;
				young_modulus.push_back(std::make_pair(avg_Emodulus_xx,avg_Emodulus_yy)) ;
			}
			
			if (tries >= maxtries)
				break ;
		}
		
	std::string outfilename("outresultfile") ;
	outfilename.append(itoa(totit++, 10)) ;
	std::fstream outfilestream  ;
	outfilestream.open(outfilename.c_str(), std::ios::out) ;
		
	for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
		outfilestream << expansion_reaction[i].first << ",   " 
		<< expansion_reaction[i].second << ",   " 
		<< expansion_stress_xx[i].first << ",   " 
		<< expansion_stress_xx[i].second << ",   " 
		<< expansion_stress_yy[i].first << ",   " 
		<< expansion_stress_yy[i].second << ",   " 
		<< young_modulus[i].first << ",   "
		<< young_modulus[i].second << ",   "
		<< apparent_extension[i].first  << ",   " 
		<< apparent_extension[i].second  << ",   " 
		<< cracked_volume[i]  << ",   " 
		<< damaged_volume[i]  << ";   " 
		<< "\n" ;

	}

/*	for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
	std::cout << expansion_reaction[i].first << "   " 
	<< expansion_reaction[i].second << "   " 
	<< expansion_stress_xx[i].first << "   " 
	<< expansion_stress_xx[i].second << "   " 
	<< expansion_stress_yy[i].first << "   " 
	<< expansion_stress_yy[i].second << "   " 
	<< apparent_extension[i]  << "   " 
	<< cracked_volume[i]  << "   " 
	<< damaged_volume[i]  << "   " 
	<< std::endl ;*/
}


std::vector<std::pair<ExpansiveZone *, EllipsoidalInclusion *> > generateExpansiveZonesHomogeneously(int n, std::vector<EllipsoidalInclusion * > & incs , FeatureTree & F)
{
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
//	double nu_incompressible = .499997 ;
	
	double E = 1*E_csh ;
	double nu = nu_csh ; //nu_incompressible ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<std::pair<ExpansiveZone *, EllipsoidalInclusion *> > ret ;
	aggregateArea = 0 ;
	double radius = 0.0000005 ;
	Vector a(double(0), 3) ;
	a[0] = 0.5 ;
	a[1] = 0.5 ;
	a[2] = 0.00 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		Point pos(((double)rand()/RAND_MAX)*(sample.width()-radius*60),((double)rand()/RAND_MAX)*(sample.height()-radius*60)) ;
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
	std::map<EllipsoidalInclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			if(dist(zonesToPlace[i]->getCenter(), incs[j]->getCenter()) < incs[j]->getMinorRadius()-radius*60)
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
	for(std::map<EllipsoidalInclusion *, int>::iterator i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
//		std::cout << aggregateArea << "  " << count << std::endl ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}

std::vector<EllipsoidalInclusion *> circleToEllipse(std::vector<Inclusion *> circle)
{
	std::vector<EllipsoidalInclusion *> ellipse ;
	double alea_radius = 1 ;
	double alea_dir = 0 ;
	Point newcenter(0,0) ;
	Point newaxis(1,0) ;
	double ellipse_area = 0 ;
	for(int i = 0 ; i < circle.size() ; i++)
	{
		alea_radius = 1 / (1 + (double)rand()/((double)RAND_MAX)) ;
		alea_dir = 0.25 - (double)rand()/((double)RAND_MAX) / 2 ;
		newaxis.setY(alea_dir) ;
		ellipse.push_back(new EllipsoidalInclusion(circle[i]->getRadius(), alea_radius*circle[i]->getRadius(), newcenter, newaxis)) ;
		if(i<100)
		{
			std::cout << ellipse[i]->getMajorRadius() << " ; " << ellipse[i]->getMinorRadius() << std::endl ;
		}
		ellipse_area += ellipse[i]->area() ;
	}
	std::cout << "\n" << ellipse_area << "\n" << std::endl ;
	return ellipse ;
}

std::vector<EllipsoidalInclusion *> sortByMajorRadius(std::vector<EllipsoidalInclusion *> unsorted)
{
	std::vector<EllipsoidalInclusion *> ret ;
	std::vector<bool> sorted ;
	for(int i = 0 ; i < unsorted.size() ; i++)
		sorted.push_back(true) ;
	int done = 0 ;
	while(done<unsorted.size())
	{
		int j = 0 ;
		double radius = 0 ;
		for(int i = 0 ; i < unsorted.size() ; i++)
		{
			if(sorted[i])
			{
				if(unsorted[i]->getMajorRadius() > radius)
				{
					radius = unsorted[i]->getMajorRadius() ;
					j = i ;
				}
			}
		}
		sorted[j] = false ;
		ret.push_back(unsorted[j]) ;
		done = ret.size() ;
	}
	return ret ;
}

std::vector<EllipsoidalInclusion *> importEllipseList(std::string ellipsefile, int nell)
{
	std::cout << "importing ellipses from file... " << ellipsefile << std::endl ; 
	std::vector<EllipsoidalInclusion *> inc ;
	std::fstream ellipsein ;
	ellipsein.open(ellipsefile.c_str(),std::ios::in) ;
	double a ;
	double b ;
	double cx ;
	double cy ;
	double ax ;
	double ay ;
	char buff [256] ;
	int nimp = 0 ;
	while(!ellipsein.eof() && nimp < nell)
	{
		nimp++ ;
		ellipsein >> buff ;
		a = atof(buff) ;
		ellipsein >> buff ;
		b = atof(buff) ;
		ellipsein >> buff ;
		cx = atof(buff) ;
		ellipsein >> buff ;
		cy = atof(buff) ;
		ellipsein >> buff ;
		ax = atof(buff) ;
		ellipsein >> buff ;
		ay = atof(buff) ;
		inc.push_back(new EllipsoidalInclusion(a,b,cx,cy,ax,ay)) ;
	}
	ellipsein.close() ;
	inc.pop_back() ;
	std::cout << inc.size() << " ellipses imported" << std::endl ;
	return inc ;
}


int main(int argc, char *argv[])
{
	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1.-nu*nu) ; m0_paste[0][1] =E_paste/(1.-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1.-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1.-nu*nu) ; m0_paste[1][2] = 0 ; 
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1.-nu*nu)*(1.-nu)/2. ; 

	Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ; 
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ; 

	FeatureTree F(&sample) ;
	featureTree = &F ;

/*  	std::vector<double> columns ;
	columns.push_back(0.) ;
	columns.push_back(0.48) ;
	columns.push_back(0.52) ;
	columns.push_back(0.) ;
	GranuloFromFile tgranulo("granulo_true.txt",columns,0.001,1) ;
	tgranulo.resize(0.0005) ;
	std::vector<EllipsoidalInclusion *> test = tgranulo.getEllipsoidalInclusion(2.3,10000,0.00112) ;
	std::vector<EllipsoidalInclusion *> inc = sortByMajorRadius(test) ;
	double area_test = 0 ;
	for(int i = 0 ; i < inc.size() ; i++)
	{
		if(i%100 == 0)		
			std::cout << inc[i]->getMajorAxis().y << "     " ;
		area_test += inc[i]->area() ;
	}
	std::cout << area_test << " (" << (100*area_test/0.0016) << "%)" << std::endl ;*/

	std::string ellipselist = "ellipse_good.txt" ;
	std::vector<EllipsoidalInclusion *> inc = importEllipseList(ellipselist,9000) ;



//	sample.setBehaviour(new StiffnessAndFracture(m0_paste, new MohrCoulomb(13500000,-8*13500000))) ;
	sample.setBehaviour(new WeibullDistributedStiffness(m0_paste,13500000)) ;
//	sample.setBehaviour(new Stiffness(m0_paste)) ;
//	std::vector<Feature *> feats ;
//	for(size_t i = 0; i < inc.size() ; i++)
//		feats.push_back(inc[i]) ;
//	inc.clear() ;
//	std::cout << width*height << std::endl ;
//	int n_agg = feats.size() ;
//	feats=placement(sample.getPrimitive(), feats, &n_agg, 640000);
//	for(size_t i = 0; i < feats.size() ; i++)
//		inc.push_back(static_cast<EllipsoidalInclusion *>(feats[i])) ;
//
//	std::string ellipsefile = "ellipse_placed_sorted_9999" ;
//	std::fstream ellipseout ;
//	ellipseout.open(ellipsefile.c_str(), std::ios::out) ;
//	for(int i = 0 ; i < inc.size() ; i++)
//	{
//		ellipseout << inc[i]->getMajorRadius() << "    " ;
//		ellipseout << inc[i]->getMinorRadius() << "    " ;
//		ellipseout << inc[i]->getCenter().x << "    " ;
//		ellipseout << inc[i]->getCenter().y << "    " ;
//		ellipseout << inc[i]->getMajorAxis().x << "    " ;
//		ellipseout << inc[i]->getMajorAxis().y << "\n" ;
//	}
//	return 0 ;	
//	StiffnessAndFracture * stiff = new StiffnessAndFracture(m0_agg, new MohrCoulomb(57000000,-8*57000000));
	WeibullDistributedStiffness * stiff = new WeibullDistributedStiffness(m0_agg,57000000) ;
//	Stiffness * stiff = new Stiffness(m0_agg) ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		if(inc[i]->getCenter().x == 0 && inc[i]->getCenter().y == 0)
		{
			std::cout << "fail to place all inclusions" << std::endl ;
			std::cout << "last inclusion placed => " << i << std::endl ;
			return 1 ;
		}
		inc[i]->setBehaviour(stiff) ;
		F.addFeature(&sample, inc[i]) ;
//		inc[i]->getMajorAxis().print() ;
		placed_area += inc[i]->area() ;
	}
	
        zones = generateExpansiveZonesHomogeneously(10000,inc,F) ;
	
	F.sample(1024) ;
	F.setOrder(LINEAR) ;
        F.generateElements() ;

	std::cout << " => " << F.getTriangles().size() << std::endl ;

// 	stepOLD() ;
        step() ;

//	dt->print() ;
 	
// 	delete dt ;*/
	
	return 0 ;
}
