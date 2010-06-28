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
#include "../physics/stiffness.h"
#include "../physics/spatially_distributed_stiffness.h"
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
#include "../utilities/itoa.h"

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

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

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

void setBC()
{
  //if (/*un truc géométrique pour choper les bons points*/)
       //l'assemblage sur lequel on fout les BC
       //featureTree->getAssembly()->
       // au choix
       //* setPoint(x,y, point.id)
       //* setPointAlong({XI, ETA}, x, point.id)
       //* setForceOn({XI, ETA}, f, point.id)
       //* setForce(fx, fy, point.id)
       //* setPeriodic(point.id, point2.id)

	triangles = featureTree->getTriangles() ;
	double bctol  = 1e-6;
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			if(triangles[k]->getBoundingPoint(c).t < 0)
			{
				featureTree->getAssembly()->setPoint( 0,0 ,triangles[k]->getBoundingPoint(c).id) ;
			}
			else
			{
				if(std::abs(triangles[k]->getBoundingPoint(c).y - .00500) < bctol 
				&& (std::abs(triangles[k]->getBoundingPoint(c).x - .01) < bctol /*|| std::abs(triangles[k]->getBoundingPoint(c).x + .01) < bctol */))
				{
					featureTree->getAssembly()->setPoint( 0,0.001 ,triangles[k]->getBoundingPoint(c).id) ;
				}
// 				else if(std::abs(triangles[k]->getBoundingPoint(c).y + .00500) < bctol 
// 				&& (std::abs(triangles[k]->getBoundingPoint(c).x - .01) < bctol || std::abs(triangles[k]->getBoundingPoint(c).x + .01) < bctol ))
// 				{
// 					featureTree->getAssembly()->setPoint( 0,-0.001 ,triangles[k]->getBoundingPoint(c).id) ;
// 				}
				else if(std::abs(triangles[k]->getBoundingPoint(c).y - .00500) < bctol)// top boundary			
				{
					featureTree->getAssembly()->setPointAlong( ETA,0.001 ,triangles[k]->getBoundingPoint(c).id) ;
				}
				else if(std::abs(triangles[k]->getBoundingPoint(c).y + .00500) < bctol)// bottom boundary			
				{
					featureTree->getAssembly()->setPointAlong( ETA,-0.001 ,triangles[k]->getBoundingPoint(c).id) ;
				}
				else if(std::abs(triangles[k]->getBoundingPoint(c).x - .0100) < bctol)// top boundary			
				{
					featureTree->getAssembly()->setPointAlong( XI,0 ,triangles[k]->getBoundingPoint(c).id) ;
				}
				else if(std::abs(triangles[k]->getBoundingPoint(c).x + .0100) < bctol)// bottom boundary			
				{
					featureTree->getAssembly()->setPointAlong( XI,0 ,triangles[k]->getBoundingPoint(c).id) ;
				}
			}
		}
	}
}

void step()
{
	
	bool cracks_did_not_touch = true;
	size_t max_growth_steps = 1;
	size_t countit = 0;	
	
	while ( (cracks_did_not_touch) && (countit < max_growth_steps) )
	{
		countit++;
		std::cout << "\r iteration " << countit << "/" << max_growth_steps << std::flush ;
//		setBC() ;
      
		int limit = 0 ;
		while(!featureTree->step(timepos) && limit < 50)//as long as we can update the features
		{
			std::cout << "." << std::flush ;
// 			timepos-= 0.0001 ;
//			setBC() ;
			triangles = featureTree->getTriangles() ;
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
	
  // Prints the crack geo to a file for each crack
		if (cracks_did_not_touch == false) // if cracks touched
		{
			std::cout << "** Cracks touched exporting file **" << endl;	
				// Print the state of the cracks to a file  
				std::string filename = "crackGeo.txt";
				fstream filestr;
				filestr.open (filename.c_str(), fstream::in | fstream::out | fstream::app);
				filestr << "Crack vertices" << std::endl;
				filestr << "x" << " " << "y" << std::endl ; 
		}
  

		timepos+= 0.0001 ;
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
			
		std::cout << "energy index :" << enr << std::endl ;
		energy.push_back(enr) ;

		if(limit < 2)
			break ;
	}
	for(size_t i = 0 ; i < energy.size() ; i++)
		std::cout << energy[i] << std::endl ;
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
			x *= 1.5 ;
// 			sigma11 *= 1.5 ;
// 			sigma22 *= 1.5 ;
// 			sigma12 *= 1.5 ;
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
	int start = 0 ;
	for (unsigned int j=0 ; j< triangles.size() ; j++ )
	{
		if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			int timePlanes = triangles[j]->timePlanes() ;
			start = ((triangles[j]->getBoundingPoints().size()-1)*(timePlanes-1))/timePlanes ;
			break ;
		}
	}

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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min)))*300., 1., 1. ) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = start+1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma11[j*triangles[j]->getBoundingPoints().size()]-sigma11_min)/(sigma11_max-sigma11_min), 1., 1.) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = start+1 ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vonMises[j*triangles[j]->getBoundingPoints().size()]-vonMises_min)/(vonMises_max-vonMises_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(angle[j*triangles[j]->getBoundingPoints().size()]-angle_min)/(angle_max-angle_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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

		double crit_max = fracCrit.max() ;
		double crit_min = fracCrit.min() ;
		glNewList(  DISPLAY_LIST_FRAC_CRIT,  GL_COMPILE ) ;
		for (unsigned int j=0 ; j< triangles.size() ; j++ )
		{
			
			if(triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR && !cracked[j]&& !triangles[j]->getBehaviour()->fractured())
			{
				double c1 ;
				double c2 ;
				double c3 ;
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(fracCrit[j*triangles[j]->getBoundingPoints().size()]-crit_min)/(crit_max-crit_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
				{
					vx = x[triangles[j]->getBoundingPoint(k).id*2];
					vy = x[triangles[j]->getBoundingPoint(k).id*2+1]; 
					
					HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(fracCrit[j*triangles[j]->getBoundingPoints().size()+k]-crit_min)/(crit_max-crit_min), 1., 1.) ;
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma22[j*triangles[j]->getBoundingPoints().size()]-sigma22_min)/(sigma22_max-sigma22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(sigma12[j*triangles[j]->getBoundingPoints().size()]-sigma12_min)/(sigma12_max-sigma12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(start)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				
				Point a = triangles[j]->inLocalCoordinates(triangles[j]->getBoundingPoint(start)) ;
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(triangles[j]->getBehaviour()->getTensor(a)[0][0]-E_min)/(E_max-E_min), 1., .2) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon11[j*triangles[j]->getBoundingPoints().size()]-epsilon11_min)/(epsilon11_max-epsilon11_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon22[j*triangles[j]->getBoundingPoints().size()]-epsilon22_min)/(epsilon22_max-epsilon22_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				
				double vx = x[triangles[j]->getBoundingPoint(start).id*2]; 
				double vy = x[triangles[j]->getBoundingPoint(start).id*2+1]; 
				
				glBegin(GL_TRIANGLE_FAN);
				HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(epsilon12[j*triangles[j]->getBoundingPoints().size()]-epsilon12_min)/(epsilon12_max-epsilon12_min), 1., 1.) ;
				glColor3f(c1, c2, c3) ;
				glVertex2f(double(triangles[j]->getBoundingPoint(start).x + vx) , double(triangles[j]->getBoundingPoint(start).y + vy) );
				
				for(size_t k = 1+start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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
				for(size_t k = start ; k < triangles[j]->getBoundingPoints().size() ; k++)
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

		for(size_t k  = 0 ; k < crack.size() ; k++)
		{
			glColor3f(1, 0, 0) ;
// 			for(unsigned int j=0 ; j< tris__.size() ; j++ )
// 			{
// 				glBegin(GL_LINE_LOOP);
// 				double vx = x[tris__[j]->first->id*2]; 
// 				double vy = x[tris__[j]->first->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->first->x/*+ vx*/) ,
// 				            double(tris__[j]->first->y/*+ vy*/) );
// 				
// 				vx = x[tris__[j]->second->id*2]; 
// 				vy = x[tris__[j]->second->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->second->x/*+ vx*/) ,
// 				            double(tris__[j]->second->y/*+ vy*/) );
// 				
// 				vx = x[tris__[j]->third->id*2]; 
// 				vy = x[tris__[j]->third->id*2+1]; 
// 				
// 				glVertex2f( double(tris__[j]->third->x/*+ vx*/) ,
// 				            double(tris__[j]->third->y/*+ vy*/) );
// 				glEnd();
// 			}
// 			
// 			glColor3f(0, 1, 1) ;
			glBegin(GL_LINES) ;
			for(size_t j=0 ; j< crack[k]->getBoundingPoints().size()-1 ; j++ )
			{
				glVertex2f( double(crack[k]->getBoundingPoint(j).x) ,
				            double(crack[k]->getBoundingPoint(j).y) );
				glVertex2f( double(crack[k]->getBoundingPoint(j+1).x) ,
				            double(crack[k]->getBoundingPoint(j+1).y) );
			}
			glEnd();
		}
		
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

std::vector<EllipsoidalInclusion *> circleToEllipse(std::vector<Inclusion *> circle, double biais)
{
	std::vector<EllipsoidalInclusion *> ellipse ;
	double alea_radius = 0.5 ;
	double alea_dir = 0 ;
	Point newcenter(0.2,0.2) ;
	Point newaxis(1,0) ;
	double ellipse_area = 0 ;
	for(int i = 0 ; i < circle.size() ; i++)
	{
//                std::cout << circle[i]->getRadius()/10 << ";" ;
//		alea_radius = 1. / 2 ;
		alea_dir = biais * (1 - 2*(double)rand()/((double)RAND_MAX)) ;
		newaxis.setY(alea_dir) ;
                newaxis /= newaxis.norm() ;
		ellipse.push_back(new EllipsoidalInclusion(newcenter,newaxis*circle[i]->getRadius(),alea_radius)) ;
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
	double cx ;
	double cy ;
	double ax ;
	double ay ;
        double bx ;
        double by ;
	char buff [256] ;
	int nimp = 0 ;
	while(!ellipsein.eof() && nimp < nell)
	{
		nimp++ ;
		ellipsein >> buff ;
		cx = atof(buff) ;
		ellipsein >> buff ;
		cy = atof(buff) ;
		ellipsein >> buff ;
		ax = atof(buff) ;
		ellipsein >> buff ;
		ay = atof(buff) ;
		ellipsein >> buff ;
		bx = atof(buff) ;
		ellipsein >> buff ;
		by = atof(buff) ;
		Point C(cx,cy) ;
		Point A(ax,ay) ;
		Point B(bx,by) ;
		inc.push_back(new EllipsoidalInclusion(C,A,B)) ;
	}
	ellipsein.close() ;
	inc.pop_back() ;
	std::cout << inc.size() << " ellipses imported" << std::endl ;
	return inc ;
}



int main(int argc, char *argv[])
{
	int biais = atof(argv[1]) ;


        Sample box(NULL, 0.04,0.04, 0.02, 0.02) ;

	FeatureTree F(&box) ;
	featureTree = &F ;

	std::vector<Inclusion *> inclusions = GranuloBolome(0.000000743, 1., BOLOME_D)(.00025, .00001, 8000, 0.00005);
        std::vector<EllipsoidalInclusion *> ellipses = circleToEllipse(inclusions,((double) biais)/100.) ;
        ellipses = sortByMajorRadius(ellipses) ;

        std::vector<Feature *> feats ;
        for(size_t i = 0 ; i < ellipses.size() ; i++) {
            feats.push_back(static_cast<Feature *>(ellipses[i])) ;
        }
        ellipses.clear() ;

        int nAgg = 1 ;
        feats=placement(box.getPrimitive(), feats, &nAgg, 10000);


/*        Grid grid(box.getPrimitive()->getRadius()*2, box.getPrimitive()->getRadius()*2, 100, box.getCenter()) ;

        std::cout << "checking for intersection... takes some time..." << std::endl ;*/
        for(size_t i = 0 ; i < 200 ; i++) {
            ellipses.push_back(static_cast<EllipsoidalInclusion *>(feats[i])) ;
 //           if(!grid.add(ellipses[i]->getPrimitive()))
 //               std::cout << i << std::endl ;
        }

/*	std::string ellipsefile = "ellipse_" ;
        double bcent = biais * 100 ;

        std::string bstring = itoa((int) biais, 10) ;
        ellipsefile.append(bstring) ;

	std::cout << ellipsefile << std::endl ;

	std::fstream ellipseout ;
	ellipseout.open(ellipsefile.c_str(), std::ios::out) ;
	for(int i = 0 ; i < ellipses.size() ; i++)
	{
		ellipseout << ellipses[i]->getCenter().x << "    " ;
		ellipseout << ellipses[i]->getCenter().y << "    " ;
		ellipseout << ellipses[i]->getMajorAxis().x << "    " ;
		ellipseout << ellipses[i]->getMajorAxis().y << "    " ;
		ellipseout << ellipses[i]->getMinorAxis().x << "    " ;
		ellipseout << ellipses[i]->getMinorAxis().y << "\n" ;
	}*/

/*	Point O(0.,0.) ;
	Point A(1.,0.) ;
	Point B(0.,0.5) ;

	Point C(2.000000001,1.) ;

	Ellipse E1(O,A,B) ;
	Ellipse E2(C,A,B) ;

	if(E1.intersects(&E2))
		std::cout << "intersection" << std::endl ;

	std::vector<Point> point = E1.getBoundingBox() ;
	for(int i = 0 ; i < point.size() ; i++)
		point[i].print() ;

	point = E2.getBoundingBox() ;
	for(int i = 0 ; i < point.size() ; i++)
		point[i].print() ;

        return 0 ;*/

/*	std::string ellipsefile = "ellipse_" ;
        std::string bstring = itoa(biais,10) ;
        ellipsefile.append(bstring) ;

        std::vector<EllipsoidalInclusion *> ellipses = importEllipseList(ellipsefile,6001) ;*/

/*	for(size_t i = 0 ; i < ellipses.size() ; i++)
		ellipses[i]->sampleBoundingSurface(40) ;

	std::cout << "checking for intersection... again" << std::endl ;
	std::vector<std::pair<size_t, size_t> > inter ;
	for(size_t i = 0 ; i < ellipses.size() ; i++)
	{
		if(i%100 == 0)
			std::cout << i << " => " ;
		for(size_t j = 0 ; j < i ; j++)
		{
			if(j%100 == 0 && i%100 == 0)
				std::cout << j << " ; " ;
//			std::cout << j << std::endl ;
			if((ellipses[i]->getPrimitive()->intersection(ellipses[j]->getPrimitive())).size() > 0)
				inter.push_back(std::make_pair(i,j)) ;
		}
		if(i%100 == 0)
			std::cout << inter.size() << " intersections found so far..." << std::endl ;
	}
	for(size_t n = 0 ; n < inter.size() ; n++)
	{
		size_t i = inter[n].first ;
		size_t j = inter[n].second ;
		std::cout << "----" << std::endl ;
		std::cout << ellipses[i]->getMajorRadius() + ellipses[i]->getCenter().x << ";";
		std::cout << - ellipses[i]->getMajorRadius() + ellipses[i]->getCenter().x << std::endl ;
		std::cout << ellipses[i]->getMajorRadius() + ellipses[i]->getCenter().y << ";" ;
		std::cout << - ellipses[i]->getMajorRadius() + ellipses[i]->getCenter().y << std::endl ;
		std::cout << ellipses[j]->getMajorRadius() + ellipses[j]->getCenter().x << ";" ;
		std::cout << - ellipses[j]->getMajorRadius() + ellipses[j]->getCenter().x << std::endl ;
		std::cout << ellipses[j]->getMajorRadius() + ellipses[j]->getCenter().y << ";" ;
		std::cout << - ellipses[j]->getMajorRadius() + ellipses[j]->getCenter().y << std::endl ;
	}
//		std::cout << inter[i].first << ";" << inter[i].second << std::endl ;

//	return 0 ;*/


        double E_agg = 70e9 ;
        double E_paste = 25e9 ;
        double nu = 0.2 ;
        Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ;
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ;

	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ;
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ;

        box.setBehaviour(new Stiffness(m0_paste)) ;
        for(size_t i = 0 ; i < ellipses.size() ; i++) {
            F.addFeature(&box,ellipses[i]) ;
            ellipses[i]->setBehaviour(new Stiffness(m0_agg)) ;
        }

        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 0.1)) ;
        F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	F.sample(500) ;
	F.setOrder(LINEAR) ;

	F.generateElements() ;

	step() ;
	
/*	glutInit(&argc, argv) ;	
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
	glutAddMenuEntry(" Frac Crit     ", ID_FRAC_CRIT);
	
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
	
// 	delete dt ;*/
	
	return 0 ;
}
