// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics.h"
#include "../physics/mohrcoulomb.h"
#include "../physics/laplacian.h"
#include "../physics/ruptureenergy.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 

#define ID_QUIT 1
#define ID_ZOOM 2
#define ID_UNZOOM 3
#define ID_NEXT10 4
#define ID_NEXT100 5
#define ID_NEXT1000 6
#define ID_NEXT 7
#define ID_NEXT_TIME 8
#define ID_REFINE 9
#define ID_AMPLIFY 10
#define ID_DEAMPLIFY 11

#define ID_DISP 12
#define ID_STRAIN_XX 13
#define ID_STRAIN_XY 14
#define ID_STRAIN_XZ 15
#define ID_STRAIN_YZ 16
#define ID_STRAIN_YY 17
#define ID_STRAIN_ZZ 18
#define ID_STRESS_XX 19
#define ID_STRESS_XY 20
#define ID_STRESS_YY 21
#define ID_STRESS_ZZ 22
#define ID_STRESS_XZ 23
#define ID_STRESS_YZ 24

#define ID_STIFNESS 25
#define ID_ELEM 26
#define ID_VON_MISES 27
#define ID_ANGLE 28
#define ID_ENRICHMENT 29

GLuint DISPLAY_LIST_DISPLACEMENT = 0 ;
GLuint DISPLAY_LIST_ELEMENTS = 0 ;
GLuint DISPLAY_LIST_STRAIN_XX = 0 ;
GLuint DISPLAY_LIST_STRAIN_YY = 0 ;
GLuint DISPLAY_LIST_STRAIN_ZZ = 0 ;
GLuint DISPLAY_LIST_STRAIN_XY = 0 ;
GLuint DISPLAY_LIST_STRAIN_XZ = 0 ;
GLuint DISPLAY_LIST_STRAIN_YZ = 0 ;
GLuint DISPLAY_LIST_STRESS_XX = 0 ;
GLuint DISPLAY_LIST_STRESS_YY = 0 ;
GLuint DISPLAY_LIST_STRESS_ZZ = 0 ;
GLuint DISPLAY_LIST_STRESS_XY = 0 ;
GLuint DISPLAY_LIST_STRESS_XZ = 0 ;
GLuint DISPLAY_LIST_STRESS_YZ = 0 ;
GLuint DISPLAY_LIST_STIFFNESS = 0 ;
GLuint DISPLAY_LIST_VON_MISES  = 0 ;
GLuint DISPLAY_LIST_ANGLE  = 0 ;
GLuint DISPLAY_LIST_ENRICHMENT = 0 ;
GLuint DISPLAY_LIST_STIFFNESS_DARK = 0 ;

using namespace Mu ;

int viewangle = 0 ;
int viewangle2 = 0 ;

FeatureTree * featureTree ;
std::vector<DelaunayTetrahedron *> tets ;
std::vector<bool> cracked ;
DelaunayTree *dt ; //(pts) ;
std::vector<Crack *> crack ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;
double z_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;
double z_min = 0 ;

GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double timepos = 0.00 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

int windowWidth = 600 ;
int windowHeight = 600 ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma33(0) ; 
Vector sigma12(0) ; 
Vector sigma13(0) ; 
Vector sigma23(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon33(0) ; 
Vector epsilon12(0) ; 
Vector epsilon13(0) ; 
Vector epsilon23(0) ; 
Vector vonMises(0) ; 
Vector stiffness(0) ; 
Vector angle(0) ; 
Vector damage(0) ; 

double nu = 0.2 ;
double E_agg = 100 ;//softest
double E_paste = 1 ;//stiff
double E_stiff = E_agg*10 ;//stiffer
double E_soft = E_agg/10; //stiffest

size_t current_list = DISPLAY_LIST_DISPLACEMENT ;
double factor = 0.3 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;


void setBC()
{
  
  tets = featureTree->getTetrahedrons() ;
  
	double bctol = 0.0001; // tolerance used to identify points where to impose bcs  
	double left = 0 ;
	double bottom = 0 ;
	double right = 400 ;
	double top = 400 ;
	for(size_t k = 0 ; k < tets.size() ;k++)
	{
		for(size_t c = 0 ;  c < tets[k]->getBoundingPoints().size() ; c++ )
		{
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).x - (left)) < bctol) 
// 				&& (std::abs(tets[k]->getBoundingPoint(c).y - (left)) < bctol) 
// 				&& (std::abs(tets[k]->getBoundingPoint(c).z - (left)) < bctol)
// 				)
// 			{
// 				featureTree->getAssembly()->setPoint( -50, 0, 0, tets[k]->getBoundingPoint(c).id) ;
// 			}
			if( (std::abs(tets[k]->getBoundingPoint(c).x - (left)) < bctol) )
			{
				featureTree->getAssembly()->setPoint( -50, 0, 0,tets[k]->getBoundingPoint(c).id) ;
			}
			else if ( (std::abs(tets[k]->getBoundingPoint(c).x - (right)) < bctol) )
			{
				featureTree->getAssembly()->setPoint( 50, 0,0, tets[k]->getBoundingPoint(c).id) ;
			}
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).x - (right)) < bctol) 
// 				&& (std::abs(tets[k]->getBoundingPoint(c).y - (right)) < bctol) 
// 				&& (std::abs(tets[k]->getBoundingPoint(c).z - (right)) < bctol)
// 				)
// 			{
// 				featureTree->getAssembly()->setPoint( 50, 0, 0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).y - (right)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA,0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).y - (left)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPointAlong( ETA,0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).z - (right)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPointAlong( ZETA,0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			if ( (std::abs(tets[k]->getBoundingPoint(c).z - (left)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPointAlong( ZETA,0, tets[k]->getBoundingPoint(c).id) ;
// 			}
		}

// 		for(size_t c = 0 ;  c < tets[k]->getBoundingPoints().size() ; c++ )
// 		{
// /*		  	if ( (std::abs(tets[k]->getBoundingPoint(c).x - (left)) < bctol) &&  (std::abs(tets[k]->getBoundingPoint(c).z - (top)) < bctol)&&  (std::abs(tets[k]->getBoundingPoint(c).y - (top)) < bctol))
// 			{
// 				featureTree->getAssembly()->setPoint( 0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			else */if ( (std::abs(tets[k]->getBoundingPoint(c).x - (left)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPoint( 10, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			else if ( (std::abs(tets[k]->getBoundingPoint(c).x - (right)) < bctol) )
// 			{
// 				featureTree->getAssembly()->setPoint( 0, tets[k]->getBoundingPoint(c).id) ;
// 			}
// 			
// 		}
	}
}

void step()
{
	
  int nsteps = 1;// number of steps between two clicks on the opengl thing


	for(size_t i = 0 ; i < nsteps ; i++)
	{
		std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;

		bool go_on = true;
		size_t tries = 0;
		while(go_on && tries < 50)
		{
			setBC() ;
			featureTree->step(timepos) ;
			go_on = featureTree->solverConverged() &&  (featureTree->meshChanged() || featureTree->enrichmentChanged());
			std::cout << "." << std::flush ;
			// 			timepos-= 0.0001 ;
			
			tries++ ;
		}
		std::cout << " " << tries << " tries." << std::endl ;
		timepos+= 0.0001 ;
		
		tets= featureTree->getDelaunayTree3D()->getTetrahedrons() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
		dt = featureTree->getDelaunayTree() ;

		std::pair<Vector, Vector > sigma_epsilon ;
		sigma_epsilon.first.resize(24*tets.size()) ;
		sigma_epsilon.second.resize(24*tets.size()) ;
		sigma_epsilon = featureTree->getStressAndStrain(tets) ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		sigma11.resize(sigma.size()/6, 0.) ;
		sigma22.resize(sigma.size()/6, 0.) ;
		sigma33.resize(sigma.size()/6, 0.) ;
		sigma12.resize(sigma.size()/6, 0.) ;
		sigma13.resize(sigma.size()/6, 0.) ;
		sigma23.resize(sigma.size()/6, 0.) ;
		
		epsilon11.resize(sigma.size()/6, 0.) ;
		epsilon22.resize(sigma.size()/6, 0.) ;
		epsilon33.resize(sigma.size()/6, 0.) ;
		epsilon12.resize(sigma.size()/6, 0.) ;
		epsilon13.resize(sigma.size()/6, 0.) ;
		epsilon23.resize(sigma.size()/6, 0.) ;
		stiffness.resize(sigma.size()/6, 0.) ;
		vonMises.resize(sigma.size()/6, 0.) ;
		angle.resize(sigma.size()/6, 0.) ;
		damage.resize(sigma.size()/6, 0.) ;
	
		std::cout << "unknowns :" << x.size() << std::endl ;
	
	
		int npoints = 4 ;
	
		double volume = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_zz = 0;
		double avg_e_xy = 0;
		double avg_e_xz = 0;
		double avg_e_yz = 0;
		double avg_s_xx = 0;
		double avg_s_zz = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double avg_s_xz = 0;
		double avg_s_yz = 0;
		double e_xx = 0 ;
		double ex_count = 0 ;
		double xavg = 0 ;
		
		for(size_t k = 0 ; k < tets.size() ; k++)
		{
			
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				
// 				for(size_t p = 0 ;p < 4 ; p++)
// 				{
// 					if(x[tets[k]->getBoundingPoint(p).id*3] > x_max)
// 						x_max = x[tets[k]->getBoundingPoint(p).id*3];
// 					if(x[tets[k]->getBoundingPoint(p).id*3] < x_min)
// 						x_min = x[tets[k]->getBoundingPoint(p).id*3];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+1] > y_max)
// 						y_max = x[tets[k]->getBoundingPoint(p).id*3+1];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+1] < y_min)
// 						y_min = x[tets[k]->getBoundingPoint(p).id*3+1];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+2] > z_max)
// 						z_max = x[tets[k]->getBoundingPoint(p).id*3+2];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+2] < z_min)
// 						z_min = x[tets[k]->getBoundingPoint(p).id*3+2];
// 					if(tets[k]->getBoundingPoint(p).x > 0.0799)
// 					{
// 						e_xx+=x[tets[k]->getBoundingPoint(p).id*3] ;
// 						ex_count++ ;
// 					}
// 				}
				volume += tets[k]->volume() ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] > E_max)
						E_max = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] < E_min)
						E_min = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					
					stiffness[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					damage[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
				}
					
				sigma11[k*npoints] = sigma[k*npoints*6];
				sigma22[k*npoints] = sigma[k*npoints*6+1];
				sigma33[k*npoints] = sigma[k*npoints*6+2];
				sigma12[k*npoints] = sigma[k*npoints*6+3];
				sigma13[k*npoints] = sigma[k*npoints*6+4];
				sigma23[k*npoints] = sigma[k*npoints*6+5];
				
				sigma11[k*npoints+1] = sigma[k*npoints*6+6];
				sigma22[k*npoints+1] = sigma[k*npoints*6+7];
				sigma33[k*npoints+1] = sigma[k*npoints*6+8];
				sigma12[k*npoints+1] = sigma[k*npoints*6+9];
				sigma13[k*npoints+1] = sigma[k*npoints*6+10];
				sigma23[k*npoints+1] = sigma[k*npoints*6+11];
				
				sigma11[k*npoints+2] = sigma[k*npoints*6+12];
				sigma22[k*npoints+2] = sigma[k*npoints*6+13];
				sigma33[k*npoints+2] = sigma[k*npoints*6+14];
				sigma12[k*npoints+2] = sigma[k*npoints*6+15];
				sigma13[k*npoints+2] = sigma[k*npoints*6+16];
				sigma23[k*npoints+2] = sigma[k*npoints*6+17];
				
				sigma11[k*npoints+3] = sigma[k*npoints*6+18];
				sigma22[k*npoints+3] = sigma[k*npoints*6+19];
				sigma33[k*npoints+3] = sigma[k*npoints*6+20];
				sigma12[k*npoints+3] = sigma[k*npoints*6+21];
				sigma13[k*npoints+3] = sigma[k*npoints*6+22];
				sigma23[k*npoints+3] = sigma[k*npoints*6+23];
				
				epsilon11[k*npoints] = epsilon[k*npoints*6];
				epsilon22[k*npoints] = epsilon[k*npoints*6+1];
				epsilon33[k*npoints] = epsilon[k*npoints*6+2];
				epsilon12[k*npoints] = epsilon[k*npoints*6+3];
				epsilon13[k*npoints] = epsilon[k*npoints*6+4];
				epsilon23[k*npoints] = epsilon[k*npoints*6+5];
				
				epsilon11[k*npoints+1] = epsilon[k*npoints*6+6];
				epsilon22[k*npoints+1] = epsilon[k*npoints*6+7];
				epsilon33[k*npoints+1] = epsilon[k*npoints*6+8];
				epsilon12[k*npoints+1] = epsilon[k*npoints*6+9];
				epsilon13[k*npoints+1] = epsilon[k*npoints*6+10];
				epsilon23[k*npoints+1] = epsilon[k*npoints*6+11];
				
				epsilon11[k*npoints+2] = epsilon[k*npoints*6+12];
				epsilon22[k*npoints+2] = epsilon[k*npoints*6+13];
				epsilon33[k*npoints+2] = epsilon[k*npoints*6+14];
				epsilon12[k*npoints+2] = epsilon[k*npoints*6+15];
				epsilon13[k*npoints+2] = epsilon[k*npoints*6+16];
				epsilon23[k*npoints+2] = epsilon[k*npoints*6+17];
				
				epsilon11[k*npoints+3] = epsilon[k*npoints*6+18];
				epsilon22[k*npoints+3] = epsilon[k*npoints*6+19];
				epsilon33[k*npoints+3] = epsilon[k*npoints*6+20];
				epsilon12[k*npoints+3] = epsilon[k*npoints*6+21];
				epsilon13[k*npoints+3] = epsilon[k*npoints*6+22];
				epsilon23[k*npoints+3] = epsilon[k*npoints*6+23];
				
				double vm0 = 0 ;
				double agl = 0 ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					vm0 = tets[k]->getState().getMaximumVonMisesStress() ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					agl = tets[k]->getState().getPrincipalAngle(tets[k]->getCenter()) ;
				for(size_t l = 0 ; l < 4 ; l++)
				{
					vonMises[k*npoints+l]  = vm0 ;
					angle[k*npoints+l]  = agl ;
				}
				
				double ar = tets[k]->volume() ;
				for(size_t l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
					avg_e_zz += (epsilon33[k*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
					avg_e_xz += (epsilon13[k*npoints+l]/npoints)*ar;
					avg_e_yz += (epsilon23[k*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
					avg_s_zz += (sigma33[k*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
					avg_s_xz += (sigma13[k*npoints+l]/npoints)*ar;
					avg_s_yz += (sigma23[k*npoints+l]/npoints)*ar;
					xavg += x[tets[k]->getBoundingPoint(l).id]*ar/npoints ;
				}
	
			}
		}
			
		xavg /= volume ;
		
		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "avg value :" << xavg << std::endl ;
		std::cout << "max sigma11 :" << sigma11.max() << std::endl ;
		std::cout << "min sigma11 :" << sigma11.min() << std::endl ;
		std::cout << "max sigma12 :" << sigma12.max() << std::endl ;
		std::cout << "min sigma12 :" << sigma12.min() << std::endl ;
		std::cout << "max sigma13 :" << sigma13.max() << std::endl ;
		std::cout << "min sigma13 :" << sigma13.min() << std::endl ;
		std::cout << "max sigma22 :" << sigma22.max() << std::endl ;
		std::cout << "min sigma22 :" << sigma22.min() << std::endl ;
		std::cout << "max sigma23 :" << sigma23.max() << std::endl ;
		std::cout << "min sigma23 :" << sigma23.min() << std::endl ;
		std::cout << "max sigma33 :" << sigma33.max() << std::endl ;
		std::cout << "min sigma33 :" << sigma33.min() << std::endl ;
		
		std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
		std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
		std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
		std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
		std::cout << "max epsilon13 :" << epsilon13.max() << std::endl ;
		std::cout << "min epsilon13 :" << epsilon13.min() << std::endl ;
		std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
		std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;
		std::cout << "max epsilon23 :" << epsilon23.max() << std::endl ;
		std::cout << "min epsilon23 :" << epsilon23.min() << std::endl ;
		std::cout << "max epsilon33 :" << epsilon33.max() << std::endl ;
		std::cout << "min epsilon33 :" << epsilon33.min() << std::endl ;
		
		std::cout << "max von Mises :" << vonMises.max() << std::endl ;
		std::cout << "min von Mises :" << vonMises.min() << std::endl ;
		
		std::cout << "average sigma11 : " << avg_s_xx/volume << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy/volume << std::endl ;
		std::cout << "average sigma33 : " << avg_s_zz/volume << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy/volume << std::endl ;
		std::cout << "average sigma13 : " << avg_s_xz/volume << std::endl ;
		std::cout << "average sigma23 : " << avg_s_yz/volume << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx/volume << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy/volume << std::endl ;
		std::cout << "average epsilon33 : " << avg_e_zz/volume << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy/volume << std::endl ;
		std::cout << "average epsilon13 : " << avg_e_xz/volume << std::endl ;
		std::cout << "average epsilon23 : " << avg_e_yz/volume << std::endl ;
		
	}

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
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB) ;
	glShadeModel (GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL) ;
	GLfloat mat_specular[] = {0.01f, 0.01f, 0.01f, 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glColorMaterial(GL_AMBIENT_AND_DIFFUSE, GL_FRONT_AND_BACK) ;
// 	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	
// 	glEnable(GL_NORMALIZE) ;
	glClearColor(0.0f,0.0f,0.0f,1.0f);                                      // Black Background
	glClearDepth(-0.f);                                                     // Depth Buffer Setup
	
	
	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	
	GLfloat light_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat light_diffuse[] = { 0.6, 0.6, 0.6, 1.0 };
	GLfloat light_specular[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_position[] = { 1, 0, 0, 0 };
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
// 	glEnable(GL_BLEND) ;
// 	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glDepthRange( 0, 1 ) ;
	
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
	case ID_STRAIN_ZZ : 
		{
			current_list = DISPLAY_LIST_STRAIN_ZZ ;
			break ;
		}
	case ID_STRAIN_XY : 
		{
			current_list = DISPLAY_LIST_STRAIN_XY ;
			break ;
		}
	case ID_STRAIN_XZ : 
		{
			current_list = DISPLAY_LIST_STRAIN_XZ ;
			break ;
		}
	case ID_STRAIN_YZ : 
		{
			current_list = DISPLAY_LIST_STRAIN_YZ ;
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
	case ID_STRESS_ZZ : 
		{
			current_list = DISPLAY_LIST_STRESS_ZZ ;
			break ;
		}
	case ID_STRESS_XY : 
		{
			current_list = DISPLAY_LIST_STRESS_XY ;
			break ;
		}
	case ID_STRESS_XZ : 
		{
			current_list = DISPLAY_LIST_STRESS_XZ ;
			break ;
		}
	case ID_STRESS_YZ : 
		{
			current_list = DISPLAY_LIST_STRESS_YZ ;
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

std::pair<double, double> centile(const Vector & v)
{
	Vector vs(v) ;
	std::sort(&vs[0], &vs[vs.size()]) ;
	return std::make_pair(vs[(vs.size()-1)*0.01], vs[(vs.size()-1)*0.99]) ;
}

void reshape(int w, int h)
{
	windowWidth = w ;
	windowHeight = h ;
	glViewport (0, 0, w, h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(45.0, 1., (double)w/(double)h, 2.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity() ;
}

void printScreen()
{
	std::valarray<unsigned int> frame(windowWidth*windowHeight) ;
	glReadPixels(0,0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, &frame[0]) ;
}

void displayTet(DelaunayTetrahedron * t, int j, std::valarray<double> & vals, double min, double max)
{
	if(t->getCenter().y > 200 && t->getCenter().x < 200 && t->getCenter().z < 200)
	{
		return ;
	}
	double c1 ;
	double c2 ;
	double c3 ;
// 	glBegin(GL_POINTS);
	double vx = x[t->first->id*3]+x[t->second->id*3]+x[t->third->id*3]+x[t->fourth->id*3]; 
	double vy = x[t->first->id*3+1]+x[t->second->id*3+1]+x[t->third->id*3+1]+x[t->fourth->id*3+1]; 
	double vz = x[t->first->id*3+2]+x[t->second->id*3+2]+x[t->third->id*3+2]+x[t->fourth->id*3+2]; 
// 	vx = 0 ; vy = 0 ; vz = 0 ;
// 	
	double v = (vals[j*4]+vals[j*4+1]+vals[j*4+2]+vals[j*4+3])*.25 ;
	HSVtoRGB( &c1, &c2, &c3, 300. - 240.*(v-min)/(max-min), 1., 1.) ;
	glColor3f(c1, c2, c3) ;
	glVertex3f(double(t->getCenter().x + vx*.25) , double(t->getCenter().y + vy*.25), double(t->getCenter().z + vz*.25) );
// 	glEnd() ;
// 	glBegin(GL_LINE_LOOP);
// 		double vx = x[t->first->id*3]; 
// 		double vy = x[t->first->id*3+1]; 
// 		double vz = x[t->first->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->first->x + vx) , double(t->first->y + vy), double(t->first->z + vz) );
// 		
// 		vx = x[t->second->id*3]; 
// 		vy = x[t->second->id*3+1]; 
// 		vz = x[t->second->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+1]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->second->x + vx) , double(t->second->y + vy),  double(t->second->z + vz) );
// 		
// 		vx = x[t->third->id*3]; 
// 		vy = x[t->third->id*3+1]; 
// 		vz = x[t->third->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+2]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->third->x + vx) , double(t->third->y + vy), double(t->third->z + vz) );
// 	glEnd() ;
// 	
// 	glBegin(GL_LINE_LOOP);
// 		vx = x[t->first->id*3]; 
// 		vy = x[t->first->id*3+1]; 
// 		vz = x[t->first->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+0]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->first->x + vx) , double(t->first->y + vy), double(t->first->z + vz) );
// 		
// 		vx = x[t->second->id*3]; 
// 		vy = x[t->second->id*3+1]; 
// 		vz = x[t->second->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+1]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->second->x + vx) , double(t->second->y + vy), double(t->second->z + vz) );
// 		
// 		vx = x[t->fourth->id*3]; 
// 		vy = x[t->fourth->id*3+1]; 
// 		vz = x[t->fourth->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+3]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->fourth->x + vx) , double(t->fourth->y + vy), double(t->fourth->z + vz) );
// 	glEnd() ;
// 	
// 	glBegin(GL_LINE_LOOP);
// 		vx = x[t->first->id*3]; 
// 		vy = x[t->first->id*3+1]; 
// 		vz = x[t->first->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+0]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->first->x + vx) , double(t->first->y + vy), double(t->first->z + vz) );
// 		
// 		vx = x[t->third->id*3]; 
// 		vy = x[t->third->id*3+1]; 
// 		vz = x[t->third->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+2]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->third->x + vx) , double(t->third->y + vy), double(t->third->z + vz) );
// 		
// 		vx = x[t->fourth->id*3]; 
// 		vy = x[t->fourth->id*3+1]; 
// 		vz = x[t->fourth->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+3]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->fourth->x + vx) , double(t->fourth->y + vy), double(t->fourth->z + vz) );
// 	glEnd() ;
// 	
// 	glBegin(GL_LINE_LOOP);
// 		vx = x[t->third->id*3]; 
// 		vy = x[t->third->id*3+1]; 
// 		vz = x[t->third->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+2]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->third->x + vx) , double(t->third->y + vy), double(t->third->z + vz) );
// 		
// 		vx = x[t->second->id*3]; 
// 		vy = x[t->second->id*3+1]; 
// 		vz = x[t->second->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+1]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->second->x + vx) , double(t->second->y + vy), double(t->second->z + vz) );
// 		
// 		vx = x[t->fourth->id*3]; 
// 		vy = x[t->fourth->id*3+1]; 
// 		vz = x[t->fourth->id*3+2]; 
// 		HSVtoRGB( &c1, &c2, &c3, 300. - 300.*(vals[j*4+3]-min)/(max-min), 1., 1.) ;
// 		glColor3f(c1, c2, c3) ;
// 		glVertex3f(double(t->fourth->x + vx) , double(t->fourth->y + vy), double(t->fourth->z + vz) );
// 	glEnd() ;
}

void Display(void)
{
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	glMatrixMode(GL_MODELVIEW);	
	glLoadIdentity() ;
	glScalef(.002, .002, .002) ;
	glRotatef( viewangle, 0, -1, 0) ;
	glRotatef(-viewangle2 , -1, 0, 0) ;
	glRotatef(20 , 0, 0, 1) ;
	glTranslatef(-200,-200, -200) ;
	glEnable(GL_DEPTH_TEST) ;
	if(!dlist)
	{
		
		DISPLAY_LIST_ENRICHMENT = glGenLists(1);
		glNewList( DISPLAY_LIST_ENRICHMENT,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && damage[j*4]<1)
			{
				displayTet(tets[j],j , damage, 0, 1) ;
			}
		}
		glEnd() ;
// 			for (unsigned int j=0 ; j< tets.size() ; j++ )
// 			{
// 				
// 				if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR  && tets[j]->getBehaviour()->fractured())
// 				{
// 					double c1 ;
// 					double c2 ;
// 					double c3 ;
// 					
// 					glBegin(GL_LINE_LOOP);
// 					double vx = x[tets[j]->first->id*3]; 
// 					double vy = x[tets[j]->first->id*3+1]; 
// 					double vz = x[tets[j]->first->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
// 					
// 					vx = x[tets[j]->second->id*3]; 
// 					vy = x[tets[j]->second->id*3+1]; 
// 					vz = x[tets[j]->second->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy),  double(tets[j]->second->z + vz) );
// 					
// 					vx = x[tets[j]->third->id*3]; 
// 					vy = x[tets[j]->third->id*3+1]; 
// 					vz = x[tets[j]->third->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
// 					glEnd() ;
// 					
// 					glBegin(GL_LINE_LOOP);
// 					vx = x[tets[j]->first->id*3]; 
// 					vy = x[tets[j]->first->id*3+1]; 
// 					vz = x[tets[j]->first->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
// 					
// 					vx = x[tets[j]->second->id*3]; 
// 					vy = x[tets[j]->second->id*3+1]; 
// 					vz = x[tets[j]->second->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy), double(tets[j]->second->z + vz) );
// 					
// 					vx = x[tets[j]->fourth->id*3]; 
// 					vy = x[tets[j]->fourth->id*3+1]; 
// 					vz = x[tets[j]->fourth->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
// 					glEnd() ;
// 					
// 					glBegin(GL_LINE_LOOP);
// 					vx = x[tets[j]->first->id*3]; 
// 					vy = x[tets[j]->first->id*3+1]; 
// 					vz = x[tets[j]->first->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
// 					
// 					vx = x[tets[j]->third->id*3]; 
// 					vy = x[tets[j]->third->id*3+1]; 
// 					vz = x[tets[j]->third->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
// 					
// 					vx = x[tets[j]->fourth->id*3]; 
// 					vy = x[tets[j]->fourth->id*3+1]; 
// 					vz = x[tets[j]->fourth->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
// 					glEnd() ;
// 					
// 					glBegin(GL_LINE_LOOP);
// 					vx = x[tets[j]->third->id*3]; 
// 					vy = x[tets[j]->third->id*3+1]; 
// 					vz = x[tets[j]->third->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
// 					
// 					vx = x[tets[j]->second->id*3]; 
// 					vy = x[tets[j]->second->id*3+1]; 
// 					vz = x[tets[j]->second->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy), double(tets[j]->second->z + vz) );
// 					
// 					vx = x[tets[j]->fourth->id*3]; 
// 					vy = x[tets[j]->fourth->id*3+1]; 
// 					vz = x[tets[j]->fourth->id*3+2]; 
// 					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
// 					glColor3f(c1, c2, c3) ;
// 					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
// 					glEnd() ;
// 				}
// 			}
		glEndList() ;
		
		DISPLAY_LIST_DISPLACEMENT = glGenLists(1);
		glNewList( DISPLAY_LIST_DISPLACEMENT,  GL_COMPILE ) ;
		double xmin = x.min() ; double xmax = x.max() ;
		for (unsigned int j=0 ; j< tets.size() ; j++ )
			{
				
				if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR  && !tets[j]->getBehaviour()->fractured())
				{
					double c1 ;
					double c2 ;
					double c3 ;
					
					glBegin(GL_LINE_LOOP);
					double vx = x[tets[j]->first->id*3]; 
					double vy = x[tets[j]->first->id*3+1]; 
					double vz = x[tets[j]->first->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
					
					vx = x[tets[j]->second->id*3]; 
					vy = x[tets[j]->second->id*3+1]; 
					vz = x[tets[j]->second->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy),  double(tets[j]->second->z + vz) );
					
					vx = x[tets[j]->third->id*3]; 
					vy = x[tets[j]->third->id*3+1]; 
					vz = x[tets[j]->third->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
					glEnd() ;
					
					glBegin(GL_LINE_LOOP);
					vx = x[tets[j]->first->id*3]; 
					vy = x[tets[j]->first->id*3+1]; 
					vz = x[tets[j]->first->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
					
					vx = x[tets[j]->second->id*3]; 
					vy = x[tets[j]->second->id*3+1]; 
					vz = x[tets[j]->second->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy), double(tets[j]->second->z + vz) );
					
					vx = x[tets[j]->fourth->id*3]; 
					vy = x[tets[j]->fourth->id*3+1]; 
					vz = x[tets[j]->fourth->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
					glEnd() ;
					
					glBegin(GL_LINE_LOOP);
					vx = x[tets[j]->first->id*3]; 
					vy = x[tets[j]->first->id*3+1]; 
					vz = x[tets[j]->first->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->first->x + vx) , double(tets[j]->first->y + vy), double(tets[j]->first->z + vz) );
					
					vx = x[tets[j]->third->id*3]; 
					vy = x[tets[j]->third->id*3+1]; 
					vz = x[tets[j]->third->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
					
					vx = x[tets[j]->fourth->id*3]; 
					vy = x[tets[j]->fourth->id*3+1]; 
					vz = x[tets[j]->fourth->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
					glEnd() ;
					
					glBegin(GL_LINE_LOOP);
					vx = x[tets[j]->third->id*3]; 
					vy = x[tets[j]->third->id*3+1]; 
					vz = x[tets[j]->third->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->third->x + vx) , double(tets[j]->third->y + vy), double(tets[j]->third->z + vz) );
					
					vx = x[tets[j]->second->id*3]; 
					vy = x[tets[j]->second->id*3+1]; 
					vz = x[tets[j]->second->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->second->x + vx) , double(tets[j]->second->y + vy), double(tets[j]->second->z + vz) );
					
					vx = x[tets[j]->fourth->id*3]; 
					vy = x[tets[j]->fourth->id*3+1]; 
					vz = x[tets[j]->fourth->id*3+2]; 
					HSVtoRGB( &c1, &c2, &c3, 300. - sqrt(((vx-x_min)*(vx-x_min) + (vy-y_min)*(vy-y_min)+ (vz-z_min)*(vz-z_min))/((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min)))*300., 1., 1. ) ;
					glColor3f(c1, c2, c3) ;
					glVertex3f(double(tets[j]->fourth->x + vx) , double(tets[j]->fourth->y + vy), double(tets[j]->fourth->z + vz) );
					glEnd() ;
				}
			}
		glEndList() ;
		
		std::pair<double, double> minmax = centile(sigma11) ;
		DISPLAY_LIST_STRAIN_XX = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_XX,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
			for (unsigned int j=0 ; j< tets.size() ; j++ )
			{
				
				if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
				{
					displayTet(tets[j],j , sigma11, minmax.first, minmax.second) ;
				}
			}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(stiffness) ;
		DISPLAY_LIST_STIFFNESS = glGenLists(1);
		glNewList(  DISPLAY_LIST_STIFFNESS,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , stiffness, minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(vonMises) ;
		DISPLAY_LIST_VON_MISES = glGenLists(1);
		glNewList(  DISPLAY_LIST_VON_MISES,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , vonMises, minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		
		minmax = centile(angle) ;
		DISPLAY_LIST_ANGLE = glGenLists(1);
		glNewList(  DISPLAY_LIST_ANGLE,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR&& !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , angle,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(sigma22) ;
		DISPLAY_LIST_STRAIN_YY = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_YY,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR&& !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , sigma22, minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(sigma33) ;
		DISPLAY_LIST_STRAIN_ZZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_ZZ,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR&& !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , sigma33, minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(sigma12) ;
		DISPLAY_LIST_STRAIN_XY = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_XY,  GL_COMPILE ) ;	
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR &&  !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , sigma12, minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(sigma13) ; 
		DISPLAY_LIST_STRAIN_XZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_XZ,  GL_COMPILE ) ;	
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR &&  !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , sigma13,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(sigma23) ; 
		DISPLAY_LIST_STRAIN_YZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRAIN_YZ,  GL_COMPILE ) ;	
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR &&  !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , sigma23,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(epsilon11) ; 
		DISPLAY_LIST_STRESS_XX = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_XX,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon11,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
				
		minmax = centile(epsilon22) ; 
		DISPLAY_LIST_STRESS_YY = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_YY,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR&& !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon22,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(epsilon33) ; 
		DISPLAY_LIST_STRESS_ZZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_ZZ,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR&& !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon33,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(epsilon12) ; 
		DISPLAY_LIST_STRESS_XY = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_XY,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon12,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(epsilon13) ; 
		DISPLAY_LIST_STRESS_XZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_XZ,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon13,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		minmax = centile(epsilon23) ; 
		DISPLAY_LIST_STRESS_YZ = glGenLists(1);
		glNewList(  DISPLAY_LIST_STRESS_YZ,  GL_COMPILE ) ;
		glBegin(GL_POINTS);
		for (unsigned int j=0 ; j< tets.size() ; j++ )
		{
			
			if(tets[j]->getBehaviour()->type != VOID_BEHAVIOUR && !tets[j]->getBehaviour()->fractured())
			{
				displayTet(tets[j],j , epsilon23,  minmax.first, minmax.second) ;
			}
		}
		glEnd() ;
		glEndList() ;
		
		dlist = true ;
		glCallList(current_list) ;
	}
	else
	{

		double c1, c2, c3 = 0;
		HSVtoRGB( &c1, &c2, &c3, 180. + 0, 1., 1.) ;
		
		glCallList(current_list) ;
		glColor3f(1, 1, 1) ;

		
	}
	
// 	glTranslatef(20,20,20) ;
// 	glBegin(GL_LINES) ;
// 	glColor3f(1, 1, 1) ;
// 	glVertex3f(20, 20, 40) ;
// 	glVertex3f(20, 20, 0) ;
// 	
// 	glVertex3f(0, 40, 40) ;
// 	glVertex3f(40, 40, 40) ;
// 	
// 	glVertex3f(40, 40, 40) ;
// 	glVertex3f(40, 0, 40) ;
// 	
// 	glVertex3f(40, 0, 40) ;
// 	glVertex3f(0, 0, 40) ;
// 	
// 	glVertex3f(0, 0, 40) ;
// 	glVertex3f(0, 40, 40) ;
// 	
// 	glVertex3f(0, 40, 0) ;
// 	glVertex3f(40, 40, 0) ;
// 	
// 	glVertex3f(40, 4, 0) ;
// 	glVertex3f(40, 0, 0) ;
// 	
// 	glVertex3f(40, 0, 0) ;
// 	glVertex3f(0, 0, 0) ;
// 	
// 	glVertex3f(0, 0, 0) ;
// 	glVertex3f(0, 40, 0) ;
// 	glEnd() ;

	glColor3f	(1, 0, 0) ;
	glFlush();
	glutSwapBuffers();
}

void processMouseActiveMotion(int x, int y) {
	viewangle++ ;
	viewangle %=360 ;
	viewangle2++ ;
	viewangle2++ ;
	viewangle2 %=360 ;
	Display() ;
}

int main(int argc, char *argv[])
{
	double nu = 0.2 ;
	double E = 1 ;

	Sample3D sample(NULL, 0.04,0.04,0.04,0.02,0.02,0.02) ;
	Sample3D samplers(NULL, 400,400,400,200,200,200) ;

	FeatureTree F(&samplers) ;
	featureTree = &F ;
	Matrix d0(3,3) ;
	d0[0][0] = 1 ;
	d0[1][1] = 1 ;
	d0[2][2] = 1 ;

	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;

	nu = 0.2 ;
	E = 4 ;
	Matrix m1(6,6) ;
	m1[0][0] = 1. - nu ; m1[0][1] = nu ; m1[0][2] = nu ;
	m1[1][0] = nu ; m1[1][1] = 1. - nu ; m1[1][2] = nu ;
	m1[2][0] = nu ; m1[2][1] = nu ; m1[2][2] = 1. - nu ;
	m1[3][3] = 0.5 - nu ;
	m1[4][4] = 0.5 - nu ;
	m1[5][5] = 0.5 - nu ;
	m1 *= E/((1.+nu)*(1.-2.*nu)) ;

	Matrix d1(3,3) ;
	d1[0][0] = 1e2 ;
	d1[1][1] = 1e2 ;
	d1[2][2] = 1e2 ;

	double itzSize = 0.000001;
	int inclusionNumber = 0 ; //48000*2 ;

	std::vector<Inclusion3D *> inclusions = GranuloBolome(5.44e-05, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);

// 	if(inclusionNumber)
// 		itzSize = inclusions[inclusions.size()-1]->getRadius() ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		delete inclusions[i] ;
// 
// 	inclusions = GranuloBolome(3.84e-05/1.8, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;

	int nAgg = 6000 ;
	feats=placement(sample.getPrimitive(), feats, &nAgg, 60000);
	
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;
	Stiffness * s = new Stiffness(m0) ;
	Stiffness * ss = new Stiffness(m1) ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		feats[i]->setCenter(feats[i]->getCenter()*10000.);
		dynamic_cast<Inclusion3D *>(feats[i])->setRadius(feats[i]->getRadius()*10000.)  ;
		F.addFeature(&samplers, feats[i]) ;
		feats[i]->setBehaviour(ss) ;
// 		feats[i]->setBehaviour(new VoidForm()) ;
// 		feats[i]->setBehaviour(new Laplacian(d1)) ;
// 		std::cout << feats[i]->getRadius() << "   " << feats[i]->getCenter().x << "   "  << feats[i]->getCenter().y << "   " << feats[i]->getCenter().z << std::endl ; 
	}
	samplers.setBehaviour(new WeibullDistributedStiffness(m0,0.1)) ;
// 	samplers.setBehaviour(new Laplacian(d0)) ;
	Vector a(6.,0) ;
// 	ExpansiveZone3D * inc = new ExpansiveZone3D(&samplers,100/*166.11322368308294*/, 200, 200, 200, m1, a) ;
	Inclusion3D * inc0 = new Inclusion3D(100/*166.11322368308294*/, 200, 200, 200) ;
// 	OctahedralInclusion * inc0 = new OctahedralInclusion(208.40029238347645, 200, 200, 200) ;
	inc0->setBehaviour(new WeibullDistributedStiffness(m1,0.4)) ;
// 	inc0->setBehaviour(new Laplacian(d1)) ;
// 	F.addFeature(&samplers, inc) ;
	F.addFeature(&samplers, inc0) ;
	F.sample(2048) ;
	F.setOrder(LINEAR) ;
	F.generateElements() ;
	step() ;


	glutInit(&argc, argv) ;	
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize(600, 600) ;
	glutReshapeFunc(reshape) ;
	glutCreateWindow("coucou !") ;
	glutMotionFunc(processMouseActiveMotion);
	int submenu = glutCreateMenu(Menu) ;
	
	glutAddMenuEntry(" Displacements ", ID_DISP);
	glutAddMenuEntry(" Stress (s) xx ", ID_STRAIN_XX);
	glutAddMenuEntry(" Stress (s) yy ", ID_STRAIN_YY);
	glutAddMenuEntry(" Stress (s) zz ", ID_STRAIN_ZZ);
	glutAddMenuEntry(" Stress (s) xy ", ID_STRAIN_XY);
	glutAddMenuEntry(" Stress (s) xz ", ID_STRAIN_XZ);
	glutAddMenuEntry(" Stress (s) yz ", ID_STRAIN_YZ);
	glutAddMenuEntry(" Strain (e) xx ", ID_STRESS_XX);
	glutAddMenuEntry(" Strain (e) yy ", ID_STRESS_YY);
	glutAddMenuEntry(" Strain (e) zz ", ID_STRESS_ZZ);
	glutAddMenuEntry(" Strain (e) xy ", ID_STRESS_XY);
	glutAddMenuEntry(" Strain (e) xz ", ID_STRESS_XZ);
	glutAddMenuEntry(" Strain (e) yz ", ID_STRESS_YZ);
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
	
	return 0 ;
}