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
#include "../physics/laplacian.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../utilities/writer/voxel_writer.h"
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

void step()
{
	
  int nsteps = 1;// number of steps between two clicks on the opengl thing
	featureTree->setMaxIterationsPerStep(50) ;

	for(size_t i = 0 ; i < nsteps ; i++)
	{

		bool go_on = featureTree->step() ;
		
		tets= featureTree->getElements3D() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;

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
				
				Vector vm0(0.,1) ;
				Vector agl(0.,3) ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					tets[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD,  tets[k]->getCenter(), vm0, false) ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					tets[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, tets[k]->getCenter(), agl, false) ;
				for(size_t l = 0 ; l < 4 ; l++)
				{
					vonMises[k*npoints+l]  = vm0[0] ;
					angle[k*npoints+l]  = agl[0] ;
				}
				
				double ar = tets[k]->volume() ;
				Vector avgsig(6) ;
				Vector avgeps(6) ;
				tets[k]->getState().getAverageField(REAL_STRESS_FIELD, avgsig);
				tets[k]->getState().getAverageField(STRAIN_FIELD, avgeps);

				avg_e_xx += avgeps[0] * ar;
				avg_e_yy += avgeps[1] * ar;
				avg_e_zz += avgeps[2] * ar;
				avg_e_xy += avgeps[3] * ar;
				avg_e_xz += avgeps[4] * ar;
				avg_e_yz += avgeps[5] * ar;
				
				avg_s_xx += avgsig[0] * ar;
				avg_s_yy += avgsig[1] * ar;
				avg_s_zz += avgsig[2] * ar;
				avg_s_xy += avgsig[3] * ar;
				avg_s_xz += avgsig[4] * ar;
				avg_s_yz += avgsig[5] * ar;
			
				for(size_t l = 0 ; l < npoints ;l++)
				{
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
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();
	
	VoxelWriter vw("sphere_stress", 200) ;
	vw.getField(featureTree, VWFT_STRESS) ;
	vw.write();
	VoxelWriter vw0("sphere_strain", 200) ;
	vw0.getField(featureTree, VWFT_STRAIN) ;
	vw0.write();
	exit(0) ;
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

void addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep, TetrahedralElement *father )
{

	if( nodes_per_side > 1 )
		nodes_per_side = 1 ;


		std::vector<std::pair<Point, Point> > sides ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 0 ), father->getBoundingPoint( 1 ) ) ) ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 1 ), father->getBoundingPoint( 2 ) ) ) ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 2 ), father->getBoundingPoint( 3 ) ) ) ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 3 ), father->getBoundingPoint( 0 ) ) ) ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 3 ), father->getBoundingPoint( 1 ) ) ) ;
		sides.push_back( std::make_pair( father->getBoundingPoint( 0 ), father->getBoundingPoint( 2 ) ) ) ;
		std::vector<size_t> positions ;

		if( nodes_per_side )
		{
			positions.push_back( 0 ) ;
			positions.push_back( 1 ) ;
			positions.push_back( 2 ) ;

			positions.push_back( 2 ) ;
			positions.push_back( 3 ) ;
			positions.push_back( 4 ) ;

			positions.push_back( 4 ) ;
			positions.push_back( 5 ) ;
			positions.push_back( 6 ) ;

			positions.push_back( 0 ) ;
			positions.push_back( 7 ) ;
			positions.push_back( 6 ) ;

			positions.push_back( 6 ) ;
			positions.push_back( 8 ) ;
			positions.push_back( 2 ) ;

			positions.push_back( 0 ) ;
			positions.push_back( 9 ) ;
			positions.push_back( 4 ) ;
		}
		else
		{
			positions.push_back( 0 ) ;
			positions.push_back( 1 ) ;

			positions.push_back( 1 ) ;
			positions.push_back( 2 ) ;

			positions.push_back( 2 ) ;
			positions.push_back( 3 ) ;

			positions.push_back( 3 ) ;
			positions.push_back( 0 ) ;

			positions.push_back( 3 ) ;
			positions.push_back( 1 ) ;

			positions.push_back( 0 ) ;
			positions.push_back( 2 ) ;
		}

		size_t nodes_per_plane = nodes_per_side * 6 + 4 ;

		std::valarray<Point *> newPoints( ( Point * )nullptr, nodes_per_plane * time_planes ) ;
		std::valarray<bool> done( false, nodes_per_plane * time_planes ) ;

		for( size_t plane = 0 ; plane < time_planes ; plane++ )
		{
			size_t current = 0 ;

			for( size_t side = 0 ; side < 6 ; side++ )
			{
				Point a( sides[side].first ) ;
				Point b( sides[side].second ) ;

				if( time_planes > 1 )
				{
					a.t = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
					b.t = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
// 					std::cout << a.t << std::endl ;
				}

				for( size_t node = 0 ; node < nodes_per_side + 2 ; node++ )
				{
					double fraction = ( double )( node ) / ( ( double )nodes_per_side + 1 ) ;
					Point proto = a * ( 1. - fraction ) + b * fraction ;

					if( !done[positions[current] + plane * nodes_per_plane] )
					{
						newPoints[positions[current] + plane * nodes_per_plane]  = new Point( proto ) ;
						done[positions[current] + plane * nodes_per_plane] = true ;
					}

					current++ ;
				}
			}
		}

		for( size_t k = 0 ; k < newPoints.size() ; k++ )
			if( !newPoints[k] )
			{
				std::cout << "ouch !" << std::endl ;

				for( size_t k = 0 ; k < newPoints.size() ; k++ )
					if( newPoints[k] )
						newPoints[k]->print() ;

				exit( 0 ) ;
			}

		father->setBoundingPoints( newPoints ) ;

}


int main(int argc, char *argv[])
{

//   Function test("1 x - y - 0.5 0.5 t * - *") ;
// 	VirtualMachine().print(test);
//              std::cout << VirtualMachine().eval(test, Point(0,0,0,-1)) 
// << std::endl ;
//              std::cout << VirtualMachine().eval(test, Point(0,1,0,-1)) 
// << std::endl ;
	
// 	Function x("x") ;
// 	Function xm("1 x 2 * 1 - abs -") ;
// 	x = 1-f_abs(x*2.-1) ;
// 	VirtualMachine().print(x) ;
// 	std::cout << VirtualMachine().eval(x, 0.) <<"  " << VirtualMachine().eval(x, 1.) <<"  " <<VirtualMachine().eval(x, 2.) <<"  " <<std::endl ;
// 	VirtualMachine().print(xm) ;
// 	std::cout << VirtualMachine().eval(xm, 0.) <<"  " << VirtualMachine().eval(xm, 1.) <<"  " <<VirtualMachine().eval(xm, 2.) <<"  " <<std::endl ;
// 	exit(0) ;
	
	
	TetrahedralElement toto(QUADRATIC) ;
	addSharedNodes(1, 1, 0, &toto) ;
// 	for(int i = 0 ; i < 10 ; i++)
// 	{
// 		for(int j = 0 ; j < 10 ; j++)
// 		{
// 			std::cout << VirtualMachine().eval(toto.getShapeFunction(i), toto.getBoundingPoint(j).x, toto.getBoundingPoint(j).y, toto.getBoundingPoint(j).z) <<"  "<<std::flush ;
// 		}
// 		std::cout <<std::endl ;
// 	}
// 	exit(0) ;
// 	VirtualMachine().print(toto.getShapeFunction(0));
// 	VirtualMachine().print(toto.getShapeFunction(1));
// 	VirtualMachine().print(toto.getShapeFunction(2));
// 	VirtualMachine().print(toto.getShapeFunction(3));
// 	toto.getBoundingPoint(0).x -=.3 ;
// 	toto.getBoundingPoint(1).x -=.3 ;
// 	toto.getBoundingPoint(2).x -=.3 ;
// 	toto.getBoundingPoint(3).x -=.3 ;
// 	Function xtrans = toto.getXTransform() ;
// 	VirtualMachine().print(xtrans);
// 	Function ytrans = toto.getYTransform() ;
// 	VirtualMachine().print(ytrans);
// 	std::cout << VirtualMachine().eval(ytrans, 1., 1., 1.) <<"  " << VirtualMachine().eval(ytrans, 2., 2., 2.) <<"  " <<VirtualMachine().eval(ytrans, 3., 3., 3.) <<"  " <<std::endl ;
// 	Function ztrans = toto.getZTransform() ;
// 	VirtualMachine().print(ztrans);
// 	
// 	
		Function position(Point(0,0,0), &toto) ;
// // 		VirtualMachine().print(position);
// // 		exit(0) ;
// // 		VirtualMachine().print(position-1.);
// // 		VirtualMachine().print(f_abs(position-1.));
		Function hat = 1.-f_abs(position-.5);
// 		VirtualMachine().print(hat);
// 		VirtualMachine().print(1.-hat);
// 		exit(0) ;
		
// 		VirtualMachine().print(toto.getShapeFunction(1));
// 		VirtualMachine().print(toto.getShapeFunction(1)*(hat-1.));
// 		for(size_t k = 6 ; k < 7 ; k++)
// 		{
// 			for(double i = -0.1 ; i < 1.1 ; i += 0.01)
// 			{
// 				for(double j = -0.1 ; j < 1.1 ; j += 0.01)
// 				{
// 					if(toto.in(Point(i, j, 0.25)))
// 						std::cout << VirtualMachine().eval(toto.getShapeFunction(k)*(hat-VirtualMachine().eval(hat,toto.getBoundingPoint(k).x, toto.getBoundingPoint(k).y, toto.getBoundingPoint(k).z)), i, j, 0.25) <<"  " << std::flush ;
// 					else
// 						std::cout << 0 <<"  " << std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// 		}
// 			for(double i = -1.1 ; i < 1.1 ; i += 0.01)
// 			{
// 				for(double j = -1.1 ; j < 1.1 ; j += 0.01)
// 				{
// // 					if(toto.in(Point(i, j)))
// 						std::cout << VirtualMachine().eval(ytrans, i, j) <<"  " << std::flush ;
// // 					else
// // 						std::cout << 0 <<"  " << std::flush ;
// 				}
// 				std::cout << std::endl ;
// 			}
// 		}
// 		exit(0) ;
// 	VirtualMachine().print(position);
// 	Function f =  toto.getShapeFunction(1)*(hat - VirtualMachine().eval(hat,toto.getBoundingPoint(0).x, toto.getBoundingPoint(0).y)) ;
// 	exit(0) ;
	double nu = 0.2 ;
	double E = 1 ;
	Sample3D samplers(nullptr, 400,400,400,200,200,200) ;

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
	E = 10 ;
	Matrix m1(6,6) ;
	m1[0][0] = 1. - nu ; m1[0][1] = nu ;      m1[0][2] = nu ;
	m1[1][0] = nu ;      m1[1][1] = 1. - nu ; m1[1][2] = nu ;
	m1[2][0] = nu ;      m1[2][1] = nu ;      m1[2][2] = 1. - nu ;
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

//	std::vector<Inclusion3D *> inclusions = GranuloBolome(5.44e-05, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);
// 	std::vector<Inclusion3D *> inclusions = ParticleSizeDistribution::get3DInclusions(0.0025, 5.44e-5, BOLOME_A, PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;
//	GranuloBolome(5.44e-05, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);
// 	if(inclusionNumber)
// 		itzSize = inclusions[inclusions.size()-1]->getRadius() ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		delete inclusions[i] ;
// 
// 	inclusions = GranuloBolome(3.84e-05/1.8, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);

// 	std::vector<Feature *> feats ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		feats.push_back(inclusions[i]) ;

// 	int nAgg = 6000 ;
	
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;
	Stiffness * s = new Stiffness(m0) ;
	Stiffness * ss = new Stiffness(m1) ;
// 	std::cout << feats.size() << std::endl ;
// 	for(size_t i = 0 ; i < feats.size() ; i++)
// 	{
// 		feats[i]->setCenter(feats[i]->getCenter()*10000.);
// 		dynamic_cast<Inclusion3D *>(feats[i])->setRadius(feats[i]->getRadius()*10000.)  ;
// 		F.addFeature(&samplers, feats[i]) ;
// 		feats[i]->setBehaviour(ss) ;
// // 		feats[i]->setBehaviour(new VoidForm()) ;
// // 		feats[i]->setBehaviour(new Laplacian(d1)) ;
// // 		std::cout << feats[i]->getRadius() << "   " << feats[i]->getCenter().x << "   "  << feats[i]->getCenter().y << "   " << feats[i]->getCenter().z << std::endl ; 
// 	}
	samplers.setBehaviour(new /*WeibullDistributed*/Stiffness(m0/*,0.1*/)) ;
// 	samplers.setBehaviour(new Laplacian(d0)) ;
	
	Vector a(0.,6) ;// a[0] = 1 ; a[1] = 1 ; a[2] = 1 ; 
// 	ExpansiveZone3D * inc = new ExpansiveZone3D(&samplers,100, 200, 200, 200, m1*4., a) ;
	
	Inclusion3D * inc = new Inclusion3D(100, 200, 200, 200) ;
// 	OctahedralInclusion * inc0 = new OctahedralInclusion(208.40029238347645, 200, 200, 200) ;
	inc->setBehaviour(new StiffnessWithImposedDeformation(m1*4.,a)) ;
// 	inc->setBehaviour(new Stiffness(m1)) ;
// 	inc0->setBehaviour(new Laplacian(d1)) ;
	
	F.addFeature(&samplers, inc) ;
	
// 	F.addFeature(&samplers, inc0) ;
	F.setSamplingNumber(atoi(argv[1])) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_RIGHT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP_LEFT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_FRONT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_FRONT)) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ZETA, FRONT, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_LEFT_BACK)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, RIGHT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, TOP)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;

	F.setOrder(QUADRATIC) ;

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
