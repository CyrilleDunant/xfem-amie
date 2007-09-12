// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "samplingcriterion.h"
#include "solvers/assembly.h"
#include "features/pore3d.h"
#include "features/inclusion3d.h"
#include "features/sample3d.h"
#include "delaunay_3d.h"
#include "filters/voxelporefilter.h"
#include "physics/void_form.h"

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <time.h> 
#define DEBUG 

#define ID_QUIT 0.5
#define ID_ZOOM 5
#define ID_UNZOOM 6
#define ID_NEXT10 7
#define ID_NEXT100 3
#define ID_NEXT1000 4
#define ID_NEXT 1
#define ID_REFINE 8
#define ID_AMPLIFY 9
#define ID_DEAMPLIFY 10



using namespace Mu ;

double factor = 0.5 ;

std::vector<double> sizered;
std::vector<DelaunayTetrahedron *> myTets ;
std::vector<HexahedralElement *> myHexs ;
std::vector<Point *> points ;
DelaunayTree_3D *dt ;
GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double max_x = 5.;
double maxs = 0;
double mins = 0;

Vector * x ;
Vector * sigma ;

std::vector<Point *> pts ;




void HSV2RGB( double *r, double *g, double *b, double h, double s, double v )
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



void computeDisplayList()
{
	glEnable(GL_DEPTH_TEST);
	glNewList(1,GL_COMPILE) ;
// 	myTets= dt->getTetrahedrons() ;
	double min = sigma->min() ;
	double max = sigma->max() ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{

// 		if(dist(*myTets[i]->getCenter(), cen) <1.8 || 
// 		   dist(*myTets[i]->getCenter(), cen1) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen2) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen3) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen4) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen5) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen6) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen7) <0.9 ||
// 		   dist(*myTets[i]->getCenter(), cen8) <0.9 
// 		  )
// 				glColor4ub (255,0,255,10);
// 			else
// 				glColor4ub (0,0,255,10); 			
			
			double r, g, b ;
			
// 		HSV2RGB(&r, &g, &b, 360.-360.*sqrt(((*sigma)[i]-min)/(max-min)), 1., 1.) ;
// 		glColor4f (r,g,b, .3);
			
		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
// 			glColor3f (0,0,0);
			
// 			if(!s.in(myTets[i]->first)  ||
// 			   !s.in(myTets[i]->second) ||
// 			   !s.in(myTets[i]->third) ||
// 			   !s.in(myTets[i]->fourth) ||
// 			   myTets[i]->Tetrahedron::volume() < 1e-6 
// 			  )
// 					glColor4ub (255,0,0,255); 		

			
// 			HSV2RGB(&r, &g, &b, 180.-180.*(myTets[i]->getBehaviour()->param[0][0]-mins)/(maxs-mins), 1., 1.) ;
// 			glColor4f (r,g,b, .3);
			
			HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glBegin(GL_LINES);
			
	
				glVertex3f(myTets[i]->first->x, myTets[i]->first->y, myTets[i]->first->z);

				glVertex3f(myTets[i]->second->x, myTets[i]->second->y, myTets[i]->second->z );
				
				glVertex3f(myTets[i]->first->x, myTets[i]->first->y, myTets[i]->first->z );
				
				glVertex3f(myTets[i]->third->x, myTets[i]->third->y, myTets[i]->third->z);

				glVertex3f(myTets[i]->first->x, myTets[i]->first->y, myTets[i]->first->z );
				
				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z );
	
				glVertex3f(myTets[i]->second->x, myTets[i]->second->y, myTets[i]->second->z );
				
				glVertex3f(myTets[i]->third->x, myTets[i]->third->y, myTets[i]->third->z);

				glVertex3f(myTets[i]->second->x, myTets[i]->second->y, myTets[i]->second->z );

				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z);
				
				glVertex3f(myTets[i]->third->x, myTets[i]->third->y, myTets[i]->third->z);

				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z);

			glEnd();

	}
	}
	
	
	
	glEndList() ;
}


void init(void)
{
	glShadeModel (GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL) ;
	GLfloat mat_specular[] = {0.01f, 0.01f, 0.01f, 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glColorMaterial(GL_AMBIENT_AND_DIFFUSE, GL_FRONT_AND_BACK) ;
// 	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	
// 	glEnable(GL_NORMALIZE) ;
	glClearColor(0.0f,0.0f,0.0f,1.0f);                                      // Black Background
	glClearDepth(1.0f);                                                     // Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);                                               // Disables Depth Testing
	
	
// 	glEnable(GL_LIGHTING) ;
// 	glEnable(GL_LIGHT0) ;
	
	GLfloat light_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat light_diffuse[] = { 0.6, 0.6, 0.6, 1.0 };
	GLfloat light_specular[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_position[] = { 4, 0, 0, 0 };
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
// 	glEnable(GL_BLEND) ;
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	computeDisplayList() ;
	
}


void display(void)
{
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity() ;
	glTranslatef(0,0,-10) ;
	gluLookAt(0, 0, 2, 0, 0, 0, 0, 1, 0) ;
	glRotatef(xangle , 0, -1, 0) ;
	glRotatef(yangle , -1, 0, 0) ;
	glRotatef(zangle , 0, 0, 1) ;
	glEnable(GL_DEPTH_TEST);
	
	GLdouble eq[] = { 0, 1, 0, max_x} ;
	
	glClipPlane(GL_CLIP_PLANE0, eq) ;
	glEnable(GL_CLIP_PLANE0) ;
	
	glCallList(1) ;
	glutSwapBuffers();
}

void reshape (int w, int h)
{
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(45.0, (GLfloat) w/(GLfloat) h, 1., 40);
}
 
void keyboard (unsigned char key, int x, int y)
   {
  	switch (key) {

  	case 'x':
	  	{
		  	xangle++ ;
	  		break;
	  	}
  	case 'y':
	  	{
		  	yangle++ ;
		  	break;
	  	}
  	case 'z':
	  	{
		  	zangle++ ;
		  	break;
	  	}
  	case '+':
	  	{
		  	max_x +=0.1 ;
		  	break;
	  	}
  	case '-':
	  	{
		  	max_x -=0.1 ;
		  	break;
	  	}
	default:
  		break;
  	}
	   
	   glutPostRedisplay() ;
  }

void mouse(int button, int state, int x, int y)
{

	switch (button) {
	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN) {
			
		}
		if (state == GLUT_UP) {
			
		}
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN) {
		
		}
		if (state == GLUT_UP) {
			
		}
		glutPostRedisplay();
		
		break;
	
	default:
		break;
		
	
	
	
	
	
// 	switch (button) {
// 	case GLUT_LEFT_BUTTON:
// 		if (state == GLUT_DOWN) {
// 			shoulder = (shoulder + .05) ;
// 			glutPostRedisplay();
// 			break;
// 		}
// 		break;
// 	case GLUT_RIGHT_BUTTON:
// 		if (state == GLUT_DOWN)
// 			shoulder = (shoulder - .05);
// 			glutPostRedisplay();
// 		break;
// 		
// 	default:
// 		break;
	}
}


int main(int argc, char *argv[])
{
	
	Assembly * K = new Assembly() ;
	
	Matrix diffusionMatrix(3,3) ;
	
	diffusionMatrix[0][0] = 0.005;
	diffusionMatrix[1][1] = 0.005;
	diffusionMatrix[2][2] = 0.005;
	
	//1 Alite
	//2 C-S-H
	//3 CH
	//4 C-S-H
	//5   "
	
	VoxelPoreFilter microstruct ;
	
	microstruct.behaviour = new Diffusion(diffusionMatrix) ;
	
	microstruct.read("/home/cyrille/xfem++/pixels20.txt") ;
	std::cout << "reading done" << std::endl ;
	
	for(size_t i = 0 ; i < microstruct.getPoints().size() ; i++)	
	{
		microstruct.getPoints()[i]->id = -1 ;
	}
	
	int index = 0 ;
	for(size_t i = 0 ; i < microstruct.getElements().size() ; i++)
	{
		
		if(microstruct.getElements()[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			for(size_t j = 0 ;j < microstruct.getElements()[i]->getBoundingPoints().size() ; j++)
			{
				if(microstruct.getElements()[i]->getBoundingPoint(j).id == -1)
					microstruct.getElements()[i]->getBoundingPoint(j).id = index++ ;
			}
			
			K->add(microstruct.getElements()[i]) ;
			microstruct.getElements()[i]->getState()->initialize() ;
		}
	}
	
	std::cout << "adding done" << std::endl ;
	
	for(size_t i = 0 ; i < microstruct.getPoints().size() ; i++)	
	{
		if(microstruct.getPoints()[i]->id != -1)
		{
			if(i%100 == 0)
				std::cout << "\rBC point " << i << "/" << microstruct.getPoints().size() << std::flush ;
	
			if(microstruct.getPoints()[i]->t < 1e-9)
				K->setPoint(0.,microstruct.getPoints()[i]->id) ;
// 				if(microstruct.getPoints()[i]->z < 1e-9 && microstruct.getPoints()[i]->t > 1e-9)
// 					K->setPoint(0,microstruct.getPoints()[i]->id) ;
			if(std::abs(microstruct.getPoints()[i]->x -7.5)< 1e-9 && microstruct.getPoints()[i]->t > 1e-9)
				K->setPoint(0.,microstruct.getPoints()[i]->id) ;
			if(microstruct.getPoints()[i]->x < 1e-9 && microstruct.getPoints()[i]->t > 1e-9)
				K->setPoint(.2,microstruct.getPoints()[i]->id) ;
		}
	}
	
	K->cgsolve() ;
	
	
	x = new Vector(K->getDisplacements()) ;

	std::cerr << " stepping through elements... " << std::flush ;
	for(size_t i = 0 ; i < microstruct.getElements().size() ;i++)
	{	
		if(i%1000 == 0)
			std::cerr << "\r stepping through elements... " << i << "/" << microstruct.getElements().size() << std::flush ;
		if(microstruct.getElements()[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			microstruct.getElements()[i]->step(.1, &K->getDisplacements()) ;
			microstruct.getElements()[i]->getBehaviour()->step(.1, microstruct.getElements()[i]->getState()) ;
		}
	}
	std::cerr << " ...done" << std::endl ;

	myTets =  microstruct.getElements() ;
	sigma = new Vector(myTets.size()*4) ;
	
	int count = 0 ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTets.size() << std::flush ;


		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			count++;
			
			for(size_t j = 4 ; j < 8 ;j++)
			{
				
				(*sigma)[i*4+j-4] = (*x)[microstruct.getElements()[i]->getBoundingPoint(j).id] ; 
			}
		}
		else
		{
			for(size_t j = 0 ; j < 4 ;j++)
			{
				(*sigma)[i*4+j] = 0 ; 
			}
		}
		
	}
	
	std::cout << " ... done" << std::endl ;
	std::cout << "sigma max = " << sigma->max() << std::endl ;
	std::cout << "sigma min = " << sigma->min() << std::endl  ;
	
	
	std::cout << " ... done" << std::endl ;
	
	
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);
	glutCreateWindow (argv[0]);
	init ();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();
	return 0;
 	
} 


