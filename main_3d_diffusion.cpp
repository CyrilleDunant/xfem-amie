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
#include "features/inclusion.h"
#include "features/sample3d.h"
#include "features/sample.h"
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
std::vector<DelaunayTriangle *> myTets ;
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
			
			HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*6]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glBegin(GL_LINES);
			
			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(0+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->first->x, myTets[i]->first->y, (*x)[myTets[i]->getBoundingPoint(0+6).id]);

			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(2+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->second->x, myTets[i]->second->y, (*x)[myTets[i]->getBoundingPoint(2+6).id] );

			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(0+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->first->x, myTets[i]->first->y, (*x)[myTets[i]->getBoundingPoint(0+6).id] );

			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(4+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->third->x, myTets[i]->third->y, (*x)[myTets[i]->getBoundingPoint(4+6).id]);

// 				glVertex3f(myTets[i]->first->x, myTets[i]->first->y, myTets[i]->first->z );
// 				
// 				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z );

			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(2+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->second->x, myTets[i]->second->y, (*x)[myTets[i]->getBoundingPoint(2+6).id] );

			HSV2RGB(&r, &g, &b, 360.*((*x)[myTets[i]->getBoundingPoint(4+6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glVertex3f(myTets[i]->third->x, myTets[i]->third->y, (*x)[myTets[i]->getBoundingPoint(4+6).id]);

// 				glVertex3f(myTets[i]->second->x, myTets[i]->second->y, myTets[i]->second->z );
// 
// 				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z);
// 				
// 				glVertex3f(myTets[i]->third->x, myTets[i]->third->y, myTets[i]->third->z);
// 
// 				glVertex3f(myTets[i]->fourth->x, myTets[i]->fourth->y, myTets[i]->fourth->z);

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
	gluLookAt(4, 0, 4, 4, 1, 0.1, 0, 0, 1) ;
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

	Matrix diffusionMatrix(2,2) ;
	
	diffusionMatrix[0][0] =10;
	diffusionMatrix[1][1] =10;

	//1 Alite
	//2 C-S-H
	//3 CH
	//4 C-S-H
	//5   "
	

	Sample * sample = new Sample(8, 2, 4,1) ;
	sample->setBehaviour(new Diffusion(diffusionMatrix)) ;

	Inclusion * inc = new Inclusion(.6, 4, 1) ;
	inc->setBehaviour(new Diffusion(diffusionMatrix*0.01)) ;
	
	FeatureTree ft(sample) ;
	ft.addFeature(sample,inc) ;

	ft.setOrder(QUADRATIC_TIME_LINEAR) ;

	ft.sample(64) ;
	
	ft.generateElements() ;

	std::vector<DelaunayTriangle *> elems = ft.getTriangles() ;

	
	for(size_t i = 0 ; i < elems.size() ; i++)
	{
		if(i% 100 == 0)
		{
			std::cout << "\r seting BC : elem " << i << "/" << elems.size() << std::flush ;
		}
		
		for(size_t j = 0 ; j < elems[i]->getBoundingPoints().size() ; j++)
		{
			ft.getAssembly()->setPoint(0., elems[i]->getBoundingPoint(j).id) ;
		}
	}

	ft.step(0.1) ;
	std::vector<std::vector<Matrix> > elementaryMatrix = elems[0]->getElementaryMatrix() ;
	
	x = &ft.getAssembly()->getDisplacements() ;

	for(size_t i = 0 ; i < elems.size() ; i++)
	{
		if(i% 100 == 0)
		{
			std::cout << "\r seting BC : elem " << i << "/" << elems.size() << std::flush ;
		}
		
		for(size_t j = 0 ;j < elems[i]->getBoundingPoints().size() ; j++)
		{
			if(elems[i]->getBoundingPoint(j).x == 8  )
			{
				ft.getAssembly()->setPoint(0,elems[i]->getBoundingPoint(j).id) ;
			}
			
			
			if(elems[i]->getBoundingPoint(j).t == -1 )
			{
				ft.getAssembly()->setPoint(0,elems[i]->getBoundingPoint(j).id) ;
			}
			else if(elems[i]->getBoundingPoint(j).x  == 0  && elems[i]->getBoundingPoint(j).t == 1)
			{
				ft.getAssembly()->setPoint(.5,elems[i]->getBoundingPoint(j).id) ;
			}

		}
		
	}

	ft.step(0.1) ;

	x = &ft.getAssembly()->getDisplacements() ;
	
	for(size_t i = 0 ; i < elems.size() ; i++)
	{
		if(i% 100 == 0)
		{
			std::cout << "\r seting BC : elem " << i << "/" << elems.size() << std::flush ;
		}
		
		for(size_t j = 0 ;j < elems[i]->getBoundingPoints().size() ; j++)
		{

			if(elems[i]->getBoundingPoint(j).x == 8  )
			{
				ft.getAssembly()->setPoint(0,elems[i]->getBoundingPoint(j).id) ;
			}
			
			if( elems[i]->getBoundingPoint(j).x == 0)
			{
				ft.getAssembly()->setPoint(.5,elems[i]->getBoundingPoint(j).id) ;
			}
			else if(elems[i]->getBoundingPoint(j).t == -1 )
			{
				ft.getAssembly()->setPoint((*x)[elems[i]->getBoundingPoint(j+6).id],elems[i]->getBoundingPoint(j).id) ;
			}

			
		}
		
	}

	
	ft.step(0.1) ;

	for(size_t timestep = 0 ; timestep < 5 ; timestep++)
	{
		x = &ft.getAssembly()->getDisplacements() ;
		
		for(size_t i = 0 ; i < elems.size() ; i++)
		{
			if(i% 100 == 0)
			{
				std::cout << "\r seting BC : elem " << i << "/" << elems.size() << std::flush ;
			}
			
			for(size_t j = 0 ;j < elems[i]->getBoundingPoints().size() ; j++)
			{

				if(elems[i]->getBoundingPoint(j).x == 8  )
				{
					ft.getAssembly()->setPoint(0,elems[i]->getBoundingPoint(j).id) ;
				}
				
				if(elems[i]->getBoundingPoint(j).x == 0 )
				{
					ft.getAssembly()->setPoint(.5,elems[i]->getBoundingPoint(j).id) ;
				}
				else if(elems[i]->getBoundingPoint(j).t == -1|| elems[i]->getBoundingPoint(j).t == 0)
				{
					ft.getAssembly()->setPoint((*x)[elems[i]->getBoundingPoint(j+6).id],elems[i]->getBoundingPoint(j).id) ;
				}
			}
			
		}
		
		
		ft.step(0.1) ;
	}
	
	x =  &ft.getAssembly()->getDisplacements() ;
	
	for(size_t i = 0 ; i < elems.size() ; i++)
	{

		for(size_t j = 0 ;j < elems[i]->getBoundingPoints().size() ; j++)
		{
			std::cout << elems[i]->getBoundingPoint(j).x << "  " << elems[i]->getBoundingPoint(j).t << "  "<< (*x)[elems[i]->getBoundingPoint(j).id] << std::endl ;
			
		}
	}



	myTets =  elems ;
	sigma = new Vector(myTets.size()*6) ;
	
	int count = 0 ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTets.size() << std::flush ;


		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			count++;

// 			for(size_t j = 0 ; j < 10 ;j++)
			for(size_t j = 6 ; j < 12 ;j++)
			{
				
				(*sigma)[i*6+j-6] = (*x)[elems[i]->getBoundingPoint(j).id] ;
			}
		}
		else
		{
			for(size_t j = 0 ; j < 6 ;j++)
			{
				(*sigma)[i*6+j] = 0 ;
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


