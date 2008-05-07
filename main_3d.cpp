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
#include "filters/voxelfilter.h"
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

std::pair<std::vector</*Virtual*/Inclusion3D * >, std::vector<Pore3D * > > generateInclusionsAndPores(size_t n, double fraction, Form * behaviour, Feature * father, FeatureTree * F)
{
// 	srandom(time(NULL)) ;
	size_t nombre_de_pores = static_cast<size_t>(round(n*fraction)) ;
	size_t nombre_d_inclusions = static_cast<size_t>(round(n*(1. - fraction))) ;
	
	std::pair<std::vector</*Virtual*/Inclusion3D * >, std::vector<Pore3D * > > ret ;
	ret.first = std::vector</*Virtual*/Inclusion3D * >() ;
	ret.second = std::vector<Pore3D * >() ;
	
	std::vector<Sphere *> cercles ;
	for(size_t j =0 ; j < n ; j++)
	{
		
		double radius = 0.5 + .9*(double)random()/((double)RAND_MAX+1) ;
		radius*=radius ;
		
		Point center = Point(
		                      (2.*(double)random()/((double)RAND_MAX+1)-1.)*(3.-2.*radius-0.05),
		                      (2.*(double)random()/((double)RAND_MAX+1)-1.)*(3.-2.*radius-0.05),
		                      (2.*(double)random()/((double)RAND_MAX+1)-1.)*(3.-2.*radius-0.05)
		                    ); 
		bool alone  = true ;
		
		for(size_t k = 0 ; k < cercles.size() ; k++ )
		{
			if (dist(center, cercles[k]->getCenter()) < (radius+cercles[k]->getRadius()+0.02))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
		{
			cercles.push_back(new Sphere(radius, center)) ;
		}
		else
			j-- ;
		
	}
	
	Feature * current = father ;
	for(size_t j =0 ; j < nombre_d_inclusions ; j++)
	{
		/*Virtual*/Inclusion3D * temp = new /*Virtual*/Inclusion3D(cercles[j]->getRadius(), cercles[j]->getCenter()) ;
		ret.first.push_back(temp) ;
		temp->setBehaviour(behaviour->getCopy()) ;
		double fact = (2.*random()/RAND_MAX-1.)*0.5+1. ;
		temp->getBehaviour(Point())->param *= fact ;
		F->addFeature(current, temp) ;
		current = temp ;
	}
	
	for(size_t j =0 ; j < nombre_de_pores ; j++)
	{
		Pore3D * temp = new Pore3D(cercles[j+nombre_d_inclusions]->getRadius(), cercles[j+nombre_d_inclusions]->getCenter()) ;
		ret.second.push_back(temp) ;
		F->addFeature(current, temp) ;
		current = temp ;
	}
	
	for(size_t k = 0 ; k < cercles.size() ; k++ )
	{
		delete cercles[k] ;
	}
	
	return ret ;
}



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
			

// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2]);
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );

// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2] );
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->third->x+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2] );

// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2] );
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );
	
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->third->x+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2]);

// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );

// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f((myTets[i]->third->x)+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2] );
				
// 				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
// 				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );

			glEnd();

	}
	}
	
	for(size_t i = 0 ; i < myHexs.size(); i++)
	{

		double r, g, b ;
		
		HSV2RGB(&r, &g, &b, 360.-360.*sqrt(((*sigma)[i]-min)/(max-min)), 1., 1.) ;
		glColor4f (r,g,b, .3);
		
		if(myHexs[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
// 			glColor3f (0,0,0);
			
// 			if(!s.in(myTets[i]->first)  ||
// 			   !s.in(myTets[i]->second) ||
// 			   !s.in(myTets[i]->third) ||
// 			   !s.in(myTets[i]->fourth) ||
// 			   myTets[i]->Tetrahedron::volume() < 1e-6 
// 			  )
// 					glColor4ub (255,0,0,255); 		
			
			
			
			glBegin(GL_LINE_LOOP);
			

				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(0).x+(*x)[myHexs[i]->getBoundingPoint(0).id*3],
			           myHexs[i]->getBoundingPoint(0).y+(*x)[myHexs[i]->getBoundingPoint(0).id*3+1],
			           myHexs[i]->getBoundingPoint(0).z+(*x)[myHexs[i]->getBoundingPoint(0).id*3+2]);
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(1).x+(*x)[myHexs[i]->getBoundingPoint(1).id*3],
			           myHexs[i]->getBoundingPoint(1).y+(*x)[myHexs[i]->getBoundingPoint(1).id*3+1],
			           myHexs[i]->getBoundingPoint(1).z+(*x)[myHexs[i]->getBoundingPoint(1).id*3+2] );

				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+2]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(3).x+(*x)[myHexs[i]->getBoundingPoint(3).id*3],
			           myHexs[i]->getBoundingPoint(3).y+(*x)[myHexs[i]->getBoundingPoint(3).id*3+1],
			           myHexs[i]->getBoundingPoint(3).z+(*x)[myHexs[i]->getBoundingPoint(3).id*3+2] );

				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(2).x+(*x)[myHexs[i]->getBoundingPoint(2).id*3],
			           myHexs[i]->getBoundingPoint(2).y+(*x)[myHexs[i]->getBoundingPoint(2).id*3+1],
			           myHexs[i]->getBoundingPoint(2).z+(*x)[myHexs[i]->getBoundingPoint(2).id*3+2] );
			glEnd();
			glBegin(GL_LINE_LOOP);
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+3]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(4).x+(*x)[myHexs[i]->getBoundingPoint(4).id*3],
			           myHexs[i]->getBoundingPoint(4).y+(*x)[myHexs[i]->getBoundingPoint(4).id*3+1],
			           myHexs[i]->getBoundingPoint(4).z+(*x)[myHexs[i]->getBoundingPoint(4).id*3+2] );

				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(5).x+(*x)[myHexs[i]->getBoundingPoint(5).id*3],
			           myHexs[i]->getBoundingPoint(5).y+(*x)[myHexs[i]->getBoundingPoint(5).id*3+1],
			           myHexs[i]->getBoundingPoint(5).z+(*x)[myHexs[i]->getBoundingPoint(5).id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+2]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(7).x+(*x)[myHexs[i]->getBoundingPoint(7).id*3],
			           myHexs[i]->getBoundingPoint(7).y+(*x)[myHexs[i]->getBoundingPoint(7).id*3+1],
			           myHexs[i]->getBoundingPoint(7).z+(*x)[myHexs[i]->getBoundingPoint(7).id*3+2] );

				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
			glVertex3f(myHexs[i]->getBoundingPoint(6).x+(*x)[myHexs[i]->getBoundingPoint(6).id*3],
			           myHexs[i]->getBoundingPoint(6).y+(*x)[myHexs[i]->getBoundingPoint(6).id*3+1],
			           myHexs[i]->getBoundingPoint(6).z+(*x)[myHexs[i]->getBoundingPoint(6).id*3+2] );
				
			glEnd();
			glBegin(GL_LINES);
			for(size_t t = 0 ; t < 4 ; t++)
			{
			
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*8+3]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myHexs[i]->getBoundingPoint(t).x+(*x)[myHexs[i]->getBoundingPoint(t).id*3],
			           myHexs[i]->getBoundingPoint(t).y+(*x)[myHexs[i]->getBoundingPoint(t).id*3+1],
			           myHexs[i]->getBoundingPoint(t).z+(*x)[myHexs[i]->getBoundingPoint(t).id*3+2] );
				
				glVertex3f(myHexs[i]->getBoundingPoint(t+4).x+(*x)[myHexs[i]->getBoundingPoint(t+4).id*3],
				           myHexs[i]->getBoundingPoint(t+4).y+(*x)[myHexs[i]->getBoundingPoint(t+4).id*3+1],
				           myHexs[i]->getBoundingPoint(t+4).z+(*x)[myHexs[i]->getBoundingPoint(t+4).id*3+2] );
			}
			
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
	
	double E_csh = 31 ;      // in GPa
	double nu_csh = 0.28 ;
	
	double E_ch = 40 ;      // in GPa
	double nu_ch = 0.3 ;
	
	double E_c3s = 135 ;      // in GPa
	double nu_c3s = 0.31 ;
	
	
	Matrix cgStressC3S(6,6) ;
	cgStressC3S[0][0] = 1. - nu_c3s ; cgStressC3S[0][1] = nu_c3s ; cgStressC3S[0][2] = nu_c3s ;
	cgStressC3S[1][0] = nu_c3s ; cgStressC3S[1][1] = 1. - nu_c3s ; cgStressC3S[1][2] = nu_c3s ;
	cgStressC3S[2][0] = nu_c3s ; cgStressC3S[2][1] = nu_c3s ; cgStressC3S[2][2] = 1. - nu_c3s ;
	cgStressC3S[3][3] = 0.5 - nu_c3s ;
	cgStressC3S[4][4] = 0.5 - nu_c3s ;
	cgStressC3S[5][5] = 0.5 - nu_c3s ;
	cgStressC3S *= E_c3s/((1.+nu_c3s)*(1.-2.*nu_c3s)) ;
	
	Matrix cgStressCSH(6,6) ;
	cgStressCSH[0][0] = 1. - nu_csh ; cgStressCSH[0][1] = nu_csh ; cgStressCSH[0][2] = nu_csh ;
	cgStressCSH[1][0] = nu_csh ; cgStressCSH[1][1] = 1. - nu_csh ; cgStressCSH[1][2] = nu_csh ;
	cgStressCSH[2][0] = nu_csh ; cgStressCSH[2][1] = nu_csh ; cgStressCSH[2][2] = 1. - nu_csh ;
	cgStressCSH[3][3] = 0.5 - nu_csh ;
	cgStressCSH[4][4] = 0.5 - nu_csh ;
	cgStressCSH[5][5] = 0.5 - nu_csh ;
	cgStressCSH *= E_csh/((1.+nu_csh)*(1.-2.*nu_csh)) ;
	
	Matrix cgStressCH(6,6) ;
	cgStressCH[0][0] = 1. - nu_ch ; cgStressCH[0][1] = nu_ch ; cgStressCH[0][2] = nu_ch ;
	cgStressCH[1][0] = nu_ch ; cgStressCH[1][1] = 1. - nu_ch ; cgStressCH[1][2] = nu_ch ;
	cgStressCH[2][0] = nu_ch ; cgStressCH[2][1] = nu_ch ; cgStressCH[2][2] = 1. - nu_ch ;
	cgStressCH[3][3] = 0.5 - nu_ch ;
	cgStressCH[4][4] = 0.5 - nu_ch ;
	cgStressCH[5][5] = 0.5 - nu_ch ;
	cgStressCH *= E_ch/((1.+nu_ch)*(1.-2.*nu_ch)) ;
	
	maxs = std::max(std::max(cgStressC3S[0][0],cgStressCSH[0][0]),cgStressCH[0][0]) ;
	mins = std::min(std::min(cgStressC3S[0][0],cgStressCSH[0][0]),cgStressCH[0][0]) ;
	
	Assembly * K = new Assembly() ;
	
	
	//1 Alite
	//2 C-S-H
	//3 CH
	//4 C-S-H
	//5   "
	
	VoxelFilter microstruct ;
	
	microstruct.behaviourMap[0] = new VoidForm() ;
	microstruct.behaviourMap[1] = new Stiffness(cgStressC3S) ;
	microstruct.behaviourMap[2] = new Stiffness(cgStressCSH) ;
	microstruct.behaviourMap[3] = new Stiffness(cgStressCH) ;
	microstruct.behaviourMap[4] = new Stiffness(cgStressCSH) ;
	microstruct.behaviourMap[5] = new Stiffness(cgStressCSH) ;
	
	microstruct.read("/home/cyrille/Documents/Data/Microstructures/size64/pixels26.txt") ;
	std::cout << "reading done" << std::endl ;
	
	
	for(size_t i = 0 ; i < microstruct.getElements().size() ; i++)
	{
		for(size_t j = 0 ;j < microstruct.getElements()[i]->getBoundingPoints().size() ; j++)
		{
			microstruct.getElements()[i]->getBoundingPoint(j).id = -1 ;
		}
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
			microstruct.getElements()[i]->getState().initialize() ;
		}
	}
	
	std::cout << "adding done" << std::endl ;
	
	std::cout << "setting Boundary Conditions" << std::endl ;
	for(size_t i = 0 ; i < microstruct.getElements().size() ; i++)	
	{
		if(microstruct.getElements()[i]->getBehaviour()->type != VOID_BEHAVIOUR)
		{
			if((std::abs(microstruct.getElements()[i]->getCenter().x) > 7 
			    || std::abs(microstruct.getElements()[i]->getCenter().x) < .5)
			   || (std::abs(microstruct.getElements()[i]->getCenter().y) > 7 
			       || std::abs(microstruct.getElements()[i]->getCenter().y) < .5)
			   || (std::abs(microstruct.getElements()[i]->getCenter().z) > 7 
			       || std::abs(microstruct.getElements()[i]->getCenter().z) < .5)
			)
				for(size_t j = 0 ; j < microstruct.getElements()[i]->getBoundingPoints().size()  ; j++)
				{
					if(microstruct.getElements()[i]->getBoundingPoint(j).x < 1e-9)
						K->setPoint(0,0,0,microstruct.getElements()[i]->getBoundingPoint(j).id) ;
					if(std::abs(microstruct.getElements()[i]->getBoundingPoint(j).x- 7.5) < 1e-9)
						K->setPoint(.2,0,0,microstruct.getElements()[i]->getBoundingPoint(j).id) ;
					if(microstruct.getElements()[i]->getBoundingPoint(j).y < 1e-9 
					|| std::abs(microstruct.getElements()[i]->getBoundingPoint(j).y-7.5) < 1e-9)
						K->setPointAlong(ETA,0,microstruct.getElements()[i]->getBoundingPoint(j).id) ;
					if(microstruct.getElements()[i]->getBoundingPoint(j).z < 1e-9 
					|| std::abs(microstruct.getElements()[i]->getBoundingPoint(j).z-7.5) < 1e-9)
						K->setPointAlong(ZETA,0,microstruct.getElements()[i]->getBoundingPoint(j).id) ;
				}
		}
	}
	std::cout << "BC done" << std::endl ;
	
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
	
	Vector avg_sigma(0.,6);
	Vector avg_epsilon(0.,6);
	int count = 0 ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTets.size() << std::flush ;


		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			Vector s = myTets[i]->getState().getStress(Point(0.25,0.25,0.25), true) ;
			Vector e = myTets[i]->getState().getStrain(Point(0.25,0.25,0.25), true) ;
			avg_sigma += s ;
			avg_epsilon += e ;
			count++;
			
			for(size_t j = 0 ; j < 4 ;j++)
			{
				
				(*sigma)[i*4+j] = s[0] ; 
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
	
	avg_sigma /= count ;
	avg_epsilon /= count ;
	std::cout << " ... done" << std::endl ;
	std::cout << "sigma max = " << sigma->max() << std::endl ;
	std::cout << "sigma min = " << sigma->min() << std::endl ;
	std::cout << "sigma avg = " << avg_sigma[0] <<", " <<
		avg_sigma[1] <<", " << 
		avg_sigma[2] <<", " << 
		avg_sigma[3] <<", " <<
		avg_sigma[4] <<", " <<
		avg_sigma[5] <<", " <<std::endl ;
	std::cout << "epsilon avg = " << avg_epsilon[0] <<", " <<
		avg_epsilon[1] <<", " << 
		avg_epsilon[2] <<", " << 
		avg_epsilon[3] <<", " <<
		avg_epsilon[4] <<", " <<
		avg_epsilon[5] <<", " <<std::endl ;
	
	
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


