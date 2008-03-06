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
std::vector<Point * > points ;
DelaunayTree_3D *dt ;
GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double max_x = 5.;

Vector * x ;
Vector * sigma ;

std::vector<Point * > pts ;

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
		
		double radius = 0.2 + .9*(double)random()/((double)RAND_MAX+1) ;
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
		temp->getBehaviour()->param *= fact ;
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
			
		HSV2RGB(&r, &g, &b, 360.-360.*sqrt(((*sigma)[i]-min)/(max-min)), 1., 1.) ;
		glColor4f (r,g,b, .3);
			
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

			
			HSV2RGB(&r, &g, &b, 360.*(myTets[i]->getBehaviour()->param[0][0]-3.8)/(11.1-3.8), 1., 1.) ;
			glColor4f (r,g,b, .3);
			
			glBegin(GL_LINES);
			
			if(*myTets[i]->first != Point() && *myTets[i]->second != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2]);
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );
			}
			
			if(*myTets[i]->first != Point() && *myTets[i]->third != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->third->x+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2] );
			}
			
			if(*myTets[i]->first != Point() && *myTets[i]->fourth != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->first->x+(*x)[myTets[i]->first->id*3],
						myTets[i]->first->y+(*x)[myTets[i]->first->id*3+1],
						myTets[i]->first->z+(*x)[myTets[i]->first->id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );
			}
			
			if(*myTets[i]->second != Point() && *myTets[i]->third != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
	
				glVertex3f(myTets[i]->third->x+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2]);
			}
			
			if(*myTets[i]->second != Point() && *myTets[i]->fourth != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+1]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->second->x+(*x)[myTets[i]->second->id*3],
						myTets[i]->second->y+(*x)[myTets[i]->second->id*3+1],
						myTets[i]->second->z+(*x)[myTets[i]->second->id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );
			}
			
			if(*myTets[i]->third != Point() && *myTets[i]->fourth != Point())
			{
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+2]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f((myTets[i]->third->x)+(*x)[myTets[i]->third->id*3],
						myTets[i]->third->y+(*x)[myTets[i]->third->id*3+1],
						myTets[i]->third->z+(*x)[myTets[i]->third->id*3+2] );
				
				HSV2RGB(&r, &g, &b, 360.*((*sigma)[i*4+3]-min)/(max-min), 1., 1.) ;
				glColor4f (r,g,b, .3);
				
				glVertex3f(myTets[i]->fourth->x+(*x)[myTets[i]->fourth->id*3],
						myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
						myTets[i]->fourth->z+(*x)[myTets[i]->fourth->id*3+2] );
			}
			glEnd();
			
// 			glBegin(GL_LINES) ;
// 			glVertex3f(myTets[i]->first->x +(*x)[myTets[i]->first->id*3],
// 			           myTets[i]->first->y +(*x)[myTets[i]->first->id*3+1],
// 			           myTets[i]->first->z +(*x)[myTets[i]->first->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(1)->x +(*x)[myTets[i]->getBoundingPoint(1)->id*3],
// 			           myTets[i]->getBoundingPoint(1)->y +(*x)[myTets[i]->getBoundingPoint(1)->id*3+1],
// 			           myTets[i]->getBoundingPoint(1)->z +(*x)[myTets[i]->getBoundingPoint(1)->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(1)->x +(*x)[myTets[i]->getBoundingPoint(1)->id*3],
// 			           myTets[i]->getBoundingPoint(1)->y +(*x)[myTets[i]->getBoundingPoint(1)->id*3+1],
// 			           myTets[i]->getBoundingPoint(1)->z +(*x)[myTets[i]->getBoundingPoint(1)->id*3+2]);
// 			glVertex3f(myTets[i]->second->x +(*x)[myTets[i]->second->id*3],
// 			           myTets[i]->second->y +(*x)[myTets[i]->second->id*3+1],
// 			           myTets[i]->second->z +(*x)[myTets[i]->second->id*3+2]);
// 
// 			glVertex3f(myTets[i]->first->x +(*x)[myTets[i]->first->id*3],
// 			           myTets[i]->first->y +(*x)[myTets[i]->first->id*3+1],
// 			           myTets[i]->first->z +(*x)[myTets[i]->first->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(9)->x +(*x)[myTets[i]->getBoundingPoint(9)->id*3],
// 			           myTets[i]->getBoundingPoint(9)->y +(*x)[myTets[i]->getBoundingPoint(9)->id*3+1],
// 			           myTets[i]->getBoundingPoint(9)->z +(*x)[myTets[i]->getBoundingPoint(9)->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(9)->x +(*x)[myTets[i]->getBoundingPoint(9)->id*3],
// 			           myTets[i]->getBoundingPoint(9)->y +(*x)[myTets[i]->getBoundingPoint(9)->id*3+1],
// 			           myTets[i]->getBoundingPoint(9)->z +(*x)[myTets[i]->getBoundingPoint(9)->id*3+2]);
// 			glVertex3f(myTets[i]->third->x +(*x)[myTets[i]->third->id*3],
// 			           myTets[i]->third->y +(*x)[myTets[i]->third->id*3+1],
// 			           myTets[i]->third->z +(*x)[myTets[i]->third->id*3+2]);
// 			
// 			glVertex3f(myTets[i]->first->x +(*x)[myTets[i]->first->id*3],
// 			           myTets[i]->first->y +(*x)[myTets[i]->first->id*3+1],
// 			           myTets[i]->first->z +(*x)[myTets[i]->first->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(7)->x +(*x)[myTets[i]->getBoundingPoint(7)->id*3],
// 			           myTets[i]->getBoundingPoint(7)->y +(*x)[myTets[i]->getBoundingPoint(7)->id*3+1],
// 			           myTets[i]->getBoundingPoint(7)->z +(*x)[myTets[i]->getBoundingPoint(7)->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(7)->x +(*x)[myTets[i]->getBoundingPoint(7)->id*3],
// 			           myTets[i]->getBoundingPoint(7)->y +(*x)[myTets[i]->getBoundingPoint(7)->id*3+1],
// 			           myTets[i]->getBoundingPoint(7)->z +(*x)[myTets[i]->getBoundingPoint(7)->id*3+2]);
// 			glVertex3f(myTets[i]->fourth->x +(*x)[myTets[i]->fourth->id*3],
// 			           myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
// 			           myTets[i]->fourth->z +(*x)[myTets[i]->fourth->id*3+2]);
// 			glEnd();
// 			glBegin(GL_LINE_LOOP) ;
// 			glVertex3f(myTets[i]->second->x +(*x)[myTets[i]->second->id*3],
// 			           myTets[i]->second->y +(*x)[myTets[i]->second->id*3+1],
// 			           myTets[i]->second->z +(*x)[myTets[i]->second->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(3)->x +(*x)[myTets[i]->getBoundingPoint(3)->id*3],
// 			           myTets[i]->getBoundingPoint(3)->y +(*x)[myTets[i]->getBoundingPoint(3)->id*3+1],
// 			           myTets[i]->getBoundingPoint(3)->z +(*x)[myTets[i]->getBoundingPoint(3)->id*3+2]);
// 			glVertex3f(myTets[i]->third->x + (*x)[myTets[i]->third->id*3],
// 			           myTets[i]->third->y + (*x)[myTets[i]->third->id*3+1],
// 			           myTets[i]->third->z + (*x)[myTets[i]->third->id*3+2] );
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(5)->x +(*x)[myTets[i]->getBoundingPoint(5)->id*3],
// 			           myTets[i]->getBoundingPoint(5)->y +(*x)[myTets[i]->getBoundingPoint(5)->id*3+1],
// 			           myTets[i]->getBoundingPoint(5)->z +(*x)[myTets[i]->getBoundingPoint(5)->id*3+2]);
// 			glVertex3f(myTets[i]->fourth->x +(*x)[myTets[i]->fourth->id*3],
// 			           myTets[i]->fourth->y+(*x)[myTets[i]->fourth->id*3+1],
// 			           myTets[i]->fourth->z +(*x)[myTets[i]->fourth->id*3+2]);
// 			if(myTets[i]->getOrder() == QUADRATIC)
// 			glVertex3f(myTets[i]->getBoundingPoint(8)->x +(*x)[myTets[i]->getBoundingPoint(8)->id*3],
// 			           myTets[i]->getBoundingPoint(8)->y +(*x)[myTets[i]->getBoundingPoint(8)->id*3+1],
// 			           myTets[i]->getBoundingPoint(8)->z +(*x)[myTets[i]->getBoundingPoint(8)->id*3+2]);
// 			glEnd();
			
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
	
	Sample3D s3d(6,6,6,0,0,0) ;
	Inclusion3D inclusion(0.5, 1, 0,0) ;
	Inclusion3D inclusion0(0.5, -1, 0,0) ;
	Inclusion3D inclusion1(0.7, 1, 0,0) ;
	Inclusion3D inclusion2(0.7, -1, 0,0) ;
	Inclusion3D inclusion3(1, 1, 0,0) ;
	Inclusion3D inclusion4(1, -1, 0,0) ;
	Inclusion3D inclusion5(1.2, 1, 0,0) ;
	Inclusion3D inclusion6(1.2, -1, 0,0) ;
	Inclusion3D inclusion7(1.5, 1, 0,0) ;
	Inclusion3D inclusion8(1.5, -1, 0,0) ;
	Inclusion3D inclusion9(1.7, 1, 0,0) ;
	Inclusion3D inclusion10(1.7, -1, 0,0) ;
	Inclusion3D inclusion11(1.9, 1, 0,0) ;
	Inclusion3D inclusion12(1.9, -1, 0,0) ;
	Inclusion3D inclusion13(2.1, 1, 0,0) ;
	Inclusion3D inclusion14(2.1, -1, 0,0) ;
	Inclusion3D inclusion15(2.4, 1, 0,0) ;
	Inclusion3D inclusion16(2.7, -1, 0,0) ;
	Inclusion3D inclusion17(2.3, -5, 0,0) ;
	Inclusion3D inclusion18(2.7, 5, 0,0) ;
	FeatureTree ft(&s3d) ;
// 	ft.addFeature(&s3d,&inclusion) ;
	
	double E_c3s = 4. ;      // in MPa
	double nu_c3s = 0.3 ;

	Matrix cgStressC3S(6,6) ;
	cgStressC3S[0][0] = 1. - nu_c3s ; cgStressC3S[0][1] = nu_c3s ; cgStressC3S[0][2] = nu_c3s ;
	cgStressC3S[1][0] = nu_c3s ; cgStressC3S[1][1] = 1. - nu_c3s ; cgStressC3S[1][2] = nu_c3s ;
	cgStressC3S[2][0] = nu_c3s ; cgStressC3S[2][1] = nu_c3s ; cgStressC3S[2][2] = 1. - nu_c3s ;
	cgStressC3S[3][3] = 0.5 - nu_c3s ;
	cgStressC3S[4][4] = 0.5 - nu_c3s ;
	cgStressC3S[5][5] = 0.5 - nu_c3s ;
	cgStressC3S *= E_c3s/((1.+nu_c3s)*(1.-2.*nu_c3s)) ;
	
	Matrix cgStress(6,6) ;
	cgStress = cgStressC3S*2. ;
	Matrix cgStress0(6,6) ;
	cgStress0 = cgStressC3S*8. ;
	
	Stiffness sc3s(cgStressC3S) ;
	Stiffness soth(cgStress) ;
	Stiffness soth0(cgStress0) ;


	std::pair<std::vector</*Virtual*/Inclusion3D * >, std::vector<Pore3D * > > features = generateInclusionsAndPores(256, 0, &soth, &s3d, &ft) ;

// 	ft.addFeature(&s3d, &inclusion15) ;
// 	ft.addFeature(&inclusion15, &inclusion17) ;
// 	ft.addFeature(&inclusion17, &inclusion16) ;
// 	ft.addFeature(&inclusion16, &inclusion18) ;
// // // 	ft.addFeature(&s3d, &inclusion16) ;
// // // 	ft.addFeature(&inclusion16, &inclusion15) ;
// 	ft.addFeature(&inclusion18, &inclusion14) ;
// 	ft.addFeature(&inclusion14, &inclusion13) ;
// 	ft.addFeature(&inclusion13, &inclusion12) ;
// 	ft.addFeature(&inclusion12, &inclusion11) ;
// 	ft.addFeature(&inclusion11, &inclusion10) ;
// 	ft.addFeature(&inclusion10, &inclusion9) ;
// 	ft.addFeature(&inclusion9, &inclusion8) ;
// 	ft.addFeature(&inclusion8, &inclusion7) ;
// 	ft.addFeature(&inclusion7, &inclusion6) ;
// 	ft.addFeature(&inclusion6, &inclusion5) ;
// 	ft.addFeature(&inclusion5, &inclusion4) ;
// 	ft.addFeature(&inclusion4, &inclusion3) ;
// 	ft.addFeature(&inclusion3, &inclusion2) ;
// 	ft.addFeature(&inclusion2, &inclusion1) ;
// 	ft.addFeature(&inclusion1, &inclusion0) ;
// // 	ft.addFeature(&inclusion0, &inclusion) ;
// 	inclusion.setBehaviour(&sc3s) ;
// 	inclusion0.setBehaviour(&sc3s) ;
// 	inclusion1.setBehaviour(&soth) ;
// 	inclusion2.setBehaviour(&soth0) ;
// 	inclusion3.setBehaviour(&soth) ;
// 	inclusion4.setBehaviour(&soth0) ;
// 	inclusion5.setBehaviour(&soth) ;
// 	inclusion6.setBehaviour(&soth0) ;
// 	inclusion7.setBehaviour(&soth) ;
// 	inclusion8.setBehaviour(&soth0) ;
// 	inclusion9.setBehaviour(&soth) ;
// 	inclusion10.setBehaviour(&soth0) ;
// 	inclusion11.setBehaviour(&soth) ;
// 	inclusion12.setBehaviour(&soth0) ;
// 	inclusion13.setBehaviour(&soth) ;
// 	inclusion14.setBehaviour(&soth0) ;
// 	inclusion15.setBehaviour(&soth) ;
// 	inclusion16.setBehaviour(&soth0) ;
// 	inclusion17.setBehaviour(&soth) ;
// 	inclusion18.setBehaviour(&soth0) ;
// 	s3d.setBehaviour(new VoidForm()) ;
	s3d.setBehaviour(&sc3s) ;
	
	ft.sample(128) ;

	ft.setOrder(LINEAR) ;

	ft.generateElements(3) ;
// 	ft.refine(1) ;	
// 	
	
// 	ft.assemble() ;
	myTets= ft.getTetrahedrons() ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		for(size_t j = 0 ; j < myTets[i]->getBoundingPoints().size() ; j++)
		{
			if(std::abs(myTets[i]->getBoundingPoint(j).x-3.) <  .001 )
			{
				ft.getAssembly()->setPoint(0.5, 0, 0, myTets[i]->getBoundingPoint(j).id) ;
				
			}
		
			if(std::abs(myTets[i]->getBoundingPoint(j).x+3.) < .001 )
			{
				ft.getAssembly()->setPoint(-0.5, 0, 0, myTets[i]->getBoundingPoint(j).id) ;
			}
		}
	}

	ft.step(0) ;
	
	x = new Vector(ft.getAssembly()->getDisplacements()) ;//Vector(ft.getAssembly()->getDisplacements()) ;
// 	for(size_t i = 0 ; i < myTets->size(); i++)
// 	{
// 		if(i%1000 == 0)
// 			std::cout << "\r stepping through elements ..." << i+1 << "/" << myTets->size() << std::flush ;
// 		myTets[i]->step(0, x) ;
// 	}
// 	std::cout << " ... done" << std::endl ;
	sigma = new Vector(myTets.size()*4) ;
	
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTets.size() << std::flush ;

// 		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
// 		{
// 			(*sigma)[i*4] = myTets[i]->getBehaviour()->param[0][0];
// 			(*sigma)[i*4+1] =  myTets[i]->getBehaviour()->param[0][0];
// 			(*sigma)[i*4+2] =  myTets[i]->getBehaviour()->param[0][0];
// 			(*sigma)[i*4+3] =  myTets[i]->getBehaviour()->param[0][0];
// 		}
		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			Vector s = myTets[i]->getState().getStress(*myTets[i]->first) ;
			(*sigma)[i*4] = sqrt(((s[0]-s[1])*(s[0]-s[1]) + (s[2]-s[1])*(s[2]-s[1]) + (s[2]-s[0])*(s[2]-s[0]))/2.);
// 			(*sigma)[i*4] = s[3];
 			if((*sigma)[i*4] > 5.)
 				(*sigma)[i*4] = 5. ;
// 			if(s[0] < -5.)
// 				(*sigma)[i*4] = -5. ;

			s = myTets[i]->getState().getStress(*myTets[i]->second) ;
			(*sigma)[i*4+1] = sqrt(((s[0]-s[1])*(s[0]-s[1]) + (s[2]-s[1])*(s[2]-s[1]) + (s[2]-s[0])*(s[2]-s[0]))/2.);
// 			(*sigma)[i*4+1] = s[3];
 			if((*sigma)[i*4+1] > 5.)
 				(*sigma)[i*4+1] = 5. ;
// 			if(s[0] < -5.)
// 				(*sigma)[i*4+1] = -5. ;
	
			s = myTets[i]->getState().getStress(*myTets[i]->third) ;
			(*sigma)[i*4+2] = sqrt(((s[0]-s[1])*(s[0]-s[1]) + (s[2]-s[1])*(s[2]-s[1]) + (s[2]-s[0])*(s[2]-s[0]))/2.);
// 			(*sigma)[i*4+2] = s[3];
 			if((*sigma)[i*4+2] > 5.)
 				(*sigma)[i*4+2] = 5. ;
// 			if(s[0] < -5.)
// 				(*sigma)[i*4+2] = -5. ;

			s = myTets[i]->getState().getStress(*myTets[i]->fourth) ;
			(*sigma)[i*4+3] = sqrt(((s[0]-s[1])*(s[0]-s[1]) + (s[2]-s[1])*(s[2]-s[1]) + (s[2]-s[0])*(s[2]-s[0]))/2.);
// 			(*sigma)[i*4+3] = s[3];
 			if((*sigma)[i*4+3] > 5.)
 				(*sigma)[i*4+3] = 5. ;
// 			if(s[0] < -5.)
// 				(*sigma)[i*4+3] = -5. ;

		
// 			(*sigma)[i*4] = sqrt((*x)[myTets[i]->first->id*3]*(*x)[myTets[i]->first->id*3] + 
// 			                     (*x)[myTets[i]->first->id*3+1]*(*x)[myTets[i]->first->id*3+1]+ 
// 			                          (*x)[myTets[i]->first->id*3+2]*(*x)[myTets[i]->first->id*3+2]);
// 			
// 			(*sigma)[i*4+1] = sqrt((*x)[myTets[i]->first->id*3]*(*x)[myTets[i]->second->id*3] + 
// 			                       (*x)[myTets[i]->first->id*3+1]*(*x)[myTets[i]->second->id*3+1]+ 
// 			                            (*x)[myTets[i]->first->id*3+2]*(*x)[myTets[i]->second->id*3+2]);
// 		
// 			(*sigma)[i*4+2] = sqrt((*x)[myTets[i]->first->id*3]*(*x)[myTets[i]->third->id*3] + 
// 			                       (*x)[myTets[i]->first->id*3+1]*(*x)[myTets[i]->third->id*3+1]+ 
// 			                            (*x)[myTets[i]->first->id*3+2]*(*x)[myTets[i]->third->id*3+2]);
// 		
// 			(*sigma)[i*4+3] = sqrt((*x)[myTets[i]->first->id*3]*(*x)[myTets[i]->fourth->id*3] + 
// 			                       (*x)[myTets[i]->first->id*3+1]*(*x)[myTets[i]->fourth->id*3+1]+ 
// 			                            (*x)[myTets[i]->first->id*3+2]*(*x)[myTets[i]->fourth->id*3+2]);
		}
		else
		{
			(*sigma)[i*4] = 0 ;
			(*sigma)[i*4+1] = 0 ;
			(*sigma)[i*4+2] = 0 ;
			(*sigma)[i*4+3] = 0 ;
		}
		
	}
	
	std::cout << " ... done" << std::endl ;
	std::cout << "sigma max = " << sigma->max() << std::endl ;
	std::cout << "sigma min = " << sigma->min() << std::endl ;
	
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


