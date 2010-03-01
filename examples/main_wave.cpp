// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../solvers/assembly.h"
#include "../features/pore3d.h"
#include "../features/inclusion3d.h"
#include "../features/inclusion.h"
#include "../features/sample3d.h"
#include "../features/sample.h"
#include "../features/vibratingcircularmembrane.h"
#include "../mesher/delaunay_3d.h"
#include "../filters/voxelporefilter.h"
#include "../physics/void_form.h"
#include "../physics/wave.h"
#include "../polynomial/vm_base.h"

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
#include <fstream>
#include <valarray>

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
std::vector<DelaunayTriangle *> myTris ;
std::vector<HexahedralElement *> myHexs ;
std::vector<Point *> points ;
FeatureTree *FT ;
GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double max_x = 5.;
double maxs = 0;
double mins = 0;

Vector * x ;
Vector * sigma ;

std::vector<Point *> pts ;

int tstep = 0 ;

std::pair<Matrix,Matrix> arnoldi(const Matrix & A)
{
	std::vector<Vector> Q(A.numCols()+1,Vector(double(0),A.numRows())) ;
	Matrix H(A.numRows(),A.numRows()) ;

	// the first Arnoldi vector is { 1, 0, 0, 0, ...} 
	Q[0][0] = 1 ;
	
	for(size_t k = 1 ; k < A.numCols()+1 ; k++)
	{
		Q[k] = A*Q[k-1] ;

		for(size_t j = 0 ; j < k ; j++)
		{
			H[j][k-1] = std::inner_product(&Q[j][0], &Q[j][Q[j].size()], &Q[k][0], double(0)) ;
			Q[k] -= Q[j]*H[j][k-1] ;
		}

		H[k][k-1] = sqrt(std::inner_product(&Q[k][0], &Q[k][Q[k].size()], &Q[k][0], double(0))) ;
		Q[k]/=H[k][k-1] ;
	}

	Matrix Q_(A.numCols(), A.numRows()) ;

	for(size_t i = 0 ; i < A.numRows() ; i++)
	{
		for(size_t j = 0 ; j < A.numCols() ; j++)
		{
			Q_[i][j] = Q[i][j] ;
		}
	}

	H.print() ;
	std::cout << "\n" << std::endl ;
	Q_.print() ;
	return std::make_pair(H,Q_) ;
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
	double min = 0 ;
	double max = 0 ;
	for(size_t i = 0 ; i < myTris.size(); i++)
	{
		for(size_t j = 6 ; j < 9 ; j++)
		{
			double v = (*x)[myTris[i]->getBoundingPoint(j).id] ;
			if(min > v) min = v ;
			if(max < v) max = v ;
		}
	}

	for(size_t i = 0 ; i < myTris.size(); i++)
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
			
// 		HSV2RGB(&r, &g, &b, 360.-180.-180.*sqrt(((*sigma)[i]-min)/(max-min)), 1., 1.) ;
// 		glColor4f (r,g,b, .3);
			
		if(myTris[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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
			
// 			HSV2RGB(&r, &g, &b, 180.-180.*((*sigma)[i*10]-min)/(max-min), 1., 1.) ;
// 			glColor4f (r,g,b, .3);
			
			glBegin(GL_LINES);
			
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->first->x*.1 + 5, myTris[i]->first->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(6).id]*.1);
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(7).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->second->x*.1 + 5, myTris[i]->second->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(7).id]*.1 );

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(6).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->first->x*.1 + 5, myTris[i]->first->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(6).id]*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(8).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->third->x*.1 + 5, myTris[i]->third->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(8).id]*.1) ;

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(8).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->third->x*.1 + 5, myTris[i]->third->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(8).id]*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTris[i]->getBoundingPoint(7).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTris[i]->second->x*.1 + 5, myTris[i]->second->y*.1 + 5, (*x)[myTris[i]->getBoundingPoint(7).id]*.1);


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
	gluPerspective(45.0, (GLfloat) w/(GLfloat) h, 1., 120);
}
 
void keyboard (unsigned char key, int x_, int y_)
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
  	case 's':
	  	{
			for(size_t t = 0 ; t < 10 ; t++)
			{
				x = &(FT->getAssembly()->getDisplacements()) ;
				
				std::set<Point *> source ;
				std::set<Point *> border ;
				std::set<std::pair<Point *, Point *> > init ;
				for(size_t i = 0 ; i < myTris.size() ; i++)
				{
					for(size_t j = 0 ; j < myTris[i]->getBoundingPoints().size() ; j++)
					{
						if(std::abs(sqrt((myTris[i]->getBoundingPoint(j).x)*(myTris[i]->getBoundingPoint(j).x)+(myTris[i]->getBoundingPoint(j).y)*(myTris[i]->getBoundingPoint(j).y)) - 50) < 1e-8)
						{
							border.insert(&myTris[i]->getBoundingPoint(j)) ;
						}
						else if(myTris[i]->getBoundingPoint(j).x == 0 
							&& myTris[i]->getBoundingPoint(j).y == 0 
							&& myTris[i]->getBoundingPoint(j).t == 1)
						{
							source.insert(&myTris[i]->getBoundingPoint(j)) ;
						}
						else if(myTris[i]->getBoundingPoint(j).t == -1 || myTris[i]->getBoundingPoint(j).t == 0)
						{
							init.insert(std::make_pair(&myTris[i]->getBoundingPoint(j),
										  &myTris[i]->getBoundingPoint(j+3))) ;
						}
					}
				}
				
				for(std::set<Point *>::iterator i = source.begin() ; i != source.end() ; ++i)
				{
					FT->getAssembly()->setPoint(10.*sin(2.*M_PI*(double)(tstep+1)/100.),(*i)->id) ;
				}
// 				for(std::set<Point *>::iterator i = border.begin() ; i != border.end() ; ++i)
// 					FT->getAssembly()->setPoint(0,(*i)->id) ;
				
				for(std::set<std::pair<Point *, Point *> >::iterator i = init.begin() ; i != init.end() ; ++i)
				{
					FT->getAssembly()->setPoint((*x)[(*i).second->id],(*i).first->id) ;
				}
				
				FT->step(0.1) ;
				x =  &(FT->getAssembly()->getDisplacements()) ;
				tstep++ ;
			}
// 		  	FT->getAssembly()->print() ;
		  	computeDisplayList() ;
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

	}
}


int main(int argc, char *argv[])
{
	
	Matrix diffusionMatrix(2,2) ;

	
	diffusionMatrix[0][0] = 50;
	diffusionMatrix[1][1] = 50;

	Sample * base = new Sample(NULL, 200, 200, 0, 0) ;
	base->setBehaviour(new VoidForm()) ;
	Inclusion * sample = new Inclusion(80, 0, 0) ;
	sample->setBehaviour(new Wave(diffusionMatrix)) ;
	FeatureTree ft(base) ;
	ft.addFeature(base, sample) ;
// 	ft.addFeature(base, new VibratingMembrane(sample,50, 0, 0, diffusionMatrix)) ;


	ft.sample(512) ;
	ft.setOrder(LINEAR_TIME_QUADRATIC) ;
	ft.generateElements(0, false) ;
	

	std::vector<DelaunayTriangle *> elems = ft.getTriangles() ;

	std::set<Point *> points ;
	for(size_t i = 0 ; i < elems.size() ; i++)
	{
		
		if(i% 100 == 0)
		{
			std::cout << "\r seting BC : elem " << i << "/" << elems.size() << std::flush ;
		}
		
		for(size_t j = 0 ; j < elems[i]->getBoundingPoints().size() ; j++)
		{
			points.insert(&elems[i]->getBoundingPoint(j)) ;
		}
	}
// 
// 	for(std::set<Point *>::iterator i = points.begin() ; i != points.end() ; ++i)
// 	{
// 		
// 		ft.getAssembly()->setPoint(0., (*i)->id) ;
// 
// 	}
// 	ft.step(0.1) ;
// 
// 	x = &ft.getAssembly()->getDisplacements() ;

	
	for(std::set<Point *>::iterator i = points.begin() ; i != points.end() ; ++i)
	{
		if(std::abs(sqrt(((*i)->x)*((*i)->x)+((*i)->y)*((*i)->y)) - 50) < 1e-8)
		{
			ft.getAssembly()->setPoint(0, (*i)->id) ;
		}
		else if((*i)->x  == 0  && (*i)->y == 0 && (*i)->t == -1)
		{
			ft.getAssembly()->setPoint(100.*sin(M_PI*(double)0/10), (*i)->id) ;
		}
		else if((*i)->x  == 0  && (*i)->y == 0 && (*i)->t == 1)
		{
			ft.getAssembly()->setPoint(100.*sin(M_PI*(double)0/10), (*i)->id) ;
		}
	}

	ft.step(0.1) ;
	FT = &ft ;

// 	ft.getAssembly()->print() ;

	for(size_t timestep = 2 ; timestep < 0 ; timestep++)
	{
		x = &ft.getAssembly()->getDisplacements() ;


		std::set<Point *> source ;
		std::set<Point *> border ;
		std::set<std::pair<Point *, Point *> > init ;
		for(size_t i = 0 ; i < elems.size() ; i++)
		{
			for(size_t j = 0 ; j < elems[i]->getBoundingPoints().size() ; j++)
			{
				if(std::abs(sqrt((elems[i]->getBoundingPoint(j).x)*(elems[i]->getBoundingPoint(j).x)+(elems[i]->getBoundingPoint(j).y)*(elems[i]->getBoundingPoint(j).y)) - 50) < 1e-8)
				{
					border.insert(&elems[i]->getBoundingPoint(j)) ;
				}
				else if(elems[i]->getBoundingPoint(j).x == 0 
				   && elems[i]->getBoundingPoint(j).y == 0 
				   && elems[i]->getBoundingPoint(j).t == 1)
				{
					source.insert(&elems[i]->getBoundingPoint(j)) ;
				}
				else if(elems[i]->getBoundingPoint(j).t == -1)
				{
					init.insert(std::make_pair(&elems[i]->getBoundingPoint(j),
					                           &elems[i]->getBoundingPoint(j+3))) ;
				}
			}
		}
		
		for(std::set<Point *>::iterator i = source.begin() ; i != source.end() ; ++i)
			ft.getAssembly()->setPoint(100.*sin(M_PI*(double)timestep/10),(*i)->id) ;
		for(std::set<Point *>::iterator i = border.begin() ; i != border.end() ; ++i)
			ft.getAssembly()->setPoint(0,(*i)->id) ;
		
		for(std::set<std::pair<Point *, Point *> >::iterator i = init.begin() ; i != init.end() ; ++i)
		{
			ft.getAssembly()->setPoint((*x)[(*i).second->id],(*i).first->id) ;
		}

		ft.step(0.1) ;
	}

	x =  &ft.getAssembly()->getDisplacements() ;

// 	for(size_t i = 0 ; i < elems.size() ; i++)
// 	{
// 
// 		for(size_t j = 0 ;j < elems[i]->getBoundingPoints().size() ; j++)
// 		{
// 			std::cout << elems[i]->getBoundingPoint(j).x << "  " << elems[i]->getBoundingPoint(j).t << "  "<< (*x)[elems[i]->getBoundingPoint(j).id] << std::endl ;
// 
// 		}
// 	}

	myTris =  elems ;

	sigma = new Vector(myTris.size()*3) ;
	
	int count = 0 ;
	for(size_t i = 0 ; i < myTris.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTris.size() << std::flush ;


		if(myTris[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			count++;

			for(size_t j = 6 ; j < 9 ;j++)
			{
				
				(*sigma)[i*3+j-6] = (*x)[elems[i]->getBoundingPoint(j).id] ;
			}
		}
		else
		{
			for(size_t j = 0 ; j < 3 ;j++)
			{
				(*sigma)[i*3+j] = 0 ;
			}
		}
		
	}

// 	(*x) *= 10 ;
	
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


