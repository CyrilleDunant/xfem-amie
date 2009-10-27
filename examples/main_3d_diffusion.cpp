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
#include "../mesher/delaunay_3d.h"
#include "../filters/voxelporefilter.h"
#include "../physics/void_form.h"

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
DelaunayTree3D *dt ;
GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double max_x = 5.;
double maxs = 0;
double mins = 0;

Vector * x ;
Vector * sigma ;

std::vector<Point *> pts ;


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
			
// 		HSV2RGB(&r, &g, &b, 360.-180.-180.*sqrt(((*sigma)[i]-min)/(max-min)), 1., 1.) ;
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
			
// 			HSV2RGB(&r, &g, &b, 180.-180.*((*sigma)[i*10]-min)/(max-min), 1., 1.) ;
// 			glColor4f (r,g,b, .3);
			
			glBegin(GL_LINES);
			
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(0+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->first->x*.1, myTets[i]->first->y*.1, myTets[i]->first->z*.1);
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(2+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->second->x*.1, myTets[i]->second->y*.1, myTets[i]->second->z*.1 );

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(0+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->first->x*.1, myTets[i]->first->y*.1, myTets[i]->first->z*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(4+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->third->x*.1, myTets[i]->third->y*.1, myTets[i]->third->z*.1) ;

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(0+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->first->x*.1, myTets[i]->first->y*.1, myTets[i]->first->z*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(6+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->fourth->x*.1, myTets[i]->fourth->y*.1, myTets[i]->fourth->z*.1);

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(2+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->second->x*.1, myTets[i]->second->y*.1, myTets[i]->second->z*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(4+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->third->x*.1, myTets[i]->third->y*.1, myTets[i]->third->z*.1);

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(2+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->second->x*.1, myTets[i]->second->y*.1, myTets[i]->second->z*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(6+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->fourth->x*.1, myTets[i]->fourth->y*.1, myTets[i]->fourth->z*.1);

			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(4+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->third->x*.1, myTets[i]->third->y*.1, myTets[i]->third->z*.1 );
			HSV2RGB(&r, &g, &b, 180.-180.*((*x)[myTets[i]->getBoundingPoint(6+10).id]-min)/(max-min), 1., 1.) ;
			glColor4f (r,g,b, .3);
			glVertex3f(myTets[i]->fourth->x*.1, myTets[i]->fourth->y*.1, myTets[i]->fourth->z*.1);


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

	}
}


int main(int argc, char *argv[])
{
	
	Matrix diffusionMatrix(3,3) ;

	diffusionMatrix[0][0] =100;
	diffusionMatrix[1][1] =100;
	diffusionMatrix[2][2] =100;


	Sample3D * sample = new Sample3D(100, 100, 100, 50, 50, 50) ;
	sample->setBehaviour(new Diffusion(diffusionMatrix)) ;
	FeatureTree ft(sample) ;
	
// 	std::ifstream file("InputFiles/1000Part/porein20.csv") ;		
// 	char comma ;
// 	int dummy ;
// 	double r ;
// 	double xx ;
// 	double y ;
// 	double z ;
// 	file >> dummy >> comma >> r >> comma >> xx >> comma >> y >> comma >> z ;
// 	Pore3D * lastPore = new Pore3D(r, xx, y, z) ;
// 	/*Pore3D*/
// 	ft.addFeature(sample,lastPore) ;
// 	while(!file.eof())
// 	{
// 		file >> dummy >> comma >> r >> comma >> xx >> comma >> y >> comma >> z ;
// // 		std::cout << dummy << ", " << r << ", " << xx << ", " << y << ", " << z << std::endl ;
// 		Pore3D * newPore = new Pore3D(r, xx, y, z) ;
// 		ft.addFeature(lastPore, newPore) ;
// 		Pore3D * originalPore = newPore;
// 		/*VirtualInclusion3D*/
// 		if(!sample->in(originalPore->getCenter() + Point(originalPore->getRadius(), 0, 0)))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx-100, y, z) ;
// 			ft.addFeature(lastPore, newPore) ;
// 			lastPore = newPore ;
// 		}
// 		if(!sample->in(originalPore->getCenter() + Point(-originalPore->getRadius(), 0, 0)))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx+100, y, z) ;
// 			ft.addFeature(lastPore, newPore) ;
// 		}
// 		if(!sample->in(originalPore->getCenter() + Point(0,originalPore->getRadius() , 0)))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx, y-100, z) ;
// 			ft.addFeature(lastPore, newPore) ;
// 		}
// 		if(!sample->in(originalPore->getCenter() + Point(0,-originalPore->getRadius() , 0)))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx, y+100, z) ;
// 			ft.addFeature(lastPore, newPore) ;
// 		}
// 		if(!sample->in(originalPore->getCenter() + Point(0, 0, originalPore->getRadius())))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx, y, z-100) ;
// 			ft.addFeature(lastPore, newPore) ;
// 		}
// 		if(!sample->in(originalPore->getCenter() + Point(0,0 ,-originalPore->getRadius() )))
// 		{
// 			Pore3D * newPore = new Pore3D(r, xx, y, z+100) ;
// 			ft.addFeature(lastPore, newPore) ;
// 		}
// 	}
// 	
	ft.addFeature(sample,new Pore3D(25, 50, 50, 50)) ;
	ft.setOrder(QUADRATIC_TIME_LINEAR) ;
	ft.sample(512) ;
	ft.generateElements(0, true) ;
	std::vector<DelaunayTetrahedron *> elems = ft.getTetrahedrons() ;

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
	std::cout << " ...done." << std::endl ;
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
		if(std::abs((*i)->t  +1) < 1e-6)
		{
			ft.getAssembly()->setPoint(0, (*i)->id) ;
		}
		else if(std::abs((*i)->x ) < 1e-6 && std::abs((*i)->t  - 1) < 1e-6)
		{
			ft.getAssembly()->setPoint(.5, (*i)->id) ;
		}
	}

	ft.step(0.1) ;

	for(size_t timestep = 0 ; timestep < 50 ; timestep++)
	{
		x = &ft.getAssembly()->getDisplacements() ;


		std::set<Point *> source ;
		std::set<std::pair<Point *, Point *> > init ;
		for(size_t i = 0 ; i < elems.size() ; i++)
		{
			for(size_t j = 0 ; j < elems[i]->getBoundingPoints().size() ; j++)
			{
				if(elems[i]->getBoundingPoint(j).x == 0 )
				{
					source.insert(&elems[i]->getBoundingPoint(j)) ;
				}
				else if(elems[i]->getBoundingPoint(j).t == -1)
				{
					init.insert(std::make_pair(&elems[i]->getBoundingPoint(j),
					                           &elems[i]->getBoundingPoint(j+10))) ;
				}
			}
		}
		
		for(std::set<Point *>::iterator i = source.begin() ; i != source.end() ; ++i)
			ft.getAssembly()->setPoint(.5,(*i)->id) ;
		
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

	myTets =  elems ;

	sigma = new Vector(myTets.size()*10) ;
	
	int count = 0 ;
	for(size_t i = 0 ; i < myTets.size(); i++)
	{
		if(i%1000 == 0)
			std::cout << "\r getting strains ..." << i+1 << "/" << myTets.size() << std::flush ;


		if(myTets[i]->getBehaviour()->type != VOID_BEHAVIOUR  )
		{
			count++;

			for(size_t j = 10 ; j < 20 ;j++)
			{
				
				(*sigma)[i*10+j-10] = (*x)[elems[i]->getBoundingPoint(j).id] ;
			}
		}
		else
		{
			for(size_t j = 0 ; j < 10 ;j++)
			{
				(*sigma)[i*10+j] = 0 ;
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


