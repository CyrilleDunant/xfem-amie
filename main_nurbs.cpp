//
// C++ Implementation: main_nurbs
//
// Description: 
//
//
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2006
//
//
//
#include<iostream>
#include<cmath>
#include "geometry/geometry_base.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

using namespace Mu ;

Nurb *N ;
std::vector<Point> *sampleNurb ;
std::valarray<Point> *insertedPoint;

void init(void)
{
	glShadeModel (GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL) ;
	float mat_specular[] = {0.01f, 0.01f, 0.01f, 1.0f};
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
	
// 	computeDisplayList() ;
	
}

void display(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
	glFlush();
	glutSwapBuffers();
	glMatrixMode(GL_PROJECTION) ;
	glLoadIdentity() ;
	glOrtho(-4.5, 4.5, -4.5, 4.5, -4.5, 4.5);
		
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity() ;

	glBegin(GL_LINE_LOOP) ;
//	for(size_t i = 0 ; i < sampledPoints->size() ; i++ )
//		glVertex2f((*sampledPoints)[i].x, (*sampledPoints)[i].y) ;
	glColor3f(1.,1.,1.) ;
	for(size_t i = 0 ; i < sampleNurb->size() ; i++ )
	{
		glVertex2f((*sampleNurb)[i].x, (*sampleNurb)[i].y) ;
		std::cout << (*sampleNurb)[i].x << ", " << (*sampleNurb)[i].y << std::endl ;
	}
	glEnd() ;
	glutSwapBuffers();
}

void reshape (int w, int h)
{
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(45.0, (GLfloat) w/(GLfloat) h, 1., 40);
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


int main(/*int argc, char *argv[]*/)
{
// 	const char* outputFileName = "nurb.out";
// 	ofstream myOut(outputFileName);
	
	size_t degree;
	std::vector<Point> controlPoint;
	std::vector<long double> weight;
	std::vector<long double> knot;

	// enter control points
	controlPoint.push_back(Point(1,0));
	controlPoint.push_back(Point(1,1));
	controlPoint.push_back(Point(0,1));
	controlPoint.push_back(Point(-1,1));
	controlPoint.push_back(Point(-1,0));
	controlPoint.push_back(Point(-1,-1));
	controlPoint.push_back(Point(0,-1));
	controlPoint.push_back(Point(1,-1));
	controlPoint.push_back(Point(1,0));

	// enter weight vector
	weight.push_back(1);
	weight.push_back(sqrt(2)/2);
	weight.push_back(1);
	weight.push_back(sqrt(2)/2);
	weight.push_back(1);
	weight.push_back(sqrt(2)/2);
	weight.push_back(1);
	weight.push_back(sqrt(2)/2);
	weight.push_back(1);

	// enter knot vector
	knot.push_back(0);
	knot.push_back(0);
	knot.push_back(0);
	knot.push_back(0.25);
	knot.push_back(0.25);
	knot.push_back(0.5);
	knot.push_back(0.5);
	knot.push_back(0.75);
	knot.push_back(0.75);
	knot.push_back(1);
	knot.push_back(1);
	knot.push_back(1);

	degree = 2;

// 	// enter control points
// 	controlPoint.push_back(Point(-7.0/2,-7.6/2));
// 	controlPoint.push_back(Point(-7.0/2,7.5/2));
// 	controlPoint.push_back(Point(7.4/2,7.5/2));
// 	controlPoint.push_back(Point(7.4/2,-7.7/2));
// 	controlPoint.push_back(Point(-4.0/2,-7.6/2));
// 	
// 
// 	// enter weight vector
// 	weight.push_back(1);
// 	weight.push_back(0.5);
// 	weight.push_back(4);
// 	weight.push_back(5);
// 	weight.push_back(1);
// 	
// 
// 	// enter knot vector
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0.5);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);

// 	// enter knot vector - testing insertion
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0.2);
// 	knot.push_back(0.4);
// 	knot.push_back(0.6);
// 	knot.push_back(0.8);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);
	
// 	degree = 3;

// 	// enter knot vector - testing insertion
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0.125);
// 	knot.push_back(0.25);
// 	knot.push_back(0.375);
// 	knot.push_back(0.5);
// 	knot.push_back(0.625);
// 	knot.push_back(0.75);
// 	knot.push_back(0.875);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);

// //	enter knot vector - testing subdivision
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0);
// 	knot.push_back(0.25);
// 	knot.push_back(0.5);
// 	knot.push_back(0.75);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);
// 	knot.push_back(1);

//  	degree = 4;

// // 	Nurb N(knot,degree);
	Nurb N(controlPoint,weight,knot,degree);
	N.pointOnNurb(0.8);
	N.getNurbPoint(0.8);
	WeightedPoint P;
	P.print();
// 	N = new Nurb(controlPoint,weight,knot);
	//sampleNurb = new std::vector<Point> (N->sampleNurb(100));
// 	/insertedPoint = new std::valarray<Point>  (N->insertedPoint(0.5));


 	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);
	glutCreateWindow (argv[0]);
	init ();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
//	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();


}
