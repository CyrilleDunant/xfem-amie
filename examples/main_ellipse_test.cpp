// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../physics/mohrcoulomb.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_2D.h"

#include <fstream>

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


using namespace Mu ;

/*FeatureTree * featureTree ;
std::vector<DelaunayTetrahedron *> triangles ;
DelaunayTree3D *dt ; //(pts) ;

double timepos = 0.1 ;
bool firstRun = true ;


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

double nu = 0.3 ;
double E = 20 ;

bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ; 



void setBC()
{
	triangles = featureTree->getTetrahedrons() ;
	
	for(size_t k = 0 ; k < triangles.size() ;k++)
	{
		for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
		{
			
			if (triangles[k]->getBoundingPoint(c).x < .00001)
			{
				featureTree->getAssembly()->setPoint(-timepos,0 ,0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).x > .03999)
			{
				featureTree->getAssembly()->setPoint( timepos,0, 0,triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).y < .00001)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).y > .03999)
			{
				featureTree->getAssembly()->setPointAlong( ETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if (triangles[k]->getBoundingPoint(c).z < .00001)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
			if(triangles[k]->getBoundingPoint(c).z > .03999)
			{
				featureTree->getAssembly()->setPointAlong( ZETA, 0, triangles[k]->getBoundingPoint(c).id) ;
			}
		}
	}
	
}


void step()
{
	
	for(size_t i = 0 ; i < 1 ; i++)
	{
		std::cout << "\r iteration " << i << "/2" << std::flush ;
		setBC() ;
		featureTree->step(timepos) ;
		
// 		timepos+= 0.01 ;
	}
	std::cout << "\r iteration " << "2/2 ... done" << std::endl ;
	x.resize(featureTree->getDelaunayTree3D()->numPoints()*3) ;
	x = featureTree->getDisplacements() ;
	dt = featureTree->getDelaunayTree3D() ;
	sigma.resize(triangles.size()*6*4) ;

	epsilon.resize(triangles.size()*6*4) ;

// 	sigma = F.strainFromDisplacements() ;
// 	epsilon = F.stressFromDisplacements() ;
	std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
	sigma.resize(sigma_epsilon.first.size()) ;
	sigma = sigma_epsilon.first ;
	epsilon.resize(sigma_epsilon.second.size()) ;
	epsilon = sigma_epsilon.second ;

	sigma11.resize(sigma.size()/6) ;
	sigma22.resize(sigma.size()/6) ;
	sigma12.resize(sigma.size()/6) ;
	epsilon11.resize(sigma.size()/6) ;
	epsilon22.resize(sigma.size()/6) ;
	epsilon12.resize(sigma.size()/6) ;
	vonMises.resize(sigma.size()/6) ;

	std::cout << "unknowns :" << x.size() << std::endl ;
	std::cout << "max value :" << x.max() << std::endl ;
	std::cout << "min value :" << x.min() << std::endl ;


	for(size_t k = 0 ; k < triangles.size() ; k++)
	{

		if(!triangles[k]->getBehaviour()->fractured())
		{
			sigma11[k*4] = sigma[k*6*4];
			sigma22[k*4] = sigma[k*6*4+1];
			sigma12[k*4] = sigma[k*6*4+3];
			sigma11[k*4+1] = sigma[k*6*4+6];
			sigma22[k*4+1] = sigma[k*6*4+7];
			sigma12[k*4+1] = sigma[k*6*4+9];
			sigma11[k*4+2] = sigma[k*6*4+12];
			sigma22[k*4+2] = sigma[k*6*4+13];
			sigma12[k*4+2] = sigma[k*6*4+15];

			epsilon11[k*4] = epsilon[k*6*4];
			epsilon22[k*4] = epsilon[k*6*4+1];
			epsilon12[k*4] = epsilon[k*6*4+3];
			epsilon11[k*3+1] = epsilon[k*6*4+6];
			epsilon22[k*4+1] = epsilon[k*6*4+7];
			epsilon12[k*4+1] = epsilon[k*6*4+8];
			epsilon11[k*4+2] = epsilon[k*6*4+12];
			epsilon22[k*4+2] = epsilon[k*6*4+13];
			epsilon12[k*4+2] = epsilon[k*6*4+15];

			Vector vm0 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->first) ;
			Vector vm1 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->second) ;
			Vector vm2 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->third) ;
			vonMises[k*4]   = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
			vonMises[k*4+1] = sqrt(((vm1[0]-vm1[1])*(vm1[0]-vm1[1]))/2.) ;
			vonMises[k*4+2] = sqrt(((vm2[0]-vm2[1])*(vm2[0]-vm2[1]))/2.) ;
			vonMises[k*4+3] = sqrt(((vm2[0]-vm2[1])*(vm2[0]-vm2[1]))/2.) ;

// 			vonMises[k*3]   = sigma11[k*3]*epsilon11[k*3]     + sigma22[k*3]*epsilon22[k*3]     + sigma12[k*3]*epsilon12[k*3];
// 			vonMises[k*3+1] = sigma11[k*3+1]*epsilon11[k*3+1] + sigma22[k*3+1]*epsilon22[k*3+1] + sigma12[k*3+1]*epsilon12[k*3+2];
// 			vonMises[k*3+2] = sigma11[k*3+2]*epsilon11[k*3+2] + sigma22[k*3+2]*epsilon22[k*3+2] + sigma12[k*3+2]*epsilon12[k*3+2];
		}
		else
		{
			sigma11[k*4] = 0;
			sigma22[k*4] = 0;
			sigma12[k*4] = 0;
			sigma11[k*4+1] = 0;
			sigma22[k*4+1] = 0;
			sigma12[k*4+1] = 0;
			sigma11[k*4+2] = 0;
			sigma22[k*4+2] = 0;
			sigma12[k*4+2] = 0;

			epsilon11[k*4] = 0;
			epsilon22[k*4] = 0;
			epsilon12[k*4] = 0;
			epsilon11[k*3+1] = 0;
			epsilon22[k*4+1] = 0;
			epsilon12[k*4+1] = 0;
			epsilon11[k*4+2] = 0;
			epsilon22[k*4+2] = 0;
			epsilon12[k*4+2] = 0;

// 			Vector vm0 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->first) ;
// 			Vector vm1 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->second) ;
// 			Vector vm2 = triangles[k]->getState().getPrincipalStresses(*triangles[k]->third) ;
			vonMises[k*4]   = 0;
			vonMises[k*4+1] = 0;
			vonMises[k*4+2] = 0;
			vonMises[k*4+3] = 0;
		}
	}
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
	
}*/



int main(int argc, char *argv[])
{

	// main ellipse test

/*
	// initialize points
	Point cone(0,0) ;
	Point ctwo(3,3) ;
	Point vatwo(-2,1) ;
	Point pone(0,5) ;
	Point ptwo(-3,7) ;

	std::cout << "first center" << std::endl  ;	
	cone.print();
	std::cout << "second center" << std::endl  ;
	ctwo.print();
	std::cout << "second main axis" << std::endl  ;
	vatwo.print();
	std::cout << "first test point" << std::endl  ;
	pone.print();
	std::cout << "second test point" << std::endl  ;
	ptwo.print();

	Circle c(1., cone) ;
	
	// initialize ellipses
	Ellipse eone(5,4,0,0,1,0) ;
	Ellipse etwo(7,1,ptwo,vatwo) ;

	// get initial axis angle
	double aone = eone.getAxisAngle() ;
	double atwo = etwo.getAxisAngle() ;
	
	std::cout << "angle axis (first and second)" << std::endl  ;
	std::cout << aone << std::endl  ;
	std::cout << atwo << std::endl  ;

	// get foci for ellipse 1
	Point fone = eone.getFocus(true) ;
	Point fonebis = eone.getFocus(false) ;

	std::cout << "first ellipse foci" << std::endl  ;
	fone.print();
	fonebis.print();

	// get initial minor axis
	Point maone = eone.getMinorAxis() ;
	Point matwo = etwo.getMinorAxis() ;

	std::cout << "first ellipse minor axis" << std::endl  ;
	maone.print();
	std::cout << "second ellipse minor axis" << std::endl  ;
	matwo.print();

	// transform ellipse 1 and verify foci have been changed
	eone.setRadius(5,2) ;
	Point foneter = eone.getFocus(true) ;
	Point foneqater = eone.getFocus(false) ;
	
	std::cout << "modified first ellipse foci" << std::endl  ;
	foneter.print();
	foneqater.print();

	// choose 10 points on eone at different angle
	std::vector<double> vtheta(10) ;
	std::vector<double> rtheta(10) ;
	vtheta[0] = 0 ;
	vtheta[1] = 1 ;
	vtheta[2] = 2 ;
	vtheta[3] = 3 ;
	vtheta[4] = 5 - 2 * M_PI ;
	vtheta[5] = 8 - 2 * M_PI ;
	vtheta[6] = 13 - 4 * M_PI ;
	vtheta[7] = 21 - 6 * M_PI ;
	vtheta[8] = 34 - 10 * M_PI ;
	vtheta[9] = 55 - 18 * M_PI ;
	double coordx = 0 ;
	double coordy = 0 ;
	for(size_t i = 0 ; i < vtheta.size() ; i++)
	{
		rtheta[i] = eone.getRadiusOnEllipse(vtheta[i]) ;
		std::cout << rtheta[i] << std::endl ;
		coordx = rtheta[i]*cos(vtheta[i]) ;
		coordy = rtheta[i]*sin(vtheta[i]) ;
		coordx = coordx * coordx / (eone.getMajorRadius() * eone.getMajorRadius()) ;
		coordy = coordy * coordy / (eone.getMinorRadius() * eone.getMinorRadius()) ;
		std::cout << coordx + coordy << std::endl ;

	}



	// rotate ellipse 2 and verify axis has been changed
	etwo.setAxis(ctwo) ;
	Point matwobis = etwo.getMinorAxis() ;

	std::cout << "modified second ellipse minor axis (rotation)" << std::endl  ;
	matwobis.print();

	// project test points on ellipses and verify consistency
	Point projone = eone.project(pone) ;
//	double vpone = (pone.x * pone.x) / (eone.getMajorRadius() * eone.getMajorRadius() ) + (pone.y * pone.y) / (eone.getMinorRadius() * eone.getMinorRadius() ) ;
	Point projtwo = eone.project(ptwo) ;
//	double vptwo = (ptwo.x * pone.x) / (eone.getMajorRadius() * eone.getMajorRadius() ) + (ptwo.y * ptwo.y) / (eone.getMinorRadius() * eone.getMinorRadius() ) ; 

	std::cout << "project first point on first ellipse" << std::endl ;
	projone.print() ;
	std::cout << "project second point on first ellipse" << std::endl ;
	projtwo.print() ;

	// test pointer-wise projection
	Point * ppone = &pone ;
	Point * pptwo = &ptwo ;
//	ppone->x = pone.x ;
//	ppone->y = pone.y ;
//	pptwo->x = ptwo.x ;
//	pptwo->y = ptwo.y ;
	eone.project(ppone) ;
	eone.project(pptwo) ;
	std::cout << "project first point on first ellipse" << std::endl ;
	ppone->print() ;
	std::cout << "project second point on first ellipse" << std::endl ;
	pptwo->print() ;


	std::cout << "sample first ellipse" << std::endl ;
	// get sampling bounding points... todo
	std::vector<Point> bound = eone.getSamplingBoundingPoints(10) ;
	for(size_t i = 0 ; i < bound.size() ; i++)
	{
		bound[i].print() ;
		coordx = bound[i].x ;
		coordy = bound[i].y ;
		coordx = coordx * coordx / (eone.getMajorRadius() * eone.getMajorRadius()) ;
		coordy = coordy * coordy / (eone.getMinorRadius() * eone.getMinorRadius()) ;
		std::cout << coordx + coordy << std::endl ;
	}
	

	// bounding box
	std::cout << "get first ellipse bounding box" << std::endl ;
	std::vector<Point> bboxone = eone.getBoundingBox() ;
	for(size_t i = 0 ; i < bboxone.size() ; i++)
	{
		bboxone[i].print() ;
	}


	std::cout << "get second ellipse bounding box" << std::endl ;
	std::vector<Point> bboxtwo = etwo.getBoundingBox() ;
	for(size_t i = 0 ; i < bboxtwo.size() ; i++)
	{
		bboxtwo[i].print() ;
	}
	
	// in?
	Point pinone(4,0) ;
	bool inone = eone.in(pinone) ;
	Point pintwo(10,5) ;
	bool intwo = eone.in(pintwo) ;
	std::cout << "is inside (4,0)" << inone << std::endl ;
	std::cout << "is inside (10,5)" << intwo << std::endl ;

	eone.sampleSurface(10) ;
	std::cout << "first ellipse bounding points" << std::endl ;
	for(size_t i = 0 ; i < eone.getBoundingPoints().size() ; i++)
	{
		eone.getBoundingPoint(i).print() ;
	}

	std::cout << "first ellipse inner sampling points" << std::endl ;
	std::cout << eone.getInPoints().size() << std::endl ;
	for(size_t i = 0 ; i < eone.getInPoints().size() ; i++)
	{
		eone.getInPoint(i).print() ;
	}

	// tranform ellipse 1 and verify major and minor axis have been inverted
	eone.setRadius(1,2) ;
	Point maonebis = eone.getMinorAxis() ;

	std::cout << "modified first ellipse minor axis (radius change)" << std::endl  ;
	maonebis.print();

	double aonebis = eone.getAxisAngle() ;



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


	std::cout << 2. * POINT_TOLERANCE << std::endl ;

	// test for intersection*/
	Point origin(0.0,0.0) ;
	Point axisX(1,0) ;
	Point axisY(0,1) ;
/*	Ellipse ethree(5,2,origin,axisX) ;
	Geometry * geom ;
	geom = &ethree ;

	Line horone(origin + axisY,axisX) ;
	Line horthree(origin + axisY * 3,axisX) ;
	Line verthree(origin + axisX * 3,axisY) ;
	Line oneone(origin, axisX + axisY) ;
	Line yaxis(origin, axisY) ;

	bool interhorone = horone.intersects(geom) ;
	bool interhorthree = horthree.intersects(geom) ;
	bool interverthree = verthree.intersects(geom) ;
	bool interoneone = oneone.intersects(geom) ;
	bool intery = yaxis.intersects(geom) ;

	std::cout << interhorone << interhorthree << interverthree << interoneone << intery << std::endl ;

	std::vector<Point> intershorone = horone.intersection(geom) ;
	std::vector<Point> intershorthree = horthree.intersection(geom) ;
	std::vector<Point> intersverthree = verthree.intersection(geom) ;
	std::vector<Point> intersoneone = oneone.intersection(geom) ;
	std::vector<Point> intersy = yaxis.intersection(geom) ;

	std::cout << "intersection between Ellipse and y=1" << std::endl ;
	for(size_t i = 0 ; i < intershorone.size() ; i++)
		intershorone[i].print() ;

	std::cout << "intersection between Ellipse and y=3" << std::endl ;
	for(size_t i = 0 ; i < intershorthree.size() ; i++)
		intershorthree[i].print() ;

	std::cout << "intersection between Ellipse and x=3" << std::endl ;
	for(size_t i = 0 ; i < intersverthree.size() ; i++)
		intersverthree[i].print() ;

	std::cout << "intersection between Ellipse and y=x" << std::endl ;
	for(size_t i = 0 ; i < intersoneone.size() ; i++)
		intersoneone[i].print() ;

	std::cout << "intersection between Ellipse and x=0" << std::endl ;
	for(size_t i = 0 ; i < intersy.size() ; i++)
		intersy[i].print() ;

*/	std::cout << "origin on (-2,-2)-(1,1) ?" << std::endl ;
	Segment segtest(origin + axisX + axisY, origin - axisX - axisY) ;
	segtest.first().print() ;
	segtest.second().print() ;
	origin.print() ;
	std::cout << isAligned(origin,segtest.first(),segtest.second()) << std::endl ;/*

	Segment seghorone(origin,origin + axisX) ;
	Segment segverfive(origin,origin + axisY * 5) ;
	Segment segseven(origin + axisX * 7,origin + axisX * 10) ;

	interhorone = seghorone.intersects(geom) ;
	double interverfive = segverfive.intersects(geom) ;
	double interseven = segseven.intersects(geom) ;

	std::cout << interhorone << interverfive << interseven << std::endl ;

	intershorone = seghorone.intersection(geom) ;
	std::vector<Point> intersverfive = segverfive.intersection(geom) ;
	std::vector<Point> intersseven = segseven.intersection(geom) ;

	std::cout << "intersection between Ellipse and (0,0)-(1,0)" << std::endl ;
	for(size_t i = 0 ; i < intershorone.size() ; i++)
		intershorone[i].print() ;

	std::cout << "intersection between Ellipse and (0,0)-(0,5)" << std::endl ;
	for(size_t i = 0 ; i < intersverfive.size() ; i++)
		intersverfive[i].print() ;

	std::cout << "intersection between Ellipse and (7,0)-(10,0)" << std::endl ;
	for(size_t i = 0 ; i < intersseven.size() ; i++)
		intersseven[i].print() ;

	Rectangle rone(4,1,4,0) ;
	Rectangle rtwo(0.5,0.5,0,0) ;
	Geometry * geom1 = &rone ;
	std::vector<Point> interrone = ethree.intersection(geom1) ;
//	std::vector<Point> interrtwo = ethree.intersection(geom2) ;
	std::cout << "intersection between Ellipse and first rectangle" << std::endl ;
	for(size_t i = 0 ; i < interrone.size() ; i++)
		interrone[i].print() ;

	ethree.sampleSurface(10) ;
	Circle cone(3,0,0) ;
	Geometry * geom2 = &cone ;
	std::vector<Point> intercone = ethree.intersection(geom2) ;
	std::cout << "intersection between Ellipse and first circle" << std::endl ;
	for(size_t i = 0 ; i < intercone.size() ; i++)
		intercone[i].print() ;


//	std::cout << "intersection between Ellipse and second rectangle" << std::endl ;
//	for(size_t i = 0 ; i < interrtwo.size() ; i++)
//		interrtwo[i].print() ;


	Ellipse ei(5,2,origin,axisX) ;
	Ellipse eii(6,1,origin,axisX) ;
	Geometry * g ;
	g = &eii ;
	ei.sampleBoundingSurface(9) ;
	std::vector<Point> interv = ei.intersection(g) ;
//	std::cout << interv.size() << std::endl ;
	Point ref = interv[1] ;
	Point last = origin ;
//	for(size_t i = 0 ; i < interv.size() ; i++)
//		interv[i].print() ;
	
	double conv = 1. ;
	int iter = 10 ;	
	while((conv > POINT_TOLERANCE) && (iter < 1000))
	{
//		std::cout << "iteration " << iter << std::endl ;
		ei.sampleBoundingSurface(iter) ;
		interv = ei.intersection(g) ;
//		std::cout << interv.size() << std::endl ;
//		ref.print() ;
		iter += 1 ;
		conv = (last - ref).norm() + (last - interv[1]).norm() ;
		last = ref ;
		ref = interv[1] ;
		std::cout << "convergence = " << conv << std::endl ;
	}

	std::cout << iter << " iterations;" << std::endl ;
	std::cout << "convergence = " << conv << std::endl ;

	double z = 3 ;
	double o = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		Ellipse testellipse(z,z,o,o) ;
		for(size_t j = 20 ; j > i ; j--)
			z = z + o ;
	}*/

}

/*
	// create sample (NULL, size, center)
	Sample3D sample(NULL, 0.04,0.04,0.04,0.02,0.02,0.02) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;

	// initialize material behaviour matrix
	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;
	
	// define ITZ size
	double itzSize = 0.00001;

	// create granulometry
	int inclusionNumber = 3200 ;
	std::vector<Inclusion3D *> inclusions = GranuloBolome(0.0000384, 1, BOLOME_D)(true, .002, .0001, inclusionNumber, itzSize);

	// adjust ITZ size to smallest inclusion radius
	if(inclusionNumber)
		itzSize = inclusions[inclusions.size()-1]->getRadius() ;

	// delete inclusions and create new inclusion distribution with modified ITZ size
	for(size_t i = 0; i < inclusions.size() ; i++)
		delete inclusions[i] ;
	inclusions = GranuloBolome(0.0000384, 1, BOLOME_D)(true, .002, .0001, inclusionNumber, itzSize);

	// puts inclusions in a feature vector
	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;

	// place inclusions in sample
	int nAgg = 3200 ;
	feats=placement(sample.getPrimitive(), feats, &nAgg, 64000);
	
	// place feats (=inclusions) in feature tree, set inclusions behaviour
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		F.addFeature(&sample, feats[i]) ;
		feats[i]->setBehaviour(new Stiffness(m0*4)) ;
		std::cout << feats[i]->getRadius() << "   " << feats[i]->getCenter().x << "   "  << feats[i]->getCenter().y << "   " << feats[i]->getCenter().z << std::endl ; 
	}

	// define MohrCoulomb rupture criteria
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;

	// define linear + rupture behaviour
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;

	// set material behaviour in sample
	sample.setBehaviour(new Stiffness(m0)) ;

	// set "mesh density" with number of points on the boudary
	F.sample(2048) ;
	// set element order (existing: linear, quadratic, +time)
	F.setOrder(LINEAR) ;
	// generate element
	F.generateElements() ;
	// compute
	step() ;
	
	return 0 ;
} */
