// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/stiffness.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/homogenization/homogenization_base.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/features.h"
#include "../physics/viscoelasticity_with_internal_variable.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

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

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress ;

Vector b(0) ;
Vector x(0) ;
Vector lastx(0) ;
Vector sigma(0) ;
Vector sigma11(0) ;
Vector sigma22(0) ;
Vector sigma12(0) ;
Vector epsilon(0) ;
Vector epsilon11(0) ;
Vector epsilon22(0) ;
Vector epsilon12(0) ;
Vector vonMises(0) ;
Vector angle(0) ;

double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

size_t current_list = DISPLAY_LIST_STRAIN_XX ;
double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

double time_step = 60*60 ;
double true_alpha = 0.5 ;

int nsteps = 0 ;


// void setBC()
// {
// 	triangles = featureTree->getState().getElements2D() ;
// 	std::vector<Point *> points ;
// 	if(firstRun)
// 	{
// 		for(size_t k = 0 ; k < triangles.size() ;k++)
// 		{
// 			for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
// 			{
// 				if(triangles[k]->getBoundingPoint(c).t >= 0)
// 				{
//
// 					if(triangles[k]->getBoundingPoint(c).x < -.01999 /*&& triangles[k]->getBoundingPoint(c).y < -0.0199*/)
// 					{
// 						featureTree->getAssembly()->setPoint( 0,0, triangles[k]->getBoundingPoint(c).id) ;
// 					}
// 					if (triangles[k]->getBoundingPoint(c).x > 0.01999 /*&& triangles[k]->getBoundingPoint(c).x > .0199*/)
// 					{
// // 						featureTree->getAssembly()->setPoint( .005,0 ,triangles[k]->getBoundingPoint(c).id) ;
//
// 						points.push_back(&triangles[k]->getBoundingPoint(c)) ;
// 					}
// 				}
// // 				else if(triangles[k]->getBoundingPoint(c).t == 0)
// // 				{
// // 					if(triangles[k]->getBoundingPoint(c).x < -.0199 /*&& triangles[k]->getBoundingPoint(c).y < -0.0199*/)
// // 					{
// // 						featureTree->getAssembly()->setPoint( 0,0, triangles[k]->getBoundingPoint(c).id) ;
// // 					}
// // 					if (triangles[k]->getBoundingPoint(c).x > 0.0199 /*&& triangles[k]->getBoundingPoint(c).x > .0199*/)
// // 					{
// // // 						featureTree->getAssembly()->setPoint( .0025,0 ,triangles[k]->getBoundingPoint(c).id) ;
// // 						featureTree->getAssembly()->setForceOn(XI, 0.5, triangles[k]->getBoundingPoint(c).id) ;
// // 					}
// // 				}
// 				else
// 				{
// 					featureTree->getAssembly()->setPoint( 0,0 ,triangles[k]->getBoundingPoint(c).id) ;
// 				}
// 			}
// 		}
// 		firstRun = false ;
// 		std::sort(points.begin(), points.end()) ;
// 		double mul = 1 ;
// 		Point * curr = points.front() ;
// 		for(size_t i = 1 ; i < points.size() ; i++)
// 		{
// 			if(points[i] != curr && mul != 1)
// 			{
// 				featureTree->getAssembly()->setForceOn(XI, 10000, curr->id) ;
// 				mul = 1 ;
// 				curr = points[i] ;
// 			}
// 			else if(mul != 1)
// 			{
// 				mul++ ;
// 			}
// 			else
// 			{
// 				featureTree->getAssembly()->setForceOn(XI, 5000, curr->id) ;
// 				mul = 1 ;
// 				curr = points[i] ;
// 			}
// 		}
// 		if(mul != 1)
// 			featureTree->getAssembly()->setForceOn(XI, 10000, points.back()->id) ;
// 		else
// 			featureTree->getAssembly()->setForceOn(XI, 5000, points.back()->id) ;
// 	}
// 	else
// 	{
// 		for(size_t k = 0 ; k < triangles.size() ;k++)
// 		{
// 			for(size_t c = 0 ;  c < triangles[k]->getBoundingPoints().size() ; c++ )
// 			{
// 				if(triangles[k]->getBoundingPoint(c).t >= 0)
// 				{
//
// 					if(triangles[k]->getBoundingPoint(c).x < -.01999 /*&& triangles[k]->getBoundingPoint(c).y < -0.0199*/)
// 					{
// 						featureTree->getAssembly()->setPoint( 0,0, triangles[k]->getBoundingPoint(c).id) ;
// 					}
// 					if (triangles[k]->getBoundingPoint(c).x > 0.01999 /*&& triangles[k]->getBoundingPoint(c).x > .0199*/)
// 					{
// // 						featureTree->getAssembly()->setPoint( .005,0 ,triangles[k]->getBoundingPoint(c).id) ;
// 						points.push_back(&triangles[k]->getBoundingPoint(c)) ;
// 					}
// 				}
// 				else
// 				{
// 					featureTree->getAssembly()->setPoint( x[triangles[k]->getBoundingPoint(c+6).id*2],x[triangles[k]->getBoundingPoint(c+6).id*2+1], triangles[k]->getBoundingPoint(c).id) ;
// 				}
// 			}
// 		}
//
// 		std::sort(points.begin(), points.end()) ;
// 		double mul = 1 ;
// 		Point * curr = points.front() ;
// 		for(size_t i = 1 ; i < points.size() ; i++)
// 		{
// 			if(points[i] != curr && mul != 1)
// 			{
// 				featureTree->getAssembly()->setForceOn(XI, 10000, curr->id) ;
// 				mul = 1 ;
// 				curr = points[i] ;
// 			}
// 			else if(mul != 1)
// 			{
// 				mul++ ;
// 			}
// 			else
// 			{
// 				featureTree->getAssembly()->setForceOn(XI, 5000, curr->id) ;
// 				curr = points[i] ;
// 				mul = 1 ;
// 			}
// 		}
// 		if(mul != 1)
// 			featureTree->getAssembly()->setForceOn(XI, 10000, points.back()->id) ;
// 		else
// 			featureTree->getAssembly()->setForceOn(XI, 5000, points.back()->id) ;
//
// 	}
//
// }

const GaussPointArray generateGaussPoints()
{
	size_t ordre = 0;
	std::valarray< std::pair<Point, double> > fin ;

	ordre = 14 ;
	fin.resize(ordre);
	fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
	fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
	fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,-0.577350269189626), 0.062969590272413) ;
	fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,-0.577350269189626), 0.066197076394253) ;
	fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
	fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
	fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.1125) ;
	fin[7] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
	fin[8] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
	fin[9] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,0.577350269189626), 0.062969590272413) ;
	fin[10] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,0.577350269189626), 0.066197076394253) ;
	fin[11] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
	fin[12] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
	fin[13] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.1125) ;

	GaussPointArray gauss(fin, LINEAR_TIME_QUADRATIC)  ;
	return gauss ;
}

Point inGlobalCoordinates(TriElement* tri, Point p)
{
	Point gp ;

	gp += tri->getBoundingPoint(1)*p.x ;
	gp += tri->getBoundingPoint(2)*p.y ;

	double t0 = tri->getBoundingPoint(3).t-tri->getBoundingPoint(0).t ;
	t0 /= 2 ;

	gp.t = 	t0*p.t ;

	if(!tri->in(gp))
		std::cerr << "blah" << std::endl ;

	return gp ;
}


void step()
{

	int nsteps = 1 ;
//	featureTree->setMaxIterationsPerStep(600) ;
	for(size_t i = 0 ; i < nsteps ; i++)
	{

		bool go_on = featureTree->step() ;
//		go_on = featureTree->step() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
		std::cerr << x.size() << "\t" << x.min() << "\t" << x.max() << std::endl ;
		triangles = featureTree->getElements2D() ;

		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;

	// 	sigma = F.strainFromDisplacements() ;
	// 	epsilon = F.stressFromDisplacements() ;
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		std::cerr << sigma.size() << "\t" << epsilon.size() << std::endl ;

		sigma11.resize(sigma.size()/3) ;
		sigma22.resize(sigma.size()/3) ;
		sigma12.resize(sigma.size()/3) ;
		epsilon11.resize(sigma.size()/3) ;
		epsilon22.resize(sigma.size()/3) ;
		epsilon12.resize(sigma.size()/3) ;
		vonMises.resize(sigma.size()/3) ;
		angle.resize(sigma.size()/3) ;

		std::cout << "unknowns :" << x.size() << std::endl ;


		int npoints = triangles[0]->getBoundingPoints().size() ;
		std::cout << npoints << std::endl ;

		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx = 0 ;
		double ex_count = 0 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;

		y_min = 0 ;
		y_max = 0 ;

		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
	/*		bool in = !triangles[k]->getEnrichmentFunctions().empty() ;*/
			bool in = false ;
			for(size_t m = 0 ; m < tris__.size() ; m++)
			{
				if(triangles[k] == tris__[m])
				{
					in = true ;
					break ;
				}
			}
			cracked.push_back(in) ;



			if(!in && !triangles[k]->getBehaviour()->fractured())
			{

				for(size_t p = 0 ;p < triangles[k]->getBoundingPoints().size() ; p++)
				{
					if(x[triangles[k]->getBoundingPoint(p).id*2] > x_max)
						x_max = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2] < x_min)
						x_min = x[triangles[k]->getBoundingPoint(p).id*2];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] > y_max)
						y_max = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(x[triangles[k]->getBoundingPoint(p).id*2+1] < y_min)
						y_min = x[triangles[k]->getBoundingPoint(p).id*2+1];
					if(triangles[k]->getBoundingPoint(p).x > 0.0799)
					{
						e_xx+=x[triangles[k]->getBoundingPoint(p).id*2] ;
						ex_count++ ;
					}
				}
				area += triangles[k]->area() ;
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(triangles[k]->getBehaviour()->param[0][0] > E_max)
						E_max = triangles[k]->getBehaviour()->param[0][0] ;
					if(triangles[k]->getBehaviour()->param[0][0] < E_min)
						E_min = triangles[k]->getBehaviour()->param[0][0] ;
				}

				sigma11[k*npoints] = sigma[k*npoints*3];
				sigma22[k*npoints] = sigma[k*npoints*3+1];
				sigma12[k*npoints] = sigma[k*npoints*3+2];
				sigma11[k*npoints+1] = sigma[k*npoints*3+3];
				sigma22[k*npoints+1] = sigma[k*npoints*3+4];
				sigma12[k*npoints+1] = sigma[k*npoints*3+5];
				sigma11[k*npoints+2] = sigma[k*npoints*3+6];
				sigma22[k*npoints+2] = sigma[k*npoints*3+7];
				sigma12[k*npoints+2] = sigma[k*npoints*3+8];

				if(npoints >3)
				{
					sigma11[k*npoints+3] = sigma[k*npoints*3+9];
					sigma22[k*npoints+3] = sigma[k*npoints*3+10];
					sigma12[k*npoints+3] = sigma[k*npoints*3+11];
					sigma11[k*npoints+4] = sigma[k*npoints*3+12];
					sigma22[k*npoints+4] = sigma[k*npoints*3+13];
					sigma12[k*npoints+4] = sigma[k*npoints*3+14];
					sigma11[k*npoints+5] = sigma[k*npoints*3+15];
					sigma22[k*npoints+5] = sigma[k*npoints*3+16];
					sigma12[k*npoints+5] = sigma[k*npoints*3+17];
				}

				if(npoints > 6)
				{
					sigma11[k*npoints+6] = sigma[k*npoints*3+18];
					sigma22[k*npoints+6] = sigma[k*npoints*3+19];
					sigma12[k*npoints+6] = sigma[k*npoints*3+20];
					sigma11[k*npoints+7] = sigma[k*npoints*3+21];
					sigma22[k*npoints+7] = sigma[k*npoints*3+22];
					sigma12[k*npoints+7] = sigma[k*npoints*3+23];
					sigma11[k*npoints+8] = sigma[k*npoints*3+24];
					sigma22[k*npoints+8] = sigma[k*npoints*3+25];
					sigma12[k*npoints+8] = sigma[k*npoints*3+26];
				}

				if(npoints > 9)
				{
					sigma11[k*npoints+9] = sigma[k*npoints*3+27];
					sigma22[k*npoints+9] = sigma[k*npoints*3+28];
					sigma12[k*npoints+9] = sigma[k*npoints*3+29];
					sigma11[k*npoints+10] = sigma[k*npoints*3+30];
					sigma22[k*npoints+10] = sigma[k*npoints*3+31];
					sigma12[k*npoints+10] = sigma[k*npoints*3+32];
					sigma11[k*npoints+11] = sigma[k*npoints*3+33];
					sigma22[k*npoints+11] = sigma[k*npoints*3+34];
					sigma12[k*npoints+11] = sigma[k*npoints*3+35];
				}

				if(npoints > 12)
				{
					sigma11[k*npoints+12] = sigma[k*npoints*3+36];
					sigma22[k*npoints+12] = sigma[k*npoints*3+37];
					sigma12[k*npoints+12] = sigma[k*npoints*3+38];
					sigma11[k*npoints+13] = sigma[k*npoints*3+39];
					sigma22[k*npoints+13] = sigma[k*npoints*3+40];
					sigma12[k*npoints+13] = sigma[k*npoints*3+41];
					sigma11[k*npoints+14] = sigma[k*npoints*3+42];
					sigma22[k*npoints+14] = sigma[k*npoints*3+43];
					sigma12[k*npoints+14] = sigma[k*npoints*3+44];
				}

				if(npoints > 15)
				{
					sigma11[k*npoints+15] = sigma[k*npoints*3+45];
					sigma22[k*npoints+15] = sigma[k*npoints*3+46];
					sigma12[k*npoints+15] = sigma[k*npoints*3+47];
					sigma11[k*npoints+16] = sigma[k*npoints*3+48];
					sigma22[k*npoints+16] = sigma[k*npoints*3+49];
					sigma12[k*npoints+16] = sigma[k*npoints*3+50];
					sigma11[k*npoints+17] = sigma[k*npoints*3+51];
					sigma22[k*npoints+17] = sigma[k*npoints*3+52];
					sigma12[k*npoints+17] = sigma[k*npoints*3+53];
				}

				epsilon11[k*npoints] = epsilon[k*npoints*3];
				epsilon22[k*npoints] = epsilon[k*npoints*3+1];
				epsilon12[k*npoints] = epsilon[k*npoints*3+2];
				epsilon11[k*npoints+1] = epsilon[k*npoints*3+3];
				epsilon22[k*npoints+1] = epsilon[k*npoints*3+4];
				epsilon12[k*npoints+1] = epsilon[k*npoints*3+5];
				epsilon11[k*npoints+2] = epsilon[k*npoints*3+6];
				epsilon22[k*npoints+2] = epsilon[k*npoints*3+7];
				epsilon12[k*npoints+2] = epsilon[k*npoints*3+8];

				if(npoints > 3)
				{
					epsilon11[k*npoints+3] = epsilon[k*npoints*3+9];
					epsilon22[k*npoints+3] = epsilon[k*npoints*3+10];
					epsilon12[k*npoints+3] = epsilon[k*npoints*3+11];
					epsilon11[k*npoints+4] = epsilon[k*npoints*3+12];
					epsilon22[k*npoints+4] = epsilon[k*npoints*3+13];
					epsilon12[k*npoints+4] = epsilon[k*npoints*3+14];
					epsilon11[k*npoints+5] = epsilon[k*npoints*3+15];
					epsilon22[k*npoints+5] = epsilon[k*npoints*3+16];
					epsilon12[k*npoints+5] = epsilon[k*npoints*3+17];
				}

				if(npoints > 6)
				{
					epsilon11[k*npoints+6] = epsilon[k*npoints*3+18];
					epsilon22[k*npoints+6] = epsilon[k*npoints*3+19];
					epsilon12[k*npoints+6] = epsilon[k*npoints*3+20];
					epsilon11[k*npoints+7] = epsilon[k*npoints*3+21];
					epsilon22[k*npoints+7] = epsilon[k*npoints*3+22];
					epsilon12[k*npoints+7] = epsilon[k*npoints*3+23];
					epsilon11[k*npoints+8] = epsilon[k*npoints*3+24];
					epsilon22[k*npoints+8] = epsilon[k*npoints*3+25];
					epsilon12[k*npoints+8] = epsilon[k*npoints*3+26];
				}

				if(npoints > 9)
				{
					epsilon11[k*npoints+9] = epsilon[k*npoints*3+27];
					epsilon22[k*npoints+9] = epsilon[k*npoints*3+28];
					epsilon12[k*npoints+9] = epsilon[k*npoints*3+29];
					epsilon11[k*npoints+10] = epsilon[k*npoints*3+30];
					epsilon22[k*npoints+10] = epsilon[k*npoints*3+31];
					epsilon12[k*npoints+10] = epsilon[k*npoints*3+32];
					epsilon11[k*npoints+11] = epsilon[k*npoints*3+33];
					epsilon22[k*npoints+11] = epsilon[k*npoints*3+34];
					epsilon12[k*npoints+11] = epsilon[k*npoints*3+35];
				}

				if(npoints > 12)
				{
					epsilon11[k*npoints+12] = epsilon[k*npoints*3+36];
					epsilon22[k*npoints+12] = epsilon[k*npoints*3+37];
					epsilon12[k*npoints+12] = epsilon[k*npoints*3+38];
					epsilon11[k*npoints+13] = epsilon[k*npoints*3+39];
					epsilon22[k*npoints+13] = epsilon[k*npoints*3+40];
					epsilon12[k*npoints+13] = epsilon[k*npoints*3+41];
					epsilon11[k*npoints+14] = epsilon[k*npoints*3+42];
					epsilon22[k*npoints+14] = epsilon[k*npoints*3+43];
					epsilon12[k*npoints+14] = epsilon[k*npoints*3+44];
				}

				if(npoints > 15)
				{
					epsilon11[k*npoints+15] = epsilon[k*npoints*3+45];
					epsilon22[k*npoints+15] = epsilon[k*npoints*3+46];
					epsilon12[k*npoints+15] = epsilon[k*npoints*3+47];
					epsilon11[k*npoints+16] = epsilon[k*npoints*3+48];
					epsilon22[k*npoints+16] = epsilon[k*npoints*3+49];
					epsilon12[k*npoints+16] = epsilon[k*npoints*3+50];
					epsilon11[k*npoints+17] = epsilon[k*npoints*3+51];
					epsilon22[k*npoints+17] = epsilon[k*npoints*3+52];
					epsilon12[k*npoints+17] = epsilon[k*npoints*3+53];
				}

				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;

					double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l))[0] ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl ;
				}

				double ar = triangles[k]->area() ;
				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ;l++)
				{
				    avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
				    avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
				    avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
				    avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
				    avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
				    avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
				}

				if(triangles[k]->getEnrichmentFunctions().size() > 0)
				{
/*					for(size_t l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k*npoints+l]/npoints)*ar;

					}*/
					nogel_area+= ar ;
				}

			}
			else
			{
				sigma11[k*npoints] = 0 ;
				sigma22[k*npoints] = 0 ;
				sigma12[k*npoints] = 0 ;
				sigma11[k*npoints+1] = 0 ;
				sigma22[k*npoints+1] = 0 ;
				sigma12[k*npoints+1] = 0 ;
				sigma11[k*npoints+2] = 0 ;
				sigma22[k*npoints+2] = 0 ;
				sigma12[k*npoints+2] = 0 ;

				if(npoints >3)
				{
					sigma11[k*npoints+3] = 0 ;
					sigma22[k*npoints+3] = 0 ;
					sigma12[k*npoints+3] = 0 ;
					sigma11[k*npoints+4] = 0 ;
					sigma22[k*npoints+4] = 0 ;
					sigma12[k*npoints+4] = 0 ;
					sigma11[k*npoints+5] = 0 ;
					sigma22[k*npoints+5] = 0 ;
					sigma12[k*npoints+5] =0 ;
				}

				epsilon11[k*npoints] = 0 ;
				epsilon22[k*npoints] = 0 ;
				epsilon12[k*npoints] = 0 ;
				epsilon11[k*npoints+1] = 0 ;
				epsilon22[k*npoints+1] = 0 ;
				epsilon12[k*npoints+1] = 0 ;
				epsilon11[k*npoints+2] = 0 ;
				epsilon22[k*npoints+2] = 0 ;
				epsilon12[k*npoints+2] = 0 ;

				if(npoints > 3)
				{
					epsilon11[k*npoints+3] = 0 ;
					epsilon22[k*npoints+3] = 0 ;
					epsilon12[k*npoints+3] =0 ;
					epsilon11[k*npoints+4] = 0 ;
					epsilon22[k*npoints+4] = 0 ;
					epsilon12[k*npoints+4] =0 ;
					epsilon11[k*npoints+5] = 0 ;
					epsilon22[k*npoints+5] =0 ;
					epsilon12[k*npoints+5] = 0 ;
				}

				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = 0 ;
				}
			}
		}


		std::cout << std::endl ;
		std::cout << "max value :" << y_max << std::endl ;
		std::cout << "min value :" << y_min << std::endl ;
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

		std::cout << "average sigma11 : " << avg_s_xx/area << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy/area << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy/area << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx/area << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy/area << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy/area << std::endl ;

/*		std::cout << "average sigma11 (no gel): " << avg_s_xx_nogel/nogel_area << std::endl ;
		std::cout << "average sigma22 (no gel): " << avg_s_yy_nogel/nogel_area << std::endl ;
		std::cout << "average sigma12 (no gel): " << avg_s_xy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon11 (no gel): " << avg_e_xx_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon22 (no gel): " << avg_e_yy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon12 (no gel): " << avg_e_xy_nogel/nogel_area << std::endl ;*/

		std::cout << "apparent extension " << e_xx/ex_count << std::endl ;
		//(1./epsilon11.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);

		double delta_r = sqrt(aggregateArea*0.12/((double)zones.size()*M_PI))/36. ;
		double reactedArea = 0 ;

		for(size_t z = 0 ; z < zones.size() ; z++)
		{
			zones[z].first->setRadius(zones[z].first->getRadius()+delta_r) ;
	// 		zones[z].first->reset() ;
			reactedArea += zones[z].first->area() ;
		}

		std::cout << "reacted Area : " << reactedArea << std::endl ;

		if (go_on)
		{
			expansion_reaction.push_back(std::make_pair(y_min, y_max)) ;
			expansion_stress.push_back(std::make_pair(avg_s_yy/area, avg_e_yy/area)) ;
		}

		for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
			std::cout << expansion_reaction[i].first << "   "
			<< expansion_reaction[i].second << "   "
			<< expansion_stress[i].first << "   "
			<< expansion_stress[i].second << "   "
			<< std::endl ;

		if (go_on)
			break ;
	}


}

int FD_implicit()
{
	bool elastic = true ;

	Matrix c(3,3) ;
	double E = 3.*1e9 ;
	double nu = 0.3 ;
	c[0][0] = E/(1-nu*nu) ; c[0][1] = E/(1-nu*nu)*nu ; c[0][2] = 0 ;
	c[1][0] = E/(1-nu*nu)*nu ; c[1][1] = E/(1-nu*nu) ; c[1][2] = 0 ;
	c[2][0] = 0 ; c[2][1] = 0 ; c[2][2] = E/(1-nu*nu)*(1.-nu)/2. ;

	Matrix e(3,3) ;
	double eta = 24*3600*60 ; // 2 month
	e = c * eta ;

	double tau = time_step ; //1 s
	double alpha = true_alpha ;
	double at1 = 1./(alpha*tau) ;
	std::cout << at1 << std::endl ;

	Matrix ce = c + (e*at1) ;

	c.print();
	e.print();
	ce.print();

	Vector u ;
	Vector v ;
	Vector du ;
	std::vector<double> u_min ;
	std::vector<double> v_min ;
	std::vector<double> time ;

	Sample mainSample(NULL, 0.02, 0.02,0,0) ;
	FeatureTree mainFT(&mainSample) ;
	mainFT.setSamplingNumber(56) ;
	mainFT.setOrder(LINEAR) ;

	Sample helpSample(NULL, 0.02,0.02,0,0) ;
	FeatureTree helpFT(&helpSample) ;
	helpFT.setSamplingNumber(56) ;
	helpFT.setOrder(LINEAR) ;

	// initialization step
	mainSample.setBehaviour(new Stiffness(ce)) ;
	if(elastic)
		mainSample.setBehaviour(new Stiffness(c)) ;

	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -1e6)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	srand(0) ;
	mainFT.step() ;
	du = mainFT.getDisplacements() ;
	std::cout << du.min() << "\t" << du.max() << std::endl ;

	if(elastic)
		return 0 ;

	u.resize(du.size(),0.) ;
	v.resize(du.size(),0.) ;

	u = du ;
	v = du*at1 ;

	helpSample.setBehaviour(new Stiffness(c)) ;
	srand(0) ;
	helpFT.step() ;

	// iteration
	for(size_t i = 0 ; i < (long) 60*60*24*365*100/time_step ; i++)
	{
		double true_time = i*time_step ;
		double true_day = true_time/(60*60*24) ;

		mainFT.resetBoundaryConditions() ;

		// boundary conditions
		mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
		mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

		// viscoelastic forces
		Vector k(v.size()) ;
		k = helpFT.getAssembly()->getMatrix()*v ;
		k *= tau ;
		for(size_t j = 0 ; j < k.size()/2 ; j++)
		{
			mainFT.getAssembly()->addForceOn(XI,-k[j*2],j);
			mainFT.getAssembly()->addForceOn(ETA,-k[j*2+1],j);
		}

		mainFT.step() ;
		Vector du(u.size()) ;
		du = mainFT.getDisplacements() ;
//		std::cout << du.min() << "\t" << du.max() << std::endl ;

		u = u + v*tau + du ;
		v = v + du*at1 ;

		std::cout << k.min() << "\t" << k.max() << std::endl ;
		std::cout << v.min() << "\t" << v.max() << std::endl ;
		std::cout << du.min() << "\t" << du.max() << std::endl ;
		std::cout << u.min() << "\t" << u.max() << std::endl ;
		if(true_day == (int) true_day)
		{
			u_min.push_back(u.min()) ;
			v_min.push_back(v.min()) ;
			time.push_back(true_day) ;
			std::cout << "day = " << true_day << std::endl ;
		}
	}

	std::fstream file ;
	file.open("arrow_"+itoa(time_step)+"_"+itoa((int) ((double) 100*alpha))+".txt", std::ios::out) ;
	for(unsigned long i = 0 ; i < u_min.size() ; i++)
		file <<  time[i] << "," << u_min[i] << std::endl ;
	file.close() ;


	return 0 ;
}

int FD_explicit()
{
	Matrix c(3,3) ;
	double E = 3.*1e9 ;
	double nu = 0.3 ;
	c[0][0] = E/(1-nu*nu) ; c[0][1] = E/(1-nu*nu)*nu ; c[0][2] = 0 ;
	c[1][0] = E/(1-nu*nu)*nu ; c[1][1] = E/(1-nu*nu) ; c[1][2] = 0 ;
	c[2][0] = 0 ; c[2][1] = 0 ; c[2][2] = E/(1-nu*nu)*(1.-nu)/2. ;

	Matrix e(3,3) ;
	double eta = 24*3600*60 ;
	e = c * eta ;

	double tau = time_step ;

	c.print();
	e.print();

	Vector u ;
	Vector ue ;
	Vector uv ;
	Vector v ;
	Vector k ;

	Sample mainSample(NULL, 0.02, 0.02,0,0) ;
	FeatureTree mainFT(&mainSample) ;
	mainFT.setSamplingNumber(56) ;
	mainFT.setOrder(LINEAR) ;

	Sample helpSample(NULL, 0.02,0.02,0,0) ;
	FeatureTree helpFT(&helpSample) ;
	helpFT.setSamplingNumber(56) ;
	helpFT.setOrder(LINEAR) ;

	mainSample.setBehaviour(new Stiffness(c)) ;
	helpSample.setBehaviour(new Stiffness(e)) ;

	// initialization ; u = 0 ; v = 0 ; sigma = 0 ;
	// first step ; sigma = sigma_1
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -1e6)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	mainFT.step() ;
	ue.resize(mainFT.getDisplacements().size()) ;
	ue = mainFT.getDisplacements() ;

	// blank step for assembly
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	helpFT.step() ;

	// get viscoelastic forces
	Vector f(ue.size()) ;
	f = helpFT.getAssembly()->getMatrix() * ue ;
	f *= (1./tau) ;
	mainFT.resetBoundaryConditions() ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	for(size_t j = 0 ; j < f.size()/2 ; j++)
	{
		mainFT.getAssembly()->addForceOn(XI, -f[j*2],  j);
		mainFT.getAssembly()->addForceOn(ETA,-f[j*2+1],j);
	}

	// solve auxiliary problem
	mainFT.step() ;
	uv.resize(mainFT.getDisplacements().size()) ;
	uv = mainFT.getDisplacements() ;
	k.resize(uv.size()) ;
	for(size_t i = 0 ; i < k.size() ; i++)
		k[i] = -ue[i] / uv[i] ;
	std::cout << k[0] << std::endl ;
//	return 0 ;

	u.resize(uv.size()) ;
	v.resize(uv.size()) ;
	std::vector<double> u_min ;
	std::vector<double> v_min ;
	std::vector<double> time ;

	std::vector<DelaunayTriangle *> mainTri = mainFT.getElements2D() ;
	std::vector<DelaunayTriangle *> helpTri = helpFT.getElements2D() ;

	mainTri[0]->getBoundingPoint(1).print() ;
//	for(size_t i = 0 ; i < helpTri.size() ; i++)
	helpTri[0]->getBoundingPoint(1).print() ;
//	return 0 ;

/*	Vector Auv(u.size()) ;
	Vector Avv(u.size()) ;
	Vector Auuv(u.size()) ;
	Vector Avuv(u.size()) ;
	for(size_t i = 0 ; i < Auv.size() ; i++)
	{
		Auv[i] = (1.-std::exp(-k[i]))*tau/k[i] ;
		Avv[i] = std::exp(-k[i]) ;
		Auuv[i] = 1. - std::exp(-k[i]) ;
		Avuv[i] = k[i] * std::exp(-k[i]) ;
		std::cout << k[i] << std::endl ;
	}*/

	for(size_t i = 0 ; i < u.size() ; i++)
	{
		u[i] = (1.-std::exp(-k[i]))*uv[i] + ue[i] ;
		v[i] = ((k[i]*std::exp(-k[i]))*uv[i] + ue[i] ) / tau ;
	}
	std::cerr << u.max() << "\t" << u.min() << std::endl ;
	time.push_back(0) ;
	u_min.push_back(u.min()) ;
	v_min.push_back(v.min()) ;

	for(size_t i = 0 ; i < (long) 60*60*24*365*100/time_step ; i++)
	{
		double true_time = i*time_step ;
		double true_day = true_time/(60*60*24) ;

		for(size_t j = 0 ; j < u.size() ; j++)
		{
			uv[j] = tau*v[j]/k[j] ;
			u[j] = u[j] + (1.-std::exp(-k[j]))*uv[j] ;
			v[j] = ((k[j]*std::exp(-k[j]))*uv[j]) / tau ;

//			u[j] = u[j] + Auv[j]*v[j] ;
//			v[j] = Avv[j]*v[j] ;
		}
		if(true_day == (int) true_day)
		{
			u_min.push_back(u.min()) ;
			v_min.push_back(v.min()) ;
			time.push_back(true_day) ;
			std::cout << "day = " << true_day << std::endl ;
		}
	}

	std::fstream file ;
	file.open("explicit_"+itoa(time_step)+".txt", std::ios::out) ;
	for(unsigned long i = 0 ; i < u_min.size() ; i++)
		file <<  time[i] << "," << u_min[i] << std::endl ;
	file.close() ;

	return 0 ;


}

int burger()
{
	double nu = 0.3 ;

	Matrix c_aggregates(3,3) ;
	double E_aggregates = 59*1e9 ;
	c_aggregates = Material::cauchyGreen(std::make_pair(E_aggregates, nu), true, SPACE_TWO_DIMENSIONAL) ;

	Matrix c_le_paste(3,3) ;
	double E_paste = 12*1e9 ;
	c_le_paste = Material::cauchyGreen(std::make_pair(E_paste, nu), true, SPACE_TWO_DIMENSIONAL) ;

	double kappa = 1.5 ;
	Matrix c_kv_paste = c_le_paste * kappa ;

	double eta = 3600*30*24 ;
	Matrix e_kv_paste = c_kv_paste * eta ;

	Matrix c_eq_paste = c_le_paste * (kappa / (1. + kappa) ) ;

	Matrix c_null = c_aggregates*0 ;

	double tau = 3600*24 ;

	double itzSize = 0.000002;
	double densite = 1.;
	int inclusionNumber = 40/*96*/ ;
	std::vector<Inclusion *> inclusions = GranuloBolome(4.79263e-07, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize) ;

	std::vector<Feature *> features ;
	for( size_t i = 0; i < inclusions.size() ; i++ )
		features.push_back( inclusions[i] ) ;
	inclusions.clear() ;

	Rectangle placeGeometry( 0.04, 0.04, 0, 0 ) ;
	int nAgg = 1 ;
	features = placement( &placeGeometry, features, &nAgg, 0, 6400 );
	double volume = 0 ;

	for( size_t i = 0 ; i < features.size() ; i++ )
		volume += features[i]->area() ;
	std::cout << "volume fraction of aggregates: " << volume/(0.04*0.04) <<std::endl ;

	Sample mainSample(NULL, 0.04, 0.04, 0, 0) ;
	FeatureTree mainFT(&mainSample) ;
	mainFT.setSamplingNumber(256) ;
	mainFT.setOrder(LINEAR) ;

	Sample helpSample(NULL, 0.04, 0.04, 0, 0) ;
	FeatureTree helpFT(&helpSample) ;
	helpFT.setSamplingNumber(256) ;
	helpFT.setOrder(LINEAR) ;

//	Sample lastSample(NULL, 0.04, 0.04, 0, 0) ;
//	FeatureTree lastFT(&lastSample) ;
//	lastFT.setSamplingNumber(256) ;
//	lastFT.setOrder(LINEAR) ;

	std::vector<Inclusion *> mainInclusions ;
	std::vector<Inclusion *> helpInclusions ;
//	std::vector<Inclusion *> lastInclusions ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		mainInclusions.push_back(new Inclusion(features[i]->getRadius(), features[i]->getCenter())) ;
		helpInclusions.push_back(new Inclusion(features[i]->getRadius(), features[i]->getCenter())) ;
//		lastInclusions.push_back(new Inclusion(features[i]->getRadius(), features[i]->getCenter())) ;
	}

	mainSample.setBehaviour(new Stiffness(c_eq_paste)) ;
	helpSample.setBehaviour(new Stiffness(c_le_paste)) ;
//	lastSample.setBehaviour(new Stiffness(e_kv_paste * (1./tau))) ;

	Stiffness * inclusionStiffness = new Stiffness(c_aggregates) ;
	Stiffness * nullStiffness = new Stiffness(c_null) ;
	for(size_t i = 0 ; i < mainInclusions.size() ; i++)
	{
		mainInclusions[i]->setBehaviour(inclusionStiffness) ;
		mainFT.addFeature(&mainSample, mainInclusions[i]) ;
		helpInclusions[i]->setBehaviour(nullStiffness) ;
		helpFT.addFeature(&helpSample, helpInclusions[i]) ;
//		lastInclusions[i]->setBehaviour(inclusionStiffness) ;
//		lastFT.addFeature(&lastSample, lastInclusions[i]) ;
	}

	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -10e6)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	mainFT.step() ;

	size_t ndof = mainFT.getDisplacements().size() ;
	Vector u(ndof) ;
	Vector v(ndof) ;
	Vector u_le_v(ndof) ;
	Vector u_le_e(ndof) ;
	Vector u_kv_v(ndof) ;
	Vector u_kv_e(ndof) ;
	Vector u_e(ndof) ;
	Vector f(ndof) ;
	Vector k(ndof) ;
	std::valarray<bool> in_paste(false, ndof/2) ;

	u_e = mainFT.getDisplacements() ;

	// blank steps to force assembly
	srand(0) ;
	helpFT.step() ;
//	srand(0) ;
//	lastFT.step() ;

	// boundary conditions
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -10e6)) ;
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;

	std::vector<DelaunayTriangle *> triangles = mainFT.getElements2D() ;
	for(size_t i = 0 ; i < triangles.size() ; i++)
	{
		if(triangles[i]->getBehaviour()->getTensor(Point(0.3,0.3,0.,0.)) != c_aggregates)
		{
			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
				in_paste[triangles[i]->getBoundingPoint(j).id] = true ;
		}
	}
	triangles.clear() ;

	// flush known displacements
	for(size_t i = 0 ; i < ndof/2 ; i++)
	{
		if(!in_paste[i])
		{
//			helpFT.getAssembly()->setPointAlong(XI, u_e[i*2+0],i) ;
//			helpFT.getAssembly()->setPointAlong(ETA,u_e[i*2+1],i) ;
		}
	}

	helpFT.step() ;
	u_le_e = helpFT.getDisplacements() ;
	u_kv_e = u_e - u_le_e ;

//	f = lastFT.getAssembly()->getMatrix()*u_kv_e ;

//	helpFT.resetBoundaryConditions() ;
//	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
//	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
//	for(size_t i = 0 ; i < ndof/2 ; i++)
//	{
//		if(in_paste[i])
//		{
//			helpFT.getAssembly()->addForceOn(XI, -f[i*2],  i);
//			helpFT.getAssembly()->addForceOn(ETA,-f[i*2+1],i);
//		}
//		else
//		{
//			helpFT.getAssembly()->setPointAlong(XI, u_e[i*2+0],i) ;
//			helpFT.getAssembly()->setPointAlong(ETA,u_e[i*2+1],i) ;
//		}
//	}
//	helpFT.step() ;
//	u_kv_v = helpFT.getDisplacements() ;
	for(size_t i = 0 ; i < ndof ; i++)
	{
//		if(in_paste[i/2])
		{
			k[i] = tau/eta ;//- u_kv_e[i] / u_kv_v[i] ;
			u_kv_v[i] = - u_kv_e[i] * k[i] ;
//			u_le_v[i] = 0 ;
		}
//		else
//		{
//			u_kv_v[i] = - u_kv_e[i] / k[i] ;
//			u_le_v[i] = - u_kv_v[i] ;
//		}
	}


	for(size_t i = 0 ; i < ndof ; i++)
	{
		u[i] = u_e[i] + (1. - std::exp(-k[i]))*u_kv_v[i] ;
		v[i] = u_e[i]/tau + k[i]/tau * std::exp(-k[i])*u_kv_v[i] ;
	}

	std::cout << u.max() << "\t" << u.min() << std::endl ;

	for(size_t j = 0 ; j < 400 ; j++)
	{
		for(size_t i = 0 ; i < ndof ; i++)
		{
//			if(in_paste[i/2])
				u_kv_v[i] = tau/k[i] * v[i] ;
//			else
//				u_kv_v[i] = 0 ;
			u[i] = u[i] + (1. - std::exp(-k[i]))*u_kv_v[i] ;
			v[i] = k[i]/tau * std::exp(-k[i])*u_kv_v[i] ;
		}
		std::cout << u.max() << "\t" << u.min() << std::endl ;
	}

	return 0 ;
}


/*	srand(0) ;
	helpFT.step() ;


	Sample mainSample(NULL, 0.02, 0.02,0,0) ;
	FeatureTree mainFT(&mainSample) ;
	mainFT.setSamplingNumber(56) ;
	mainFT.setOrder(LINEAR) ;

	Sample helpSample(NULL, 0.02,0.02,0,0) ;
	FeatureTree helpFT(&helpSample) ;
	helpFT.setSamplingNumber(56) ;
	helpFT.setOrder(LINEAR) ;

	mainSample.setBehaviour(new Stiffness(c)) ;
	helpSample.setBehaviour(new Stiffness(e)) ;

	// initialization ; u = 0 ; v = 0 ; sigma = 0 ;
	// first step ; sigma = sigma_1
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -1e6)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	mainFT.step() ;
	ue.resize(mainFT.getDisplacements().size()) ;
	ue = mainFT.getDisplacements() ;

	// blank step for assembly
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	helpFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	srand(0) ;
	helpFT.step() ;

	// get viscoelastic forces
	Vector f(ue.size()) ;
	f = helpFT.getAssembly()->getMatrix() * ue ;
	f *= (1./tau) ;
	mainFT.resetBoundaryConditions() ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	for(size_t j = 0 ; j < f.size()/2 ; j++)
	{
		mainFT.getAssembly()->addForceOn(XI, -f[j*2],  j);
		mainFT.getAssembly()->addForceOn(ETA,-f[j*2+1],j);
	}

	// solve auxiliary problem
	mainFT.step() ;
	uv.resize(mainFT.getDisplacements().size()) ;
	uv = mainFT.getDisplacements() ;
	k.resize(uv.size()) ;
	for(size_t i = 0 ; i < k.size() ; i++)
		k[i] = -ue[i] / uv[i] ;
	std::cout << k[0] << std::endl ;
//	return 0 ;

	u.resize(uv.size()) ;
	v.resize(uv.size()) ;
	std::vector<double> u_min ;
	std::vector<double> v_min ;
	std::vector<double> time ;

	std::vector<DelaunayTriangle *> mainTri = mainFT.getElements2D() ;
	std::vector<DelaunayTriangle *> helpTri = helpFT.getElements2D() ;

	mainTri[0]->getBoundingPoint(1).print() ;
//	for(size_t i = 0 ; i < helpTri.size() ; i++)
	helpTri[0]->getBoundingPoint(1).print() ;
//	return 0 ;

/*	Vector Auv(u.size()) ;
	Vector Avv(u.size()) ;
	Vector Auuv(u.size()) ;
	Vector Avuv(u.size()) ;
	for(size_t i = 0 ; i < Auv.size() ; i++)
	{
		Auv[i] = (1.-std::exp(-k[i]))*tau/k[i] ;
		Avv[i] = std::exp(-k[i]) ;
		Auuv[i] = 1. - std::exp(-k[i]) ;
		Avuv[i] = k[i] * std::exp(-k[i]) ;
		std::cout << k[i] << std::endl ;
	}

	for(size_t i = 0 ; i < u.size() ; i++)
	{
		u[i] = (1.-std::exp(-k[i]))*uv[i] + ue[i] ;
		v[i] = ((k[i]*std::exp(-k[i]))*uv[i] + ue[i] ) / tau ;
	}
	std::cerr << u.max() << "\t" << u.min() << std::endl ;
	time.push_back(0) ;
	u_min.push_back(u.min()) ;
	v_min.push_back(v.min()) ;

	for(size_t i = 0 ; i < (long) 60*60*24*365*100/time_step ; i++)
	{
		double true_time = i*time_step ;
		double true_day = true_time/(60*60*24) ;

		for(size_t j = 0 ; j < u.size() ; j++)
		{
			uv[j] = tau*v[j]/k[j] ;
			u[j] = u[j] + (1.-std::exp(-k[j]))*uv[j] ;
			v[j] = ((k[j]*std::exp(-k[j]))*uv[j]) / tau ;

//			u[j] = u[j] + Auv[j]*v[j] ;
//			v[j] = Avv[j]*v[j] ;
		}
		if(true_day == (int) true_day)
		{
			u_min.push_back(u.min()) ;
			v_min.push_back(v.min()) ;
			time.push_back(true_day) ;
			std::cout << "day = " << true_day << std::endl ;
		}
	}

	std::fstream file ;
	file.open("explicit_"+itoa(time_step)+".txt", std::ios::out) ;
	for(unsigned long i = 0 ; i < u_min.size() ; i++)
		file <<  time[i] << "," << u_min[i] << std::endl ;
	file.close() ;

	return 0 ;


}*/

int STFE()
{
	Matrix c(3,3) ;
	double E = 3.*1e9 ;
	double nu = 0.3 ;
	c[0][0] = E/(1-nu*nu) ; c[0][1] = E/(1-nu*nu)*nu ; c[0][2] = 0 ;
	c[1][0] = E/(1-nu*nu)*nu ; c[1][1] = E/(1-nu*nu) ; c[1][2] = 0 ;
	c[2][0] = 0 ; c[2][1] = 0 ; c[2][2] = E/(1-nu*nu)*(1.-nu)/2. ;

	Matrix e(3,3) ;
	double eta = 24*3600*60 ;
	e = c * eta ;

	double tau = time_step ;

	c.print();
	e.print();
	
	Vector u ;

	Sample mainSample(NULL, 0.02, 0.02,0,0) ;
	FeatureTree mainFT(&mainSample) ;
	mainFT.setSamplingNumber(56) ;
	mainFT.setOrder(LINEAR) ;
	mainFT.setDeltaTime(tau) ;

//	mainSample.setBehaviour(new KelvinVoight(c,e)) ;
 	mainSample.setBehaviour(new Stiffness(c)) ;
	
//	mainFT.generateElements() ;
	std::vector<DelaunayTriangle *> tri = mainFT.getElements2D() ;
	std::set<std::pair<std::pair<Point *, Point *>, DelaunayTriangle *> > pointList ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{	  
/*		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(0),&tri[i]->getBoundingPoint(6)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(1),&tri[i]->getBoundingPoint(7)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(2),&tri[i]->getBoundingPoint(8)), tri[i])) ;*/
/*		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(3),&tri[i]->getBoundingPoint(6)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(4),&tri[i]->getBoundingPoint(7)), tri[i])) ;
		pointList.insert(std::make_pair(std::make_pair(&tri[i]->getBoundingPoint(5),&tri[i]->getBoundingPoint(8)), tri[i])) ;*/
	}
	
	std::set<std::pair<DofDefinedBoundaryCondition *, size_t> > pointBC ;
	for(auto i = pointList.begin() ; i != pointList.end() ; i++)
	{
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_XI, i->second, i->first.first->id, 0),i->first.second->id*2)) ;
	  pointBC.insert(std::make_pair(new DofDefinedBoundaryCondition(SET_ALONG_ETA,i->second, i->first.first->id, 0),i->first.second->id*2+1)) ;
	}
	
/*	for(auto i = pointBC.begin() ; i != pointBC.end() ; i++)
	  mainFT.addBoundaryCondition(i->first) ;*/
	

	BoundingBoxDefinedBoundaryCondition * stressAfter = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, -1e6) ;
	BoundingBoxDefinedBoundaryCondition * stressNow = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_NOW, -1e6) ;
	
	mainFT.addBoundaryCondition(stressAfter) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
// 	mainFT.addBoundaryCondition(stressNow) ;
// 	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_NOW)) ;
// 	mainFT.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_NOW)) ;

	mainFT.step() ;
	u.resize(mainFT.getDisplacements().size()) ;
	u = mainFT.getDisplacements() ;
	
//	mainFT.getAssembly()->print() ;
//	exit(0) ;
	
//	MultiTriangleWriter writera("visco.txt","visco_0",&mainFT) ;
//	writera.write() ;

//	MultiTriangleWriter writerb("viscoa.txt","visco_1",&mainFT,1) ;
//	writerb.write() ;
	
//	MultiTriangleWriter writerc("visco.txt","visco_2",&mainFT) ;
	//writera.reset(&mainFT,2) ;
//	writerc.write() ;

	TriangleWriter writer("toto",&mainFT) ;

	std::cout << u.max() << std::endl ;

	std::vector<double> u_min ;
	std::vector<double> v_min ;
	std::vector<double> time ;

	u_min.push_back(u.max()) ;
	time.push_back(0) ;
	
	std::vector<double> u_before ;
	std::vector<double> u_now ;
	std::vector<double> u_after ;
	
	double u_max_before = 0 ;
	double u_max_now = 0 ;
	double u_max_after = 0 ;

	

	for(size_t n = 0 ; n < nsteps/*(long) 60*60*24*365*2/time_step*/ ; n++)
	{
		double true_time = n*tau ;
		double true_day = true_time/(60*60*24) ;

		for(auto j = pointBC.begin() ; j != pointBC.end() ; j++)
		    j->first->setData(u[j->second]) ;
		
//		stressNow->setData(0.001) ;

		mainFT.step() ;
		u = mainFT.getDisplacements() ;

/*		writera.reset(&mainFT,0) ;
		writera.write() ;
		writerb.reset(&mainFT,1) ;
		writerb.write() ;*/
/*		writerc.reset(&mainFT,2) ;
		writerc.write() ;*/


		u_max_before = 0 ;
		u_max_now = 0 ;
		u_max_after = 0 ;

//		tri = mainFT.getElements2D() ;
		for(size_t i = 0 ; i < tri.size() ; i++)
		{
			for(size_t j = 0 ; j < 3 ; j++)
			{
				size_t id_before = tri[i]->getBoundingPoint(j).id ;
				size_t id_now = tri[i]->getBoundingPoint(j+3).id ;
//				size_t id_after = tri[i]->getBoundingPoint(j+6).id ;

//				if(u[id_before*2+0] > u_max_before)
//					u_max_before = u[id_before*2+0] ;
// 				if(u[id_before*2+1] > u_max_before)
// 					u_max_before = u[id_before*2+1] ;

//				if(u[id_now*2+0] > u_max_now)
//					u_max_now = u[id_now*2+0] ;
				if(u[id_now*2+1] < u_max_now)
					u_max_now = u[id_now*2+1] ;

//				if(u[id_after*2+0] > u_max_after)
//					u_max_after = u[id_after*2+0] ;
/*				if(u[id_after*2+1] > u_max_after)
					u_max_after = u[id_after*2+1] ;*/
			
			}
		}
	
		u_before.push_back(u.min()) ;
		u_now.push_back(u_max_now) ;
// 		u_after.push_back(u_max_after) ;

	}

	std::fstream file ;
	file.open("space-time.txt", std::ios::out) ;
	for(unsigned long i = 0 ; i < u_before.size() ; i++)
	{
		file <<  i+1 << "," << u_before[i] /*<< "," << u_now[i] << "," << u_after[i]*/ << std::endl ;
	}
	file.close() ;

	return 0 ;


}

int analytical()
{
	Matrix c(3,3) ;
	double E = 3.*1e9 ;
	double nu = 0.3 ;
	c[0][0] = E/(1-nu*nu) ; c[0][1] = E/(1-nu*nu)*nu ; c[0][2] = 0 ;
	c[1][0] = E/(1-nu*nu)*nu ; c[1][1] = E/(1-nu*nu) ; c[1][2] = 0 ;
	c[2][0] = 0 ; c[2][1] = 0 ; c[2][2] = E/(1-nu*nu)*(1.-nu)/2. ;

	Matrix e(3,3) ;
	double eta = 24*3600*60 ;
	e = c * eta ;

	double tau = 60*60 ;

	Vector sigma(3) ;
	sigma[0] = -1e6 ;

	Vector epsilon_elastic(3) ;
	Vector epsilon_visco(3) ;

	Matrix s = c ;
	invert3x3Matrix(s) ;

	double k = tau / eta ;
	std::cout << k << std::endl ;

	epsilon_elastic = s*sigma ;
	epsilon_visco = -epsilon_elastic / k ;

	double uv = epsilon_visco[0]*0.02 ;
	double ue = epsilon_elastic[0]*0.02 ;

	double u = (1-std::exp(-k))*uv + ue ;
	double v = (k/tau * exp(-k))*uv + ue/tau ;

	std::cout << u << "," << v << std::endl ;

	for(unsigned long i = 0 ; i < 24*365*10 ; i++)
	{
		uv = tau*v/k ;
		u = u +(1-exp(-k))*uv ;
		v = k/tau*exp(-k)*uv ;
		if(i%24 == 0)
			std::cout << u << "," << v << std::endl ;
	}

	return 0 ;

}

int main(int argc, char *argv[])
{
      
  
  
  
  time_step = atof(argv[1]) ;
	if(argc > 2)
		nsteps = atof(argv[2]) ;
	else
		nsteps = 0 ;
//	true_alpha = atof(argv[1]) ;
	return STFE() ;

	Matrix m0_agg(3,3) ;
	m0_agg[0][0] = E_agg/(1-nu*nu) ; m0_agg[0][1] =E_agg/(1-nu*nu)*nu ; m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg/(1-nu*nu)*nu ; m0_agg[1][1] = E_agg/(1-nu*nu) ; m0_agg[1][2] = 0 ;
	m0_agg[2][0] = 0 ; m0_agg[2][1] = 0 ; m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ;

	Matrix m0_paste(3,3) ;
	m0_paste[0][0] = E_paste/(1-nu*nu) ; m0_paste[0][1] =E_paste/(1-nu*nu)*nu ; m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste/(1-nu*nu)*nu ; m0_paste[1][1] = E_paste/(1-nu*nu) ; m0_paste[1][2] = 0 ;
	m0_paste[2][0] = 0 ; m0_paste[2][1] = 0 ; m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ;

	Matrix eta(3,3) ;
	double e = atof(argv[1]) ;
	std::cout << e << std::endl ;
	eta = m0_paste * e ;

	Sample sample(NULL, 0.02, 0.02,0,0) ;
//	Sample3D sample(NULL, 2000, 2000, 2000, 0,0,0) ;

	FeatureTree F(&sample) ;
	featureTree = &F ;

	double itzSize = 0.00005;
	int inclusionNumber = 16 ;
	double masseInitiale = .00000743;
	double densite = 1.;

	std::vector<Inclusion *> inclusions = GranuloBolome(4.79263e-07, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);

	std::vector<Feature *> feats ;
	for(size_t i = 0; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;

	feats=placement(sample.getPrimitive(), feats, &inclusionNumber, 0, 6400);

	inclusions.clear() ;
	for(size_t i = 0 ; i < feats.size() ; i++)
		inclusions.push_back(dynamic_cast<Inclusion*>(feats[i])) ;

	Matrix m0(3,3) ;
//	KelvinVoight* agg = new KelvinVoight(m0_agg, m0_agg*e) ;
	Stiffness* agg = new Stiffness(Material::cauchyGreen(std::make_pair(59e9,0.3), true, SPACE_TWO_DIMENSIONAL)) ;
//	AggregateBehaviour * agg = new AggregateBehaviour() ;

	for(size_t i = 0 ; i < inclusions.size() ; i++)
	{
		inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
		inclusions[i]->setBehaviour(agg) ;
//		F.addFeature(&sample,inclusions[i]) ;
	}
//	sample.setBehaviour(new KelvinVoight(m0_paste, m0_paste*e)) ;
//	sample.setBehaviour(new Stiffness(m0_paste));
	sample.setBehaviour(new IncrementalKelvinVoight(m0_paste, m0_paste*e,0.0001));

/*	Vector k_csh(8);
	Vector g_csh(8);
	for(size_t i = 0 ; i < k_csh.size() ; i++)
		{
			k_csh[i]= 1200. ;
		}
	for(size_t i = 0 ; i < g_csh.size() ; i++)
		{
			g_csh[i]= 880. ;
		}*/
//	sample.setBehaviour(new ViscoElasticity(0.005,0.005,k_csh,g_csh)) ;
//	sample.setBehaviour(new Stiffness(Material::cauchyGreen(std::make_pair(12e9,0.3), true, SPACE_TWO_DIMENSIONAL))) ;
//	sample.setBehaviour(new PasteBehaviour()) ;

	F.setSamplingNumber(128) ;
	F.setDeltaTime(3600*24) ;
	F.setMaxIterationsPerStep(1000) ;
	F.setOrder(LINEAR) ;
//	F.setOrder(LINEAR_TIME_QUADRATIC) ;

	Function t("t") ;
	Function nt = (t+3600*12)*0.00001 ;
//	Function top = nt*(0.00001) ;
	double ntt = VirtualMachine().eval(nt, Point(0,0,0,3600*12)) ;
	ntt = 1e6 ;
	std::cerr << ntt << std::endl ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP, 0)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP, nt)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA,TOP, VirtualMachine().eval(nt, Point(0,0,0,0.025)))) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,BEFORE)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BEFORE)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM_AFTER)) ;

	step() ;

	MultiTriangleWriter writer("visco_head","visco",featureTree) ;
	writer.getField(TWFT_STRAIN) ;
	writer.write() ;

	F.resetBoundaryConditions();
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP,ntt)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT)) ;

	step() ;

	writer.reset(featureTree) ;
	writer.getField(TWFT_STRAIN) ;
	writer.append() ;

	for(size_t i = 0 ; i < 10 ; i++)
	{
//	    nt = nt + 5 ;
	    F.resetBoundaryConditions() ;
//	    F.addBoundaryCondition(new TimeContinuityBoundaryCondition());
	    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;
	    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT)) ;


	    step() ;
	writer.reset(featureTree) ;
	writer.getField(TWFT_STRAIN) ;
	writer.append() ;

	}
/*	F.resetBoundaryConditions() ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP, -5e6)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP,-1e-6)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM)) ;
	step() ;*/

// 	delete dt ;

	return 0 ;
}

