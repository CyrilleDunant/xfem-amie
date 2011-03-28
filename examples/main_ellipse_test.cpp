// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/random.h"
#include "../utilities/placement.h"
#include "../physics/stiffness.h"
#include "../utilities/writer/voxel_writer.h"
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
double percent = 0.10 ;
double placed_area = 0 ;

double stress = 15e6 ;

Sample sample(0.04,0.04,0.,0.) ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<Inclusion *> inclusions ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;

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
Vector angle(0) ; 

double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

int totit = 5 ;

double factor = 200 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

double shape ;
double orientation ;
double spread ;

struct HomogenizationInformation
{
	Matrix C_hom ;
	Matrix C_paste ;
	Matrix C_aggregate ;
	Matrix A_paste ;
	Matrix A_aggregate ;
	
	HomogenizationInformation()
	{
		C_hom = Matrix(3,3) ;
		C_paste = Matrix(3,3) ;
		C_aggregate = Matrix(3,3) ;
		A_paste = Matrix(3,3) ;
		A_aggregate = Matrix(3,3) ;
	}
	
} ;

bool isInInclusions(Point p)
{
	for(size_t i = 0 ; i < inclusions.size() ; i++)
	{
		if(inclusions[i]->in(p))
			return true ;	
	}
	return false ;
}

Matrix getStiffness(Vector sigma_1, Vector sigma_2, Vector epsilon_1, Vector epsilon_2)
{
	Matrix K(4,4) ;
	K[0][0] = epsilon_1[0] ;
	K[0][1] = epsilon_1[1] ;
	K[2][2] = epsilon_1[0] ;
	K[2][3] = epsilon_1[1] ;
	K[1][0] = epsilon_2[0] ;
	K[1][1] = epsilon_2[1] ;
	K[3][2] = epsilon_2[0] ;
	K[3][3] = epsilon_2[1] ;
	Matrix Ki = inverse4x4Matrix(K) ;
	
	Vector s(4) ;
	s[0] = sigma_1[0] ;
	s[2] = sigma_1[1] ;
	s[1] = sigma_2[0] ;
	s[3] = sigma_2[1] ;
	Vector c = Ki*s ;
	
	Matrix stiffness(3,3) ;
	stiffness[0][0] = c[0] ;
	stiffness[0][1] = c[2] ;
	stiffness[1][0] = c[1] ;
	stiffness[1][1] = c[3] ;
	stiffness[2][2] = sigma_1[3]/epsilon_1[3] ;
	
	return stiffness ;
}

HomogenizationInformation getHomogenizationInformation()
{
	HomogenizationInformation hom ;

	Vector sigma_hom_1(3) ;
	Vector sigma_paste_1(3) ;
	Vector sigma_aggregate_1(3) ;
	Vector sigma_hom_2(3) ;
	Vector sigma_paste_2(3) ;
	Vector sigma_aggregate_2(3) ;
	
	Vector epsilon_hom_1(3) ;
	Vector epsilon_paste_1(3) ;
	Vector epsilon_aggregate_1(3) ;
	Vector epsilon_hom_2(3) ;
	Vector epsilon_paste_2(3) ;
	Vector epsilon_aggregate_2(3) ;
	
	double area_hom ;
	double area_paste ;
	double area_aggregate ;

	Vector a(3) ;	
	for(size_t i = 0 ; i < zones.size() ; i++)
		zones[i].first->setExpansion(a) ;
//	featureTree->forceEnrichmentChange() ;
	
	featureTree->resetBoundaryConditions() ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 0.000001)) ;
			
	featureTree->elasticStep() ;
//	featureTree->elasticStep() ;

	std::cout << "first step done" << std::endl ;
			
	triangles = featureTree->getElements2D() ;
	int npoints = triangles[0]->getBoundingPoints().size() ;
	std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;

	for(size_t i = 0 ; i < triangles.size() ; i++)
	{
		double area = triangles[i]->area() ;
		Point c = triangles[i]->getCenter() ;

		double sgm_11 = 0 ;
		double sgm_22 = 0 ;
		double sgm_12 = 0 ;

		double eps_11 = 0 ;
		double eps_22 = 0 ;
		double eps_12 = 0 ;
		
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3];
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3+3];
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3+6];
		
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+1];
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+4];
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+7];

		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+2];
		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+5];
		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+8];

		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3];
		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3+3];
		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3+6];
		
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+1];
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+4];
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+7];

		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+2];
		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+5];
		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+8];
		
		area_hom += area ;
		sigma_hom_1[0] += sgm_11 ;
		sigma_hom_1[1] += sgm_22 ;
		sigma_hom_1[2] += sgm_12 ;
			
		epsilon_hom_1[0] += eps_11 ;
		epsilon_hom_1[1] += eps_22 ;
		epsilon_hom_1[2] += eps_12 ;
		
		if(isInInclusions(c))
		{
			area_aggregate += area ;

			sigma_aggregate_1[0] += sgm_11 ;
			sigma_aggregate_1[1] += sgm_22 ;
			sigma_aggregate_1[2] += sgm_12 ;
			
			epsilon_aggregate_1[0] += eps_11 ;
			epsilon_aggregate_1[1] += eps_22 ;
			epsilon_aggregate_1[2] += eps_12 ;		
		}
		else
		{
			area_paste += area ;

			sigma_paste_1[0] += sgm_11 ;
			sigma_paste_1[1] += sgm_22 ;
			sigma_paste_1[2] += sgm_12 ;
			
			epsilon_paste_1[0] += eps_11 ;
			epsilon_paste_1[1] += eps_22 ;
			epsilon_paste_1[2] += eps_12 ;		
		}
	}	

	std::cout << area_aggregate << std::endl ;

	sigma_hom_1 /= area_hom ;
	sigma_paste_1 /= area_paste ;
	sigma_aggregate_1 /= area_aggregate ;

	epsilon_hom_1 /= area_hom ;
	epsilon_paste_1 /= area_paste ;
	epsilon_aggregate_1 /= area_aggregate ;

	for(size_t i = 0 ; i < zones.size() ; i++)
		zones[i].first->setExpansion(a) ;

	featureTree->resetBoundaryConditions() ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP, 0.000001)) ;
			
	featureTree->elasticStep() ;
//	featureTree->elasticStep() ;

	triangles = featureTree->getElements2D() ;
	sigma_epsilon = featureTree->getStressAndStrain() ;

	for(size_t i = 0 ; i < triangles.size() ; i++)
	{
		double area = triangles[i]->area() ;
		Point c = triangles[i]->getCenter() ;

		double sgm_11 = 0 ;
		double sgm_22 = 0 ;
		double sgm_12 = 0 ;

		double eps_11 = 0 ;
		double eps_22 = 0 ;
		double eps_12 = 0 ;
		
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3];
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3+3];
		sgm_11 += area/npoints * sigma_epsilon.first[i*npoints*3+6];
		
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+1];
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+4];
		sgm_22 += area/npoints * sigma_epsilon.first[i*npoints*3+7];

		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+2];
		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+5];
		sgm_12 += area/npoints * sigma_epsilon.first[i*npoints*3+8];

		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3];
		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3+3];
		eps_11 += area/npoints * sigma_epsilon.second[i*npoints*3+6];
		
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+1];
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+4];
		eps_22 += area/npoints * sigma_epsilon.second[i*npoints*3+7];

		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+2];
		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+5];
		eps_12 += area/npoints * sigma_epsilon.second[i*npoints*3+8];
		
		sigma_hom_2[0] += sgm_11 ;
		sigma_hom_2[1] += sgm_22 ;
		sigma_hom_2[2] += sgm_12 ;
			
		epsilon_hom_2[0] += eps_11 ;
		epsilon_hom_2[1] += eps_22 ;
		epsilon_hom_2[2] += eps_12 ;
		
		if(isInInclusions(c))
		{
			sigma_aggregate_2[0] += sgm_11 ;
			sigma_aggregate_2[1] += sgm_22 ;
			sigma_aggregate_2[2] += sgm_12 ;
			
			epsilon_aggregate_2[0] += eps_11 ;
			epsilon_aggregate_2[1] += eps_22 ;
			epsilon_aggregate_2[2] += eps_12 ;		
		}
		else
		{
			sigma_paste_2[0] += sgm_11 ;
			sigma_paste_2[1] += sgm_22 ;
			sigma_paste_2[2] += sgm_12 ;
			
			epsilon_paste_2[0] += eps_11 ;
			epsilon_paste_2[1] += eps_22 ;
			epsilon_paste_2[2] += eps_12 ;		
		}
	}	

	sigma_hom_2 /= area_hom ;
	sigma_paste_2 /= area_paste ;
	sigma_aggregate_2 /= area_aggregate ;

	epsilon_hom_2 /= area_hom ;
	epsilon_paste_2 /= area_paste ;
	epsilon_aggregate_2 /= area_aggregate ;
	
	std::cout << epsilon_aggregate_1[0] << " " << epsilon_aggregate_1[1] << std::endl ;
	std::cout << sigma_aggregate_1[0] << " " << sigma_aggregate_1[1] << std::endl ;

	std::cout << epsilon_aggregate_2[0] << " " << epsilon_aggregate_2[1] << std::endl ;
	std::cout << sigma_aggregate_2[0] << " " << sigma_aggregate_2[1] << std::endl ;

	hom.C_hom = getStiffness(sigma_hom_1,sigma_hom_2,epsilon_hom_1,epsilon_hom_2) ;
	hom.C_paste = getStiffness(sigma_paste_1,sigma_paste_2,epsilon_paste_1,epsilon_paste_2) ;
	hom.C_aggregate = getStiffness(sigma_aggregate_1,sigma_aggregate_2,epsilon_aggregate_1,epsilon_aggregate_2) ;
	
	hom.A_paste = getStiffness(epsilon_hom_1,epsilon_hom_2,epsilon_paste_1,epsilon_paste_2) ;
	hom.A_aggregate = getStiffness(epsilon_hom_1,epsilon_hom_2,epsilon_aggregate_1,epsilon_aggregate_2) ;
	
	return hom ;
}



void step()
{

	int nsteps = 30;
	int nstepstot = 30;
	int maxtries = 2000 ;
	int tries = 0 ;
	
	featureTree->setDeltaTime(0.0001);
	featureTree->setMaxIterationsPerStep(maxtries) ;
	for(size_t i = 0 ; i < nsteps ; i++)
	{
		std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;


		featureTree->resetBoundaryConditions() ;
		featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
		featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM)) ;
		
		Vector a(3) ;
		a[0] = 0.5 ;
		a[1] = 0.5 ;
		for(size_t j = 0 ; j < zones.size() ; j++)
			zones[j].first->setExpansion(a) ;

		bool go_on = featureTree->step() ;

		if(featureTree->solverConverged())
		{
			cracked_volume.push_back(featureTree->crackedVolume) ;
			damaged_volume.push_back(featureTree->damagedVolume) ;
		}
	
		triangles = featureTree->getElements2D() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;
		sigma.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		epsilon.resize(triangles.size()*triangles[0]->getBoundingPoints().size()*3) ;
		
		std::pair<Vector, Vector > sigma_epsilon = featureTree->getStressAndStrain() ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		
		sigma11.resize(sigma.size()/3) ;
		sigma22.resize(sigma.size()/3) ;
		sigma12.resize(sigma.size()/3) ;
		epsilon11.resize(sigma.size()/3) ;
		epsilon22.resize(sigma.size()/3) ;
		epsilon12.resize(sigma.size()/3) ;
		vonMises.resize(sigma.size()/3) ;
		angle.resize(sigma.size()/3) ;
		
		std::cout << "unknowns :" << x.size() << std::endl ;
		
		cracked.clear() ;
		
		int npoints = triangles[0]->getBoundingPoints().size() ;
		
		double area = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_xy = 0;
		double avg_s_xx = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double e_xx_max = 0 ;
		double e_xx_min = 0 ;
                double e_yy_max = 0 ;
                double e_yy_min = 0 ;
                double ex_count = 0 ;
		double avg_e_xx_nogel = 0;
		double avg_e_yy_nogel = 0;
		double avg_e_xy_nogel = 0;
		double avg_s_xx_nogel = 0;
		double avg_s_yy_nogel = 0;
		double avg_s_xy_nogel = 0;
		double nogel_area = 0 ;
		
		for(size_t k = 0 ; k < triangles.size() ; k++)
		{
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
					if(triangles[k]->getBoundingPoint(p).x > sample.width()*.4999)
					{
						if(e_xx_max < x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_max=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
					if(triangles[k]->getBoundingPoint(p).x < -sample.width()*.4999)
					{
						if(e_xx_min > x[triangles[k]->getBoundingPoint(p).id*2])
							e_xx_min=x[triangles[k]->getBoundingPoint(p).id*2] ;
// 						ex_count++ ;
					}
                                        if(triangles[k]->getBoundingPoint(p).y > sample.height()*.4999)
                                        {
                                                if(e_yy_max < x[triangles[k]->getBoundingPoint(p).id*2+1])
                                                        e_yy_max=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
// 						ex_count++ ;
                                        }
                                        if(triangles[k]->getBoundingPoint(p).y < -sample.height()*.4999)
                                        {
                                                if(e_yy_min > x[triangles[k]->getBoundingPoint(p).id*2+1])
                                                        e_yy_min=x[triangles[k]->getBoundingPoint(p).id*2+1] ;
// 						ex_count++ ;
                                        }
                                }
				area += triangles[k]->area() ;
				if(triangles[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] > E_max)
						E_max = triangles[k]->getBehaviour()->param[0][0] ;
					if(!triangles[k]->getBehaviour()->param.isNull() && triangles[k]->getBehaviour()->param[0][0] < E_min)
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
				
				for(size_t l = 0 ; l < triangles[k]->getBoundingPoints().size() ; l++)
				{
					Vector vm0 = triangles[k]->getState().getPrincipalStresses(triangles[k]->getBoundingPoint(l)) ;
					vonMises[k*triangles[k]->getBoundingPoints().size()+l]  = sqrt(((vm0[0]-vm0[1])*(vm0[0]-vm0[1]))/2.) ;
	
					double agl = triangles[k]->getState().getPrincipalAngle(triangles[k]->getBoundingPoint(l))[0] ;
					agl = (vm0[0] <= 0 && vm0[1] <= 0)*180. ;
					angle[k*triangles[k]->getBoundingPoints().size()+l]  = agl ;
				}
				
				double ar = triangles[k]->area() ;
				for(size_t l = 0 ; l < npoints ;l++)
				{
					avg_e_xx += (epsilon11[k*npoints+l]/npoints)*ar;
					avg_e_yy += (epsilon22[k*npoints+l]/npoints)*ar;
					avg_e_xy += (epsilon12[k*npoints+l]/npoints)*ar;
					avg_s_xx += (sigma11[k*npoints+l]/npoints)*ar;
					avg_s_yy += (sigma22[k*npoints+l]/npoints)*ar;
					avg_s_xy += (sigma12[k*npoints+l]/npoints)*ar;
				}
				
				if(triangles[k]->getEnrichmentFunctions().size() == 0)
				{
					for(size_t l = 0 ; l < npoints ;l++)
					{
						avg_e_xx_nogel += (epsilon11[k*npoints+l]/npoints)*ar;
						avg_e_yy_nogel += (epsilon22[k*npoints+l]/npoints)*ar;
						avg_e_xy_nogel += (epsilon12[k*npoints+l]/npoints)*ar;
						avg_s_xx_nogel += (sigma11[k*npoints+l]/npoints)*ar;
						avg_s_yy_nogel += (sigma22[k*npoints+l]/npoints)*ar;
						avg_s_xy_nogel += (sigma12[k*npoints+l]/npoints)*ar;
						
					}
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
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
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
		
		std::cout << "average sigma11 (no gel): " << avg_s_xx_nogel/nogel_area << std::endl ;
		std::cout << "average sigma22 (no gel): " << avg_s_yy_nogel/nogel_area << std::endl ;
		std::cout << "average sigma12 (no gel): " << avg_s_xy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon11 (no gel): " << avg_e_xx_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon22 (no gel): " << avg_e_yy_nogel/nogel_area << std::endl ;
		std::cout << "average epsilon12 (no gel): " << avg_e_xy_nogel/nogel_area << std::endl ;
		
                std::cout << "apparent extension X " << e_xx_max-e_xx_min << std::endl ;
                std::cout << "apparent extension Y " << e_yy_max-e_yy_min << std::endl ;

		std::cout << tries << std::endl ;

		if (go_on && i < nstepstot-1)
		{
			TriangleWriter writer("itz_"+itoa(i),featureTree) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;
			
/*			HomogenizationInformation inf = getHomogenizationInformation() ;
			std::fstream hom ;
			hom.open("homogenization.txt",std::ios::out|std::ios::app) ;
			hom << inf.C_hom[0][0] << "\t" << inf.C_hom[1][1] << "\t"
				<< inf.C_paste[0][0] << "\t" << inf.C_paste[1][1] << "\t"
				<< inf.C_aggregate[0][0] << "\t" << inf.C_aggregate[1][1] << "\t"
				<< inf.A_paste[0][0] << "\t" << inf.A_paste[1][1] << "\t"
				<< inf.A_aggregate[0][0] << "\t" << inf.A_aggregate[1][1] << std::endl ;
			hom.close() ;*/
			
		
			double delta_r = sqrt(aggregateArea*0.03/((double)zones.size()*M_PI))/(double)nstepstot ;
			double reactedArea = 0 ;
			
			Inclusion * current = NULL ;
			if(!zones.empty())
				current = zones[0].second ;
			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;
			for(size_t z = 0 ; z < zones.size() ; z++)
			{
				
				zones[z].first->setRadius(zones[z].first->getRadius()+delta_r) ;	
		// 		zones[z].first->reset() ;
				if(zones[z].second == current)
				{
					current_area += zones[z].first->area() ;
					current_number++ ;
				}
				else
				{
					if(current_area/zones[z-1].second->area() > 0.03)
					{
						stopped_reaction++ ;
						for(size_t m = 0 ; m < current_number ; m++)
						{
							reactedArea -= zones[z-1-m].first->area() ;
							zones[z-1-m].first->setRadius(zones[z].first->getRadius()-delta_r) ;
							reactedArea += zones[z-1-m].first->area() ;
						}
					}
					current_area = zones[z].first->area() ;
					current_number = 1 ;
					current = zones[z].second ;
				}
				reactedArea += zones[z].first->area() ;
			}
			
			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in "<< stopped_reaction << " aggs."<< std::endl ;

		
			if(go_on)
			{
				expansion_reaction.push_back(std::make_pair(reactedArea/placed_area, avg_e_xx/area)) ;
				expansion_stress_xx.push_back(std::make_pair((avg_e_xx_nogel)/(nogel_area), (avg_s_xx_nogel)/(nogel_area))) ;
				expansion_stress_yy.push_back(std::make_pair((avg_e_yy_nogel)/(nogel_area), (avg_s_yy_nogel)/(nogel_area))) ;
				apparent_extension.push_back(std::make_pair(e_xx_max-e_xx_min, e_yy_max-e_yy_min)) ;
			}
			
			if (!go_on)
				break ;

			for(size_t i = 0 ; i < expansion_reaction.size() ; i++)
			{
				std::cout << expansion_reaction[i].first << "   " 
				<< expansion_reaction[i].second << "   " 
				<< expansion_stress_xx[i].first << "   " 
				<< expansion_stress_xx[i].second << "   " 
				<< expansion_stress_yy[i].first << "   " 
				<< expansion_stress_yy[i].second << "   " 
				<< apparent_extension[i].first  << "   " 
				<< apparent_extension[i].second  << "   " 				 
				<< cracked_volume[i]  << "   " 
				<< damaged_volume[i]  << "   " 
				<< std::endl ;
			}
			
		}

	}
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZones(int n, Inclusion* inc , FeatureTree & F)
{
	zones.clear() ;
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	double nu_incompressible = .499997 ;
	double radiusFraction = 10 ;
	double radius = 0.00001 ;
	
	double E = percent*E_csh ;
	double nu = nu_csh ; //nu_incompressible ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;
	
	Vector a(double(0), 3) ;
	a[0] = 0.5 ;
	a[1] = 0.5 ;
	a[2] = 0.00 ;

	Point c = inc->getCenter() ;
	double r = inc->getRadius() ;
	
	RandomNumber rnd ;
	
	while(zones.size() < n)
	{
		Point x(rnd.uniform(-r,r), rnd.uniform(-r,r)) ;
		x += c ;
		if(inc->in(x))
		{
			bool alone = true ;
			for(size_t i = 0 ; i < zones.size() ; i++)
			{
				Point cz = zones[i].first->getCenter() ;
				if (squareDist(x, cz) < (radius*radiusFraction+radius*radiusFraction)*(radius*radiusFraction+radius*radiusFraction))
				{
					alone = false ;
					break ;
				}
			}
			if(alone)
			{
				ExpansiveZone* z = new ExpansiveZone(NULL, radius, x.x, x.y, m0, a) ;
				F.addFeature(inc, z) ;
				zones.push_back(std::make_pair(z, inc)) ;
			}
		}
	}
	return zones ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZonesHomogeneously(int n, std::vector<Inclusion * > & incs , FeatureTree & F)
{
	double E_csh = 31e9 ;
	double nu_csh = .28 ;
	double nu_incompressible = .499997 ;
	double radiusFraction = 10 ;
	double radius = 0.00001 ;
	
	double E = percent*E_csh ;
	double nu = nu_csh ; //nu_incompressible ;
	
	Matrix m0(3,3) ;
	m0[0][0] = E/(1.-nu*nu) ; m0[0][1] =E/(1.-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E/(1.-nu*nu)*nu ; m0[1][1] = E/(1.-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;
	
	Vector a(double(0), 3) ;
	a[0] = 0.5 ;
	a[1] = 0.5 ;
	a[2] = 0.00 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		Point pos(((double)rand()/RAND_MAX-.5)*(sample.width()-radius*radiusFraction),((double)rand()/RAND_MAX-.5)*(sample.height()-radius*radiusFraction)) ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*radiusFraction+radius*radiusFraction)*(radius*radiusFraction+radius*radiusFraction))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(NULL, radius, pos.x, pos.y, m0, a)) ;
		else
			i-- ;
	}
	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			if(dist(zonesToPlace[i]->getCenter(), incs[j]->getCenter()) < incs[j]->getRadius()-radius*radiusFraction /*&& incs[j]->getRadius() <= 0.008 && incs[j]->getRadius() > 0.004*/ && sample.in(zonesToPlace[i]->getCenter()))
			{
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
	}
	
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
		std::cout << aggregateArea << "  " << count << std::endl ;
	}
	
	std::cout << "initial Reacted Area = " << M_PI*radius*radius*ret.size() << " in "<< ret.size() << " zones"<< std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
}


int main(int argc, char *argv[])
{
	srand(0) ;

	FeatureTree F(&sample) ;
	featureTree = &F ;

//	for(int i = 0 ; i < 10 ; i++)
//		inclusions.push_back(new Inclusion(0.009/(i/2+1),0.,0.)) ;

	inclusions.push_back(new Inclusion(0.001, 0., 0.)) ;
	inclusions[0]->sample(40) ;

	Inclusion * cement = new Inclusion(0.0012, 0., 0.) ;
//	cement->sample(150) ;
		
	double a = 0 ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		a += inclusions[i]->area() ;
	
	std::cout << a << std::endl ;
	
/*	std::vector<Feature* > feats ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		feats.push_back(inclusions[i]) ;
	inclusions.clear() ;
	
	int n = 0 ; 
	Sample* box = new Sample(NULL, 0.035, 0.035, 0., 0.);
	feats = placement(dynamic_cast<Rectangle*>(box), feats, &n) ;
	
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		inclusions.push_back(dynamic_cast<Inclusion* >(feats[i])) ;
		inclusions[i]->getCenter().print() ;
	}*/
	
	sample.setBehaviour(new PasteBehaviour()) ;
	cement->setBehaviour(new PasteBehaviour()) ;
	F.addFeature(&sample, cement) ;
	
	AggregateBehaviour* agg = new AggregateBehaviour() ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
	{
		inclusions[i]->setBehaviour(agg) ;
		F.addFeature(cement, inclusions[i]) ;
	}
	
	zones = generateExpansiveZones(5, inclusions[0] , F) ;
	
	F.setSamplingNumber(500) ;

	step() ;
	
/*	TriangleWriter writer("test_elastic",featureTree) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;*/
	
	HomogenizationInformation inf = getHomogenizationInformation() ;
	inf.C_hom.print() ;
	inf.C_paste.print() ;
	inf.C_aggregate.print() ;
	inf.A_paste.print() ;
	inf.A_aggregate.print() ;
	
	return 0 ;

/*	sample.setBehaviour(new ElasticOnlyPasteBehaviour()) ;
	cem->setBehaviour(new ElasticOnlyPasteBehaviour()) ;
	inc->setBehaviour(new ElasticOnlyAggregateBehaviour()) ;

	F.addFeature(&sample, inc) ;
//	F.addFeature(cem, inc) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT,-0.00005)) ;
	
	step() ;	
	
	TriangleWriter trg("test",featureTree) ;
	trg.getField(TWFT_STRAIN_AND_STRESS) ;
	trg.getField(TWFT_STIFFNESS) ;
	trg.write() ;
	return 0 ;*/
}
