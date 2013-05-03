// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
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
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../utilities/writer/voxel_writer.h"
#include "../features/expansiveZone3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
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

GLuint DISPLAY_LIST_DISPLACEMENT = 0 ;
GLuint DISPLAY_LIST_ELEMENTS = 0 ;
GLuint DISPLAY_LIST_STRAIN_XX = 0 ;
GLuint DISPLAY_LIST_STRAIN_YY = 0 ;
GLuint DISPLAY_LIST_STRAIN_ZZ = 0 ;
GLuint DISPLAY_LIST_STRAIN_XY = 0 ;
GLuint DISPLAY_LIST_STRAIN_XZ = 0 ;
GLuint DISPLAY_LIST_STRAIN_YZ = 0 ;
GLuint DISPLAY_LIST_STRESS_XX = 0 ;
GLuint DISPLAY_LIST_STRESS_YY = 0 ;
GLuint DISPLAY_LIST_STRESS_ZZ = 0 ;
GLuint DISPLAY_LIST_STRESS_XY = 0 ;
GLuint DISPLAY_LIST_STRESS_XZ = 0 ;
GLuint DISPLAY_LIST_STRESS_YZ = 0 ;
GLuint DISPLAY_LIST_STIFFNESS = 0 ;
GLuint DISPLAY_LIST_VON_MISES  = 0 ;
GLuint DISPLAY_LIST_ANGLE  = 0 ;
GLuint DISPLAY_LIST_ENRICHMENT = 0 ;
GLuint DISPLAY_LIST_STIFFNESS_DARK = 0 ;

using namespace Mu ;

int viewangle = 0 ;
int viewangle2 = 0 ;

FeatureTree * featureTree ;
std::vector<DelaunayTetrahedron *> tets ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;
double z_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;
double z_min = 0 ;

GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double timepos = 0.00 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

int windowWidth = 600 ;
int windowHeight = 600 ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma33(0) ; 
Vector sigma12(0) ; 
Vector sigma13(0) ; 
Vector sigma23(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon33(0) ; 
Vector epsilon12(0) ; 
Vector epsilon13(0) ; 
Vector epsilon23(0) ; 
Vector vonMises(0) ; 
Vector stiffness(0) ; 
Vector angle(0) ; 
Vector damage(0) ; 

double nu = 0.2 ;
double E_agg = 100 ;//softest
double E_paste = 1 ;//stiff
double E_stiff = E_agg*10 ;//stiffer
double E_soft = E_agg/10; //stiffest

size_t current_list = DISPLAY_LIST_DISPLACEMENT ;
double factor = 0.3 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

void step()
{
	
  int nsteps = 1;// number of steps between two clicks on the opengl thing
	featureTree->setMaxIterationsPerStep(50) ;

	for(size_t i = 0 ; i < nsteps ; i++)
	{

		bool go_on = featureTree->step() ;
		
		tets= featureTree->getElements3D() ;
		x.resize(featureTree->getDisplacements().size()) ;
		x = featureTree->getDisplacements() ;

		std::pair<Vector, Vector > sigma_epsilon ;
		sigma_epsilon.first.resize(24*tets.size()) ;
		sigma_epsilon.second.resize(24*tets.size()) ;
		sigma_epsilon = featureTree->getStressAndStrain(tets) ;
		sigma.resize(sigma_epsilon.first.size()) ;
		sigma = sigma_epsilon.first ;
		epsilon.resize(sigma_epsilon.second.size()) ;
		epsilon = sigma_epsilon.second ;
		sigma11.resize(sigma.size()/6, 0.) ;
		sigma22.resize(sigma.size()/6, 0.) ;
		sigma33.resize(sigma.size()/6, 0.) ;
		sigma12.resize(sigma.size()/6, 0.) ;
		sigma13.resize(sigma.size()/6, 0.) ;
		sigma23.resize(sigma.size()/6, 0.) ;
		
		epsilon11.resize(sigma.size()/6, 0.) ;
		epsilon22.resize(sigma.size()/6, 0.) ;
		epsilon33.resize(sigma.size()/6, 0.) ;
		epsilon12.resize(sigma.size()/6, 0.) ;
		epsilon13.resize(sigma.size()/6, 0.) ;
		epsilon23.resize(sigma.size()/6, 0.) ;
		stiffness.resize(sigma.size()/6, 0.) ;
		vonMises.resize(sigma.size()/6, 0.) ;
		angle.resize(sigma.size()/6, 0.) ;
		damage.resize(sigma.size()/6, 0.) ;
	
		std::cout << "unknowns :" << x.size() << std::endl ;
	
	
		int npoints = 4 ;
	
		double volume = 0 ;
		double avg_e_xx = 0;
		double avg_e_yy = 0;
		double avg_e_zz = 0;
		double avg_e_xy = 0;
		double avg_e_xz = 0;
		double avg_e_yz = 0;
		double avg_s_xx = 0;
		double avg_s_zz = 0;
		double avg_s_yy = 0;
		double avg_s_xy = 0;
		double avg_s_xz = 0;
		double avg_s_yz = 0;
		double e_xx = 0 ;
		double ex_count = 0 ;
		double xavg = 0 ;
		
		for(size_t k = 0 ; k < tets.size() ; k++)
		{
			
			if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				
// 				for(size_t p = 0 ;p < 4 ; p++)
// 				{
// 					if(x[tets[k]->getBoundingPoint(p).id*3] > x_max)
// 						x_max = x[tets[k]->getBoundingPoint(p).id*3];
// 					if(x[tets[k]->getBoundingPoint(p).id*3] < x_min)
// 						x_min = x[tets[k]->getBoundingPoint(p).id*3];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+1] > y_max)
// 						y_max = x[tets[k]->getBoundingPoint(p).id*3+1];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+1] < y_min)
// 						y_min = x[tets[k]->getBoundingPoint(p).id*3+1];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+2] > z_max)
// 						z_max = x[tets[k]->getBoundingPoint(p).id*3+2];
// 					if(x[tets[k]->getBoundingPoint(p).id*3+2] < z_min)
// 						z_min = x[tets[k]->getBoundingPoint(p).id*3+2];
// 					if(tets[k]->getBoundingPoint(p).x > 0.0799)
// 					{
// 						e_xx+=x[tets[k]->getBoundingPoint(p).id*3] ;
// 						ex_count++ ;
// 					}
// 				}
				volume += tets[k]->volume() ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] > E_max)
						E_max = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					if(tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] < E_min)
						E_min = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					
					stiffness[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					stiffness[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0] ;
					damage[k*npoints] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+1] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+2] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
					damage[k*npoints+3] = tets[k]->getBehaviour()->getTensor(Point(.25, .25, .25))[0][0]/tets[k]->getBehaviour()->param[0][0] ;
				}
					
				sigma11[k*npoints] = sigma[k*npoints*6];
				sigma22[k*npoints] = sigma[k*npoints*6+1];
				sigma33[k*npoints] = sigma[k*npoints*6+2];
				sigma12[k*npoints] = sigma[k*npoints*6+3];
				sigma13[k*npoints] = sigma[k*npoints*6+4];
				sigma23[k*npoints] = sigma[k*npoints*6+5];
				
				sigma11[k*npoints+1] = sigma[k*npoints*6+6];
				sigma22[k*npoints+1] = sigma[k*npoints*6+7];
				sigma33[k*npoints+1] = sigma[k*npoints*6+8];
				sigma12[k*npoints+1] = sigma[k*npoints*6+9];
				sigma13[k*npoints+1] = sigma[k*npoints*6+10];
				sigma23[k*npoints+1] = sigma[k*npoints*6+11];
				
				sigma11[k*npoints+2] = sigma[k*npoints*6+12];
				sigma22[k*npoints+2] = sigma[k*npoints*6+13];
				sigma33[k*npoints+2] = sigma[k*npoints*6+14];
				sigma12[k*npoints+2] = sigma[k*npoints*6+15];
				sigma13[k*npoints+2] = sigma[k*npoints*6+16];
				sigma23[k*npoints+2] = sigma[k*npoints*6+17];
				
				sigma11[k*npoints+3] = sigma[k*npoints*6+18];
				sigma22[k*npoints+3] = sigma[k*npoints*6+19];
				sigma33[k*npoints+3] = sigma[k*npoints*6+20];
				sigma12[k*npoints+3] = sigma[k*npoints*6+21];
				sigma13[k*npoints+3] = sigma[k*npoints*6+22];
				sigma23[k*npoints+3] = sigma[k*npoints*6+23];
				
				epsilon11[k*npoints] = epsilon[k*npoints*6];
				epsilon22[k*npoints] = epsilon[k*npoints*6+1];
				epsilon33[k*npoints] = epsilon[k*npoints*6+2];
				epsilon12[k*npoints] = epsilon[k*npoints*6+3];
				epsilon13[k*npoints] = epsilon[k*npoints*6+4];
				epsilon23[k*npoints] = epsilon[k*npoints*6+5];
				
				epsilon11[k*npoints+1] = epsilon[k*npoints*6+6];
				epsilon22[k*npoints+1] = epsilon[k*npoints*6+7];
				epsilon33[k*npoints+1] = epsilon[k*npoints*6+8];
				epsilon12[k*npoints+1] = epsilon[k*npoints*6+9];
				epsilon13[k*npoints+1] = epsilon[k*npoints*6+10];
				epsilon23[k*npoints+1] = epsilon[k*npoints*6+11];
				
				epsilon11[k*npoints+2] = epsilon[k*npoints*6+12];
				epsilon22[k*npoints+2] = epsilon[k*npoints*6+13];
				epsilon33[k*npoints+2] = epsilon[k*npoints*6+14];
				epsilon12[k*npoints+2] = epsilon[k*npoints*6+15];
				epsilon13[k*npoints+2] = epsilon[k*npoints*6+16];
				epsilon23[k*npoints+2] = epsilon[k*npoints*6+17];
				
				epsilon11[k*npoints+3] = epsilon[k*npoints*6+18];
				epsilon22[k*npoints+3] = epsilon[k*npoints*6+19];
				epsilon33[k*npoints+3] = epsilon[k*npoints*6+20];
				epsilon12[k*npoints+3] = epsilon[k*npoints*6+21];
				epsilon13[k*npoints+3] = epsilon[k*npoints*6+22];
				epsilon23[k*npoints+3] = epsilon[k*npoints*6+23];
				
				Vector vm0(0.,1) ;
				Vector agl(0.,3) ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					tets[k]->getState().getField( PRINCIPAL_REAL_STRESS_FIELD,  tets[k]->getCenter(), vm0, false) ;
				if(tets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
					tets[k]->getState().getField( PRINCIPAL_ANGLE_FIELD, tets[k]->getCenter(), agl, false) ;
				for(size_t l = 0 ; l < 4 ; l++)
				{
					vonMises[k*npoints+l]  = vm0[0] ;
					angle[k*npoints+l]  = agl[0] ;
				}
				
				double ar = tets[k]->volume() ;
				Vector avgsig(6) ;
				Vector avgeps(6) ;
				tets[k]->getState().getAverageField(REAL_STRESS_FIELD, avgsig);
				tets[k]->getState().getAverageField(STRAIN_FIELD, avgeps);

				avg_e_xx += avgeps[0] * ar;
				avg_e_yy += avgeps[1] * ar;
				avg_e_zz += avgeps[2] * ar;
				avg_e_xy += avgeps[3] * ar;
				avg_e_xz += avgeps[4] * ar;
				avg_e_yz += avgeps[5] * ar;
				
				avg_s_xx += avgsig[0] * ar;
				avg_s_yy += avgsig[1] * ar;
				avg_s_zz += avgsig[2] * ar;
				avg_s_xy += avgsig[3] * ar;
				avg_s_xz += avgsig[4] * ar;
				avg_s_yz += avgsig[5] * ar;
			
				for(size_t l = 0 ; l < npoints ;l++)
				{
					xavg += x[tets[k]->getBoundingPoint(l).id]*ar/npoints ;
				}
	
			}
		}
			
		xavg /= volume ;
		
		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "avg value :" << xavg << std::endl ;
		std::cout << "max sigma11 :" << sigma11.max() << std::endl ;
		std::cout << "min sigma11 :" << sigma11.min() << std::endl ;
		std::cout << "max sigma12 :" << sigma12.max() << std::endl ;
		std::cout << "min sigma12 :" << sigma12.min() << std::endl ;
		std::cout << "max sigma13 :" << sigma13.max() << std::endl ;
		std::cout << "min sigma13 :" << sigma13.min() << std::endl ;
		std::cout << "max sigma22 :" << sigma22.max() << std::endl ;
		std::cout << "min sigma22 :" << sigma22.min() << std::endl ;
		std::cout << "max sigma23 :" << sigma23.max() << std::endl ;
		std::cout << "min sigma23 :" << sigma23.min() << std::endl ;
		std::cout << "max sigma33 :" << sigma33.max() << std::endl ;
		std::cout << "min sigma33 :" << sigma33.min() << std::endl ;
		
		std::cout << "max epsilon11 :" << epsilon11.max() << std::endl ;
		std::cout << "min epsilon11 :" << epsilon11.min() << std::endl ;
		std::cout << "max epsilon12 :" << epsilon12.max() << std::endl ;
		std::cout << "min epsilon12 :" << epsilon12.min() << std::endl ;
		std::cout << "max epsilon13 :" << epsilon13.max() << std::endl ;
		std::cout << "min epsilon13 :" << epsilon13.min() << std::endl ;
		std::cout << "max epsilon22 :" << epsilon22.max() << std::endl ;
		std::cout << "min epsilon22 :" << epsilon22.min() << std::endl ;
		std::cout << "max epsilon23 :" << epsilon23.max() << std::endl ;
		std::cout << "min epsilon23 :" << epsilon23.min() << std::endl ;
		std::cout << "max epsilon33 :" << epsilon33.max() << std::endl ;
		std::cout << "min epsilon33 :" << epsilon33.min() << std::endl ;
		
		std::cout << "max von Mises :" << vonMises.max() << std::endl ;
		std::cout << "min von Mises :" << vonMises.min() << std::endl ;
		
		std::cout << "average sigma11 : " << avg_s_xx/volume << std::endl ;
		std::cout << "average sigma22 : " << avg_s_yy/volume << std::endl ;
		std::cout << "average sigma33 : " << avg_s_zz/volume << std::endl ;
		std::cout << "average sigma12 : " << avg_s_xy/volume << std::endl ;
		std::cout << "average sigma13 : " << avg_s_xz/volume << std::endl ;
		std::cout << "average sigma23 : " << avg_s_yz/volume << std::endl ;
		std::cout << "average epsilon11 : " << avg_e_xx/volume << std::endl ;
		std::cout << "average epsilon22 : " << avg_e_yy/volume << std::endl ;
		std::cout << "average epsilon33 : " << avg_e_zz/volume << std::endl ;
		std::cout << "average epsilon12 : " << avg_e_xy/volume << std::endl ;
		std::cout << "average epsilon13 : " << avg_e_xz/volume << std::endl ;
		std::cout << "average epsilon23 : " << avg_e_yz/volume << std::endl ;
		
	}
// 	VoxelWriter vw1("sphere_stiffness", 100) ;
// 	vw1.getField(featureTree, VWFT_STIFFNESS) ;
// 	vw1.write();
	
	VoxelWriter vw("sphere_stress", 100) ;
	vw.getField(featureTree, VWFT_STRESS) ;
	vw.write();
	VoxelWriter vw0("sphere_strain", 100) ;
	vw0.getField(featureTree, VWFT_STRAIN) ;
	vw0.write();
	exit(0) ;
}


std::pair<double, double> centile(const Vector & v)
{
	Vector vs(v) ;
	std::sort(&vs[0], &vs[vs.size()]) ;
	return std::make_pair(vs[(vs.size()-1)*0.01], vs[(vs.size()-1)*0.99]) ;
}

void printScreen()
{
	std::valarray<unsigned int> frame(windowWidth*windowHeight) ;
	glReadPixels(0,0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, &frame[0]) ;
}

int main(int argc, char *argv[])
{
			
// 	Function myFunc("x 1 +") ; myFunc /= 2 ;
// 		
// 	VirtualMachine().print(myFunc); std::cout << VirtualMachine().eval(myFunc, 2) << std::endl ;
// 	
// 	myFunc = Function("x 1 +") ; myFunc *= 2 ;
// 		
// 	VirtualMachine().print(myFunc); std::cout << VirtualMachine().eval(myFunc, 2) << std::endl ;
// 	
// 	myFunc = Function("x 1 +") ; myFunc -= 2 ;
// 		
// 	VirtualMachine().print(myFunc); std::cout << VirtualMachine().eval(myFunc, 2) << std::endl ;
// 	
// 	myFunc = Function("x 1 +") ; myFunc += 2 ;
// 		
// 	VirtualMachine().print(myFunc); std::cout << VirtualMachine().eval(myFunc, 2) << std::endl ;
// 	
// 	exit(0) ;
	double nu = 0.2 ;
	double E = 1 ;
	Sample3D samplers(nullptr, 400,400,400,200,200,200) ;

	FeatureTree F(&samplers) ;
	featureTree = &F ;
	Matrix d0(3,3) ;
	d0[0][0] = 1 ;
	d0[1][1] = 1 ;
	d0[2][2] = 1 ;

	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;

	nu = 0.2 ;
	E = 10 ;
	Matrix m1(6,6) ;
	m1[0][0] = 1. - nu ; m1[0][1] = nu ;      m1[0][2] = nu ;
	m1[1][0] = nu ;      m1[1][1] = 1. - nu ; m1[1][2] = nu ;
	m1[2][0] = nu ;      m1[2][1] = nu ;      m1[2][2] = 1. - nu ;
	m1[3][3] = 0.5 - nu ;
	m1[4][4] = 0.5 - nu ;
	m1[5][5] = 0.5 - nu ;
	m1 *= E/((1.+nu)*(1.-2.*nu)) ;

	Matrix d1(3,3) ;
	d1[0][0] = 1e2 ;
	d1[1][1] = 1e2 ;
	d1[2][2] = 1e2 ;

	double itzSize = 0.000001;
	int inclusionNumber = 0 ; //48000*2 ;

//	std::vector<Inclusion3D *> inclusions = GranuloBolome(5.44e-05, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);
// 	std::vector<Inclusion3D *> inclusions = ParticleSizeDistribution::get3DInclusions(0.0025, 5.44e-5, BOLOME_A, PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;
//	GranuloBolome(5.44e-05, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);
// 	if(inclusionNumber)
// 		itzSize = inclusions[inclusions.size()-1]->getRadius() ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		delete inclusions[i] ;
// 
// 	inclusions = GranuloBolome(3.84e-05/1.8, 1, BOLOME_A)(true, .0025, .0001, inclusionNumber, itzSize);

// 	std::vector<Feature *> feats ;
// 	for(size_t i = 0; i < inclusions.size() ; i++)
// 		feats.push_back(inclusions[i]) ;

// 	int nAgg = 6000 ;
	
	MohrCoulomb * mc = new MohrCoulomb(30, -60) ;
	StiffnessAndFracture * sf = new StiffnessAndFracture(m0*0.5, mc) ;
	Stiffness * s = new Stiffness(m0) ;
	Stiffness * ss = new Stiffness(m1) ;
// 	std::cout << feats.size() << std::endl ;
// 	for(size_t i = 0 ; i < feats.size() ; i++)
// 	{
// 		feats[i]->setCenter(feats[i]->getCenter()*10000.);
// 		dynamic_cast<Inclusion3D *>(feats[i])->setRadius(feats[i]->getRadius()*10000.)  ;
// 		F.addFeature(&samplers, feats[i]) ;
// 		feats[i]->setBehaviour(ss) ;
// // 		feats[i]->setBehaviour(new VoidForm()) ;
// // 		feats[i]->setBehaviour(new Laplacian(d1)) ;
// // 		std::cout << feats[i]->getRadius() << "   " << feats[i]->getCenter().x << "   "  << feats[i]->getCenter().y << "   " << feats[i]->getCenter().z << std::endl ; 
// 	}
	samplers.setBehaviour(new /*WeibullDistributed*/Stiffness(m0/*,0.1*/)) ;
// 	samplers.setBehaviour(new Laplacian(d0)) ;
	Vector a(0.,6) ;// a[0] = 1 ; a[1] = 1 ; a[2] = 1 ; 
	ExpansiveZone3D inc(&samplers,100, 200, 200, 200, m1*4, a) ;
// 	Inclusion3D inc(100, 200, 200, 200) ;
// 	OctahedralInclusion * inc0 = new OctahedralInclusion(208.40029238347645, 200, 200, 200) ;
// 	inc->setBehaviour(new StiffnessWithImposedDeformation(m1*4.,a)) ;
// 	inc.setBehaviour(new Stiffness(m1*4)) ;
// 	inc0->setBehaviour(new Laplacian(d1)) ;
	
	F.addFeature(&samplers, &inc) ;
// 	F.addFeature(&samplers, inc0) ;
	F.setSamplingNumber(atoi(argv[1])) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_RIGHT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_RIGHT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, TOP_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, TOP_LEFT_BACK)) ;
// 	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_FRONT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_FRONT)) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ZETA, FRONT, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP, 1.)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_LEFT_BACK)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BOTTOM_LEFT_BACK)) ;
	
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, RIGHT)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, TOP)) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;
	F.setOrder(QUADRATIC) ;

	step() ;
	
// 	delete dt ;
	
	return 0 ;
}
