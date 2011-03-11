// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
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
#include "../physics/diffusion.h"
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
#include "../utilities/placement.h"
#include "../physics/stiffness.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#include "../utilities/writer/voxel_writer.h"
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

double prop ;

double x_max = 0 ;
double y_max = 0 ;
double z_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;
double z_min = 0 ;

GLint xangle = 0;
GLint yangle = 0;
GLint zangle = 0;

double timepos = 0.1 ;

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

Vector epsilon_bcx(6) ;
Vector sigma_bcx(6) ;
Vector epsilon_bcy(6) ;
Vector sigma_bcy(6) ;
Vector epsilon_bcz(6) ;
Vector sigma_bcz(6) ;

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

typedef enum
{
	S1,
	S2024,
	S3200,
	O1,
} BenchmarkMicrostructure ;

typedef enum
{
	DIFFUSION,
	ELASTICITY,
} BenchmarkPhenomenon ;

BenchmarkMicrostructure getMicrostructure(std::string micro)
{
	if(micro == std::string("S1"))
		return S1 ;
	if(micro == std::string("S2024"))
		return S2024 ;
	if(micro == std::string("S3200"))
		return S3200 ;
	if(micro == std::string("O1"))
		return O1 ;
	return S1 ;
}

BenchmarkPhenomenon getPhenomenon(std::string pheno)
{
	if(pheno == std::string("diffusion"))
		return DIFFUSION ;
	if(pheno == std::string("elasticity"))
		return ELASTICITY ;
	return DIFFUSION ;
}

double getLength(BenchmarkMicrostructure micro)
{
	switch(micro)
	{
	case(S1) :
		return 0.15 ;
	case(S2024) :
		return 0.15 ;
	case(S3200) :
		return 400. ;
	case(O1) :
		return 1. ;
	}
	return 1. ;
}

int main(int argc, char *argv[])
{
	std::cout << "usage: benchmark <microstructure> <phenomenon> <order> <scale> <properties> <sampling>" << std::endl ;
	std::cout << "all fields are required" << std::endl ;
	std::cout << "<microstructure>\tstring among <S1, S2024, S3200, O1>" << std::endl ;
	std::cout << "<phenomenon>\t\tstring among <diffusion, elasticity>" << std::endl ;
	std::cout << "<order>\t\t\tinteger : order of the elements (1 for linear, 2 for quadratic)" << std::endl ;
	std::cout << "<scale>\t\t\tdouble : scale factor" << std::endl ;
	std::cout << "<properties>\t\tdouble : value of the Young's Modulus or Diffusion coefficient of the inclusion(s)" << std::endl ;
	std::cout << "<sampling>\t\tinteger : number of points at the surface of the REV" << std::endl ;

	if(argc != 7)
		return 1 ;

	BenchmarkMicrostructure micro = getMicrostructure(std::string(argv[1])) ;
	BenchmarkPhenomenon pheno = getPhenomenon(std::string(argv[2])) ;
	int order = atoi(argv[3]) ;
	if(order < 1 || order > 2) { order = 1 ; }
	double scale = atof(argv[4]) ;
	prop = atof(argv[5]) ;
	int sampling = atoi(argv[6]) ;
	double length = getLength(micro) ;

	double size = scale*length ;
	double halfSize = size/2 ;

	Sample3D sample(NULL, size, size, size, halfSize, halfSize, halfSize) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;

	Form* behaviour = NULL ;

	switch(pheno)
	{
	case DIFFUSION:
	{
		Matrix d0(3,3) ;
		double lambda = 1 ;
		d0[0][0] = lambda ;
		d0[1][1] = lambda ;
		d0[2][2] = lambda ;
		sample.setBehaviour(new Laplacian(d0)) ;
		
		Matrix d1(3,3) ;
		lambda = prop ;
		d1[0][0] = lambda ;
		d1[1][1] = lambda ;
		d1[2][2] = lambda ;
		behaviour = new Laplacian(d1) ;
		break ;
	}
	case ELASTICITY:
	{
		
		double nu = 0.2 ;
		double E = 1 ;	
		Matrix m0(6,6) ;
		m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
		m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
		m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
		m0[3][3] = 0.5 - nu ;
		m0[4][4] = 0.5 - nu ;
		m0[5][5] = 0.5 - nu ;
		m0 *= E/((1.+nu)*(1.-2.*nu)) ;
		sample.setBehaviour(new Stiffness(m0)) ;
		
		E = prop ;
		std::cout << prop << std::endl ;
		Matrix m1(6,6) ;
		m1[0][0] = 1. - nu ; m1[0][1] = nu ; m1[0][2] = nu ;
		m1[1][0] = nu ; m1[1][1] = 1. - nu ; m1[1][2] = nu ;
		m1[2][0] = nu ; m1[2][1] = nu ; m1[2][2] = 1. - nu ;
		m1[3][3] = 0.5 - nu ;
		m1[4][4] = 0.5 - nu ;
		m1[5][5] = 0.5 - nu ;
		m1 *= E/((1.+nu)*(1.-2.*nu)) ;
		behaviour = new Stiffness(m1) ;
		break ;
	}
	}

	std::string str_micro = "S1" ;
	if(micro != O1)
	{
		std::vector<Inclusion3D * > inclusions ;
		if(micro == S1)
		{
			inclusions.push_back(new Inclusion3D(0.0623*scale, sample.getCenter().x, sample.getCenter().y, sample.getCenter().z)) ;
			std::cout << inclusions[0]->volume() << std::endl ;
			std::cout << sample.volume() << std::endl ;
			std::cout << inclusions[0]->volume()/sample.volume() << std::endl ;
		}
		else
		{
			str_micro = "S2024" ;
		 	int n = 2024 ;
			std::string file = "sphere_2024.txt" ;
		 	std::vector<std::string> columns ;
		 	columns.push_back("center_x") ;
		 	columns.push_back("center_y") ;
		 	columns.push_back("center_z") ;
		 	columns.push_back("radius") ;
			if(micro == S3200)
			{
				str_micro = "S3200" ;
				n = 3200 ;
				file = "sphere_3200.txt" ;
				columns.clear() ;
			 	columns.push_back("radius") ;
			 	columns.push_back("center_x") ;
			 	columns.push_back("center_y") ;
			 	columns.push_back("center_z") ;
			}
			GranuloFromFile spheres(file, columns) ;
			inclusions = spheres.getInclusion3D(n, scale) ;
		}
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			inclusions[i]->setBehaviour(behaviour) ;
			F.addFeature(&sample, inclusions[i]) ;
		}
	
	}
	else
	{
		str_micro = "O1" ;
		OctahedralInclusion* oct = new OctahedralInclusion(NULL, 0.4182554*2.*scale, sample.getCenter().x, sample.getCenter().y, sample.getCenter().z) ;
		oct->setBehaviour(behaviour) ;
		F.addFeature(&sample, oct) ;
	}

	F.setSamplingNumber(sampling) ;
	F.setMaxIterationsPerStep(2);
	F.setDeltaTime(0.001);
	F.setElementGenerationMethod(0,true) ;
	
	if(order == 2)
		F.setOrder(QUADRATIC) ;
	else
		F.setOrder(LINEAR) ;
	
	Function pos("x") ;
	Function grad = pos*0.01/length ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 0.01*scale)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, TOP, grad)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, BOTTOM, grad)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, BACK, grad)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, FRONT, grad)) ;

	F.step() ;
	
	tets= F.getElements3D() ;
	x.resize(F.getDisplacements().size()) ;
	x = featureTree->getDisplacements() ;

	switch(pheno)
	{
	case DIFFUSION:
	{
		Vector gradient11(4*tets.size()) ;
		Vector gradient22(4*tets.size()) ;
		Vector gradient33(4*tets.size()) ;
		Vector flux11(4*tets.size()) ;
		Vector flux22(4*tets.size()) ;
		Vector flux33(4*tets.size()) ;
		{
			Vector gradient(12*tets.size()) ;
			Vector flux(12*tets.size()) ;
			{
				std::pair<Vector,Vector> gradient_flux ;
				gradient_flux.first.resize(12*tets.size()) ;
				gradient_flux.second.resize(12*tets.size()) ;
				gradient_flux = featureTree->getGradientAndFlux(tets) ;		
				gradient = gradient_flux.first ;
				flux = gradient_flux.second ;
			}

			std::cout << "get gradient and flux..." << std::endl ;		
			for(size_t i = 0 ; i < tets.size() ; i++)
			{
				for(size_t j = 0 ; j < 4 ; j++)
				{
					gradient11[4*i+j] = gradient[4*3*i+3*j+0] ;
					gradient22[4*i+j] = gradient[4*3*i+3*j+1] ;
					gradient33[4*i+j] = gradient[4*3*i+3*j+2] ;
					flux11[4*i+j] = flux[4*3*i+3*j+0] ;
					flux22[4*i+j] = flux[4*3*i+3*j+1] ;
					flux33[4*i+j] = flux[4*3*i+3*j+2] ;
				}
			}
		}

		std::cout << "averaging gradient and flux..." << std::endl ;		
		Vector average_gradient(3) ;
		Vector average_flux(3) ;
		double total_volume = 0. ;
		for(size_t i = 0 ; i < tets.size() ; i++)
		{
			double volume = tets[i]->volume() ;
			total_volume += volume ;
			for(size_t j = 0 ; j < 4 ; j++)
			{
				average_gradient[0] += gradient11[4*i+j]*volume/4 ;
				average_gradient[1] += gradient22[4*i+j]*volume/4 ;
				average_gradient[2] += gradient33[4*i+j]*volume/4 ;
				average_flux[0] += flux11[4*i+j]*volume/4 ;
				average_flux[1] += flux22[4*i+j]*volume/4 ;
				average_flux[2] += flux33[4*i+j]*volume/4 ;
			}
		}
		
		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "max flux11 :" << flux11.max() << std::endl ;
		std::cout << "min flux11 :" << flux11.min() << std::endl ;
		std::cout << "max flux22 :" << flux22.max() << std::endl ;
		std::cout << "min flux22 :" << flux22.min() << std::endl ;
		std::cout << "max flux33 :" << flux33.max() << std::endl ;
		std::cout << "min flux33 :" << flux33.min() << std::endl ;
		
		std::cout << "max gradient11 :" << gradient11.max() << std::endl ;
		std::cout << "min gradient11 :" << gradient11.min() << std::endl ;
		std::cout << "max gradient22 :" << gradient22.max() << std::endl ;
		std::cout << "min gradient22 :" << gradient22.min() << std::endl ;
		std::cout << "max gradient33 :" << gradient33.max() << std::endl ;
		std::cout << "min gradient33 :" << gradient33.min() << std::endl ;
		
		std::cout << "average flux11 : " << average_flux[0]/total_volume << std::endl ;
		std::cout << "average flux22 : " << average_flux[1]/total_volume << std::endl ;
		std::cout << "average flux33 : " << average_flux[2]/total_volume << std::endl ;
		std::cout << "average gradient11 : " << average_gradient[0]/total_volume << std::endl ;
		std::cout << "average gradient22 : " << average_gradient[1]/total_volume << std::endl ;
		std::cout << "average gradient33 : " << average_gradient[2]/total_volume << std::endl ;
		
		
		std::string filebench("benchmark.txt") ;
		std::fstream out ;
		out.open(filebench.c_str(), std::ios::out|std::ios::app) ;
		out << "DIFFUSION\t" << str_micro << "\t" << "D_inc = " << prop << "\t" 
			<< "dof = " << x.size() << "\t"
			<< "D11 = " << -average_flux[0]/average_gradient[0] << std::endl ;
		out.close() ;

		VoxelWriter vw("fem_s1", 200) ;
		vw.getField(featureTree, VWFT_GRADIENT) ;
		vw.write();
		
		break ;
	}
	case ELASTICITY:
	{
		Vector strain11(4*tets.size()) ;
		Vector strain22(4*tets.size()) ;
		Vector strain33(4*tets.size()) ;
		Vector stress11(4*tets.size()) ;
		Vector stress22(4*tets.size()) ;
		Vector stress33(4*tets.size()) ;
		{
			Vector stress(24*tets.size()) ;
			Vector strain(24*tets.size()) ;
			{
				std::pair<Vector,Vector> stress_strain ;
				stress_strain.first.resize(24*tets.size()) ;
				stress_strain.second.resize(24*tets.size()) ;
				stress_strain = featureTree->getStressAndStrain(tets) ;		
				stress = stress_strain.first ;
				strain = stress_strain.second ;
			}
		
			std::cout << "get stress and strain..." << std::endl ;		
			for(size_t i = 0 ; i < tets.size() ; i++)
			{
				for(size_t j = 0 ; j < 4 ; j++)
				{
					stress11[4*i+j] = stress[4*6*i+6*j+0] ;
					stress22[4*i+j] = stress[4*6*i+6*j+1] ;
					stress33[4*i+j] = stress[4*6*i+6*j+2] ;
					strain11[4*i+j] = strain[4*6*i+6*j+0] ;
					strain22[4*i+j] = strain[4*6*i+6*j+1] ;
					strain33[4*i+j] = strain[4*6*i+6*j+2] ;
				}
			}
		}
		
		std::cout << "averaging stress and strain..." << std::endl ;		
		Vector average_stress(3) ;
		Vector average_strain(3) ;
		double total_volume = 0. ;
		for(size_t i = 0 ; i < tets.size() ; i++)
		{
			double volume = tets[i]->volume() ;
			total_volume += volume ;
			for(size_t j = 0 ; j < 4 ; j++)
			{
				average_stress[0] += stress11[4*i+j]*volume/4 ;
				average_stress[1] += stress22[4*i+j]*volume/4 ;
				average_stress[2] += stress33[4*i+j]*volume/4 ;
				average_strain[0] += strain11[4*i+j]*volume/4 ;
				average_strain[1] += strain22[4*i+j]*volume/4 ;
				average_strain[2] += strain33[4*i+j]*volume/4 ;
			}
		}
		
		std::cout << std::endl ;
		std::cout << "max value :" << x_max << std::endl ;
		std::cout << "min value :" << x_min << std::endl ;
		std::cout << "max stress11 :" << stress11.max() << std::endl ;
		std::cout << "min stress11 :" << stress11.min() << std::endl ;
		std::cout << "max stress22 :" << stress22.max() << std::endl ;
		std::cout << "min stress22 :" << stress22.min() << std::endl ;
		std::cout << "max stress33 :" << stress33.max() << std::endl ;
		std::cout << "min stress33 :" << stress33.min() << std::endl ;
		
		std::cout << "max strain11 :" << strain11.max() << std::endl ;
		std::cout << "min strain11 :" << strain11.min() << std::endl ;
		std::cout << "max strain22 :" << strain22.max() << std::endl ;
		std::cout << "min strain22 :" << strain22.min() << std::endl ;
		std::cout << "max strain33 :" << strain33.max() << std::endl ;
		std::cout << "min strain33 :" << strain33.min() << std::endl ;
		
		std::cout << "average stress11 : " << average_stress[0]/total_volume << std::endl ;
		std::cout << "average stress22 : " << average_stress[1]/total_volume << std::endl ;
		std::cout << "average stress33 : " << average_stress[2]/total_volume << std::endl ;
		std::cout << "average strain11 : " << average_strain[0]/total_volume << std::endl ;
		std::cout << "average strain22 : " << average_strain[1]/total_volume << std::endl ;
		std::cout << "average strain33 : " << average_strain[2]/total_volume << std::endl ;
		
		Matrix K(2,2) ;
		K[0][0] = average_strain[0] ;
		K[0][1] = average_strain[1]*2 ;
		K[1][0] = average_strain[1] ;
		K[1][1] = average_strain[0]+average_strain[1] ;
		invert2x2Matrix(K) ;
		K.print() ;

		Vector s(2) ;
		s[0] = average_stress[0] ;
		s[1] = average_stress[1] ;
		
		Vector c = K*s ;
		
		std::string filebench("benchmark.txt") ;
		std::fstream out ;
		out.open(filebench.c_str(), std::ios::out|std::ios::app) ;
		out << "ELASTICITY\t" << str_micro ;
			if(order==2)
				std::cout << "QUAD" ;
			std::cout << "\t" << "E_inc = " << prop << "\t" 
			<< "dof = " << x.size() << "\t"
			<< "C1111 = " << c[0] << "\t"
			<< "C1122 = " << c[1] << std::endl ;
		out.close() ;
		
		break ;
	}
	}
	
	return 0 ;
}
