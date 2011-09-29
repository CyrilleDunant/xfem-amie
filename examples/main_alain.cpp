
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
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/itoa.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../physics/void_form.h"
#include "../physics/homogenization/phase.h"
#include "../physics/homogenization/composite.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../utilities/writer/triangle_writer.h"


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

int nz ;
double fraction ;
int sampling ;
std::vector<Point> centers ;

Matrix getStiffnessTensorAndExpansionStress(std::string hom)
{
	std::cerr << hom << std::endl ;

	Inclusion * agg = new Inclusion(0,0,0) ;
	agg->setBehaviour(new ElasticOnlyAggregateBehaviour());

	Inclusion * gel = new Inclusion(0,0,0) ;
	gel->setBehaviour(new GelBehaviour());

	Phase matrix(agg) ;
	Phase inclusion(gel) ;
	
	Matrix homogenization(sampling,3) ;
	
	MatrixInclusionComposite * scheme = NULL ;

	for(size_t i = 0 ; i < sampling ; i++)
	{
		inclusion.volume = fraction * (double) (i+1) / (double) sampling ;
		matrix.volume = 1. - inclusion.volume ;
	
		if(hom == "voigt")
			scheme = new VoigtMatrixInclusionComposite(matrix, inclusion) ;
		if(hom == "reuss")
			scheme = new ReussMatrixInclusionComposite(matrix, inclusion) ;
		if(hom == "mori-tanaka")
			scheme = new MoriTanakaMatrixInclusionComposite(matrix, inclusion) ;
		if(hom == "inverse-mori-tanaka")
			scheme = new InverseMoriTanakaMatrixInclusionComposite(matrix, inclusion) ;
		if(hom == "dilute")
			scheme = new DiluteMatrixInclusionComposite(matrix, inclusion) ;
		if(hom == "self-consistent")
			scheme = new BiphasicSelfConsistentComposite(matrix, inclusion) ;
			
		std::cout << scheme->matrix.volume << std::endl ;
			
		if(scheme != NULL)
		{
			scheme->apply() ;
			homogenization[i][0] = scheme->C[0][0] ;
			homogenization[i][1] = scheme->C[0][1] ;
			homogenization[i][2] = scheme->beta[0] ;
			
			delete scheme ;
		}
		else
			homogenization.resize(0,0) ;
	}
	
	return homogenization ;
}

Vector getStiffnessTensor(bool random)
{
	ElasticOnlyAggregateBehaviour * aggregate = new ElasticOnlyAggregateBehaviour() ;
	GelBehaviour * gel = new GelBehaviour(22e9,0.28,0) ;

	Sample box(0.05, 0.05, 0.025, 0.025) ;
	box.setBehaviour(aggregate) ;

	std::vector<Inclusion *> inclusions ;
	double interval = 0.05/(nz+1) ;
	double radius = std::sqrt(fraction*(0.05*0.05)/(nz*nz*3.141592)) ;

	for(size_t i = 0 ; i < nz ; i++)
	{
		double cx = interval*(i+1) ;
		for(size_t j = 0 ; j < nz ; j++)
		{
			double cy = interval*(j+1) ;
			inclusions.push_back(new Inclusion(radius, cx, cy)) ;
		}
	}
	std::cerr << inclusions[0]->area()/(0.05*0.05) << std::endl ;
	
	if(random)
	{
		std::vector<Feature *> features ;
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			features.push_back(dynamic_cast<Feature *>(inclusions[i])) ;
		}
		inclusions.clear() ;
		int granulats = 1 ;
		features = placement(dynamic_cast<Rectangle *>(&box), features, &granulats) ;
		
		for(size_t i = 0 ; i < features.size() ; i++)
		{
			inclusions.push_back(dynamic_cast<Inclusion *>(features[i])) ;
			inclusions[i]->getCenter().print() ;
		}
		features.clear() ;
	}

	centers.clear() ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		centers.push_back(inclusions[i]->getCenter()) ;

	for(size_t i = 0 ; i < inclusions.size() ; i++)
		inclusions[i]->setBehaviour(gel) ;

	FeatureTree F(&box) ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		F.addFeature(&box, inclusions[i]) ;

	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA,TOP, -1e6));

	F.step() ;
	Vector disp(F.getDisplacements().size()) ;
	disp = F.getDisplacements() ;

	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	double sxx = 0. ;
	double syy = 0. ;
	double sxy = 0. ;
	double exx = 0. ;
	double eyy = 0. ;
	double exy = 0. ;
	double area = 0. ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		double a = tri[i]->area() ;
		Point p ;
		for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
			p += tri[i]->getBoundingPoint(j) ;
		p /= tri[i]->getBoundingPoints().size() ;
		Vector stress = tri[i]->getState().getStress(p,false) ;
		sxx += stress[0]*a ;
		syy += stress[1]*a ;
		sxy += stress[2]*a ;
		Vector strain = tri[i]->getState().getStrain(p,false) ;
		exx += strain[0]*a ;
		eyy += strain[1]*a ;
		exy += strain[2]*a ;
		area += a ;
	}

	Vector sigma(2) ;
	sigma[0] = sxx/area ;
	sigma[1] = syy/area ;

	Matrix epsilon(2,2) ;
	epsilon[0][0] = exx ;
	epsilon[0][1] = - eyy ;
	epsilon[1][0] = - eyy ;
	epsilon[1][1] = exx ;
	epsilon /= (exx*exx - eyy*eyy) ;
	epsilon *= area ;
	
	return epsilon * sigma ;
}

Vector getExpansionStress(bool random)
{
	ElasticOnlyAggregateBehaviour * aggregate = new ElasticOnlyAggregateBehaviour() ;
	GelBehaviour * gel = new GelBehaviour() ;

	Sample box(0.05, 0.05, 0.025, 0.025) ;
	box.setBehaviour(aggregate) ;

	std::vector<Inclusion *> inclusions ;
	double interval = 0.05/(nz+1) ;
	double radius = std::sqrt(fraction*(0.05*0.05)/(nz*nz*3.141592)) ;
	std::cout << radius << std::endl ;

	for(size_t i = 0 ; i < nz ; i++)
	{
		double cx = interval*(i+1) ;
		for(size_t j = 0 ; j < nz ; j++)
		{
			double cy = interval*(j+1) ;
			inclusions.push_back(new Inclusion(radius, cx, cy)) ;
		}
	}
	std::cerr << inclusions[0]->area()/(0.05*0.05) << std::endl ;

/*	if(random)
	{
		std::vector<Feature *> features ;
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			features.push_back(dynamic_cast<Feature *>(inclusions[i])) ;
		}
		inclusions.clear() ;
		int granulats = 1 ;
		features = placement(dynamic_cast<Rectangle *>(&box), features, &granulats) ;
		
		for(size_t i = 0 ; i < granulats ; i++)
		{
			inclusions.push_back(dynamic_cast<Inclusion *>(features[i])) ;
			inclusions[i]->getCenter().print() ;
		}
		features.clear() ;
	}*/

	for(size_t i = 0 ; i < inclusions.size() ; i++)
	{
		inclusions[i]->setBehaviour(gel) ;
		Point c = centers[i] ;
		dynamic_cast<Circle *>(inclusions[i])->setCenter(c) ;
	}

	FeatureTree F(&box) ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		F.addFeature(&box, inclusions[i]) ;

	F.setSamplingNumber(sampling) ;
	F.setOrder(LINEAR) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,RIGHT));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,TOP));

	F.step() ;
	Vector disp(F.getDisplacements().size()) ;
	disp = F.getDisplacements() ;

	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	double sxx = 0. ;
	double syy = 0. ;
	double sxy = 0. ;
	double area = 0. ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		double a = tri[i]->area() ;
		Point p ;
		for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
			p += tri[i]->getBoundingPoint(j) ;
		p /= tri[i]->getBoundingPoints().size() ;
		Vector stress = tri[i]->getState().getStress(p,false) ;
		sxx += stress[0]*a ;
		syy += stress[1]*a ;
		sxy += stress[2]*a ;
		area += a ;
	}

	Vector sigma(2) ;
	sigma[0] = sxx/area ;
	sigma[1] = syy/area ;
	return sigma ;

}

int main(int argc, char *argv[])
{
	Function f1 = "0.5 0.5 x * - 0.5 0.5 t * - *" ;
	Function f2 = "0.5 0.5 x * + 0.5 0.5 t * - *" ;
	Function f3 = "0.5 0.5 x * - 0.5 0.5 t * + *" ;
	Function f4 = "0.5 0.5 x * + 0.5 0.5 t * + *" ;
	
	VirtualMachine vm ;

	Function fx1 = "0.5 0.5 t * - -0.5 *" ;
	Function fx2 = "0.5 0.5 t * - 0.5 *" ;
	Function fx3 = "0.5 0.5 t * + -0.5 *" ;
	Function fx4 = "0.5 0.5 t * + 0.5 *" ;

	Matrix A(4,4) ;
	Matrix B(4,4) ;
	
	A[0][0] =  vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	A[0][1] =  vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	A[0][2] =  vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	A[0][3] =  vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;
	
	A[1][0] =  vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	A[1][1] =  vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	A[1][2] =  vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	A[1][3] =  vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;

	A[2][0] =  vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	A[2][1] =  vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	A[2][2] =  vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	A[2][3] =  vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;

	A[3][0] =  vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	A[3][1] =  vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	A[3][2] =  vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	A[3][3] =  vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;

// 	 std::endl ;

	B[0][0] =  vm.deval(f1, XI, 0,0,0,0)*vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[0][1] =  vm.deval(f1, XI, 0,0,0,0)*vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[0][2] =  vm.deval(f1, XI, 0,0,0,0)*vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[0][3] =  vm.deval(f1, XI, 0,0,0,0)*vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0) ;
	
	B[1][0] =  vm.deval(f2, XI, 0,0,0,0)*vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[1][1] =  vm.deval(f2, XI, 0,0,0,0)*vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[1][2] =  vm.deval(f2, XI, 0,0,0,0)*vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[1][3] =  vm.deval(f2, XI, 0,0,0,0)*vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0) ;

	B[2][0] =  vm.deval(f3, XI, 0,0,0,0)*vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[2][1] =  vm.deval(f3, XI, 0,0,0,0)*vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[2][2] =  vm.deval(f3, XI, 0,0,0,0)*vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[2][3] =  vm.deval(f3, XI, 0,0,0,0)*vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0) ;

	B[3][0] =  vm.deval(f4, XI, 0,0,0,0)*vm.ddeval(f1, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[3][1] =  vm.deval(f4, XI, 0,0,0,0)*vm.ddeval(f2, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[3][2] =  vm.deval(f4, XI, 0,0,0,0)*vm.ddeval(f3, TIME_VARIABLE,XI, 0,0,0,0)  ;
	B[3][3] =  vm.deval(f4, XI, 0,0,0,0)*vm.ddeval(f4, TIME_VARIABLE,XI, 0,0,0,0) ;
	
	((Matrix)(A+B)).print() ; 
	
	Matrix C(4,4) ;
	C[0][0] = vm.deval(f1, XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	C[0][1] =  vm.deval(f1, XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	C[0][2] =  vm.deval(f1, XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	C[0][3] =  vm.deval(f1, XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;
	
	C[1][0] =  vm.deval(f2, XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	C[1][1] =  vm.deval(f2, XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	C[1][2] =  vm.deval(f2, XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	C[1][3] =  vm.deval(f2, XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;

	C[2][0] =  vm.deval(f3, XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	C[2][1] =  vm.deval(f3, XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	C[2][2] =  vm.deval(f3, XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	C[2][3] =  vm.deval(f3, XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;

	C[3][0] =  vm.deval(f4, XI, 0,0,0,0)*vm.deval(f1, XI, 0,0,0,0)  ;
	C[3][1] =  vm.deval(f4, XI, 0,0,0,0)*vm.deval(f2, XI, 0,0,0,0)  ;
	C[3][2] =  vm.deval(f4, XI, 0,0,0,0)*vm.deval(f3, XI, 0,0,0,0)  ;
	C[3][3] =  vm.deval(f4, XI, 0,0,0,0)*vm.deval(f4, XI, 0,0,0,0) ;
	
	((Matrix)(C+B+B.transpose())).print() ; 
	
	
	
	
	
	return 0 ;
	
	

/*	std::string type(argv[1]) ;
	if(type == "--regular")
	{
		std::cout << "first argument is number of zones along a side" << std::endl ;
		std::cout << "second argument is percentage of area covered by all the zones" << std::endl ;
		std::cout << "third argument is sampling number" << std::endl ;

		nz = atoi(argv[2]) ;
		fraction = atof(argv[3]) ;
		sampling = atoi(argv[4]) ;

		Vector stiffness = getStiffnessTensor(false) ;
		Vector expansion = getExpansionStress(false) ;
	
		std::fstream out ;
		out.open("grid_regular", std::ios::out | std::ios::app) ;
		out << "finite-elements-" << sampling << "\t" << fraction << "\t" << stiffness[0] << "\t" << stiffness[1] << "\t" << -expansion[0] << std::endl ;
		out.close() ;

		return 0 ;
	}
	
	if(type == "--random")
	{
		std::cout << "first argument is number of zones along a side" << std::endl ;
		std::cout << "second argument is percentage of area covered by all the zones" << std::endl ;
		std::cout << "third argument is sampling number" << std::endl ;

		nz = atoi(argv[2]) ;
		fraction = atof(argv[3]) ;
		sampling = atoi(argv[4]) ;

		Vector stiffness = getStiffnessTensor(true) ;
		Vector expansion = getExpansionStress(true) ;
	
		std::fstream out ;
		out.open("grid_random", std::ios::out | std::ios::app) ;
		out << "finite-elements-" << sampling << "\t" << fraction << "\t" << stiffness[0] << "\t" << stiffness[1] << "\t" << -expansion[0] << std::endl ;
		out.close() ;

		return 0 ;
	}

	return 0 ;*/

}
