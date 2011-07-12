
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


int main(int argc, char *argv[])
{
	Inclusion * ag = new Inclusion(0,0,0) ;
	ag->setBehaviour(new ElasticOnlyAggregateBehaviour(/*59e9, 0.3, SPACE_THREE_DIMENSIONAL*/));

	Inclusion * gl = new Inclusion(0,0,0) ;
	gl->setBehaviour(new GelBehaviour(/*22e9, 0.28, 0.5, SPACE_THREE_DIMENSIONAL*/));

	Phase matrix(ag) ;
	Phase inclusion(gl) ;

	std::vector<double> mean ;
	std::vector<double> inv_mean ;
	std::vector<double> inter ;
	std::vector<double> inv_inter ;
	std::vector<double> sc ;
	std::vector<double> dil ;

	Matrix m(3,3) ;
	Vector v(3) ;
	for(size_t i = 1 ; i < 11 ; i++)
	{
		matrix.volume = 1. - (double) i/1000. ;
		inclusion.volume = 1. ;

		VoigtMatrixInclusionComposite voigt(matrix, inclusion) ;
		voigt.apply() ;
		m = voigt.C ;
		mean.push_back(m[0][0]) ;
		mean.push_back(m[0][1]) ;
		invert3x3Matrix(m) ;
		v = m * voigt.beta ;
		mean.push_back(voigt.beta[0]) ;


		matrix.volume = 1. - (double) i/1000. ;
		inclusion.volume = 1. ;

		ReussMatrixInclusionComposite reuss(matrix, inclusion) ;
		reuss.apply() ;
		m = reuss.C ;
		inv_mean.push_back(m[0][0]) ;
		inv_mean.push_back(m[0][1]) ;
		invert3x3Matrix(m) ;
		v = m * reuss.beta ;
		inv_mean.push_back(reuss.beta[0]) ;


/*		matrix.volume = 1. - (double) i/100. ;
		inclusion.volume = 1. ;

		MoriTanakaMatrixInclusionComposite mori(matrix, inclusion) ;
		mori.apply() ;
		m = mori.C ;
//		invert3x3Matrix(m) ;
		v = m * mori.beta ;

		inter.push_back(m[0][0]) ;
		inter.push_back(m[0][1]) ;

		matrix.volume = 1. - (double) i/100. ;
		inclusion.volume = 1. ;

		InverseMoriTanakaMatrixInclusionComposite inv_mori(matrix, inclusion) ;
		inv_mori.apply() ;
		m = inv_mori.C ;
//		invert3x3Matrix(m) ;
		v = m * inv_mori.beta ;

		inv_inter.push_back(m[0][0]);
		inv_inter.push_back(m[0][1]) ;

		matrix.volume = 1. - (double) i/100. ;
		inclusion.volume = 1. ;

		BiphasicSelfConsistentComposite self(matrix, inclusion) ;
		self.apply() ;
		m = self.C ;
//		invert3x3Matrix(m) ;
		v = m * self.beta ;

		sc.push_back(m[0][0]);
		sc.push_back(m[0][1]) ;

		matrix.volume = 1. - (double) i/100. ;
		inclusion.volume = 1. ;

		DiluteMatrixInclusionComposite dilute(matrix, inclusion) ;
		dilute.apply() ;
		m = dilute.C ;
//		invert3x3Matrix(m) ;
		v = m * dilute.beta ;

		dil.push_back(m[0][0]);
		dil.push_back(m[0][1]) ;*/


	}

	for(size_t i = 0 ; i < mean.size()/3 ; i++)
		std::cout << (double) (i+1)/1000 << "\t" << mean[3*i] << "\t" << mean[3*i+1] << "\t" << mean[3*i+2] << std::endl ;

	for(size_t i = 0 ; i < mean.size()/3 ; i++)
		std::cout << (double) (i+1)/1000 << "\t" << inv_mean[3*i] << "\t" << inv_mean[3*i+1] << "\t" << inv_mean[3*i+2] << std::endl ;

/*	for(size_t i = 0 ; i < mean.size()/2 ; i++)
		std::cout << (double) (i+1)/100 << "\t" << inter[2*i] << "\t" << inter[2*i+1] << std::endl ;

	for(size_t i = 0 ; i < mean.size()/2 ; i++)
		std::cout << (double) (i+1)/100 << "\t" << inv_inter[2*i] << "\t" << inv_inter[2*i+1] << std::endl ;

	for(size_t i = 0 ; i < sc.size()/2 ; i++)
		std::cout << (double) (i+1)/100 << "\t" << sc[2*i] << "\t" << sc[2*i+1] << std::endl ;

	for(size_t i = 0 ; i < mean.size()/2 ; i++)
		std::cout << (double) (i+1)/100 << "\t" << dil[2*i] << "\t" << dil[2*i+1] << std::endl ;*/

//	return 0 ;

	std::cout << "first argument is number of zones along a side" << std::endl ;
	std::cout << "second argument is percentage of area covered by all the zones" << std::endl ;
	std::cout << "third argument is sampling number" << std::endl ;

	int nz = atoi(argv[1]) ;
	double area = atof(argv[2]) ;

	ElasticOnlyAggregateBehaviour * aggregate = new ElasticOnlyAggregateBehaviour() ;
	GelBehaviour * gel = new GelBehaviour() ;

	Sample box(0.05, 0.05, 0.025, 0.025) ;
	box.setBehaviour(aggregate) ;

	std::vector<Inclusion *> inclusions ;
	double interval = 0.05/(nz+1) ;
	double radius = std::sqrt(area*(0.05*0.05)/(nz*nz*3.141592)) ;

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

	for(size_t i = 0 ; i < inclusions.size() ; i++)
		inclusions[i]->setBehaviour(gel) ;

	FeatureTree F(&box) ;
	for(size_t i = 0 ; i < inclusions.size() ; i++)
		F.addFeature(&box, inclusions[i]) ;

	F.setSamplingNumber(atoi(argv[3])) ;
	F.setOrder(LINEAR) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,LEFT));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,BOTTOM));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI,RIGHT));
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA,TOP));

	F.step() ;
	Vector disp(F.getDisplacements().size()) ;
	disp = F.getDisplacements() ;

	std::vector<DelaunayTriangle *> tri = F.getElements2D() ;
	double sxx = 0. ;
	double syy = 0. ;
	double sxy = 0. ;
	area = 0. ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		double a = tri[i]->area() ;
		Point p ;
		for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
			p += tri[i]->getBoundingPoint(j) ;
		p /= tri[i]->getBoundingPoints().size() ;
		Vector stress = tri[i]->getState().getStress(p,false) ;
		sxx -= stress[0]*a ;
		syy -= stress[1]*a ;
		sxy -= stress[2]*a ;
		area += a ;
	}

	std::cout << sxx/area << "\t" << syy /area<< std::endl ;

	FeatureTree * f = &F ;
	TriangleWriter writer("expansive_zones_"+itoa(nz)+"x"+itoa(nz)+"_"+itoa((int)((double) 100*area))+"_"+argv[3], f) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.getField(TWFT_STRAIN) ;
	writer.getField(TWFT_STRESS) ;
	writer.write() ;

	return 0 ;

}
