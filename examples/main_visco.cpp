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
#include "../physics/maxwell.h"
#include "../physics/stiffness.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/generalized_spacetime_viscoelasticity.h"
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
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <sys/time.h>
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

Sample box(nullptr, 0.1, 0.4,0.0,0.0) ;
Rectangle rbox(0.07,0.07,0.0,0.0) ;
double tau = 0.1;

double s2d(double s)
{
	return s/(24*60*60) ;
}

double getMinDisplacement(const Vector & x)
{
	double ret = x[0] ;
	for(size_t i = 0 ; i < x.size()/6 ; i++)
	{
		if(x[i*6] < ret)
			ret = x[i*6] ;
		if(x[i*6+1] < ret)
			ret = x[i*6+1] ;
	}
	return ret ;
}

int main(int argc, char *argv[])
{
	tau = atof(argv[1]) ;
//	ElasticOnlyPasteBehaviour * paste_e = new ElasticOnlyPasteBehaviour() ;
  
	FeatureTree F(&box) ;
	F.setSamplingNumber(2) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	
// 	Sample Gbox(nullptr, 0.1,0.4,0.,0.) ;
// 	FeatureTree G(&Gbox) ;
// 	G.setSamplingNumber(2) ;
// 	G.setOrder(LINEAR) ;
// 	G.setDeltaTime(tau) ;
// 
// 	Sample Hbox(nullptr, 0.1,0.4,0.,0.) ;
// 	FeatureTree H(&Hbox) ;
// 	H.setSamplingNumber(2) ;
// 	H.setOrder(LINEAR) ;
// 	H.setDeltaTime(tau) ;


	Matrix ckv = Material::cauchyGreen(30e9,0.3,true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS) ;
	Matrix ekv = ckv*300 ;
	Matrix cmx = Material::cauchyGreen(14e9,0.3,true,SPACE_TWO_DIMENSIONAL, PLANE_STRESS) ;
	Matrix emx = cmx*5000 ;
	
	GeneralizedSpaceTimeViscoelasticity * burger =  new GeneralizedSpaceTimeViscoelasticity(BURGER,ckv,ekv,cmx,emx) ;
	
	double density = 2300 ;
	Matrix stiffness = burger->param ;
	Matrix viscosity = burger->eta ;
	
	double tausquare = tau*tau ;
	
	box.setBehaviour(burger) ;
	
// 	box.setBehaviour( new MassAndViscosityAndStiffnessByBlock(stiffness, viscosity*(1.5/tau), density*(1./tausquare) ,3)) ;	
// 	Gbox.setBehaviour( new MassAndViscosityAndStiffnessByBlock(stiffness*0., viscosity*(-2./tau), density*(-2./tausquare) ,3)) ;	
// 	Hbox.setBehaviour( new MassAndViscosityAndStiffnessByBlock(stiffness*0., viscosity*(0.5/tau), density*(1./tausquare),3)) ;


	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,5)) ;
	F.addBoundaryCondition(new TimeContinuityBoundaryCondition()) ;
	srand(0) ;
	F.step() ;
// 	G.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT, 0,0)) ;
// 	G.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM, 0,1)) ;
// 	srand(0) ;
// 	G.step() ;
// 	H.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT, 0,0)) ;
// 	H.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM, 0,1)) ;
// 	srand(0) ;
// 	H.step() ;

// 	Vector uprev = F.getDisplacements() ;
// 	Vector unow = uprev ;
// 	Vector forces = unow ;
// 	
// 	forces = (Vector) (G.getAssembly()->getMatrix() * unow) + (Vector) (H.getAssembly()->getMatrix() * unow) ;
// 	forces *= -1 ;
// 	
// 	GlobalForceBoundaryCondition * bc = new GlobalForceBoundaryCondition(forces) ;
// 	F.addBoundaryCondition(bc) ;
	F.step() ;
 	size_t i = 0 ;
	BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -7e6) ;
	F.addBoundaryCondition(stress) ;

	std::string toto = "fit_bengougam_fem3_" ;
	toto.append(argv[1]) ;
	std::fstream out(toto.c_str(), std::ios::out) ;
	
	double time = 0. ;
	F.step() ;
	
	while(time < 3500)
	{
		i++ ;
		F.setDeltaTime(tau) ;
// 		uprev = unow ;
// 		unow = F.getDisplacements() ;
		std::cout << time << /*"\t" << unow.min()/(0.4*-7e-6) <<*/ "\t" << F.getAverageField(STRAIN_FIELD)[1]/(-7e-6) << std::endl ;
// 		forces = (Vector) (G.getAssembly()->getMatrix() * unow) + (Vector) (H.getAssembly()->getMatrix() * uprev) ;
// 		forces *= -1 ;
// 		bc->setDataVector(forces) ;
		F.step() ;
//		F.getAssembly()->print() ;
		time += tau ;
	}
	
		
	return 0 ;
}

