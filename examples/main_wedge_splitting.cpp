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
#include "../physics/stiffness_and_fracture.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/damagemodels/fiberbasedisotropiclineardamage.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
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

// see Cecot 2001, p98, for sample dimensions

double length = 0.2 ; //5.5 ;
double width = 0.06 ;
double depth = 0.03 ;//2.5 ;
double nnotch = 0.085 ;
double nwidth = 0.0025 ;
Sample box(nullptr, length, length,0.0,0.0) ;
Sample top(nullptr, width, depth, 0.0, length*0.5 - depth*0.5) ;
Sample notch(nullptr,nwidth, nnotch, 0.0, length*0.5 - depth - nnotch*0.5) ;

int main(int argc, char *argv[])
{
	FeatureTree F(&box) ;
	F.setSamplingNumber(atof(argv[1])) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(atof(argv[2])) ;
	double speed = atof(argv[3]) ;
//	F.setOrder(LINEAR) ;
	
	
	Rectangle refinement( 0.005, length, 0.,0.) ;
	F.addRefinementZone(&refinement);
	
	Matrix c = (new PasteBehaviour())->param ;
	
	box.setBehaviour( new ViscoelasticityAndFracture(GENERALIZED_KELVIN_VOIGT, c, c*0.3, c*0.3*10, new SpaceTimeNonLocalMohrCoulomb(0.001, -0.008, 15e9), new SpaceTimeFiberBasedIsotropicLinearDamage() ) ) ;
//	box.setBehaviour( new Viscoelasticity(PURE_ELASTICITY, c) ) ;
//	box.setBehaviour( new StiffnessAndFracture( c, new NonLocalMohrCoulomb( 0.001, -0.008, 15e9) ) ) ;
//	box.setBehaviour( new Stiffness( c ) ) ;
	box.getBehaviour()->getFractureCriterion()->setMaterialCharacteristicRadius(0.003);	
	top.setBehaviour( new VoidForm() ) ;
	notch.setBehaviour( new VoidForm() ) ;

	ElasticOnlyAggregateBehaviour hop ;
	Viscoelasticity * agg = new Viscoelasticity( PURE_ELASTICITY, hop.param, 1 ) ;
	
	/*std::vector<Inclusion *> inclusions = */ParticleSizeDistribution::get2DConcrete( &F, agg, 400, 0.008, 0.000001) ;
 	F.addFeature(&box, &top) ;
 	F.addFeature(&box, &notch) ;
	
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0 )) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2 )) ;
 
	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 1 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 3 ) ) ;

	BoundingBoxAndRestrictionDefinedBoundaryCondition * disp = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, width*0.4, width*0.6, length*0.4, length*0.6, 0. ) ;
	F.addBoundaryCondition(disp) ;

	
	F.step() ;
 	Vector x = F.getAverageField(STRAIN_FIELD) ;
 	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
//	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	size_t i = 0 ;

	std::fstream out ;
	std::string tata = "wedge_" ;
	tata.append(argv[1]) ;
	tata.append("_") ;
	tata.append(argv[2]) ;
	tata.append("_") ;
	tata.append(argv[3]) ;
	tata.append(".txt") ;
	out.open(tata.c_str(), std::ios::out) ;
	
	while(i < 300)
	{
		i++ ;
		disp->setData( speed*i ) ;

//		F.setDeltaTime(1.) ;
		F.step() ;
		
		std::string tati = "wedge" ;
		tati.append("_") ;
		tati.append(itoa(i)) ;
		std::cout << tati << std::endl ;
		TriangleWriter writer(tati, &F, 1) ;
// 		writer.getField(TWFT_STRAIN) ;
// 		writer.getField(TWFT_STRESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.write() ;
			
 		x = F.getAverageField(STRAIN_FIELD, -1, 1) ;
 		y = F.getAverageField(REAL_STRESS_FIELD, -1, 1) ;
		out << F.getCurrentTime() << "\t" << speed*i << "\t" << x[0] << "\t" << y[0]*length*length<< "\t" << F.averageDamage << std::endl ;

	}
	
//	F.getAssembly()->print() ;
	
	return 0 ;
}

