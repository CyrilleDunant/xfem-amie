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
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
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

Sample box(nullptr, 0.08, 0.08,0.,0.) ;


int main(int argc, char *argv[])
{
	double timestep = 1. ;//atof(argv[1]) ;
	FeatureTree F(&box) ;
	F.setSamplingNumber(64) ;

	ViscoElasticOnlyPasteBehaviour dpaste ;
	ViscoelasticityAndFracture paste( PURE_ELASTICITY, dpaste.param, new SpaceTimeNonLocalMaximumStrain( 0.0001, 0.0001*dpaste.param[0][0]), new SpaceTimeFiberBasedIsotropicLinearDamage( 0.1, atof(argv[1]) ) ) ;
	paste.getFractureCriterion()->setMaterialCharacteristicRadius(0.0005) ;
	box.setBehaviour(&paste);

	ViscoElasticOnlyAggregateBehaviour dagg ;
	ParticleSizeDistribution::get2DConcrete(&F, new Viscoelasticity(PURE_ELASTICITY, dagg.param), 500, 0.006, 0.0004, BOLOME_A, CIRCLE, 1., M_PI, 50000, 0.4) ;
	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;
	F.setMinDeltaTime(atof(argv[1])) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
	F.step() ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0.0001, 1) ;
	F.addBoundaryCondition(disp) ;

 	
	std::fstream out ;
	std::string name = "tensile_" ;
	name.append(argv[1]) ;
	
	out.open(name.c_str(), std::ios::out) ;

	double time = timestep ;	
	Vector stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	Vector strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	Vector rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	Vector displ = F.getDisplacements() ;
	out << std::setprecision(16) << 0. << "\t" << displ.max() << "\t" << displ.min() << "\t" <<  stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
		"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
		"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	std::string nametrg = name ;
	nametrg.append("_trg_0") ;
	TriangleWriter writer(nametrg, &F, 1) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(REAL_STRESS_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.getField(TWFT_DAMAGE) ;
	writer.write() ;

	bool goOn = false ;
	int i = 0 ;
	
	while(!goOn)
	{
		i++ ;
		goOn = F.step() ;

		stress = F.getAverageField(REAL_STRESS_FIELD,-1,-1) ;
		strain = F.getAverageField(STRAIN_FIELD,-1,-1) ;
		rate = F.getAverageField(STRAIN_RATE_FIELD,-1,-1) ;
		displ = F.getDisplacements() ;
		out << std::setprecision(16) << (F.getElements2D()[0])->getBoundingPoint(0).t << "\t" << displ.max() << "\t" << displ.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;
		
		std::string nametrg = name ;
		nametrg.append("_trg_") ;
		nametrg.append(std::to_string((int) i)) ;
		TriangleWriter writer(nametrg, &F, -1) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.getField(TWFT_DAMAGE) ;
		writer.write() ;
	}

	stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	displ = F.getDisplacements() ;
	out << std::setprecision(16) << 1. << "\t" << displ.max() << "\t" << displ.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	std::string nametrg2 = name ;
	nametrg2.append("_trg_") ;
	nametrg2.append(std::to_string((int) i+1)) ;
	TriangleWriter writer2(nametrg2, &F, 1) ;
	writer2.getField(STRAIN_FIELD) ;
	writer2.getField(REAL_STRESS_FIELD) ;
	writer2.getField(TWFT_STIFFNESS) ;
	writer2.getField(TWFT_DAMAGE) ;
	writer2.write() ;

		
	return 0 ;
}
